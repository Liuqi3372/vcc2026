import argparse
import fcntl
import json
import time
from pathlib import Path

import h5py
import numba as nb
import numpy as np
import pandas as pd
from anndata._io.specs import read_elem

ROOT = Path(__file__).resolve().parents[1]
PUB = Path('/ssd1/PubData/arc_cviplus_panel/by_cell_line')


@nb.njit(cache=True)
def accumulate(values, indices, pointers, groups, batches, gene_map, nt,
               sums, cpms, batch_n, batch_lib):
    total_projected = 0.0
    for r in range(len(groups)):
        group = groups[r]
        if group < 0:
            continue
        left, right = pointers[r], pointers[r + 1]
        depth = 0.0
        for k in range(left, right):
            v = values[k]
            if not np.isfinite(v) or v < 0 or v != np.floor(v):
                raise ValueError('Source must contain finite nonnegative integer counts')
            depth += v
        if depth <= 0:
            raise ValueError('Zero-depth source cell')
        projected = 0.0
        for k in range(left, right):
            j = gene_map[indices[k]]
            if j >= 0:
                v = values[k]
                sums[group, j] += v
                cpms[group, j] += 1000000.0 * v / depth
                projected += v
        total_projected += projected
        if group < nt:
            batch_n[group, batches[r]] += 1
            batch_lib[group, batches[r]] += projected
    return total_projected


def prepare(source):
    started = time.time()
    data = ROOT / 'data'
    nt_targets = pd.read_csv(data / 'pert_counts.csv').target_gene.astype(str).tolist()
    h1 = pd.read_csv(ROOT / 'raw/pert_counts_Training.csv').target_gene.astype(str).tolist()
    targets = np.asarray(list(dict.fromkeys(nt_targets + [t for t in h1 if t != 'non-targeting'])))
    is_k562 = source == 'K562_GWPS_CPM'
    unit = 'Replogle_K562_gwps' if is_k562 else 'XAtlas_' + source
    path = PUB / unit / '01_preprocess/adata_processed.h5ad'
    output = data / (source + '_full_statistics.npz')
    if output.exists():
        print('Already prepared', output, flush=True)
        return
    with h5py.File(path, 'r') as f:
        obs, var = read_elem(f['obs']), read_elem(f['var'])
        genes0 = var.gene_name.astype(str).to_numpy() if 'gene_name' in var else var.index.astype(str).to_numpy()
        if is_k562:
            genes = np.asarray(list(dict.fromkeys(genes0)), dtype=str)

            batch_labels = obs.index.to_series().astype(str).str.rsplit('-', n=1).str[-1]
            if not batch_labels.str.fullmatch(r'\d+').all():
                raise ValueError('Cannot reconstruct K562 gem_group from barcode suffix')
        else:
            official_genes = pd.read_csv(data / 'gene_names.csv').gene_name.astype(str).tolist()
            h1_genes = pd.read_csv(ROOT / 'raw/gene_names.csv')
            col = 'gene_name' if 'gene_name' in h1_genes else h1_genes.columns[0]
            genes = np.asarray(list(dict.fromkeys(official_genes + h1_genes[col].astype(str).tolist())))
            batch_labels = obs['sample'].astype(str)
        batches, batch_names = pd.factorize(batch_labels, sort=False)
        batches = batches.astype(np.int32)
        nt, nbat = len(targets), len(batch_names)
        lookup = {g: i for i, g in enumerate(genes)}
        gene_map = np.asarray([lookup.get(g, -1) for g in genes0], dtype=np.int32)
        tl = {t: i for i, t in enumerate(targets)}
        groups = np.asarray([tl.get(t, -1) for t in obs.perturbation.astype(str)], dtype=np.int32)
        control = obs.is_control.to_numpy(dtype=bool)

        selected_control = control.copy()
        if not is_k562:
            for b in range(nbat):
                rows = np.flatnonzero(control & (batches == b))
                selected_control[rows[250:]] = False
        groups[control] = -1
        groups[selected_control] = nt + batches[selected_control]
        n = np.bincount(groups[groups >= 0], minlength=nt + nbat)
        sums = np.zeros((nt + nbat, len(genes)), dtype=np.float64)
        cpms = np.zeros_like(sums)
        bn = np.zeros((nt, nbat), dtype=np.int64)
        bl = np.zeros((nt, nbat), dtype=np.float64)
        x = f['layers/counts']
        if x.attrs['encoding-type'] != 'csr_matrix':
            raise ValueError('Expected CSR raw count layer')
        ptr = x['indptr'][:].astype(np.int64)
        projected_total = 0.0
        print(source, 'source cells', len(obs), 'selected', int(n.sum()), 'genes', len(genes), 'batches', nbat, flush=True)

        for left in range(0, len(obs), 4096):
            right = min(left + 4096, len(obs))
            g = groups[left:right]
            if not (g >= 0).any():
                continue
            selected_rows = np.flatnonzero(g >= 0) + left
            runs = []


            gap = max(x['data'].chunks[0], x['indices'].chunks[0])
            for row in selected_rows:
                if runs and ptr[row] - ptr[runs[-1][1]] <= gap:
                    runs[-1][1] = int(row) + 1
                else:
                    runs.append([int(row), int(row) + 1])
            for first, last in runs:
                lo, hi = int(ptr[first]), int(ptr[last])
                values = x['data'][lo:hi]
                indices = x['indices'][lo:hi]
                projected_total += accumulate(values, indices, ptr[first:last + 1] - lo,
                                              groups[first:last], batches[first:last], gene_map, nt,
                                              sums, cpms, bn, bl)
            if left % (4096 * 32) == 0 or right == len(obs):
                print(source, right, '/', len(obs), 'seconds', round(time.time() - started), flush=True)
        np.testing.assert_allclose(sums.sum(), projected_total, rtol=1e-12)
        np.testing.assert_array_equal(bn.sum(axis=1), n[:nt])
        control_counts, control_n = sums[nt:], n[nt:]
        if control_n.sum() == 0:
            raise ValueError('No controls')
        global_cp = control_counts.sum(axis=0) / control_counts.sum()
        global_cm = cpms[nt:].sum(axis=0) / control_n.sum()
        ctrl_cp = np.tile(global_cp, (nbat, 1))
        ctrl_cm = np.tile(global_cm, (nbat, 1))
        valid = control_n > 0
        ctrl_cp[valid] = control_counts[valid] / control_counts[valid].sum(axis=1, keepdims=True)
        ctrl_cm[valid] = cpms[nt:][valid] / control_n[valid, None]
        cw = np.divide(bn, n[:nt, None], out=np.zeros_like(bl), where=n[:nt, None] > 0)
        lw = np.divide(bl, bl.sum(axis=1, keepdims=True), out=np.zeros_like(bl), where=bl.sum(axis=1, keepdims=True) > 0)
        means = np.divide(cpms[:nt], n[:nt, None], out=np.zeros_like(sums[:nt]), where=n[:nt, None] > 0)
        result = dict(source=np.asarray(source + '_local_QC'), targets=targets, genes=genes,
                      n_cells=n[:nt], target_count_sums=sums[:nt].astype(np.float32),
                      target_mean_cpm=means.astype(np.float32),
                      matched_control_probability=(lw @ ctrl_cp).astype(np.float32),
                      matched_control_mean_cpm=(cw @ ctrl_cm).astype(np.float32),
                      global_control_probability=global_cp.astype(np.float32),
                      global_control_mean_cpm=global_cm.astype(np.float32),
                      measured_genes=np.isin(np.arange(len(genes)), gene_map[gene_map >= 0]))
        assert np.isfinite(means).all()
        np.testing.assert_allclose(result['matched_control_probability'][n[:nt] > 0].sum(axis=1), 1, rtol=1e-5)
        tmp = output.with_suffix('.partial.npz')
        np.savez_compressed(tmp, **result)
        tmp.rename(output)
        qc = read_elem(f['uns/cell_qc'])
        report = dict(source=source, input=str(path), layer='counts', source_cells=len(obs),
                      selected_cells=int(n.sum()), control_cells=int(control_n.sum()),
                      source_gene_count=len(genes0), output_gene_count=len(genes),
                      measured_gene_count=int(result['measured_genes'].sum()),
                      targets=len(targets), targets_with_20_cells=int((n[:nt] >= 20).sum()),
                      batches_without_controls=int((~valid).sum()),
                      batch_definition='barcode suffix (gem_group)' if is_k562 else 'sample',
                      control_selection='all post-QC controls' if is_k562 else 'first 250 post-QC rows per sample',
                      cpm_denominator='row sum over retained source genes',
                      differences='Post-QC local source; K562 has a filtered gene panel; no remote X-Atlas cell_integer_id ordering.',
                      original_local_qc={k:v.item() if isinstance(v,np.generic) else v for k,v in qc.items()},
                      elapsed_seconds=time.time()-started)
        (ROOT / 'logs' / (source + '_preparation.json')).write_text(json.dumps(report, indent=2))
        print('SAVED', output, 'seconds', round(time.time()-started), flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--source', choices=['K562_GWPS_CPM', 'HCT116', 'HEK293T'], required=True)
    source = p.parse_args().source
    with (ROOT / 'logs' / (source + '.lock')).open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        prepare(source)
