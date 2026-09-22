import json
import sys
from pathlib import Path
import h5py
import numpy as np
import pandas as pd
from anndata._io.specs import read_elem

r = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(r / 'src'))
from predict import SOURCES, axis_csv
from model import load_xatlas, aligned_effect, aligned_cd4

genes = axis_csv(r/'data/gene_names.csv', 'gene_name')
targets = axis_csv(r/'data/pert_counts.csv', 'target_gene')
assert (len(targets), len(genes)) == (300, 18533)
report = {'contexts': {}, 'sources': {}, 'n_targets': len(targets), 'n_genes': len(genes)}
coverage = pd.DataFrame({'target_gene': targets})
for name in SOURCES:
    source = load_xatlas(r/'data'/name)
    for a in [source.probability, source.control_probability, source.mean_cpm, source.control_mean_cpm]:
        assert np.isfinite(a).all() and (a >= 0).all(), name
    effect, mask = aligned_effect(source, targets, genes, space='log2fc')
    assert np.isfinite(effect).all()
    coverage[source.name] = mask.sum(axis=1)
    report['sources'][source.name] = {'targets_with_20_cells': int(mask.any(axis=1).sum()),
                                   'covered_official_genes': int(mask.any(axis=0).sum()),
                                   'source_targets': len(source.targets), 'source_genes': len(source.genes)}
    print('CHECKED SOURCE', name, report['sources'][source.name], flush=True)
effect, mask = aligned_cd4(r/'data/CD4_DE_statistics.npz',targets,genes,center_scope='source')
assert np.isfinite(effect).all()
coverage['CD4_DE'] = mask.sum(axis=1)
report['sources']['CD4_DE']={'targets_passing_publisher_qc':int(mask.any(axis=1).sum()),'covered_official_genes':int(mask.any(axis=0).sum())}
for c in 'ABC':
    with h5py.File(r/f'data/context_{c}.h5ad') as f:
        obs,var=read_elem(f['obs']),read_elem(f['var'])
        np.testing.assert_array_equal(var.index.astype(str),genes)
        assert len(obs)>=1600 and obs.index.is_unique
        assert (obs.context.astype(str)==c).all()
        assert (obs.target_gene.astype(str)=='non-targeting').all()
        x=f['X'];assert x.attrs['encoding-type']=='csr_matrix'
        ptr=x['indptr'][:];depths=[]
        for left in range(0,len(obs),1024):
            right=min(left+1024,len(obs));lo,hi=int(ptr[left]),int(ptr[right]);v=x['data'][lo:hi]
            idx=x['indices'][lo:hi]
            assert np.isfinite(v).all() and (v>=0).all() and (v==np.floor(v)).all()
            assert (idx>=0).all() and (idx<len(genes)).all()

            cs=np.r_[0,np.cumsum(v,dtype=np.float64)]
            depths.extend(np.diff(cs[ptr[left:right+1]-lo]))
        depths=np.array(depths);assert (depths>0).all()
        report['contexts'][c]={'control_cells':len(obs),'depth_min':int(depths.min()),'depth_max':int(depths.max())}
        print('CHECKED CONTEXT',c,report['contexts'][c],flush=True)
pairs=pd.read_csv(r/'data/official_pairs.csv')
assert pairs.target.isin(targets).all() and pairs.neighbor.isin(genes).all()
assert pairs.distance.between(0,5000).all()
report['promoter_pairs']=len(pairs)
report['targets_without_any_source']=coverage.loc[(coverage.iloc[:,1:]==0).all(axis=1),'target_gene'].tolist()
coverage.to_csv(r/'data/source_coverage.csv',index=False)
(r/'logs/input_validation.json').write_text(json.dumps(report,indent=2))
print('ALL INPUT CHECKS PASSED',json.dumps(report),flush=True)
