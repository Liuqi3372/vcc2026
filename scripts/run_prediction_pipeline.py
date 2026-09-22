from pathlib import Path
import hashlib
import json
import os
import subprocess
import sys
import time
import shutil
import tempfile

r=Path(__file__).resolve().parents[1]
status=r/'logs/prediction_pipeline_status.json'
start=time.time()
def state(stage,**kwargs):
    value=dict(stage=stage,pid=os.getpid(),started_at=start,updated_at=time.time(),**kwargs)
    tmp=status.with_suffix('.tmp');tmp.write_text(json.dumps(value,indent=2));tmp.replace(status)
    print(json.dumps(value),flush=True)
def run(stage,args,env=None):
    state(stage)
    with (r/'logs'/(stage+'.log')).open('a') as log:
        subprocess.run(args,cwd=r,stdout=log,stderr=subprocess.STDOUT,check=True,env=env)
try:
    state('waiting_for_source_statistics')
    while True:
        s=json.loads((r/'logs/local_sources_status.json').read_text())
        if s['state']=='failed':raise RuntimeError(str(s))
        if s['state']=='complete':break
        time.sleep(20)
    prior_audit=r/'logs/input_validation.json'
    if not prior_audit.exists() or any(p.stat().st_mtime > prior_audit.stat().st_mtime for p in (r/'data').iterdir() if p.is_file()):
        run('audit_inputs',[sys.executable,'-u',str(r/'scripts/audit_inputs.py')])
    out=r/'output';out.mkdir(exist_ok=True)
    pred=out/'prediction_ABC_local.h5ad'


    compact=pred
    packed=out/'prediction_ABC_local.vcc'
    if not pred.exists():
        run('predict',[sys.executable,'-u',str(r/'src/predict.py'),'--data-dir',str(r/'data'),'--output',str(pred)])
    if not compact.exists():
        if shutil.disk_usage('/dev/shm').free < 60*1024**3:
            raise RuntimeError('Need 60 GiB temporary space for compaction')
        with tempfile.TemporaryDirectory(prefix='xiexiu_compact_',dir='/dev/shm') as tmp:
            fast=Path(tmp)/compact.name
            run('compact',[sys.executable,'-u',str(r/'src/compact.py'),str(pred),str(fast),'--threads','4'])
            state('saving_compact')
            pending=compact.with_suffix('.copying.h5ad')
            if pending.exists():raise FileExistsError(pending)
            shutil.copyfile(fast,pending)
            def digest(path):
                h=hashlib.sha256()
                with path.open('rb') as f:
                    for b in iter(lambda:f.read(8*1024**2),b''):h.update(b)
                return h.hexdigest()
            assert digest(fast)==digest(pending), 'Compact copy checksum mismatch'
            pending.rename(compact)
    if not packed.exists():
        with tempfile.TemporaryDirectory(prefix='xiexiu_pack_',dir='/dev/shm') as tmp:
            run('pack',[str(r/'vcc-env/bin/python'),'-u',str(r/'src/pack.py'),str(compact),'--data-dir',str(r/'data'),'--output',str(packed),'--scratch-dir',tmp],env={**os.environ,'TMPDIR':tmp})
    state('verifying_final_artifacts')
    import h5py
    from anndata._io.specs import read_elem
    import numpy as np
    import pandas as pd
    with h5py.File(compact) as f:
        obs,var=read_elem(f['obs']),read_elem(f['var'])
        assert (len(obs),len(var))==(360000,18533)
        assert obs.index.is_unique
        counts=obs.groupby(['context','target_gene'],observed=True).size()
        assert len(counts)==900 and (counts==400).all()
        assert set(obs.context.astype(str))==set('ABC')
        np.testing.assert_array_equal(var.index.astype(str),pd.read_csv(r/'data/gene_names.csv').gene_name.astype(str))
        assert int(f['X/indptr'][-1])==len(f['X/data'])
    files=[]
    for path in dict.fromkeys([pred,compact,packed]):
        h=hashlib.sha256()
        with path.open('rb') as f:
            for b in iter(lambda:f.read(8*1024**2),b''):h.update(b)
        files.append({'path':str(path),'bytes':path.stat().st_size,'sha256':h.hexdigest()})
    report={'state':'complete','contexts':['A','B','C'],'cells':360000,'genes':18533,'targets_per_context':300,'cells_per_target':400,'files':files,'elapsed_seconds':time.time()-start,'variant':'Local QC K562/HCT116/HEK293T; original H1/CD4/promoter preparation','official_vcc_prep_passed':True,'optional_compaction':'skipped; official package built directly from original CSR prediction'}
    (r/'logs/final_validation.json').write_text(json.dumps(report,indent=2))
    state('complete',report=str(r/'logs/final_validation.json'))
except Exception as e:
    state('failed',error=str(e));raise
