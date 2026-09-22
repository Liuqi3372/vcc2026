from pathlib import Path
import sys,hashlib
root=Path(__file__).resolve().parents[1];sys.path.insert(0,str(root/'src'))
import prepare
def local_only(name,directory):
 s=prepare.ASSETS[name];p=directory/s['filename']
 if not p.is_file():raise FileNotFoundError('Offline preparation missing '+str(p))
 if p.stat().st_size!=s['bytes']:raise ValueError('Source size mismatch '+str(p))
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8*1024**2),b''):h.update(b)
 if h.hexdigest()!=s['sha256']:raise ValueError('Source hash mismatch '+str(p))
 print('Verified',p,flush=True);return p
prepare.download=local_only
for name,func,out in [('cd4',prepare.prepare_cd4,'CD4_DE_statistics.npz'),('promoters',prepare.prepare_promoters,'official_pairs.csv')]:
 if not (root/'data'/out).exists():
  print('Preparing',name,flush=True);func(root/'data',root/'raw');print('Saved',out,flush=True)
