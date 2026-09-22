from pathlib import Path
import subprocess,sys,os,json,time
r=Path(__file__).resolve().parents[1];state=r/'logs/local_sources_status.json'
try:
 for source in ['K562_GWPS_CPM','HCT116','HEK293T']:
  state.write_text(json.dumps({'state':'preparing','source':source,'updated_at':time.time()}))
  with (r/'logs'/('prepare_'+source+'.log')).open('a') as log:
   subprocess.run([sys.executable,'-u',str(r/'scripts/prepare_local_atlas.py'),'--source',source],stdout=log,stderr=subprocess.STDOUT,check=True)
 state.write_text(json.dumps({'state':'complete','updated_at':time.time()}))
except Exception as e:
 state.write_text(json.dumps({'state':'failed','error':str(e),'updated_at':time.time()}));raise
