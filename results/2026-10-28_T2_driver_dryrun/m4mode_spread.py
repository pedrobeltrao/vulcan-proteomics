# M4-mode (non-deterministic CUDA index_add_) run-to-run spread of our localizer's E8 fields on dev D2, using M4's PSMs/stores.
import sys, os, json
sys.dont_write_bytecode=True
ROOT='$VULCAN_ROOT'
sys.path.insert(0, ROOT); sys.path.insert(0, ROOT+'/bench')
import yaml
from t2 import datasets as DS, analysis as AN
import e8 as M4E8
from vulcan.ptm import localize as LZ
from vulcan.index.build import L1Index
cfg=yaml.safe_load(open(ROOT+'/config/t2.yaml'))
lp=LZ.LocParams(**cfg['localization'])
ds=DS.dataset('D2'); idx='/tmp/vulcan_m4/index/d2_phospho_prod'; l1=L1Index(idx); truths=DS.truth_for(ds)
KEYS=[('all_thr_0.01','n_sty_eligible'),('all_thr_0.01','true_flr'),('all_thr_0.05','n_sty_eligible'),('all_thr_0.05','true_flr'),
      ('recall_at_true_flr_0.01',None),('recall_at_true_flr_0.05',None),('n_sty_at_flrhat_0.01',None)]
ref=json.load(open(ROOT+'/results/2026-10-27_M4_ptm/e8.json'))['datasets']
out={}
for acc,d in ds['accessions'].items():
    runs=[DS.run_name(r) for r in d['raws']]
    psm,_=AN.psm_list(runs,'/tmp/vulcan_m4/psm_d2_phospho',idx,l1)
    T=truths[acc]['truth']; tof=[T.get((rn,s)) if len(ph)==1 else None for rn,s,ph in zip(psm['run'],psm['seq'],psm['phos'])]
    sp,_=AN.load_spectra(['/tmp/vulcan_m4/store/%s.vspec.arrow'%r for r in runs], list(zip(psm['run'],psm['spec_idx'])))
    for dec in 'AL':
        vals=[]
        for rep in range(5):
            s,_,_=M4E8.run_ours(psm,sp,'STY'+dec,lp)
            e=M4E8.evaluate(s,psm,tof,dec)['mono']
            vals.append([e[a][b] if b else e[a] for a,b in KEYS])
        sd,_,_=AN.run_ours_deterministic(psm,sp,'STY'+dec,lp); e=M4E8.evaluate(sd,psm,tof,dec)['mono']
        det=[e[a][b] if b else e[a] for a,b in KEYS]
        r=ref[acc]['ours_'+dec]['eval']['mono']; m4=[r[a][b] if b else r[a] for a,b in KEYS]
        out[f'{acc}.{dec}']={'fields':[f'{a}.{b}' if b else a for a,b in KEYS],
            'm4_mode_min':[min(v[i] for v in vals) for i in range(len(KEYS))],'m4_mode_max':[max(v[i] for v in vals) for i in range(len(KEYS))],
            'm4_committed':m4,'deterministic':det}
        print(acc,dec,json.dumps(out[f'{acc}.{dec}']),flush=True)
json.dump(out,open(sys.argv[1],'w'),indent=1)
