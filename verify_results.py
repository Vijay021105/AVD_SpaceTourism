"""Check generated mission output consistency. Does not turn open design gates into passes."""
from pathlib import Path
import argparse,json,csv,hashlib
import numpy as np


def verify(root):
    data=root/'data';s=json.loads((root/'mission_summary.json').read_text());cfg=json.loads((root/'effective_configuration.json').read_text())
    def read(n):
        with (data/n).open() as f:return list(csv.DictReader(f))
    flight=read('ascent.csv');ev=read('events.csv');inj=read('injection_state.csv')[0]
    assert all(float(r['mass_kg'])>0 for r in flight)
    assert all(float(r['q_Pa'])>=0 for r in flight)
    assert all(float(r['h_m'])>=0 for r in flight),'Ascent crosses surface'
    assert all(float(flight[i+1]['t_s'])>=float(flight[i]['t_s']) for i in range(len(flight)-1))
    assert abs(float(flight[-1]['rx_m'])-float(inj['rx_m']))<1e-8,'Rendezvous was reset to another injection position'
    for e in ev:
        if 'separation' in e['event']:
            j=int(e['event'][1])-1;s0=s['stages'][j]
            assert abs(float(e['mass_before_kg'])-float(e['mass_after_kg'])-s0['dry_kg']-s0['reserve_kg'])<1e-4
    inputs=root.parent/'inputs'
    if not inputs.is_dir():inputs=Path(cfg['input_dir'])
    for n,h in s['input_hashes'].items():assert hashlib.sha256((inputs/n).read_bytes()).hexdigest()==h
    for name in s['figures']:assert (root/'figures'/name).stat().st_size>1000
    assert s['rendezvous']['two_period_wait_s']==2*s['rendezvous']['nominal_period_s']
    assert all(np.isfinite(float(r['relative_distance_m'])) for r in read('rendezvous.csv'))
    for p in (data/'cfd').glob('*_residual.csv'):
        with p.open() as f:rr=list(csv.DictReader(f))
        assert all(float(r['residual'])>=0 for r in rr)
    print(f'Output integrity passed: {len(flight)} ascent samples, {len(s["figures"])} figures; {sum(c["status"]=="OPEN" for c in s["checks"])} design gates remain open.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('output',type=Path,nargs='?',default=Path(__file__).parent/'results');verify(p.parse_args().output)
