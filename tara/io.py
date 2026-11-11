"""Strict CSV input contracts, SI constants and reproducible output helpers."""
from pathlib import Path
import csv
import hashlib
import json
import numpy as np

MU = 3.986004418e14
RE = 6378137.0
J2 = 1.08262668e-3
OMEGA = 7.292115e-5
G0 = 9.80665


def table(path):
    with Path(path).open(encoding='utf-8-sig', newline='') as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise ValueError(f'Empty CSV: {path}')
    for row in rows:
        for key, value in row.items():
            if key is None or value is None:
                raise ValueError(f'Malformed row: {path}')
            try:
                row[key] = float(value)
            except ValueError:
                row[key] = value
    return rows


def write_csv(path, rows):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    rows = list(rows)
    if not rows:
        raise ValueError(f'Refusing headerless empty CSV: {path}')
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with path.open('w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields); w.writeheader()
        for row in rows:
            w.writerow({k: (v.item() if isinstance(v, np.generic) else v) for k, v in row.items()})


def json_write(path, data):
    def default(x):
        if isinstance(x, np.ndarray): return x.tolist()
        if isinstance(x, np.generic): return x.item()
        if isinstance(x, Path): return str(x)
        raise TypeError(type(x).__name__)
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(data, indent=2, default=default, allow_nan=False)+'\n')


SCALARS = ('mission', 'geometry', 'aerodynamics', 'structures', 'thermal', 'dynamics', 'guidance')
TABLES = ('stages', 'engines', 'fluids', 'materials', 'wind', 'human_rating', 'payload')


def load_inputs(folder):
    folder = Path(folder)
    expected = set(SCALARS + TABLES)
    actual = {p.stem for p in folder.glob('*.csv')}
    if actual != expected:
        raise ValueError(f'Input CSV set differs: missing={expected-actual}, unexpected={actual-expected}')
    c = {'input_dir': str(folder.resolve()), 'input_hashes': {}}
    for name in sorted(expected):
        p = folder/f'{name}.csv'; rows = table(p)
        c['input_hashes'][p.name] = hashlib.sha256(p.read_bytes()).hexdigest()
        if name in SCALARS:
            keys = [r['parameter'] for r in rows]
            if len(set(keys)) != len(keys): raise ValueError(f'Duplicate parameter in {p}')
            c[name] = {r['parameter']:r['value'] for r in rows}
        else:
            c[name] = rows
    c['stages'].sort(key=lambda r:r['stage'])
    if [int(s['stage']) for s in c['stages']] != [1,2,3]: raise ValueError('Exactly stages 1,2,3 required')
    if [int(s['engine_count']) for s in c['stages']] != [10,2,1]: raise ValueError('Baseline is fixed 10/2/1')
    engines = {r['engine_id']:r for r in c['engines']}
    if c['stages'][1]['engine_id'] != 'M1D_standard': raise ValueError('S2 must use standard Merlin 1D')
    if not 0 < c['mission']['diameter_m'] <= 4: raise ValueError('External diameter must be <=4 m')
    for s in c['stages']:
        if not 0 < s['epsilon'] < .3 or not 0 <= s['reserve_fraction'] < .2: raise ValueError('Invalid stage fractions')
        if s['engine_id'] not in engines: raise ValueError('Unknown engine')
        s['engine'] = engines[s['engine_id']]
        for k in ('isp_vac_s','thrust_vac_N','engine_mass_kg'):
            if s['engine'][k] <= 0: raise ValueError(f'Engine {k} must be positive')
        e=s['engine']
        if e['thrust_sl_N'] > e['thrust_vac_N']: raise ValueError('Thrust must increase toward vacuum')
        s['mdot_full_kgps'] = s['engine_count']*e['thrust_vac_N']/(G0*e['isp_vac_s'])
    if abs(sum(p['mass_kg'] for p in c['payload'])-c['mission']['payload_kg'])>1e-6:
        raise ValueError('Payload component ledger must sum to mission payload')
    wind = np.array([r['altitude_m'] for r in c['wind']])
    if np.any(np.diff(wind)<=0): raise ValueError('Wind altitude must increase strictly')
    for group in SCALARS:
        for key,v in c[group].items():
            if isinstance(v,float) and not np.isfinite(v): raise ValueError(f'Nonfinite {group}.{key}')
    return c
