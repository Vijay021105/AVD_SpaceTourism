"""Constrained three-stage mass allocation; structural-factor sensitivity."""
import numpy as np
from scipy.optimize import minimize
from .io import G0, write_csv


def ledger(c, prop, eps=None, floors=None):
    eps=np.array([s['epsilon'] for s in c['stages']]) if eps is None else np.asarray(eps)
    fixed=np.array([s['fixed_dry_allowance_kg'] for s in c['stages']])
    dry=np.asarray(prop)*eps/(1-eps)+fixed
    if floors is not None: dry=np.maximum(dry, floors)
    upper=c['mission']['payload_kg']; rows=[]
    for j in (2,1,0):
        s=c['stages'][j];reserve=prop[j]*s['reserve_fraction']
        mi=upper+dry[j]+prop[j];mf=upper+dry[j]+reserve
        thrust=s['engine_count']*s['engine']['thrust_vac_N']
        if j==0: thrust=s['engine_count']*s['engine']['thrust_sl_N']
        rows.append(dict(stage=j+1,engine_id=s['engine_id'],engine_count=int(s['engine_count']),
            propellant_kg=prop[j],usable_propellant_kg=prop[j]-reserve,reserve_kg=reserve,dry_kg=dry[j],
            upper_stack_kg=upper,ignition_mass_kg=mi,cutoff_mass_kg=mf,epsilon_budget=eps[j],
            epsilon_effective=dry[j]/(prop[j]+dry[j]),dv_ideal_mps=G0*s['isp_sizing_s']*np.log(mi/mf),
            ignition_tw=thrust/(mi*G0),full_flow_kgps=s['mdot_full_kgps'],
            minimum_burn_s=(prop[j]-reserve)/s['mdot_full_kgps']))
        upper=mi
    return list(reversed(rows))


def optimize(c, eps=None, floors=None):
    target=c['mission']['total_dv_budget_mps']-c['mission']['spacecraft_dv_allocation_mps']
    limits=[c['mission'][k] for k in ['liftoff_tw_min','stage2_tw_min','stage3_tw_min']]
    history=[]
    def rows(x):return ledger(c,np.exp(x)*1e5,eps,floors)
    def constraint(x):
        r=rows(x)
        return np.r_[(sum(s['dv_ideal_mps'] for s in r)-target)/1000,
                     [r[j]['ignition_tw']-limits[j] for j in range(3)]]
    def callback(x):
        r=rows(x);history.append(dict(iteration=len(history),liftoff_mass_kg=r[0]['ignition_mass_kg'],
            dv_ideal_mps=sum(s['dv_ideal_mps'] for s in r),min_scaled_constraint=min(constraint(x))))
    x0=np.log(np.array([s['initial_prop_kg'] for s in c['stages']])/1e5)
    sol=minimize(lambda x:rows(x)[0]['ignition_mass_kg']/1e5,x0,method='SLSQP',
        bounds=[(np.log(.4),np.log(14)),(np.log(.2),np.log(5)),(np.log(.12),np.log(1.5))],
        constraints=[dict(type='ineq',fun=constraint)],callback=callback,
        options=dict(maxiter=220,ftol=1e-10))
    r=rows(sol.x)
    return r,dict(success=bool(sol.success and min(constraint(sol.x))>-1e-5),message=str(sol.message),
        target_launcher_dv_mps=target,achieved_dv_mps=sum(s['dv_ideal_mps'] for s in r),
        minimum_scaled_constraint=float(min(constraint(sol.x))),objective='minimum liftoff mass at selected structural factors',
        conditional_on='assumed structural factors; engine performance; dry floors; no global optimum proof'),history


def run(c,out):
    rows,status,hist=optimize(c)
    write_csv(out/'data/sizing_iterations.csv',hist)
    sensitivity=[]
    for case,key in [('optimistic','epsilon_low'),('selected','epsilon'),('conservative','epsilon_high')]:
        r,st,_=optimize(c,eps=[s[key] for s in c['stages']])
        sensitivity.append(dict(case=case,feasible=st['success'],mass_kg=r[0]['ignition_mass_kg'],
            dv_mps=st['achieved_dv_mps'],tw1=r[0]['ignition_tw'],tw2=r[1]['ignition_tw'],tw3=r[2]['ignition_tw']))
    write_csv(out/'data/structural_fraction_trade.csv',sensitivity)
    write_csv(out/'data/stages.csv',rows)
    for s,r in zip(c['stages'],rows):s.update(r)
    return rows,status,hist,sensitivity
