"""Finite, balanced reaction paths and mechanical screening for ideal rocks.

Amounts are mol, masses kg, temperature K and pressure Pa internally. This is
not a phase-equilibrium solver or a fracture-propagation simulation.
"""
from __future__ import annotations
import math

ATOMIC = {'H':1.00794e-3,'O':15.9994e-3,'Mg':24.305e-3,'Fe':55.845e-3,
          'Si':28.0855e-3,'Ca':40.078e-3,'Al':26.9815385e-3,'Na':22.989769e-3,'K':39.0983e-3}
# Nominal reference densities for volume screening, not high-P/T equations of state.
MINERALS = {
 'forsterite':({'Mg':2,'Si':1,'O':4},3270),
 'fayalite':({'Fe':2,'Si':1,'O':4},4390),
 'enstatite':({'Mg':1,'Si':1,'O':3},3200),
 'ferrosilite':({'Fe':1,'Si':1,'O':3},4000),
 'lizardite':({'Mg':3,'Si':2,'O':9,'H':4},2550),
 'brucite':({'Mg':1,'O':2,'H':2},2390),
 'magnetite':({'Fe':3,'O':4},5180),
 'hematite':({'Fe':2,'O':3},5260),
 'quartz':({'Si':1,'O':2},2650),
 'diopside':({'Ca':1,'Mg':1,'Si':2,'O':6},3280),
 'anorthite':({'Ca':1,'Al':2,'Si':2,'O':8},2760),
 'albite':({'Na':1,'Al':1,'Si':3,'O':8},2620),
 'microcline':({'K':1,'Al':1,'Si':3,'O':8},2560),
 'annite':({'K':1,'Fe':3,'Al':1,'Si':3,'O':12,'H':2},3300),
 'water':({'H':2,'O':1},None),'hydrogen':({'H':2},None)}
REACTIONS = {
 'forsterite':{'forsterite':-2,'water':-3,'lizardite':1,'brucite':1},
 'fayalite':{'fayalite':-3,'water':-2,'magnetite':2,'quartz':3,'hydrogen':2},
 'enstatite':{'enstatite':-3,'water':-2,'lizardite':1,'quartz':1},
 'ferrosilite':{'ferrosilite':-3,'water':-1,'magnetite':1,'quartz':3,'hydrogen':1},
 'magnetite':{'magnetite':-2,'water':-1,'hematite':3,'hydrogen':1}}
R = 8.314462618
YEAR = 365.25*86400


def molar_mass(name):
    return sum(ATOMIC[e]*n for e,n in MINERALS[name][0].items())


def totals(amounts):
    return {e:sum(n*MINERALS[s][0].get(e,0) for s,n in amounts.items()) for e in ATOMIC}


def solid_volume(amounts):
    return sum(n*molar_mass(s)/MINERALS[s][1] for s,n in amounts.items() if MINERALS[s][1])


def validate_reactions():
    for name,stoich in REACTIONS.items():
        if max(abs(x) for x in totals(stoich).values())>1e-12:
            raise ValueError('Unbalanced reaction: '+name)


def rate_per_year(c):
    t=c['temperature_c']+273.15;tr=c['reference_temperature_c']+273.15
    # Delta V activation may have either sign. Do not confuse it with reaction volume.
    dv=c['activation_volume_cm3_mol']*1e-6
    exponent=-c['activation_energy_kj_mol']*1000/R*(1/t-1/tr)-dv*(c['pore_pressure_mpa']-c['reference_pressure_mpa'])*1e6/(R*t)
    return math.log(2)/c['half_time_reference_years']*math.exp(exponent)


def mechanics(initial,final,c):
    v0=solid_volume(initial);v1=solid_volume(final);vb=v0/(1-c['porosity'])
    strain=(v1-v0)/vb
    # Only the declared pore fraction can accommodate precipitation without load.
    pore_allowance=c['porosity']*c['available_pore_fraction']
    constrained=max(strain-pore_allowance,0)
    trial=c['effective_bulk_modulus_gpa']*1000*constrained
    effective=c['minimum_total_stress_mpa']-c['biot_coefficient']*c['pore_pressure_mpa']
    tensile=max(trial-effective,0)
    ki=c['crack_geometry_factor']*tensile*math.sqrt(math.pi*c['flaw_half_length_m'])
    return {'initial_solid_volume_m3':v0,'final_solid_volume_m3':v1,
            'solid_volume_change_percent':100*(v1/v0-1),'bulk_reaction_strain':strain,
            'free_pore_allowance_strain':pore_allowance,'trial_expansion_stress_mpa':trial,
            'minimum_effective_stress_mpa':effective,'screening_tensile_stress_mpa':tensile,
            'stress_intensity_mpa_sqrt_m':ki,'fracture_index':ki/c['fracture_toughness_mpa_sqrt_m'],
            'cracking_screen_positive':ki>=c['fracture_toughness_mpa_sqrt_m'],
            'mechanical_scope':'Constrained-volume screening with fixed nominal mineral densities; no force equilibrium, crack aperture, propagation or permeability prediction.'}


def calculate(rock,c):
    validate_reactions()
    if set(rock['mass_fractions'])-set(MINERALS):raise ValueError('Unknown mineral')
    if abs(sum(rock['mass_fractions'].values())-1)>1e-10 or min(rock['mass_fractions'].values())<0:raise ValueError('Nonnegative rock mass fractions must sum to one')
    if any(not math.isfinite(float(v)) for v in c.values()):raise ValueError('Nonfinite input')
    if not 25<=c['temperature_c']<=300:raise ValueError('Scenario temperature must be 25–300 C; no stability domain is implied')
    if not 0<=c['accessible_fraction']<=1 or not 0<c['porosity']<1:raise ValueError('Invalid fractions')
    if not 0<=c['available_pore_fraction']<=1:raise ValueError('Invalid pore availability')
    for k in ['rock_mass_kg','half_time_reference_years','effective_bulk_modulus_gpa','flaw_half_length_m','fracture_toughness_mpa_sqrt_m']:
        if c[k]<=0:raise ValueError(k+' must be positive')
    if min(c['water_rock_mass_ratio'],c['duration_years'],c['pore_pressure_mpa'])<0:raise ValueError('Negative water, time or pressure')
    initial={s:c['rock_mass_kg']*w/molar_mass(s) for s,w in rock['mass_fractions'].items()}
    initial['water']=c['rock_mass_kg']*c['water_rock_mass_ratio']/molar_mass('water');initial['hydrogen']=0.
    k=rate_per_year(c);fraction=c['accessible_fraction']*(-math.expm1(-k*c['duration_years']))
    extents={s:initial.get(s,0)*fraction/(-REACTIONS[s][s]) for s in rock['reaction_paths']}
    demands={}
    for s,xi in extents.items():
        for species,nu in REACTIONS[s].items():
            if nu<0:demands[species]=demands.get(species,0)-nu*xi
    scale=min([1.]+[initial.get(s,0)/n for s,n in demands.items() if n>0])
    final=dict(initial)
    for s,xi in extents.items():
        extents[s]=xi*scale
        for species,nu in REACTIONS[s].items():final[species]=final.get(species,0)+nu*xi*scale
    if min(final.values()) < -1e-9:raise AssertionError('Negative phase amount')
    before,after=totals(initial),totals(final)
    error=max(abs(after[e]-before[e])/max(abs(before[e]),1) for e in before)
    if error>1e-10:raise AssertionError('Element conservation failed')
    phases=[{'mineral':s,'initial_mol':initial.get(s,0),'final_mol':n,'final_mass_kg':n*molar_mass(s),'final_solid_volume_m3':n*molar_mass(s)/MINERALS[s][1]} for s,n in sorted(final.items()) if MINERALS[s][1]]
    return {'rock':rock['name'],'scope':'Specified reaction-path assemblage; not a predicted equilibrium assemblage or calibrated geological rate.',
            'conditions':dict(c),'rate_per_year':k,'accessible_conversion_fraction':-math.expm1(-k*c['duration_years']),
            'whole_initial_reactive_mineral_fraction':fraction*scale,'water_limitation_multiplier':scale,
            'hydrogen_mol':final['hydrogen'],'hydrogen_kg':final['hydrogen']*molar_mass('hydrogen'),
            'water_consumed_kg':(initial['water']-final['water'])*molar_mass('water'),
            'remaining_water_kg':final['water']*molar_mass('water'),'phases':phases,'reaction_extents_mol':extents,
            'max_relative_element_balance_error':error, 'mechanics':mechanics(initial,final,c)}


def radiolysis(power_w_into_water,years,yield_molecules_per_100ev=.45):
    """Gross yield for supplied deposited power; G=.45 is a low-LET scenario.

    This does not infer deposited power from bulk U/Th/K, nor account for
    recombination, transport, oxidants or temperature-dependent track chemistry.
    """
    if min(power_w_into_water,years,yield_molecules_per_100ev)<0:raise ValueError('Negative radiolysis input')
    return power_w_into_water*years*YEAR*yield_molecules_per_100ev/(100*1.602176634e-19*6.02214076e23)


def reaction_gibbs(delta_g_reference_j_mol,delta_h_reference_j_mol,delta_v_m3_mol,t_k,p_pa,q=1.,t_ref=298.15,p_ref=1e5):
    """Constant delta H/S/V approximation; caller supplies matched reaction data."""
    if min(t_k,q,t_ref,p_pa)<=0:raise ValueError('Positive T, P and dimensionless reaction quotient required')
    ds=(delta_h_reference_j_mol-delta_g_reference_j_mol)/t_ref
    return delta_h_reference_j_mol-t_k*ds+delta_v_m3_mol*(p_pa-p_ref)+R*t_k*math.log(q)
