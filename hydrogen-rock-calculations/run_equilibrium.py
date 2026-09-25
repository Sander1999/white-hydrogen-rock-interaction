#!/usr/bin/env python3
"""Restricted low-pressure PHREEQC pure-endmember equilibrium demonstration.

Supply a vetted external SUPCRTBL-derived database; no database is downloaded
or substituted automatically. Not a high-pressure or brine assemblage model.
"""
import argparse,csv,hashlib,json
from pathlib import Path
from rock_model.model import MINERALS,ATOMIC,molar_mass,totals
ROOT=Path(__file__).resolve().parent
PHASES=['FORSTERITE','FAYALITE','ENSTATITE','FERROSILITE','LIZARDITE','BRUCITE',
        'MAGNETITE','HEMATITE','QUARTZ','TALC','GREENALITE','MINNESOTAITE',
        'DIOPSIDE','ANORTHITE','ALBITE','MICROCLINE','ANNITE','CLINOCHLORE(ORDERED)']
FORMULAS={s.upper():v[0] for s,v in MINERALS.items() if v[1]}
FORMULAS.update(TALC={'Mg':3,'Si':4,'O':12,'H':2},GREENALITE={'Fe':3,'Si':2,'O':9,'H':4},
                MINNESOTAITE={'Fe':3,'Si':4,'O':12,'H':2})
FORMULAS['CLINOCHLORE(ORDERED)']={'Mg':5,'Al':2,'Si':3,'O':18,'H':8}
FORMULAS['ENSTATITE']={'Mg':2,'Si':2,'O':6}
FORMULAS['FERROSILITE']={'Fe':2,'Si':2,'O':6}
ELEMENTS=['Mg','Fe','Si','Ca','Al','Na','K']


def solve(rock,database,temperature=80.,mass=.1,library=None):
    from phreeqpy.iphreeqc.phreeqc_dll import IPhreeqc
    if not 25<=temperature<=95:raise ValueError('This demonstration is restricted to 25–95 C and 1 atm; high-pressure databases require separate validation')
    if mass<=0:raise ValueError('Positive rock mass required')
    initial={s:mass*w/sum(ATOMIC[e]*v for e,v in FORMULAS[s.upper()].items()) for s,w in rock['mass_fractions'].items()}
    phase_lines='\n'.join(f' {s} 0 {initial.get(s.lower(),0):.17g}' for s in PHASES)
    text=f'''TITLE Ideal {rock['name']} pure-endmember equilibrium at 1 atm
SOLUTION 1
 temp {temperature}
 pressure 1
 pH 7
 pe 4
 units mol/kgw
 Na 0.001
 Cl 0.001
 -water 1
EQUILIBRIUM_PHASES 1
{phase_lines}
 H2(g) 0 0
SELECTED_OUTPUT
 -reset false
 -high_precision true
 -pH true
 -pe true
 -water true
 -totals {' '.join(ELEMENTS)}
 -equilibrium_phases {' '.join(PHASES)} H2(g)
 -molalities H2
END
'''
    expected=json.loads((ROOT/'inputs/equilibrium_database.json').read_text())['normalized_text_sha256']
    if hashlib.sha256(database.read_text().encode()).hexdigest()!=expected:
        raise ValueError('Database differs from the tested phase stoichiometry and constants; validate it before updating inputs/equilibrium_database.json')
    q=IPhreeqc(library);q.load_database(str(database));q.run_string(text)
    array=q.get_selected_output_array();last=dict(zip(array[0],array[-1]));first=dict(zip(array[0],array[1]))
    water=last['mass_H2O'];balance={}
    for element in ELEMENTS:
        target=sum(n*FORMULAS[s.upper()].get(element,0) for s,n in initial.items())+first[f'{element}(mol/kgw)']*first['mass_H2O']
        actual=last[f'{element}(mol/kgw)']*water+sum(last[s]*FORMULAS[s].get(element,0) for s in PHASES)
        balance[element]=abs(actual-target)/max(target,1e-9)
    if max(balance.values())>1e-6:raise AssertionError('Mineral/aqueous cation balance failed: '+str(balance))
    if min(last[s] for s in PHASES)<-1e-9:raise AssertionError('Negative equilibrium phase amount')
    r={'rock':rock['name'],'temperature_c':temperature,'pressure_atm':1.,'rock_mass_kg':mass,'initial_water_kg':1.,
       'initial_na_cl_molal':.001,'pH':last['pH'],'pe':last['pe'],
       'solid_moles':{s:last[s] for s in PHASES},'aqueous_water_mass_kg':water,
       'H2_gas_mol':last['H2(g)'],'H2_aqueous_mol':last['m_H2(mol/kgw)']*water,
       'initial_H2_aqueous_mol':first['m_H2(mol/kgw)']*first['mass_H2O'],
       'cation_relative_balance_errors':balance,'database_sha256':hashlib.sha256(database.read_bytes()).hexdigest(),
       'scope':'Equilibrium within the listed pure phases and external database at 1 atm. No kinetic duration, Fe-Mg solid solutions or high-pressure mineral corrections. H2(g) can precipitate at unit fugacity from zero initial gas; no imposed H2 supply.'}
    return r,text


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--database',type=Path,required=True)
    p.add_argument('--library',help='Optional absolute IPhreeqc shared-library path')
    p.add_argument('--input',type=Path,default=ROOT/'inputs/scenarios.json')
    p.add_argument('--temperature',type=float,default=80);p.add_argument('--rock-mass',type=float,default=.1)
    p.add_argument('--output',type=Path,default=ROOT/'results/equilibrium');a=p.parse_args()
    spec=json.loads(a.input.read_text());a.output.mkdir(parents=True,exist_ok=True);results=[]
    for i,rock in enumerate(spec['rocks']):
        r,text=solve(rock,a.database,a.temperature,a.rock_mass,a.library);results.append(r)
        (a.output/f'case_{i+1}.pqi').write_text(text)
    (a.output/'assemblages.json').write_text(json.dumps(results,indent=2)+'\n')
    rows=[{'rock':r['rock'],'mineral':s,'equilibrium_mol':n} for r in results for s,n in r['solid_moles'].items()]
    with (a.output/'assemblages.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    print(json.dumps([{'rock':r['rock'],'pH':r['pH'],'H2_mol':r['H2_gas_mol']+r['H2_aqueous_mol'],'max_cation_error':max(r['cation_relative_balance_errors'].values())} for r in results],indent=2))


if __name__=='__main__':main()
