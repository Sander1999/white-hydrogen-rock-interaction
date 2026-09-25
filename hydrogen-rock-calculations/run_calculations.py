#!/usr/bin/env python3
"""Calculate conditional reaction paths, assemblages, H2 and cracking indices."""
import argparse,csv,json
from pathlib import Path
from rock_model.model import calculate, radiolysis
ROOT=Path(__file__).resolve().parent


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input',type=Path,default=ROOT/'inputs/scenarios.json')
    p.add_argument('--output',type=Path,default=ROOT/'results')
    p.add_argument('--plot',action='store_true')
    a=p.parse_args();spec=json.loads(a.input.read_text());a.output.mkdir(parents=True,exist_ok=True)
    results=[calculate(rock,spec['conditions']) for rock in spec['rocks']]
    (a.output/'assemblages.json').write_text(json.dumps(results,indent=2)+'\n')
    rows=[]
    for r in results:
        for phase in r['phases']:rows.append({'rock':r['rock'],**phase})
    def write(name,rows):
        with (a.output/name).open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    write('assemblages.csv',rows)
    summary=[{'rock':r['rock'],'hydrogen_g_per_kg_rock':1000*r['hydrogen_kg']/r['conditions']['rock_mass_kg'],
              'water_consumed_kg_per_kg_rock':r['water_consumed_kg']/r['conditions']['rock_mass_kg'],
              'solid_volume_change_percent':r['mechanics']['solid_volume_change_percent'],
              'fracture_index':r['mechanics']['fracture_index'],
              'element_balance_error':r['max_relative_element_balance_error']} for r in results]
    write('reference_summary.csv',summary)
    sweep=[]
    for rock in spec['rocks']:
        for t in [80.,150.,250.]:
            for pressure in [10.,20.,30.]:
                for dv in [-10.,0.,10.]:
                    c=spec['conditions']|{'temperature_c':t,'pore_pressure_mpa':pressure,'activation_volume_cm3_mol':dv}
                    r=calculate(rock,c)
                    sweep.append({'rock':rock['name'],'temperature_c':t,'pore_pressure_mpa':pressure,'activation_volume_cm3_mol':dv,
                                  'hydrogen_g_per_kg_rock':1000*r['hydrogen_kg']/c['rock_mass_kg'],
                                  'fracture_index':r['mechanics']['fracture_index'],'rate_per_year':r['rate_per_year'],
                                  'scope':'Conditional kinetics only; mineral stability is not selected by this sweep'})
    write('temperature_pressure_scenarios.csv',sweep)
    radio=[{'power_absorbed_by_water_w_per_kg_rock':power,'time_years':spec['conditions']['duration_years'],
            'gross_hydrogen_mol_per_kg_rock':radiolysis(power,spec['conditions']['duration_years']),
            'G_molecules_per_100eV':.45,'scope':'Illustrative deposited power and low-LET yield; not a granite production forecast'} for power in [1e-13,1e-12,1e-11]]
    write('radiolysis_scenarios.csv',radio)
    if a.plot:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig,axes=plt.subplots(1,2,figsize=(11,4.2),layout='constrained')
        for rock in spec['rocks'][:-1]:
            rr=[r for r in sweep if r['rock']==rock['name'] and r['pore_pressure_mpa']==20 and r['activation_volume_cm3_mol']==0]
            axes[0].plot([r['temperature_c'] for r in rr],[r['hydrogen_g_per_kg_rock'] for r in rr],marker='o',label=rock['name'].replace(' analogue',''))
        axes[0].set(xlabel='Temperature (°C)',ylabel='H₂ after 1000 years (g/kg initial rock)')
        axes[0].legend(fontsize=8)
        for dv in [-10,0,10]:
            rr=[r for r in sweep if r['rock']==spec['rocks'][0]['name'] and r['temperature_c']==80 and r['activation_volume_cm3_mol']==dv]
            axes[1].plot([r['pore_pressure_mpa'] for r in rr],[r['hydrogen_g_per_kg_rock'] for r in rr],marker='o',label=f'Activation volume {dv:+g} cm³/mol')
        axes[1].set(xlabel='Pore pressure (MPa)',ylabel='Dunite H₂ (g/kg initial rock)');axes[1].legend(fontsize=8)
        for ax in axes:ax.grid(alpha=.2)
        fig.suptitle('Prescribed reaction paths and illustrative rates\nThese curves do not predict mineral stability or an optimum temperature',fontsize=11)
        fig.savefig(a.output/'temperature_pressure.png',dpi=180);plt.close(fig)
    print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
