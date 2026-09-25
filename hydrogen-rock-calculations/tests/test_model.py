import json,math,unittest
from pathlib import Path
from rock_model.model import calculate,validate_reactions,molar_mass,rate_per_year,radiolysis,reaction_gibbs
SPEC=json.loads((Path(__file__).resolve().parents[1]/'inputs/scenarios.json').read_text())
class Tests(unittest.TestCase):
 def test_stoichiometry(self):validate_reactions()
 def test_all_rock_element_balances(self):
  for rock in SPEC['rocks']:
   r=calculate(rock,SPEC['conditions']);self.assertLess(r['max_relative_element_balance_error'],1e-12)
 def test_fayalite_analytical_yield(self):
  c=SPEC['conditions']|{'accessible_fraction':1.,'duration_years':1e8}
  r=calculate({'name':'pure','mass_fractions':{'fayalite':1.},'reaction_paths':['fayalite']},c)
  self.assertAlmostEqual(r['hydrogen_mol'],(1/molar_mass('fayalite'))*2/3)
 def test_magnesium_endmember_makes_no_hydrogen(self):
  r=calculate({'name':'pure','mass_fractions':{'forsterite':1.},'reaction_paths':['forsterite']},SPEC['conditions'])
  self.assertEqual(r['hydrogen_mol'],0);self.assertGreater(r['mechanics']['solid_volume_change_percent'],0)
 def test_water_limit(self):
  r=calculate(SPEC['rocks'][0],SPEC['conditions']|{'water_rock_mass_ratio':1e-5})
  self.assertAlmostEqual(r['remaining_water_kg'],0,places=12);self.assertLess(r['water_limitation_multiplier'],1)
 def test_pressure_rate_sign(self):
  c=SPEC['conditions'];a=rate_per_year(c)
  self.assertLess(rate_per_year(c|{'pore_pressure_mpa':30.,'activation_volume_cm3_mol':10.}),a)
  self.assertGreater(rate_per_year(c|{'pore_pressure_mpa':30.,'activation_volume_cm3_mol':-10.}),a)
 def test_pressure_effective_stress(self):
  a=calculate(SPEC['rocks'][0],SPEC['conditions']);b=calculate(SPEC['rocks'][0],SPEC['conditions']|{'pore_pressure_mpa':30.})
  self.assertEqual(a['hydrogen_mol'],b['hydrogen_mol']);self.assertAlmostEqual(a['mechanics']['minimum_effective_stress_mpa']-b['mechanics']['minimum_effective_stress_mpa'],10.)
 def test_zero_access_and_zero_time(self):
  for override in [{'accessible_fraction':0.},{'duration_years':0.}]:
   r=calculate(SPEC['rocks'][0],SPEC['conditions']|override);self.assertEqual(r['hydrogen_mol'],0)
 def test_gibbs_pressure_units(self):
  a=reaction_gibbs(-10000,-20000,1e-5,298.15,1e5)
  b=reaction_gibbs(-10000,-20000,1e-5,298.15,1e5+1e6)
  self.assertAlmostEqual(b-a,10.)
 def test_radiolysis_energy_units(self):
  self.assertEqual(radiolysis(0,1),0)
  self.assertAlmostEqual(radiolysis(1,2),2*radiolysis(1,1))
if __name__=='__main__':unittest.main()
