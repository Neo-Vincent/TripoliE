# TripoliE
This is a project for post-processing of Tripoli-4 result file.
It containts 5 modules:
* BaseTripoli
* atom_mass
* postProcess
* TripoliE
* TripoliHandle

## BaseTripoli
This is a base class for other tripoli object.
It is used to parse the tripoli result file, and get the input information.
Now, it can almost get all infor from Tripoli Result file.

### Example
```python
>> a=BaseTripoli("MOX_L_T4.tri.txt")
>> a.maketable()
>> a.parsing_input()
>> a.water
[{'AL-NAT': 1.10277e-06,
  'B-NAT': 2.35598e-05,
  'CR-NAT': 2.24991e-05,
  'FE-NAT': 2.09013e-05,
  'H1': 0.0471346,
  'MN-NAT': 1.94976e-07,
  'MO-NAT': 1.89291e-06,
  'NI-NAT': 5.32188e-05,
  'O-NAT': 0.0235673,
  'ZR-NAT': 0.000418372,
  'density': 705.0087593020846,
  'name': 'WATER',
  'num': 1,
  'volume': 20}]
>> a.fuel[0]
{'NP237': 6.15876006672e-08,
 'O-NAT': 0.081406,
 'PU238': 1.39560431963e-05,
 'PU239': 0.000379067102804,
 'PU240': 0.000285451760087,
 'PU241': 0.000110469214146,
 'PU242': 5.75255886075e-05,
 'U235': 4.69701906177e-05,
 'U236': 2.55791248657e-06,
 'U238': 0.023236823207,
 'density': 11704.586234243618,
 'enrich': 0.11539736780507578,
 'name': 'COMBUS10',
 'num': 1,
 'volume': 10}
```
density unit: kg/m3
enrich unit %
And it contains:
* a.boron_concentration   
a list obj, boron concentration for each water composition
* a.concentration         
a list obj, element is dict, all compositions' concetration
* a.geoCom                
a list obj, 
* a.isotopes              
a set, all isotopes of fuel
* a.max_sigma             
a num, the max_sigma for all type score
* a.reaction              
a set, all reaction types
* a.realVol               
a dict: vol_num:real Vol of this volume
* a.volume                
a set: all volume num
* a.uranium_concentration 
a list: enrich for each fuel composition
* a.score                 
a 3-level dict of all score info:a.score[vol_num][reaction type][isotope]->score of this reaction
* a.spectrum              
a 4-level dict: a.score[vol_num][reaction type][isotope]->the spectrum of this reaction
* a.total                 
a 5-level dict: all info in this obj

## atom_mass
a module for getting the atom mass
### Example
```python
>> a=atom_mass()
>> a["H"]
1.00794075382579     the nautral atom mass
>> a["H-NAT"]
1.00794075382579     the nautral atom mass equivalent to a["H"]
>> a.getKg("H1")
1.673532811125249e-27  mass in kg
>> a.amu.H1
1.007825032            amu mass of H1
>> a.kg.U235
3.902996103592073e-25    mass in kg of U235
```

## postProcess
this module is for make some table of a result file

## TripoliE
this is a Tripoli-4 Evolution module
It can compute a nest burnup concentration of fuel and luanch a new calcule.

## TripoliHandle
this obj is used for a multiProcess of Tripoli evolution
