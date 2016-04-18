from distutils.core import setup
setup(name="TripoliE",version="2.0.1",
      decription="A Tripoli post-process/evolution package",
      author="Neo Wang",
      py_modules=[
                  "TripoliE/__init__",
                  "TripoliE/atom_mass",
                  "TripoliE/BaseTripoli",
                  "TripoliE/lancejob",
                  "TripoliE/TripoliE",
                  "TripoliE/TripoliHandle",
                  "TripoliE/postProcess"
                  ],
      packages=['TripoliE'],
      package_dir={'TripoliE': 'TripoliE'},
      package_data={'TripoliE': ['abundance','mass_rmd.mas95']},)
