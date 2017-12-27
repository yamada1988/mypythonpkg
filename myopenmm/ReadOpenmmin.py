
def read(self, filename):
    with open(filename, 'rt') as f:
        kwds = ['integrator', 'dt', 'timestep', 'temperature', 'pressure', 'barostat', 'frition_const',
                'hydrogen_mass', 'recstep', 'writestep',  ]    
