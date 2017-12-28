import sys

class InpInfo:
    def __init__(self):
        self.InpDict = {'integrator': 'Langevin', 
                        'temperature': 300.0e0, 
                        'pressure': 1.0e0, 
                        'forcefield': None, 
                        'fileformat': None, 
                        'friction_const': 1, 
                        'dt': 0.0020, 
                        'nvtstep' : 5000, 
                        'nptstep' :50000, 
                        'mdstep' : 50000, 
                        'eqrecstep' : 500, 
                        'mdrecstep' : 500, 
                        'emflag' : True, 
                        'nvtflag' : True, 
                        'nptflag' : True, 
                        'nvtrecflag' : False, 
                        'nptrecflag' : False, 
                        'mdrecflag' : False, 
                        'simlogflag' : True, 
                        'verbose' : True, 
                        'reporterformat' : None, 
                        'constraints' : None, 
                        'nonbondedmethod' : 'PME', 
                        'nonbondedcutoff' : 1.20e0, 
                        'heavyhydrogen' : False, 
                        'removeCMMotion' : False, 
                        'platform' : None, 
                        'precision' : 'double'}

    def read(self, file):
        try:
            f = open(file, 'rt')
        except IOError:
            sys.exit('Input file {0} not founded.'.format(file))

        with open(file, 'rt') as f:
            InpDict = {line.split('=')[0].strip() : line.split('=')[1].strip() \
                           for line in f if "#" not in line}

        InpDict.update(self.InpDict)
        print(InpDict)
 
        self.integrator = InpDict['integrator']
        self.temperature = float(InpDict['temperature'])
        self.pressure = float(InpDict['pressure'])
        self.forcefield = InpDict['forcefield']
        self.fileformat = InpDict['fileformat']
        self.friction_const = float(InpDict['friction_const'])
        self.dt = float(InpDict['dt'])
        self.nvtstep = int(InpDict['nvtstep'])
        self.nptstep = int(InpDict['nptstep'])
        self.mdstep = int(InpDict['mdstep'])
        self.eqrecstep = int(InpDict['eqrecstep'])
        self.mdrecstep = int(InpDict['mdrecstep'])
        self.emflag = InpDict['emflag']
        self.nvtflag = InpDict['nvtflag']
        self.nptflag = InpDict['nptflag']
        self.mdrecflag = InpDict['mdrecflag']
        self.simlogflag = InpDict['simlogflag']
        self.verbos = InpDict['verbose']
        self.reporterformat = InpDict['reporterformat']
        self.constraints = InpDict['constraints']
        self.nonbondedmethod = InpDict['nonbondedmethod']
        self.nonbondedcutoff = InpDict['nonbondedcutoff']
        self.heavyhydrogen = InpDict['heavyhydrogen']
        self.removeCMMotion = InpDict['removeCMMotion'] 
        self.platform = InpDict['platform']
        self.precision = InpDict['precision']

