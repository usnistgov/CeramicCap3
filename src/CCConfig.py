#%%
import configparser

class CCC():
    def __init__(self):
        self.capdict={'1PF':1e-12,'10PF':1e-11,'100PF':1e-10,'1000PF':1e-9,'1NF':1e-9,'10NF':1e-8,'100NF':1e-7,'1000NF':1e-6,'1UF':1e-6,'10UF':1e-5}
        self.recapdir={1e-11:'10pF',1e-10:'100pF',1e-9:'1nF',1e-8:'10nF',1e-7:'100nF', 1e-6:'1uF',1e-5:'10uF'}
        self.cfgpath=r"c:\python\CeramicCap3\data\config.ini"
        self.cp = configparser.ConfigParser()
        self.meas={}
        self.read()

    def read(self):
        self.cp.read(self.cfgpath)
        k1='MEAS'
        for k2 in list(self.cp[k1].keys()):
            self.meas[k2.upper()]=self.cp[k1][k2]
        self.flist = [float(i) for i in self.meas['FREQLIST'].split(',')]
        self.fstart = float(self.meas['FSTART'])
        self.dvfrac = float(self.meas['DVFRAC'])
        self.nrMeas = int(self.meas['MEASPERFREQ'])
        self.fsamp     = int(self.meas['FSAMP'])
        self.forcefmax = self.cp['MEAS'].getboolean('forcefmax', fallback=False)
        self.datadir    = self.meas['DATADIR']
        self.logdir     = self.meas['LOGDIR']
        self.rawdatadir = self.meas['RAWDATADIR']
        self.saverawdata = self.cp['MEAS'].getboolean('saverawdata', fallback=False)
        self.eta3_limit = float(self.meas.get('ETA3LIMIT', '0.01'))
        self.max_step_frac = float(self.meas.get('MAXSTEPFRAC', '0.50'))
        self.decay = float(self.meas.get('DECAY', '0.85'))
        self.nellipse = int(self.meas.get('NELLIPSE', '8'))
        self.fixg2 = self.cp['MEAS'].getboolean('fixg2', fallback=False)
        self.sat_threshold = float(self.meas.get('SATTHRESHOLD', '10.0'))
        self.max_nhars = int(self.meas.get('MAXNHARS', '10'))
        self.version = self.meas.get('VERSION', 'unknown')

        k1='CONFIG'
        for k2 in list(self.cp[k1].keys()):
            k2u = k2.upper()
            if k2u=='C1':
                self.C1 = self.capdict[self.cp[k1][k2].upper()]
            elif k2u=='C2':
                self.C2 = self.capdict[self.cp[k1][k2].upper()]
            elif k2u=='SNC1':
                self.SN1 = self.cp[k1][k2].upper()
            elif k2u=='SNC2':
                self.SN2 = self.cp[k1][k2].upper()


    def save(self):
        with open(self.cfgpath,'w') as configfile:
            self.cp.write(configfile)

    def gain2(self,button, checked):
        if checked:
            k1='MEAS'
            k2='gain2'
            self.cp[k1][k2]=button.text()
            self.save()

    def setGains(self, g2):
        self.cp['MEAS']['gain2'] = str(int(g2))
        self.save()



# %%
if __name__=="__main__":
    config =CCC()
    for k,v in config.meas.items():
        print(k,v)
    print(f"{config.C1=} {config.C2=}")
    print(f"{config.SN1=} {config.SN2=}")
    print(f"{config.flist}")
    print(f"{config.datadir}")


# %%
