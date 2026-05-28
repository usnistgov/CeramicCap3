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
        self.fsamp  = int(self.meas['FSAMP'])
        self.nsamp  = int(self.meas['NSAMP'])
        self.datadir    = self.meas['DATADIR']
        self.logdir     = self.meas['LOGDIR']
        self.rawdatadir = self.meas['RAWDATADIR']
        self.saverawdata = self.cp['MEAS'].getboolean('saverawdata', fallback=False)
        self.nwarmup = int(self.meas.get('NWARMUP', '2'))
        self.fixg2 = self.cp['MEAS'].getboolean('fixg2', fallback=False)
        self.sat_threshold = float(self.meas.get('SATTHRESHOLD', '10.0'))
        self.switching = self.cp['MEAS'].getboolean('switching', fallback=True)
        self.max_nhars = int(self.meas.get('MAXNHARS', '10'))
        self.version = self.meas.get('VERSION', 'unknown')

        k1='CONFIG'
        for k2 in list(self.cp[k1].keys()):
            k2u = k2.upper()
            if k2u=='TOPLEFT':
                self.C31 = self.capdict[self.cp[k1][k2].upper()]
            elif k2u=='BOTTOMLEFT':
                self.C32 = self.capdict[self.cp[k1][k2].upper()]
            elif k2u=='TOPRIGHT':
                 self.C41 = self.capdict[self.cp[k1][k2].upper()]
            elif k2u=='BOTTOMRIGHT':
                self.C42 = self.capdict[self.cp[k1][k2].upper()]
            elif k2u=='SNTOPLEFT':
                self.SN31 = self.cp[k1][k2].upper()
            elif k2u=='SNBOTTOMLEFT':
                self.SN32 = self.cp[k1][k2].upper()
            elif k2u=='SNTOPRIGHT':
                self.SN41 = self.cp[k1][k2].upper()
            elif k2u=='SNBOTTOMRIGHT':
                self.SN42 = self.cp[k1][k2].upper()


    def save(self):
        with open(self.cfgpath,'w') as configfile:
            self.cp.write(configfile)

    def gain1(self,button, checked):
        if checked:
            k1='MEAS'
            k2='gain1'
            self.cp[k1][k2]=button.text()
            self.save()

    def gain2(self,button, checked):
        if checked:
            k1='MEAS'
            k2='gain2'
            self.cp[k1][k2]=button.text()
            self.save()

    def setGains(self, g1, g2):
        self.cp['MEAS']['gain1'] = str(int(g1))
        self.cp['MEAS']['gain2'] = str(int(g2))
        self.save()



# %%
if __name__=="__main__":
    config =CCC()
    for k,v in config.meas.items():
        print(k,v)
    print(f"{config.C31=} {config.C32=} {config.C41=} {config.C42=}")
    print(f"{config.SN31=} {config.SN32=} {config.SN41=} {config.SN42=}")
    print(config.recapdir[config.C31])
    print(f"{config.flist}")
    print(f"{config.datadir}")


# %%
