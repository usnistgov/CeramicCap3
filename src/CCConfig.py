#%%
import configparser
import os
 
class CCC():
    def __init__(self):
        self.capdict={'100PF':1e-10,'1000PF':1e-9,'1NF':1e-9,'10NF':1e-8,'100NF':1e-7,'1000NF':1e-6,'1UF':1e-6,'10UF':1e-5}        
        self.recapdir={1e-10:'100pF',1e-9:'1nF',1e-8:'10nF',1e-7:'100nF', 1e-6:'1uF',1e-5:'10uF'}
        self.cwd =  os.getcwd()
        self.cfgpath=r"c:\python\CeramicCap3\data\config.ini"
        self.cp = configparser.ConfigParser()
        self.done=False
        self.meas={}
        self.read()

    def read(self):
        self.cp.read(self.cfgpath)
        k1='MEAS'
        for k2 in list(self.cp[k1].keys()):
            self.meas[k2.upper()]=self.cp[k1][k2]
        self.finelist = [float(i) for i in self.meas['FINELIST'].split(',')]
        self.coarselist = [float(i) for i in self.meas['COARSELIST'].split(',')]
        if bool(self.meas['USEFINE'])==True:
            self.flist = self.finelist
        else:
            self.flist = self.coarselist
        self.fstart=float(self.meas['FSTART'])
        self.nrMeas=int(self.meas['MEASPERFREQ'])
        self.datadir =self.meas['DATADIR']

        k1='CONFIG'
        for k2 in list(self.cp[k1].keys()):
            if k2.upper()=='TOPLEFT':
                    self.C31 = self.capdict[self.cp[k1][k2].upper()]
            if k2.upper()=='BOTTOMLEFT':
                    self.C32 = self.capdict[self.cp[k1][k2].upper()]
            if k2.upper()=='TOPRIGHT':
                    self.C41 = self.capdict[self.cp[k1][k2].upper()]
            if k2.upper()=='BOTTOMRIGHT':
                    self.C42 = self.capdict[self.cp[k1][k2].upper()]
        for k2 in list(self.cp[k1].keys()):
            if k2.upper()=='SNTOPLEFT':
                    self.SN31 = self.cp[k1][k2].upper()
            if k2.upper()=='SNBOTTOMLEFT':
                    self.SN32 = self.cp[k1][k2].upper()
            if k2.upper()=='SNTOPRIGHT':
                    self.SN41 = self.cp[k1][k2].upper()
            if k2.upper()=='SNBOTTOMRIGHT':
                    self.SN42 = self.cp[k1][k2].upper()

            
        #for k2 in list(self.cp[k1].keys()):
            #self.config[k2.upper()]=self.cp[k1][k2]
        self.done=True


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

    def setNrMeasBeforeVadj(self,value):
            k1='MEAS'
            k2='NrMeasBeforeVadj'
            self.cp[k1][k2]=f"{value:.0f}"
            self.save()

    def setNrMeasBeforefadj(self,value):
            k1='MEAS'
            k2='NrMeasBeforefadj'
            self.cp[k1][k2]=f"{value:.0f}"
            self.save()



# %%
if __name__=="__main__":
    config =CCC()
    for k,v in config.meas.items():
        print(k,v)
    print(f"{config.C31=} {config.C32=} {config.C41=} {config.C42=}")
    print(f"{config.SN31=} {config.SN32=} {config.SN41=} {config.SN42=}")
    print(config.recapdir[config.C31])
    print(f"{config.finelist}")
    print(f"{config.coarselist}")
    print(f"{config.flist}")
    print(f"{config.datadir}")


# %%
