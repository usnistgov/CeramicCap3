#%%
import configparser
import os
 
class CCC():
    def __init__(self):
        super().__init__()
        self.cwd =  os.getcwd()
        self.cfgpath=r"c:\python\CeramicCap3\data\config.ini"
        self.cp = configparser.ConfigParser()
        self.done=False
        self.meas={}
        self.config={}
        self.read()

    def read(self):
        self.cp.read(self.cfgpath)
        k1='MEAS'
        for k2 in list(self.cp[k1].keys()):
            self.meas[k2.upper()]=self.cp[k1][k2]
        k1='CONFIG'
        for k2 in list(self.cp[k1].keys()):
            self.config[k2.upper()]=self.cp[k1][k2]
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


# %%
