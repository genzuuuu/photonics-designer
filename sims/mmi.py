import matplotlib.pyplot as plt
import numpy as np
import pickle

#import lumapi



#import lumapi on MPI computers
import sys, os
sys.path.append("C:\\Program Files\\Lumerical\\v211\\api\\python\\") 
sys.path.append(os.path.dirname(__file__)) #Current directory

#import importlib.util
#spec_win = importlib.util.spec_from_file_location('lumapi', 'C:\\Program Files\\Lumerical\\v211\\api\\python\\lumapi.py')
#lumapi = importlib.util.module_from_spec(spec_win) #windows
#spec_win.loader.exec_module(lumapi)
import lumapi

import gdsfactory as gf
from gdsfactory.generic_tech import LAYER_STACK, get_generic_pdk

import gplugins.lumerical as sim
import sqlite3
import re
import math
import spicy


class mmi1x2:


    component = None

    def __init__(self,**kwargs):
        
        #defaults    
        self.parameters = dict(
            background_material = "sio2",
            port_margin = 1.5,
            port_extension = 5.0,
            mesh_accuracy = 1.5,
            wavelength_start = 1.4,
            wavelength_stop = 1.6,
            wavelength_points = 5,
            xmargin = 1,
            ymargin = 1,
            zmargin = 1,
            overwrite=True,
        ) 

        self.component = None
        self.num_port = 3
        self.sparam = None
        self.IL_SR = None

        self.mmiid = None
        self.Width_MMI = 3.8 
        self.Length_MMI = 12.8  
        self.Gap_MMI = 0.25 
        self.Taper_Length = 10 
        self.Taper_Width = 1.4  

        self.center_wavelength = None
        self.start_bandwidth = None
        self.stop_bandwidth = None
        self.mean_IL = None
        self.mean_sr = None
        self.IL_center = None
        self.SR_center = None
        self.file_path = None

        #check if an unacceptable parameter is passed in
        for settings in kwargs:
            if settings not in self.parameters:
                raise ValueError(
                    f"Invalid setting: {settings})"
                )   

        #update parameters
        self.parameters.update(kwargs)
        self.__dict__.update(self.parameters)
        print(self.parameters)
        #TODO
        #search database for design 
    
    def draw_gds(self):

        #gf.config.rich_output()
        PDK = get_generic_pdk()
        PDK.activate()

        #change such that it takes in some of the parameter values
        self.component = gf.components.cells["mmi1x2"]()

    def run(self, fname, dirpath=None):
        #ok checks
        #simulate gds to get s_parameters
        s = lumapi.FDTD(hide=True)

        a = sim.write_sparameters_lumerical(self.component, run=True, session=s, dirpath=dirpath,  **self.parameters)
        
        #check specs related to s_parameters
        
        with open(fname,'wb') as file: pickle.dump(a,file)
        
        # get S parameters from the simulation
        self.sparam = s.getsweepresult("s-parameter sweep", "S parameters")
        #print(self.sparam)
        #check if design is ok


    def splitting_ratio_insertion_loss(self):
        # num_simulation: number of simulations
        # return: insertion_loss, splitting_ratio
        num_simulation = self.parameters['wavelength_points']
        insertion_loss = []
        splitting_ratio = []
        for x in range(num_simulation):     
            print(self.sparam['S21'][x])
            T2_temp = abs(self.sparam['S21'][x])*abs(self.sparam['S21'][x])
            T3_temp = abs(self.sparam['S31'][x])*abs(self.sparam['S31'][x])
            insertion_loss.append([self.sparam['lambda'][x][0],10*math.log10(T3_temp+T2_temp)])
            splitting_ratio.append([self.sparam['lambda'][x][0],-10*math.log10(max(T2_temp,T3_temp)/min(T2_temp,T3_temp))])

            data_1 = {"insertion loss":[], "splitting ratio":[]}
            data_1.update({"insertion loss": insertion_loss, "splitting ratio": splitting_ratio})

            # output: {'insertion loss': [[wavelength,...], [wavelength,...]...], 'splitting ratio':[[wavelength,...], [wavelength,...]...]}
        print(data_1)
        self.IL_SR = data_1
        
    def insert_into_database(self, file_path, database_path): 
        num_simulation = self.parameters['wavelength_points']
        conn = sqlite3.connect(database_path)
        print("[INFO] : Successful connection!")

        # reads the biggest number of MMIID 
        cur = conn.cursor()
        sql_sel_max_query = '''SELECT MAX(MMIID) FROM MMI'''
        cur.execute(sql_sel_max_query)
        MMIID = cur.fetchall()[0][0]
        if (MMIID == None): 
            MMIID = 0
        else: 
            MMIID += 1
        self.mmiid = MMIID
        
        # find center wavelength
        temp = 1 / ( 2 * (self.IL_SR['insertion loss'][0][1]+self.IL_SR['splitting ratio'][0][1] ) )
        # print(temp)
        for i in range(num_simulation):
            if (  temp <= (1 / ( 2 * (self.IL_SR['insertion loss'][i][1]+self.IL_SR['splitting ratio'][i][1] ) ))   ):
                self.center_wavelength = self.IL_SR['insertion loss'][i][0]
                self.IL_center = self.IL_SR['insertion loss'][i][1]
                self.SR_center = self.IL_SR['splitting ratio'][i][1]
                temp = 1 / ( 2 * (self.IL_SR['insertion loss'][i][1]+self.IL_SR['splitting ratio'][i][1] ) )

        # find bandwidth
        for i in range(num_simulation):
            if (self.IL_SR['insertion loss'][i][1]>-0.5 and self.IL_SR['splitting ratio'][i][1]>-0.25):
                if (self.stop_bandwidth is None):
                    self.stop_bandwidth = self.IL_SR['insertion loss'][i][0]
                self.start_bandwidth = self.IL_SR['insertion loss'][i][0]

        # find mean IL & mean SR

        data_subset_IL = []
        data_subset_SR = []
        if (self.start_bandwidth is None or self.stop_bandwidth is None):
            self.mean_IL = None
            self.mean_sr = None
        else: 
            for i in range(num_simulation):
                if (self.IL_SR['insertion loss'][i][0] <= self.stop_bandwidth and self.IL_SR['insertion loss'][i][0] >= self.start_bandwidth ):
                    data_subset_IL.append(self.IL_SR['insertion loss'][i][1])
                    data_subset_SR.append(self.IL_SR['splitting ratio'][i][1])
            self.mean_IL = sum(data_subset_IL)/len(data_subset_IL)
            self.mean_SR = sum(data_subset_SR)/len(data_subset_SR)    

        self.file_path = file_path
        # insert .dat file path along with MMI specs into MMI table into a new row
        sql_insert_data_query = '''INSERT INTO MMI(MMIID, WidthMMI, LengthMMI, GapMMI, LengthTaper, WidthTaper, CenterWavelength, StartBandwidth, StopBandwidth, MeanIL, MeanSR, ILCenter, SRCenter,  FilePath)
        VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?);'''
        cur = conn.cursor()
        cur.execute(sql_insert_data_query, (MMIID, self.Width_MMI, self.Length_MMI, self.Gap_MMI, self.Taper_Length, self.Taper_Width,self.center_wavelength, self.start_bandwidth, self.stop_bandwidth, self.mean_IL, self.mean_SR, self.IL_center, self.SR_center,  file_path))
        conn.commit()
        print("[INFO] : ", file_path, "is in the database.") 
        print("[INFO] : This is entry number:", self.mmiid) 
        
    def search_database(self, database_path, center_wavelength, start_band, stop_band):
        # user has to provide a desired center wavelength to search the database
        conn = sqlite3.connect(database_path)
        print("[INFO] : Successful connection!")

        # reads the biggest number of MMIID 
        cur = conn.cursor()
        sql_insert_file_query = ''' SELECT * FROM (SELECT * FROM MMI WHERE %s>StartBandwidth and %s<StopBandwidth  ORDER BY ABS(CenterWavelength - %s) LIMIT 3) AS subquery_table  ORDER BY ILCenter+SRcenter ASC LIMIT 1 ; '''
        cur.execute(sql_insert_file_query, (center_wavelength))
        row = cur.fetchall()
        print("[INFO] : Successful Query!")
        print("MMIID, WidthMMI, LengthMMI, GapMMI, LengthTaper, WidthTaper, CenterWavelength, StartBandwidth, StopBandwidth, MeanIL, MeanSR, ILCenter, SRCenter,  FilePath")
        print(row)

    def alter_database_entry(self, database_path):
        conn = sqlite3.connect(database_path)
        print("[INFO] : Successful connection!")
        # reads the biggest number of MMIID 
        cur = conn.cursor()
        sql_edit_query = '''UPDATE MMI SET WidthMMI= %s, LengthMMI= %s, GapMMI= %s, LengthTaper= %s, WidthTaper= %s, CenterWavelength= %s, StartBandwidth= %s, StopBandwidth= %s, MeanIL= %s, MeanSR= %s, ILCenter= %s, SRCenter= %s,  FilePath= %s WHERE ID = %s; '''
        cur.execute(sql_edit_query, (self.Width_MMI, self.Length_MMI, self.Gap_MMI, self.Taper_Length, self.Taper_Width,self.center_wavelength, self.start_bandwidth, self.stop_bandwidth, self.mean_IL, self.mean_SR, self.IL_center, self.SR_center,  self.file_path, self.mmiid))
        conn.commit()
        print("[INFO] : Entry modified")


def fitness_function(input_param):
    # input_param[0] = length mmi, input_param[1] = gap mmi, input_param[2] = number of wavelength points
    
    # draw gds
    PDK = get_generic_pdk()
    PDK.activate()

    #change such that it takes in some of the parameter values
    component = gf.components.cells["mmi1x2"](length_mmi=input_param[0], gap_mmi=input_param[1],width_mmi=2.5, length_taper=10.0, width_taper=1.0)



    # run simulation
    #simulate gds to get s_parameters
    s = lumapi.FDTD(hide=True)

    a = sim.write_sparameters_lumerical(component, run=True, session=s)
    
    #check specs related to s_parameters
    

    # get S parameters from the simulation
    sparam = s.getsweepresult("s-parameter sweep", "S parameters")
    print(sparam)

    # calculate mean insertion loss
    insertion_loss = []
    for x in range(500):     
        T2_temp = abs(sparam['S21'][x])*abs(sparam['S21'][x])
        T3_temp = abs(sparam['S31'][x])*abs(sparam['S31'][x])
        insertion_loss.append(10*math.log10(T3_temp+T2_temp))
    print(insertion_loss)

    mean_IL = sum(insertion_loss)/len(insertion_loss)

    return mean_IL

##
if __name__ == '__main__':
    #running of an example
    

    for mesh in range(12):
        marg = 6
        
        c = mmi1x2(xmargin=marg, ymargin=marg, zmargin=marg, wavelength_points=5, mesh_accuracy=mesh+1, overwrite=False)#, Length_MMI = 12.8 , Gap_MMI = 0.25)
        c.draw_gds()
        c.run(fname=f"mesh{mesh+1}margin{marg}.pk1")
        c.splitting_ratio_insertion_loss()

    # trying optimization function
    #res = spicy.optimize.minimize(fitness_function, (5.5,0.25,100),method='COBYLA')
    #print(res)
