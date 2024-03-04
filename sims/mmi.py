import matplotlib.pyplot as plt
import numpy as np
import pickle

#import lumapi on MPI computers
import sys, os
sys.path.append("C:\\Program Files\\Lumerical\\v211\\api\\python\\") 
sys.path.append(os.path.dirname(__file__)) #Current directory
import lumapi

import gdsfactory as gf
from gdsfactory.generic_tech import LAYER_STACK, get_generic_pdk
from gplugins.common.utils.get_sparameters_path import (
    get_sparameters_path_lumerical as get_sparameters_path,
)
from gdsfactory.pdk import get_layer_stack

import gplugins.lumerical as sim
from new_write_sparameters import write_sparameters_lumerical as WL

import sqlite3
import re
import math
import spicy


class mmi1x2:


    component = None

    def __init__(self, Width_MMI = 3.8, Length_MMI=12.8, mmiid=None, Gap_MMI=0.25, Taper_Length=10, Taper_Width=1.4, count=0,**kwargs):
        
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
            dirpath = None,
        ) 

        self.component = None
        self.layer_stack=get_layer_stack()
        self.num_port = 3
        self.sparam = None
        self.IL_SR = None
        
        self.mmiid = mmiid
        self.Width_MMI = Width_MMI
        self.Length_MMI = Length_MMI
        self.Gap_MMI = Gap_MMI
        self.Taper_Length = Taper_Length
        self.Taper_Width = Taper_Width
        self.count=count

        self.center_wavelength = None
        self.start_bandwidth = None
        self.stop_bandwidth = None
        self.mean_IL = None
        self.mean_SR = None
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
        self.component = gf.components.cells["mmi1x2"](width_mmi=self.Width_MMI, length_mmi=self.Length_MMI, gap_mmi=self.Gap_MMI,
                                                       length_taper=self.Taper_Length, width_taper=self.Taper_Width)

    def run(self):
        #ok checks
        #simulate gds to get s_parameters
        s = lumapi.FDTD(hide=True)

        a = WL(self.component, run=True, session=s,  count=self.count, **self.parameters)
        
        # get S parameters from the simulation
        self.sparam = s.getsweepresult("s-parameter sweep", "S parameters")
        
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
        self.IL_SR = data_1
        
    def insert_into_database(self, database_path): 
        #file_path = get_sparameters_path(
        #    component=self.component,
        #    layer_stack=self.layer_stack,
        #    **self.parameters,
        #)
        file_path= None
        num_simulation = self.parameters['wavelength_points']
        conn = sqlite3.connect(database_path)
        print("[INFO] : Successful connection!")

        # reads the biggest number of MMIID 
        cur = conn.cursor()
        sql_sel_max_query = '''SELECT MAX(MMIID) FROM MMI1x2'''
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
            self.mean_SR = None
        else: 
            for i in range(num_simulation):
                if (self.IL_SR['insertion loss'][i][0] <= self.stop_bandwidth and self.IL_SR['insertion loss'][i][0] >= self.start_bandwidth ):
                    data_subset_IL.append(self.IL_SR['insertion loss'][i][1])
                    data_subset_SR.append(self.IL_SR['splitting ratio'][i][1])
            self.mean_IL = sum(data_subset_IL)/len(data_subset_IL)
            self.mean_SR = sum(data_subset_SR)/len(data_subset_SR)    

        self.file_path = file_path
        # insert .dat file path along with MMI specs into MMI table into a new row
        sql_insert_data_query = '''INSERT INTO MMI1x2(MMIID, WidthMMI, LengthMMI, GapMMI, LengthTaper, WidthTaper, CenterWavelength, StartBandwidth, StopBandwidth, MeanIL, MeanSR, ILCenter, SRCenter,  FilePath)
        VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?);'''
        cur = conn.cursor()
        cur.execute(sql_insert_data_query, (MMIID, self.Width_MMI, self.Length_MMI, self.Gap_MMI, self.Taper_Length, self.Taper_Width,self.center_wavelength, self.start_bandwidth, self.stop_bandwidth, self.mean_IL, self.mean_SR, self.IL_center, self.SR_center,  self.
                                            file_path))
        conn.commit()
        print("[INFO] : ", file_path, "is in the database.") 
        print("[INFO] : This is entry number:", self.mmiid) 
        
    def search_database(self, database_path, center_wavelength, start_band, stop_band):
        # user has to provide a desired center wavelength to search the database
        conn = sqlite3.connect(database_path)
        print("[INFO] : Successful connection!")

        # reads the biggest number of MMIID 
        cur = conn.cursor()
        sql_insert_file_query = ''' SELECT * FROM (SELECT * FROM MMI1x2 WHERE %s>StartBandwidth and %s<StopBandwidth  ORDER BY ABS(CenterWavelength - %s) LIMIT 3) AS subquery_table  ORDER BY ILCenter+SRcenter ASC LIMIT 1 ; '''
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
        sql_edit_query = '''UPDATE MMI1x2 SET WidthMMI= %s, LengthMMI= %s, GapMMI= %s, LengthTaper= %s, WidthTaper= %s, CenterWavelength= %s, StartBandwidth= %s, StopBandwidth= %s, MeanIL= %s, MeanSR= %s, ILCenter= %s, SRCenter= %s,  FilePath= %s WHERE ID = %s; '''
        cur.execute(sql_edit_query, (self.Width_MMI, self.Length_MMI, self.Gap_MMI, self.Taper_Length, self.Taper_Width,self.center_wavelength, self.start_bandwidth, self.stop_bandwidth, self.mean_IL, self.mean_SR, self.IL_center, self.SR_center,  self.file_path, self.mmiid))
        conn.commit()
        print("[INFO] : Entry modified")
        
    #runs full simulation and inserts into database
    def runall(self, database_path):
        self.draw_gds()
        self.run()
        self.splitting_ratio_insertion_loss()
        self.insert_into_database(database_path=database_path)


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
    #filepath = ""
    databasepath = "./MMIDB.db"
    dir_path = ""
    
    #running an example
    c = mmi1x2(xmargin=1, ymargin=1, zmargin=1, 
               wavelength_points=5, mesh_accuracy=2, 
               overwrite=True, dirpath=dir_path, Width_MMI=3, 
               Length_MMI=7, Gap_MMI=0.2, Taper_Length=8, 
               Taper_Width=9,count=0)
    
    c.runall(database_path=databasepath)

    
    # trying optimization function
    #res = spicy.optimize.minimize(fitness_function, (5.5,0.25,100),method='COBYLA')
    #print(res)
