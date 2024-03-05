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
import scipy


class mmi1x2:
    def __init__(self, db, center_wavelength, bandwidth, Width_MMI = 2.5, Length_MMI=None, 
                 mmiid=None, Gap_MMI=None, Taper_Length=10, Taper_Width=1, count=0,**kwargs):
        
        #defaults    
        self.parameters = dict(
            background_material = "sio2",
            port_margin = 1.5,
            port_extension = 5.0,
            mesh_accuracy = 2,
            wavelength_start = center_wavelength- bandwidth/2,
            wavelength_stop = center_wavelength+ bandwidth/2,
            wavelength_points = 3,
            xmargin = 1,
            ymargin = 1,
            zmargin = 1,
            overwrite=False,
            dirpath = None,
        ) 

        self.dbpath = db
        self.component = None
        self.layer_stack=get_layer_stack()
        #self.num_port = 3
        
        #MMI Params
        self.mmiid = mmiid
        self.Width_MMI = Width_MMI
        self.Length_MMI = Length_MMI
        self.Gap_MMI = Gap_MMI
        self.Taper_Length = Taper_Length
        self.Taper_Width = Taper_Width
        self.count=count

        #outputs & targets
        self.sparam = None
        self.IL_SR = None
        self.target_CW = None
        self.converged = False
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
        
        #search database for design
        #entry = self.search_database(center_wavelength=self.target_CW)
        #if entry != None: converged = True #insert database entry into mmi values

    
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
        self.count += 1
        
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
        
    def insert_into_database(self): 
        #file_path = get_sparameters_path(
        #    component=self.component,
        #    layer_stack=self.layer_stack,
        #    **self.parameters,
        #)
        file_path= None
        num_simulation = self.parameters['wavelength_points']
        conn = sqlite3.connect(self.dbpath)
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
        
    def search_database(self, center_wavelength):
        # user has to provide a desired center wavelength to search the database
        conn = sqlite3.connect(self.dbpath)
        print("[INFO] : Successful connection!")

        # reads the biggest number of MMIID 
        cur = conn.cursor()
        sql_insert_file_query = ''' SELECT * FROM (SELECT * FROM MMI1x2 WHERE %s>StartBandwidth and %s<StopBandwidth  ORDER BY ABS(CenterWavelength - %s) LIMIT 3) AS subquery_table  ORDER BY ILCenter+SRcenter ASC LIMIT 1 ; '''
        cur.execute(sql_insert_file_query, (center_wavelength))
        row = cur.fetchall()
        print("[INFO] : Successful Query!")
        print("MMIID, WidthMMI, LengthMMI, GapMMI, LengthTaper, WidthTaper, CenterWavelength, StartBandwidth, StopBandwidth, MeanIL, MeanSR, ILCenter, SRCenter,  FilePath")
        print(row)

    def alter_database_entry(self):
        conn = sqlite3.connect(self.dbpath)
        print("[INFO] : Successful connection!")
        # reads the biggest number of MMIID 
        cur = conn.cursor()
        sql_edit_query = '''UPDATE MMI1x2 SET WidthMMI= %s, LengthMMI= %s, GapMMI= %s, LengthTaper= %s, WidthTaper= %s, CenterWavelength= %s, StartBandwidth= %s, StopBandwidth= %s, MeanIL= %s, MeanSR= %s, ILCenter= %s, SRCenter= %s,  FilePath= %s WHERE ID = %s; '''
        cur.execute(sql_edit_query, (self.Width_MMI, self.Length_MMI, self.Gap_MMI, self.Taper_Length, self.Taper_Width,self.center_wavelength, self.start_bandwidth, self.stop_bandwidth, self.mean_IL, self.mean_SR, self.IL_center, self.SR_center,  self.file_path, self.mmiid))
        conn.commit()
        print("[INFO] : Entry modified")
        
    #runs full simulation and inserts into database
    def runall(self):
        self.draw_gds()
        self.run()
        self.splitting_ratio_insertion_loss()
        self.insert_into_database()


    def fitness_function(self, input_param):
        # input_param[0] = length mmi, input_param[1] = gap mmi, input_param[2] = number of wavelength points
        print(input_param)
        #reset parameters
        self.IL_SR = None
        self.mean_IL = None
        self.mean_SR = None
        self.IL_center = None
        self.SR_center = None
        
        # draw gds
        self.Length_MMI = input_param[0]
        self.Gap_MMI = input_param[1]
        
        self.draw_gds()
        self.run()
        self.splitting_ratio_insertion_loss() #this may output nothing, make sure that error doesn't propogate
        self.insert_into_database()


        print(self.mean_IL)
        if (self.mean_IL == None): return 100000
        if (self.mean_IL < -0.5): return 100000
        return self.mean_IL

#Test code
if __name__ == '__main__':
    #running of an example
    #filepath = ""
    db = "./MMIDB.db"
    #dir_path = ""
    
    #running an example
    c = mmi1x2(db=db, center_wavelength=1.5, bandwidth=0.05,xmargin=1, ymargin=1, zmargin=1)
    
    #c.runall()

    # trying optimization function
    res = scipy.optimize.minimize(c.fitness_function, (5.5,0.25),method='COBYLA', options={"maxiter": 100})
    print(res)
