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

# pyswarms testing
import pyswarms as ps
from gdsfactory.config import PATH
from functools import partial

import gplugins.lumerical as sim
from new_write_sparameters import write_sparameters_lumerical as WL
from sims.component_opt import (
    ScipyOptMin as scipyminopt,
    particleswarm as swarmopt,
)

import sqlite3
import re
import math
import scipy
import time

wrk_dir = PATH.cwd / "extra"
wrk_dir.mkdir(exist_ok=True)

class mmi1x2:

    def __init__(self, db, wavelength_start,wavelength_stop,center_wavelength=None, bandwidth=None, start_bandwidth=None,stop_bandwidth=None, Width_MMI = 2.5, Length_MMI=None, 
                 mmiid=None, Gap_MMI=None, Taper_Length=10, Taper_Width=1, count=0,**kwargs):
        
        """
        Default set of parameters which define the simulation environment
        """

        self.parameters = dict(
            background_material = "sio2",
            port_margin = 1.5,
            port_extension = 5.0,
            mesh_accuracy = 1,
            #wavelength_start = center_wavelength- bandwidth/2,
            #wavelength_stop = center_wavelength+ bandwidth/2,
            wavelength_start = wavelength_start,
            wavelength_stop = wavelength_stop,
            wavelength_points = 3,
            xmargin = 1,
            ymargin = 1,
            zmargin = 1,
            overwrite=False,
            dirpath = None,
        ) 

        
        self.dbpath = db #Database filepath for components
        self.component = None #Component object for the MMI
        self.layer_stack=get_layer_stack()
        
        """
        Parameters for the MMI where mmiid defines the indentification of the MMI within the database
        """
        self.MMIparams = dict(
            mmiid = mmiid,
            Width_MMI = Width_MMI,
            Length_MMI = Length_MMI,
            Gap_MMI = Gap_MMI,
            Taper_Length = Taper_Length,
            Taper_Width = Taper_Width,
        )
        
        self.count=count
        
        #outputs & targets
        self.sparam = None
        self.IL_SR = None
        self.target_CW = None
        self.converged = False
        #self.center_wavelength = center_wavelength
        #self.start_bandwidth = center_wavelength - bandwidth/2
        #self.stop_bandwidth = center_wavelength + bandwidth/2
        self.start_bandwidth = start_bandwidth
        self.stop_bandwidth = stop_bandwidth
        self.bandwidth = stop_bandwidth - start_bandwidth
        self.center_wavelength = (wavelength_stop + wavelength_start)/2
        self.mean_IL = None
        self.mean_SR = None
        self.IL_center = None
        self.SR_center = None
        self.filepath = None

        """Sanity check to ensure simulation settings are legal according to gdsfactory settings"""
        for settings in kwargs:
            if settings not in self.parameters:
                raise ValueError(
                    f"Invalid setting: {settings})"
                )   

        #update simulation parameters and print it
        self.parameters.update(kwargs)
        self.__dict__.update(self.parameters)
        print(self.parameters)
        
        #TODO search database for design
        """
        if self.search_database(): 
            print("Found useful design!")
            print(self.MMIparams) 
            exit() #Do something else generally

            #verify design
        """   

    """Handles the initialization of the gds from gdsfactory default components, using MMIparams"""
    def draw_gds(self):

        #gf.config.rich_output()
        PDK = get_generic_pdk()
        PDK.activate()

        #initialize the component
        self.component = gf.components.cells["mmi1x2"](width_mmi=self.MMIparams["Width_MMI"], length_mmi=self.MMIparams["Length_MMI"], gap_mmi=self.MMIparams["Gap_MMI"],
                                                       length_taper=self.MMIparams["Taper_Length"], width_taper=self.MMIparams["Taper_Width"])

    def run(self):
        #Sets lumerical session
        s = lumapi.FDTD(hide=True)

        #run FDTD simulations through a gdsfacatory helper function
        a, self.filepath = WL(self.component, run=True, session=s,  count=self.count, **self.parameters)
        self.count += 1 #Do not remove as of right now, count allows for multiple simulations to run
        
        # get S parameters from the simulation
        self.sparam = s.getsweepresult("s-parameter sweep", "S parameters")
        

    """Calculates the splitting ratio and insertion loss from sparam matrix, to be run after the simulation has been run"""
    def splitting_ratio_insertion_loss(self):
        # num_simulation: number of simulations
        # return: insertion_loss, splitting_ratio
        num_simulation = self.parameters['wavelength_points']
        insertion_loss = []
        splitting_ratio = []
        for x in range(num_simulation):     
            T2_temp = abs(self.sparam['S21'][x])*abs(self.sparam['S21'][x])
            T3_temp = abs(self.sparam['S31'][x])*abs(self.sparam['S31'][x])
            insertion_loss.append([self.sparam['lambda'][x][0],10*math.log10(T3_temp+T2_temp)])
            splitting_ratio.append([self.sparam['lambda'][x][0],-10*math.log10(max(T2_temp,T3_temp)/min(T2_temp,T3_temp))])

            data_1 = {"insertion loss":[], "splitting ratio":[]}
            data_1.update({"insertion loss": insertion_loss, "splitting ratio": splitting_ratio})

            # output: {'insertion loss': [[wavelength,...], [wavelength,...]...], 'splitting ratio':[[wavelength,...], [wavelength,...]...]}
        self.IL_SR = data_1

    """Update class with values calculated from splitting_ratio_insertion_loss values"""
    def updateClass(self):
        num_simulation = self.parameters['wavelength_points']
        # find center wavelength
        temp = 1 / ( 2 * (self.IL_SR['insertion loss'][0][1]+self.IL_SR['splitting ratio'][0][1] ) )
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
            if (not data_subset_IL):
                self.mean_IL = None
            else:
                self.mean_IL = sum(data_subset_IL)/len(data_subset_IL)
            
            if (not data_subset_SR):
                self.mean_SR = None
            else:
                self.mean_SR = sum(data_subset_SR)/len(data_subset_SR)    

    """Insert component intro database"""
    def insert_into_database(self): 
        conn = sqlite3.connect(self.dbpath)
        print("[INFO] : Successful connection!")

        # reads the biggest number of DeviceID 
        cur = conn.cursor()
        sql_sel_max_query = '''SELECT MAX(DeviceID) FROM DevicesTable'''
        cur.execute(sql_sel_max_query)
        MMIID = cur.fetchall()[0][0]

        if (MMIID == None): 
            MMIID = 0
        else: 
            MMIID += 1
        self.MMIparams["mmiid"] = MMIID

        # insert .dat file path along with MMI specs into MMI table into a new row
        sql_insert_data_query = '''INSERT INTO DevicesTable(DeviceID, Type, Parameter1, Parameter2, Parameter3, Parameter4,Parameter5, CenterWavelength, StartBandwidth, StopBandwidth, MeanIL, MeanSR, ILCenter, SRCenter,  FilePath, function_call)
        VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?);'''
        cur = conn.cursor()
        cur.execute(sql_insert_data_query, (MMIID, "MMI 1x2", f'''Width_MMI = {self.MMIparams["Width_MMI"]} ''', f'''Length_MMI = {self.MMIparams["Length_MMI"]}''', 
                                            f'''Gap_MMI = {self.MMIparams["Gap_MMI"]}''', f'''Taper_Length = {self.MMIparams["Taper_Length"]}''', 
                                            f'''Taper_Width = {self.MMIparams["Taper_Width"]}''', self.center_wavelength*1000, 
                                            self.start_bandwidth*1000, self.stop_bandwidth*1000, self.mean_IL, self.mean_SR, self.IL_center, self.SR_center,                                        
                                            self.filepath, f'''c = mmi1x2(db="sims/Devices-simulation.db", start_bandwidth={self.start_bandwidth}, stop_bandwidth={self.stop_bandwidth},wavelength_start={self.wavelength_start}, wavelength_stop={self.wavelength_stop}, Width_MMI={self.MMIparams["Width_MMI"]}, Length_MMI={self.MMIparams["Length_MMI"]}, Gap_MMI={self.MMIparams["Gap_MMI"]}, Taper_Length={self.MMIparams["Taper_Length"]}, Taper_Width={self.MMIparams["Taper_Width"]})''' ))
        conn.commit()
        conn.close()
        print("[INFO] : This MMI 1x2 is in the database.") 
        print("[INFO] : This is entry number:", self.MMIparams["mmiid"]) 
        
    """Search database for component given initial parameters"""
    def search_database(self):
        # user has to provide a desired center wavelength to search the database
        conn = sqlite3.connect(self.dbpath)
        print("[INFO] : Successful connection!")

        # reads the biggest number of DeviceID 
        cur = conn.cursor()
        sql_insert_file_query = f''' SELECT * FROM (SELECT * FROM DevicesTable WHERE {self.start_bandwidth*1000}>=StartBandwidth and {self.stop_bandwidth*1000}<=StopBandwidth and Type="MMI 1x2"  ORDER BY ABS(CenterWavelength - {self.center_wavelength*1000}) LIMIT 3) AS subquery_table  ORDER BY ILCenter ASC LIMIT 1 ; '''
        cur.execute(sql_insert_file_query)
        row = cur.fetchall()
        print("[INFO] : Successful Query!")
        print("DeviceID, Type, Parameter1, Parameter2, Parameter3, Parameter4,Parameter5, CenterWavelength, StartBandwidth, StopBandwidth, MeanIL, MeanSR, ILCenter, SRCenter,  FilePath, function_call")
        print(row)
        conn.close()

    """Update existing database entry with component"""
    def alter_database_entry(self):
        conn = sqlite3.connect(self.dbpath)
        print("[INFO] : Successful connection!")
        # reads the biggest number of DeviceID 
        cur = conn.cursor()
        sql_edit_query = '''UPDATE DevicesTable SET Parameter1= ?, Parameter2= ?, Parameter3= ?, Parameter4= ?, Parameter5= ?, CenterWavelength= ?, StartBandwidth= ?, StopBandwidth= ?, MeanIL= ?, MeanSR= ?, ILCenter= ?, SRCenter= ?,  FilePath= ? WHERE DeviceID = ?; '''
        cur.execute(sql_edit_query, (f'''Width_MMI = {self.MMIparams["Width_MMI"]} ''', f'''Length_MMI = {self.MMIparams["Length_MMI"]}''', 
                                     f'''Gap_MMI = {self.MMIparams["Gap_MMI"]}''', f'''Taper_Length = {self.MMIparams["Taper_Length"]}''', 
                                     f'''Taper_Width = {self.MMIparams["Taper_Width"]}''',
                                     self.center_wavelength*1000, self.start_bandwidth*1000, self.stop_bandwidth*1000, self.mean_IL, self.mean_SR, self.IL_center, self.SR_center,  self.filepath, self.MMIparams["mmiid"]))
        conn.commit()
        print("[INFO] : Entry modified")
        conn.close()

    """Delete existing database entry with same mmiid"""
    def delete_database_entry(self):
        conn = sqlite3.connect(self.dbpath)
        print("[INFO] : Successful connection!")
        # reads the DeviceID of the current mmi1x2 object 
        cur = conn.cursor()
        
        if self.MMIparams["mmiid"] is not None:
            sql_delete_query = f'''DELETE FROM DevicesTable WHERE DeviceID={self.MMIparams["mmiid"]}; '''
            cur.execute(sql_delete_query)
            conn.commit()
            self.MMIparams["mmiid"] -= 1
            print("[INFO] : Entry deleted")
        conn.close()

    """Run entire workflow once and insert final component into the database"""
    def runall(self):
        self.draw_gds()
        self.run()
        self.splitting_ratio_insertion_loss()
        self.updateClass()
        self.insert_into_database()

    """Iteratively simulates to find an MMI with the lowest insertion loss value"""
    def search_space(self):
        guess_length = (4*3.89*(self.MMIparams["Width_MMI"]**2))/(12*self.center_wavelength)
        print(guess_length)
        final = scipyminopt(self, guess=(guess_length,0.25))
        print(final)

        self.MMIparams["Length_MMI"] = final["x"][0]
        self.MMIparams["Gap_MMI"] = final["x"][1]
        self.mesh_accuracy=7

        self.runall()


    # input_param[0] = length mmi, input_param[1] = gap mmi, input_param[2] = number of wavelength points
    def fitness_function_scipy(self, input_param):
        # input_param[0] = length mmi, input_param[1] = gap mmi, input_param[2] = number of wavelength points
        print(input_param)

        #reset parameters
        self.IL_SR = None
        self.mean_IL = None
        self.mean_SR = None
        self.IL_center = None
        self.SR_center = None
        
        # draw gds
        self.MMIparams["Length_MMI"] = input_param[0]
        self.MMIparams["Gap_MMI"] = input_param[1]
        
        self.draw_gds()
        self.run()
        self.splitting_ratio_insertion_loss() #this may output nothing, make sure that error doesn't propogate
        self.updateClass()

        print(self.mean_IL)

        #temporary fix to None and small values
        if (self.mean_IL == None): return 1e6
        if (self.mean_IL > 0): return 1e6
        return abs(self.mean_IL)

    # ## Define the trainable function for the PSO optimization
    def fitness_function_swarm(self, x):
        """Training step, or `trainable`, function for Ray Tune to run simulations and return results."""
        loss_arr = []

        for xi in x:
            print(xi)
            #reset params
            self.IL_SR = None
            self.mean_IL = None
            self.mean_SR = None
            self.IL_center = None
            self.SR_center = None

            # Component to optimize
            self.MMIparams["Length_MMI"] = xi[0]
            self.MMIparams["Gap_MMI"] = xi[1]

            self.draw_gds()
            self.run()
            self.splitting_ratio_insertion_loss() #this may output nothing, make sure that error doesn't propogate
            self.updateClass()

            print(f"mean IL{self.mean_IL}")

            #temporary fix to None and highly negative solutions
            if (self.mean_IL == None): loss_x = 1e6
            elif (self.mean_IL > 0): loss_x = 1e6
            else: loss_x = abs(self.mean_IL)
            loss_arr.append(abs(loss_x)) #trims off negative numbers

        return np.asarray(loss_arr)
    

#Test code
if __name__ == '__main__':
    db = "sims/Devices-simulation.db" 
    
    #running an optimization example
    #initialize class
    c = mmi1x2(db=db, wavelength_start = 1.53, wavelength_stop=1.565,start_bandwidth = 1.53, stop_bandwidth = 1.565, xmargin=1, ymargin=1, zmargin=1, Width_MMI=3.8, Length_MMI=12.8, Gap_MMI=0.25, Taper_Length=10.0, Taper_Width=1.4)

    #c.insert_into_database()
    #c.delete_database_entry()
    #print(c.parameters())

    #c.search_space() #search_space testing
    #scipyminopt(c)

    #running pyswarm example
    #c = mmi1x2(db=db, center_wavelength=1.5, bandwidth=0.05,xmargin=1, ymargin=1, zmargin=1)
    #swarmopt(c)
