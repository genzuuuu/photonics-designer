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

#pyswarms testing
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

class genericsim:
    def __init__(self, ParamName, ObjectiveFunction, db, center_wavelength, bandwidth, SimParams, ComponentParams, Width_MMI = 2.5, Length_MMI=None, 
                 mmiid=None, Gap_MMI=None, Taper_Length=10, Taper_Width=1, count=0, **kwargs):
        
        #defaults    
        self.parameters = SimParams

        self.dbpath = db
        self.component = None
        self.layer_stack=get_layer_stack()
        
        #MMI Params
        self.ComponentParams = ComponentParams
        self.Name = ParamName
        self.count=count
        self.ObjectiveFunction = ObjectiveFunction
        
        #outputs & targets -> may need to change for default device
        self.sparam = None
        self.IL_SR = None
        self.target_CW = None
        self.converged = False
        self.center_wavelength = center_wavelength
        self.start_bandwidth = center_wavelength - bandwidth/2
        self.stop_bandwidth = center_wavelength + bandwidth/2
        self.mean_IL = None
        self.mean_SR = None
        self.IL_center = None
        self.SR_center = None
        self.file_path = None
        
        #search database for design
        if self.search_database(): 
            print("Found useful design!")
            print(self.ComponentParams) 
            exit() #Do something else generally

            #verify design
            

    def draw_gds(self):

        #gf.config.rich_output()
        PDK = get_generic_pdk()
        PDK.activate()

        #initialize the component
        self.component = gf.components.cells[self.Name](self.ComponentParams)

    def run(self):
        #Set lumerical session
        s = lumapi.FDTD(hide=True)

        #run FDTD simulations
        a = WL(self.component, run=True, session=s,  count=self.count, **self.parameters)
        self.count += 1
        
        # get S parameters from the simulation
        self.sparam = s.getsweepresult("s-parameter sweep", "S parameters") 

        
    def insert_into_database(self): 
        #file_path = get_sparameters_path(
        #    component=self.component,
        #    layer_stack=self.layer_stack,
        #    **self.parameters,
        #)
        file_path= None
        conn = sqlite3.connect(self.dbpath)
        print("[INFO] : Successful connection!")

        # reads the biggest number of MMIID 
        cur = conn.cursor()

        #TODO create table for other components
        
    def search_database(self):
        # user has to provide a desired center wavelength to search the database
        conn = sqlite3.connect(self.dbpath)
        print("[INFO] : Successful connection!")

        # reads the biggest number of MMIID 
        cur = conn.cursor()

        #TODO check the correct table

    def alter_database_entry(self):
        conn = sqlite3.connect(self.dbpath)
        print("[INFO] : Successful connection!")
        # reads the biggest number of MMIID 
        cur = conn.cursor()

        #TODO modift correct table entry
        
    #runs full simulation and inserts into database
    def runall(self):
        self.draw_gds()
        self.run()
        self.ObjectiveFunction()
        #self.insert_into_database()

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
        #TODO select correct parameters to optimize
        
        self.draw_gds()
        self.run()
        self.ObjectiveFunction()

        print(self.mean_IL)

        #temporary fix to None and small values
        if (self.mean_IL == None): return 1e6
        if (self.mean_IL < -0.5): return 1e6
        return self.mean_IL

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
            #TODO select correct parameters to optimize
            
            self.draw_gds()
            self.run()
            self.ObjectiveFunction()

            print(f"mean IL{self.mean_IL}")

            #temporary fix to None and highly negative solutions
            if (self.mean_IL == None): loss_x = 1e6
            elif (self.mean_IL < -0.5): loss_x = 1e6
            else: loss_x = self.mean_IL
            loss_arr.append(abs(loss_x)) #trims off negative numbers

        return np.asarray(loss_arr)

#Test code
if __name__ == '__main__':
    #running of an example
    #filepath = ""
    db = "./MMIDB.db"
    #dir_path = ""

