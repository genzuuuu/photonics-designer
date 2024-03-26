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
    def __init__(self, ParamName, ObjectiveFunction, SimParams, db, center_wavelength, bandwidth, 
                 ComponentParams, count=0, **kwargs):
        
        #defaults
        #TODO check for correctness of simparams  
        self.parameters = SimParams

        self.dbpath = db
        self.component = None
        self.layer_stack=get_layer_stack()
        
        #Component Params
        self.ComponentParams = ComponentParams
        self.Name = ParamName
        self.count=count
        self.ObjectiveFunction = ObjectiveFunction
        
        #outputs & targets -> may need to change for default device
        self.sparam = None
        self.target_CW = None
        self.converged = False
        self.center_wavelength = center_wavelength
        self.start_bandwidth = center_wavelength - bandwidth/2
        self.stop_bandwidth = center_wavelength + bandwidth/2
        self.filepath = None

        
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
        outputs = self.ObjectiveFunction()
        
        #TODO parse outputs to be inserted into database
        self.insert_into_database()

    #searches for an accurate design and verifies it afterwards
    def search_space(self):
        final = scipyminopt(self)
        print(final)

        #self.MMIparams["Length_MMI"] = final["x"][0]
        #self.MMIparams["Gap_MMI"] = final["x"][1]
        self.mesh_accuracy=7

        self.runall()

    def fitness_function_scipy(self, input_param):
        #TODO
        print(input_param)
        
        self.ComponentParams["Length_MMI"] = input_param[0]
        self.ComponentParams["Gao_MMI"] = input_param[1]

        self.draw_gds()
        self.run()
        output = scipyminopt()
        print(output)

        if (output == None): return 1e6
        if (output < -0.5): return 1e6
        return output

#Test code
if __name__ == '__main__':
    db = "./MMIDB.db"
    center_wavelength = 1.5
    bandwidth=0.05
    
    #running of 1x2MMI example

    #Convert these parameters to JSON's and fill it out here
    MMIparams = dict(
        mmiid=None,
        Width_MMI=None,
        Length_MMI=None,
        Gap_MMI=None,
        Taper_Length=None,
        Taper_Width=None,
    )

    parameters = dict(
        background_material = "sio2",
        port_margin = 1.5,
        port_extension = 5.0,
        mesh_accuracy = 1,
        wavelength_start = center_wavelength- bandwidth/2,
        wavelength_stop = center_wavelength+ bandwidth/2,
        wavelength_points = 3,
        xmargin = 1,
        ymargin = 1,
        zmargin = 1,
        overwrite=False,
        dirpath = None,
    ) 

    #TODO add objective function
    #g = genericsim(ParamName="mmi1x2", ComponentParams=MMIparams, ObjectiveFunction=, SimParams=parameters, db = db, 
    #               center_wavelength=1.5, bandwidth=0.05, xmargin=1,ymargin=1,zmargin=1)

