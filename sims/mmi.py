#import license
import lumapi
import matplotlib.pyplot as plt
import numpy as np
import pickle

import gdsfactory as gf
from gdsfactory.generic_tech import LAYER_STACK, get_generic_pdk

import gplugins.lumerical as sim

#remoteArgs = { "hostname": license.hostname,"port": license.port }

class mmi1x2:
    component = None
    parameters = None

    def __init__(self,parameters: dict):
        self.parameters = parameters
        self.component = None
        #search database for design
    
    def draw_gds(self):

        gf.config.rich_output()
        PDK = get_generic_pdk()
        PDK.activate()

        #change such that it takes in some of the parameter values
        self.component = gf.components.cells["mmi1x2"]()

    def run(self):
        #ok checks
        
        #simulate gds to get s_parameters
        s = lumapi.FDTD(hide=True, remoteArgs=remoteArgs)

        a = sim.write_sparameters_lumerical(self.component, run=True, session=s, wavelength_points=15, wavelength_start=1.1, wavelength_stop=1.3)
        #check specs related to s_parameters

        with open('sim2.pk1','wb') as file: pickle.dump(a,file)

        #check if design is ok