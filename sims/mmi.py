#import license
import lumapi
import matplotlib.pyplot as plt
import numpy as np
import pickle

import gdsfactory as gf
from gdsfactory.generic_tech import LAYER_STACK, get_generic_pdk

import gplugins.lumerical as sim

import re
import math

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

        a = sim.write_sparameters_lumerical(self.component, run=True, session=s, wavelength_points=5, wavelength_start=1.4, wavelength_stop=1.6)
        #check specs related to s_parameters

        with open('sim2.pk1','wb') as file: pickle.dump(a,file)

        #check if design is ok
    

    def process_line(line):
    # Check if line starts with special characters
    if re.match(r'^[\[\(]', line):
        return None
    else:
        # Split the line into numbers
        numbers = [float(num) for num in line.split()]
        return numbers

    def splitting_ratio_insertion_loss(filename, num_simulation):
        # filename: file path to the .dat file 
        # num_simulation: number of simulations
        # return: insertion_loss, splitting_ratio
        data = []
        with open(filename, 'r') as file:
            for line in file:
                line = line.strip()
                processed_line = process_line(line)
                if processed_line is not None:
                    data.append(processed_line)
        
        # Split the original array into chunks of 6 arrays each
        chunks = [data[i:i + num_simulation] for i in range(0, len(data), num_simulation)]

        # Assign each chunk to a new variable
        for i, chunk in enumerate(chunks, 1):
            globals()[f"result_{i}"] = chunk

        # Print the newly created variables
        for i in range(1, 10):
            print(f"result_{i}: {globals()[f'result_{i}']}") # S11, S12 ... S33

        T2 = []
        T3 = []
        insertion_loss = []
        splitting_ratio = []

        for x in range(6):
        T2_temp = result_4[x][1]
        T3_temp = result_7[x][1]
        T2.append(T2_temp)
        T3.append(T3_temp)
        insertion_loss.append([result_4[x][0],10*math.log10(T3_temp+T2_temp)])
        splitting_ratio.append([result_4[x][0],-10*math.log10(max(T2_temp,T3_temp)/min(T2_temp,T3_temp))])

        data_1 = {"insertion loss":[], "splitting ratio":[]}
        data_1.update({"insertion loss": insertion_loss, "splitting ratio": splitting_ratio})

        # output: {'insertion loss': [[wavelength,...], [wavelength,...]...], 'splitting ratio':[[wavelength,...], [wavelength,...]...]}
        return data_1