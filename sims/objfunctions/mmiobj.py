import math

#Takes in parameters as input
def splitting_ratio_insertion_loss(self, parameters, sparam):
    # num_simulation: number of simulations
    # return: insertion_loss, splitting_ratio
    num_simulation = parameters['wavelength_points']
    insertion_loss = []
    splitting_ratio = []
    for x in range(num_simulation):     
        T2_temp = abs(sparam['S21'][x])*abs(sparam['S21'][x])
        T3_temp = abs(sparam['S31'][x])*abs(sparam['S31'][x])
        insertion_loss.append([sparam['lambda'][x][0],10*math.log10(T3_temp+T2_temp)])
        splitting_ratio.append([sparam['lambda'][x][0],-10*math.log10(max(T2_temp,T3_temp)/min(T2_temp,T3_temp))])

        data_1 = {"insertion loss":[], "splitting ratio":[]}
        data_1.update({"insertion loss": insertion_loss, "splitting ratio": splitting_ratio})

        # output: {'insertion loss': [[wavelength,...], [wavelength,...]...], 'splitting ratio':[[wavelength,...], [wavelength,...]...]}
    IL_SR = data_1

    temp = 1 / ( 2 * (IL_SR['insertion loss'][0][1]+IL_SR['splitting ratio'][0][1] ) )
    for i in range(num_simulation):
        if (  temp <= (1 / ( 2 * (IL_SR['insertion loss'][i][1]+IL_SR['splitting ratio'][i][1] ) ))   ):
            center_wavelength = IL_SR['insertion loss'][i][0]
            IL_center = IL_SR['insertion loss'][i][1]
            SR_center = IL_SR['splitting ratio'][i][1]
            temp = 1 / ( 2 * (IL_SR['insertion loss'][i][1]+IL_SR['splitting ratio'][i][1] ) )

    # find bandwidth
    for i in range(num_simulation):
        if (IL_SR['insertion loss'][i][1]>-0.5 and IL_SR['splitting ratio'][i][1]>-0.25):
            if (parameters["stop_bandwidth"] is None):
                parameters["stop_bandwidth"] = IL_SR['insertion loss'][i][0]
            parameters["start_bandwidth"] = IL_SR['insertion loss'][i][0]

    # find mean IL & mean SR
    data_subset_IL = []
    data_subset_SR = []
    if (parameters["start_bandwidth"] is None or parameters["stop_bandwidth"] is None):
        mean_IL = None
        mean_SR = None
    else: 
        for i in range(num_simulation):
            if (IL_SR['insertion loss'][i][0] <= parameters["stop_bandwidth"] and IL_SR['insertion loss'][i][0] >= parameters["start_bandwidth"] ):
                data_subset_IL.append(IL_SR['insertion loss'][i][1])
                data_subset_SR.append(IL_SR['splitting ratio'][i][1])
        if (not data_subset_IL):
            mean_IL = None
        else:
            mean_IL = sum(data_subset_IL)/len(data_subset_IL)
        
        if (not data_subset_SR):
            mean_SR = None
        else:
            mean_SR = sum(data_subset_SR)/len(data_subset_SR)    

    return mean_IL

        