#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 29 13:02:14 2020

@author: butenko
"""

import numpy as np
from pandas import read_csv

# Get VTA planes from neuron array
# we need only the index of the central node

# works only for axonal arrays, not VTA arrays
def get_VTA_planes(study_number,axon_param,Axon_Model_Type):

    [ranvier_nodes, para1_nodes, para2_nodes, inter_nodes, ranvier_length, para1_length, para2_length, inter_length, deltax, fiberD]=(axon_param[:])  
    if Axon_Model_Type == 'McIntyre2002': 
        N_segm=((ranvier_nodes-1)+inter_nodes+para1_nodes+para2_nodes)/(ranvier_nodes-1)
        #N_segm=int((N_Ranv-1)*n_comp+1)        #overall number of points on Axon incl. internodal
    elif Axon_Model_Type == 'Reilly2016': 
        N_segm=2
       # N_segm=int((N_Ranv-1)*n_comp+1)        #overall number of points on Axon incl. internodal    
    
    N_Ranvier=ranvier_nodes

    i_central_node=int(N_Ranvier/2)*N_segm        # use odd number of nodes for VTA

    Vert_get=read_csv('Field_solutions/Activation/Neuron_model_results.csv', delimiter=' ', header=None)    # get only physiologically correct neuron models 
    Vert=Vert_get.values
    Vert=np.round(Vert,8)     

    VTA_planes_points=[]
    VTA_planes_activation=[]
    loc_counter=0
    
    for i in range(Vert.shape[0]):
        if loc_counter==i_central_node:
            VTA_planes_points.append(Vert[i,:3])                    #might be too slow if a lot of axons
            VTA_planes_activation.append(Vert[i,:])
        
        loc_counter+=1
        if loc_counter==(N_Ranvier-1)*N_segm+1:
            loc_counter=0
            
    VTA_planes_array=np.array(VTA_planes_points)
    VTA_planes_activation=np.array(VTA_planes_activation)
    
    np.savetxt('VTA_seeding_planes.csv', VTA_planes_array, delimiter=" ")
    np.savetxt('VTA_planes_activation.csv', VTA_planes_activation, delimiter=" ")   #you can plot single planes just by taking certain slices (defined by the number of axons in plane), for example, for the first study the second plane  is np.savetxt('VTA_2nd_seeding_plane.csv', VTA_planes_activation[130:260], delimiter=" ")

    return True