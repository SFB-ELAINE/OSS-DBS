# -*- coding: utf-8 -*-
"""
Created on Fri Nov  2 16:54:19 2018

@author: butenko
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Jun 26 15:10:10 2018

'''Data and method are adapted from Peterson2011'''

@author: konstantin
"""


import numpy as np
from pandas import read_csv

# This script takes potential in time on an axon and checks whether it is activated using modified driving force method (Peterson2011)
# Derived from "Predicting Myelinated Axon Activation Using Spatial Characteristics of the Extracellular Field" by Peterson, Izad and Tyler
# It should be parallelized in the future

def check_activation_by_MDF(Neuron_model,N_models_ext,n_Ranvier,fib_diam,t_step_extraction,pw,Ampl_scale,VTA_res): 

    #if not(isinstance(n_Ranvier,list)):
    #    n_Ranvier=[n_Ranvier]
    #    fib_diam=[fib_diam]
    
    last_point=0
    last_index=0

    Vert_get=read_csv('Neuron_model_arrays/Vert_of_Neural_model_NEURON.csv', delimiter=' ', header=None)    # get only physiologically correct neuron models 
    Vert=Vert_get.values
    Vert=np.round(Vert,8)
    
    for i in range(len(n_Ranvier)):

        if Neuron_model=="Reilly2016":
            n_segments=n_Ranvier[i]*2-1
            N_comp_in_between=2 #jump over one internodal
        elif "McIntyre2002":               
            from Axon_files.axon import Axon    
            param_ax={
            'centered':True,
            'diameter':fib_diam[i]
            }
            ax=Axon(param_ax)
            axon_dict=Axon.get_axonparams(ax)
            
            paranodes1=axon_dict["para1_nodes"]*(n_Ranvier[i]-1)/(21-1)
            paranodes2=axon_dict["para2_nodes"]*(n_Ranvier[i]-1)/(21-1)
            if axon_dict["fiberD"]>3.0:
                axoninter=(n_Ranvier[i]-1)*6
                N_comp_in_between=11
            else:
                axoninter=(n_Ranvier[i]-1)*3  
                N_comp_in_between=8
            n_segments=int(n_Ranvier[i]+paranodes1+paranodes2+axoninter)
    

        N_models=int(N_models_ext[i])
        n_Ranvier[i]=int(n_Ranvier[i])
        
        Nodes_status=np.zeros((N_models*n_segments,4),float)    #Nodes_status will contain info whether the placed(!) axon was activated
        
        
        Nodes_status[:,:3]=Vert[last_point:N_models*n_segments+last_point,:]
        last_point+=N_models*n_segments
            
        weights=np.genfromtxt('table_S1_new.csv',delimiter=' ')   # 5.7 µm, 90 µs
        weights_S3=np.genfromtxt('table_new.csv',delimiter=' ')    #5.7 µm, 90 µs, till 1V
          
        MDF_phi=np.zeros((N_models,n_Ranvier[i]),float)     #(number_of axons, number_of_nodes_of_Ranvier)
        Phi_nodes=np.zeros((N_models,n_Ranvier[i]),float)    
        
    
        for i_axon in range(N_models):      #goes through axons
                
            nodes=[]
            for point_inx in np.arange(0+last_index,n_segments+last_index,N_comp_in_between):         #goes over Ranviers

                nodes_point_in_time=np.load('Points_in_time/Signal_t_conv'+str(point_inx)+'.npy')
                nodes.append(abs(nodes_point_in_time[t_step_extraction])*(1000)*Ampl_scale)    #convert to mV    

            nodes=np.asarray(nodes)
            nodes = nodes.ravel()
            last_index=point_inx+1
            Phi_nodes[i_axon,:]=nodes
            
            
            j_node_max=(n_Ranvier[i]-1)      #number of Ranvier nodes and (-1 because indexing starts at 0)
            
            for i_node in range(n_Ranvier[i]):  
                #arr_around_nodes=np.arange(i_node-int(n_Ranvier[i]/2),i_node+int(n_Ranvier[i]/2)+1)     #10 nodes from left and right (gives there exact indicies)
                arr_around_nodes=np.arange(i_node-5,i_node+5+1)         #only 10 closest
                for j in arr_around_nodes:
                    if j>0 and j<j_node_max:   
                        #print(j)
                        MDF_phi[i_axon,i_node]=MDF_phi[i_axon,i_node]+weights[abs(j-i_node)]*(abs(Phi_nodes[i_axon,j-1]-2*Phi_nodes[i_axon,j]+Phi_nodes[i_axon,j+1]))

        N_axons_activ=0
        VTA=0.0

        for i_axon in range(N_models):
            for i_node in range(n_Ranvier[i]):  

                Ve1_index=int(Phi_nodes[i_axon,i_node]/25.0)
                Ve2_index=int(Phi_nodes[i_axon,i_node]/25.0)+1
                
                if Ve1_index>21:
                    Ve1_index=21
                    
                if Ve2_index>21:
                    Ve2_index=21   
                    Ve1_index=20

                MDF_threshold1=weights_S3[Ve1_index,3]
                if Phi_nodes[i_axon,i_node]%25==0.0:
                    MDF_threshold=MDF_threshold1                    
                else:
                    MDF_threshold2=weights_S3[Ve2_index,3]                
                    MDF_threshold=MDF_threshold1*(1-(Phi_nodes[i_axon,i_node]%25)/25.0)+MDF_threshold2*((Phi_nodes[i_axon,i_node]%25)/25.0)
                
                
                if MDF_phi[i_axon,i_node]>=MDF_threshold/1000.0:       #both in mV?
                    VTA+=VTA_res**3
                    Nodes_status[n_segments*i_axon:(n_segments*i_axon+n_segments),3]=1.0
                    N_axons_activ=N_axons_activ+1
                    #print ('Axon {} was activated'.format(i_axon))
                    #print ('MDF and MDF_threshold: {}     {} '.format(MDF_phi[i_axon,i_node], MDF_threshold))   
                    #print("==================")
                    break
        
        if len(n_Ranvier)==1:
            np.savetxt('Field_solutions/Activation/Neuron_model_results.csv', Nodes_status, delimiter=" ")
        else:
            np.savetxt('Field_solutions/Activation/Neuron_model_results_MDF_population_'+str(i)+'.csv', Nodes_status, delimiter=" ") 
        
        print("Number of axons activated: ",N_axons_activ)
   
    return True
