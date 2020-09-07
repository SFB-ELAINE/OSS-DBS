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
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import numpy as np
from pandas import read_csv

# This script takes potential in time on an axon and checks whether it is activated using activating function method (Butson2006)
# Derived from "Role of electrode design on the volume of tissue activated during deep brain stimulation" by Butson and McIntyre


def check_activation_by_activ_function(Neuron_model,N_models_ext,n_Ranvier,fib_diam,t_step_extraction,pw,Ampl_scale,VTA_res): 

    #if not(isinstance(n_Ranvier,list)):
    #    n_Ranvier=[n_Ranvier]
    #    fib_diam=[fib_diam]
    

    last_point=0
    last_index=0

    Vert_get=read_csv('Neuron_model_arrays/Vert_of_Neural_model_NEURON.csv', delimiter=' ', header=None)    # get only physiologically correct neuron models 
    Vert=Vert_get.values
    Vert=np.round(Vert,8)
    
    for i in range(len(n_Ranvier)):
        
        #print(n_Ranvier[i])
        #print(int(N_models_ext[i]))

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
            
        Butson_threshold=22.7*np.e**(-90.0*1.0/50.9)+4.3        # from the paper
        Phi_nodes=np.zeros((N_models,n_Ranvier[i]),float)    

        N_axons_activ=0
        VTA=0.0        
    
        for i_axon in range(N_models):      #goes through axons
                
            nodes=[]
            for point_inx in np.arange(0+last_index,n_segments+last_index,N_comp_in_between):         #goes over Ranviers

                nodes_point_in_time=np.load('Points_in_time/Signal_t_conv'+str(point_inx)+'.npy')
                nodes.append(nodes_point_in_time[t_step_extraction]*(1000)*Ampl_scale)    #convert to mV    

            nodes=np.asarray(nodes)
            nodes = nodes.ravel()
            last_index=point_inx+1
            Phi_nodes[i_axon,:]=nodes
                        
            for i_node in range(1,n_Ranvier[i]-1):  
                
                #if i_node==10:
                #    print(abs(Phi_nodes[i_axon,i_node-1]-2*Phi_nodes[i_axon,i_node]+Phi_nodes[i_axon,i_node+1]))
                
                if Phi_nodes[i_axon,i_node-1]-2*Phi_nodes[i_axon,i_node]+Phi_nodes[i_axon,i_node+1]>Butson_threshold:
                    N_axons_activ+=1
                    VTA+=VTA_res**3
                    Nodes_status[n_segments*i_axon:(n_segments*i_axon+n_segments),3]=1.0
                    break

            if i_axon==100:
                Second_der=np.zeros(n_Ranvier[i]-2,float)
                for i_node in range(1,n_Ranvier[i]-1):
                    Second_der[i_node-1]=Phi_nodes[i_axon,i_node-1]-2*Phi_nodes[i_axon,i_node]+Phi_nodes[i_axon,i_node+1]
                    
                plt.figure(1111221133)
                plt.plot(np.arange(1,n_Ranvier[i]-1),Second_der)
                #plt.xlim(0.000,d["T"]*5)
                plt.grid(True)
                plt.xlabel('n_Ranvier')
                plt.ylabel('Second Derivative')
                #plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
                plt.savefig('Images/Activ_function.png', format='png', dpi=750)


                    
        if len(n_Ranvier)==1:
            np.savetxt('Field_solutions/Activation/Neuron_model_results.csv', Nodes_status, delimiter=" ")
        else:
            np.savetxt('Field_solutions/Activation/Neuron_model_results_Butson_population_'+str(i)+'.csv', Nodes_status, delimiter=" ") 
        
        print("Number of axons activated: ",N_axons_activ)
        #print("VTA: ",VTA)
    
    return True
