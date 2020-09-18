#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 14:09:42 2020

@author: butenko
"""

import h5py
import numpy as np
from dolfin import *
import matplotlib.pyplot as plt
from pandas import read_csv

from dolfin import *
import mshr

#This script allows to use VTA arrays of points instead of axons in OSS-DBS

def resave_as_verts(array_name):    #in mm, in MRI space

    mesh = Mesh("Meshes/Mesh_unref.xml")
    [__,__,__,x_min,y_min,z_min,__,__,__,MRI_voxel_size_x,MRI_voxel_size_y,MRI_voxel_size_z]=np.genfromtxt('MRI_DTI_derived_data/MRI_misc.csv', delimiter=' ')
        
    if array_name[-3:]=='.h5':
        hf = h5py.File(array_name, 'r')
        lst=list(hf.keys())
        result_total=[]
        arrays_shapes=[]
        
        for i in lst:
            a=hf.get(i)
            a=np.array(a)
            
            for j in range(a.shape[0]):
                pnt=Point(a[j,0]-x_min,a[j,1]-y_min,a[j,2]-z_min)
                if not(mesh.bounding_box_tree().compute_first_entity_collision(pnt)<mesh.num_cells()*10): 
                    a[j,:]=-100000000.0
                    
            a=a[~np.all(a==-100000000.0,axis=1)] 
            arrays_shapes.append(a.shape[0])        #save to use later to recognize the array. Make sure you know the order!
            
            result_total.append(a)  
        
        Array_coord=np.concatenate(result_total)
        hf.close()
      
        # shift to the positive octant space
        Array_coord[:,0]=Array_coord[:,0]-x_min
        Array_coord[:,1]=Array_coord[:,1]-y_min
        Array_coord[:,2]=Array_coord[:,2]-z_min    
        
        EPN_cut=arrays_shapes[0]
        np.savetxt('Neuron_model_arrays/Vert_EPN.csv', Array_coord[:EPN_cut,:], delimiter=" ") 
    
    np.savetxt('Neuron_model_arrays/Vert_of_Neural_model_NEURON.csv', Array_coord, delimiter=" ")  

    return arrays_shapes

#here we need to add function that will check that these vertices do not intersect with CSF, encap and so on

# do not use parallelization here


# can I scale E-field amplitude and potential instead???
#d,FR_vector_signal,Xs_signal_norm,t_vector
def ifft_on_VTA_array(name_sol,d,FREQ_vector_signal,Xs_signal_normalized,t_vect,T,i_start_octv,arrays_shape):    #in mm, in MRI space

    num_segments=sum(arrays_shape)

    Max_signal_for_point=np.zeros(num_segments,float)

    hf = h5py.File(name_sol[:-4]+'.h5', 'r')
    solution_sort_octv = hf.get('dataset_1')
    solution_sort_octv = np.array(solution_sort_octv)
    hf.close() 

    Fr_corresp_ar = np.genfromtxt('Stim_Signal/Fr_corresp_array'+str(d["trunc_param"]*1.0)+'.csv', delimiter=' ')
    #Fr_octave_vect = np.genfromtxt('Stim_Signal/Fr_octave_vector_'+str(d["trunc_param"]*1.0)+'.csv', delimiter=' ')
    FR_vec_sign_octv = np.genfromtxt('Stim_Signal/FR_vector_signal_octaves'+str(d["trunc_param"]*1.0)+'.csv', delimiter=' ')
        
    Fr_corresp_ar=np.round(Fr_corresp_ar,6)
    #Fr_octave_vect=np.round(Fr_octave_vect,6)
    N_freq_octv=(FR_vec_sign_octv.shape[0])
        
    
    for i_point in range(num_segments):
        Xs_Tr=np.vectorize(complex)(solution_sort_octv[(i_point)*N_freq_octv:(i_point*N_freq_octv+N_freq_octv),3],solution_sort_octv[i_point*N_freq_octv:(i_point*N_freq_octv+N_freq_octv),4])         #real and im parts for the first point in VTA
        Xs_Tr_full_real=np.zeros(FREQ_vector_signal.shape[0],float)
        Xs_Tr_full_imag=np.zeros(FREQ_vector_signal.shape[0],float)
        stepper=0
        for i_inx in range(Xs_Tr.shape[0]):
            if i_inx>=i_start_octv:
                rslt=np.where(Fr_corresp_ar[:,0]==np.round(FR_vec_sign_octv[i_inx],6))
                step_octv=rslt[0].shape[0]   #size of the freq. pack in the octave
                
                Xs_Tr_full_real[stepper:stepper+step_octv]=(Xs_Tr[i_inx].real) 
                Xs_Tr_full_imag[stepper:stepper+step_octv]=(Xs_Tr[i_inx].imag) 
                stepper=stepper+step_octv
            else:
                Xs_Tr_full_real[stepper]=(Xs_Tr[i_inx].real) 
                Xs_Tr_full_imag[stepper]=(Xs_Tr[i_inx].imag) 
                stepper=stepper+1
                
        if i_point==0:
            np.savetxt('Field_solutions/Xs_Tr_full_real.csv', Xs_Tr_full_real, delimiter=" ")
            np.savetxt('Field_solutions/Xs_Tr_full_imag.csv', Xs_Tr_full_imag, delimiter=" ")
            
        Xs_Tr_full_complex=np.vectorize(complex)(Xs_Tr_full_real,Xs_Tr_full_imag)
    
        Xs_conv=Xs_signal_normalized*Xs_Tr_full_complex

        if np.mod(t_vect.shape[0], 2):  # if the FT vector is odd
            fv_conj = np.conjugate(Xs_conv[-1:0:-1])
        else:  # if the FT vector is even
            fv_conj = np.conjugate(Xs_conv[-2:0:-1])
        
        Y = np.concatenate((Xs_conv, fv_conj), axis=0)
        
        Signal_t_conv=np.fft.ifft(Y).real
        
        Max_signal_for_point[i_point]=abs(max(Signal_t_conv[:], key=abs))
        
        if i_point==1:
            plt.figure(11111234)
            plt.plot(t_vect,Signal_t_conv)
            plt.xlim(0.000,T*5)
            plt.grid(True)
            plt.xlabel('t, sec')
            plt.ylabel('Potential, V')
            plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
            plt.savefig('Images/Signal_convoluted_1st_point.png', format='png', dpi=500)
       
        #np.save('Points_in_time/Signal_t_conv'+str(i_point), Signal_t_conv.real)
        
    return(Max_signal_for_point)
        
def get_VTA_from_E_ampl(Max_signal_for_point,arrays_shape,VTA_threshold,VTA_res):

    
    VTA_Vertices = np.genfromtxt('Neuron_model_arrays/Complete_neuron_models.csv', delimiter=' ')
    
    VTA_affected=np.zeros((VTA_Vertices.shape[0],4),float)
    VTA_affected[:,:3]=VTA_Vertices
    
    VTA=0.0         # here the VTA is just a sum of VTA voxels with E-field>E_threshold
    Axons_activated=0
    #slice_num=0 #we have 8 slices in total
    
    for i in range(VTA_Vertices.shape[0]):
        if abs(Max_signal_for_point[i])>=VTA_threshold:
            VTA_affected[i,3]=1.0
            VTA+=VTA_res**3
            Axons_activated+=1

    np.savetxt('Field_solutions/VTA_affected.csv', VTA_affected, delimiter=" ")
                
    print("VTA from VTA voxels: ",VTA)
    #print("N axons activated: ",Axons_activated)
    
    return(VTA)
                
        
#goes over each plane, separates it in two sides and computes volume from rotation of discs
#it can also check how many fibers intersected the volume (we assume 1 VTA region here)
    
#this function is valid only for the VTA chapter!!! Use it manually
def get_vta_arrays_as_discs(study_number,seeding_point_MRI_coord,create_VTA_mesh=False):

    if create_VTA_mesh==True:
        radius_7th=float(input("Please, enter the estimated radius of the VTA in the horizontal plane (assess visually the VTA_planes_activation): "))
    
    Vert_get=read_csv('VTA_planes_activation.csv', delimiter=' ', header=None)    # get only physiologically correct neuron models 
    Vert=Vert_get.values
    Vert=np.round(Vert,6)

    #we have 7 planes, but for the 7th we manually fit a sphere and added the shapes in Salome
    N_planes=6
    
    if study_number==1:
        z_axons_one_side=13
    else:
        z_axons_one_side=17      
    
    VTA_res=0.5
    r_encap=0.935 # this is the electrode radius + 0.3 encap

    VTA_volumes=np.zeros(13,float)  #from both sides around the electrode (first the left sides, 13th is reserved for the 7th plane)

    if create_VTA_mesh==True:
        total_disc=0

    # This part of the code is written in haste andis far from the optimal, please do not adopt it
    # Basically, we go over central nodes (first along Z-axis, then X-axis) and check the extent of activation
    #print(study_number)
    
    if study_number==2:
        
        maximum_r=4.0
        #first left planes    
        for i in range(N_planes):           #first halves
            for z_shift in range(z_axons_one_side):                
                if z_shift<4:
                    x_axons_one_side=8
                elif z_shift>=4:
                    x_axons_one_side=7
                    
                for j in np.arange(x_axons_one_side):
                   ix_axon=z_shift+j*z_axons_one_side
                   if Vert[ix_axon+i*249,3]==1.0:
                       if z_shift>=4:
                           VTA_volumes[i]=VTA_volumes[i]+np.pi*(((4-j*0.5)-r_encap)**2)*0.5           # disc volume, starts from 4 mm distance
                       else:
                           VTA_volumes[i]=VTA_volumes[i]+np.pi*(((4-j*0.5))**2)*0.5           # cylinder below, starts from 4 mm distance

                       #print(j)

                       if create_VTA_mesh==True:                        
                           r_disc=maximum_r-j*0.5
                           
            
                           point_bottom=Point(seeding_point_MRI_coord[0],seeding_point_MRI_coord[1],seeding_point_MRI_coord[2]-maximum_r+z_shift*0.5-0.25)
                           point_top=Point(seeding_point_MRI_coord[0],seeding_point_MRI_coord[1],seeding_point_MRI_coord[2]-maximum_r+z_shift*0.5+0.25)
            
                           if total_disc==0:
                               total_disc=mshr.Cylinder(point_top,point_bottom,r_disc,r_disc)-mshr.Cylinder(point_top,point_bottom,r_encap,r_encap)                  
                           else:                            
                               total_disc=total_disc+mshr.Cylinder(point_top,point_bottom,r_disc,r_disc)-mshr.Cylinder(point_top,point_bottom,r_encap,r_encap)
    
                       break
                

            for z_shift in range(z_axons_one_side):
                if z_shift<4:
                    x_axons_one_side=8
                elif z_shift>=4:
                    x_axons_one_side=7
                    
                for j in range(16,16-x_axons_one_side,-1):
                   if j==9:                     
                       ix_axon=z_shift+j*z_axons_one_side-27 
                   else:
                       ix_axon=z_shift+(j-1)*z_axons_one_side+4-27
                       
                   if Vert[ix_axon+i*249,3]==1.0:
                       if z_shift>=4:
                           VTA_volumes[i+6]=VTA_volumes[i+6]+np.pi*((((j-8)*0.5)-r_encap)**2)*0.5           # disc volume, starts from 4 mm distance
                       else:
                           VTA_volumes[i+6]=VTA_volumes[i+6]+np.pi*((((j-8)*0.5))**2)*0.5           # cylinder below, starts from 4 mm distance

                       #print(j)

                       if create_VTA_mesh==True:                        
                           # r_disc=((j-maximum_r+1.0)*0.5)
                           r_disc=((j-maximum_r+1.0)*0.5-0.5)
            
                           point_bottom=Point(seeding_point_MRI_coord[0],seeding_point_MRI_coord[1],seeding_point_MRI_coord[2]-maximum_r+z_shift*0.5-0.25)
                           point_top=Point(seeding_point_MRI_coord[0],seeding_point_MRI_coord[1],seeding_point_MRI_coord[2]-maximum_r+z_shift*0.5+0.25)
            
                           if total_disc==0:
                               total_disc=mshr.Cylinder(point_top,point_bottom,r_disc,r_disc)-mshr.Cylinder(point_top,point_bottom,r_encap,r_encap)                  
                           else:                            
                               total_disc=total_disc+mshr.Cylinder(point_top,point_bottom,r_disc,r_disc)-mshr.Cylinder(point_top,point_bottom,r_encap,r_encap)

                       break
         
            print("_______________________")
      
    else: # for the first and the third studies, there are no neurons below the electrode

        if study_number==3:
            x_axons_one_side=7
            total_x=13
            total=7*2*z_axons_one_side
            maximum_r=4.0
        else:
            x_axons_one_side=5
            total_x=9
            total=5*2*z_axons_one_side
            maximum_r=3.0
             
        #first left planes    
        for i in range(N_planes):           #first halves
            for z_shift in range(z_axons_one_side):                                    
                for j in np.arange(x_axons_one_side):
                    ix_axon=z_shift+j*z_axons_one_side
                    if Vert[ix_axon+i*total,3]==1.0:
                        VTA_volumes[i]=VTA_volumes[i]+np.pi*(((maximum_r-j*0.5)-r_encap)**2)*0.5           # disc volume, starts from 4 mm distance

                        if create_VTA_mesh==True:                        
                            r_disc=maximum_r-j*0.5
                            
                            point_bottom=Point(seeding_point_MRI_coord[0],seeding_point_MRI_coord[1],seeding_point_MRI_coord[2]-maximum_r+z_shift*0.5-0.25)
                            point_top=Point(seeding_point_MRI_coord[0],seeding_point_MRI_coord[1],seeding_point_MRI_coord[2]-maximum_r+z_shift*0.5+0.25)
                            
                            if total_disc==0:
                                total_disc=mshr.Cylinder(point_top,point_bottom,r_disc,r_disc)-mshr.Cylinder(point_top,point_bottom,r_encap,r_encap)                  
                            else:                            
                                total_disc=total_disc+mshr.Cylinder(point_top,point_bottom,r_disc,r_disc)-mshr.Cylinder(point_top,point_bottom,r_encap,r_encap)
                        break

            for z_shift in range(z_axons_one_side):       
                for j in range(total_x,total_x-x_axons_one_side,-1):
                    ix_axon=z_shift+(j)*z_axons_one_side
                      
                    if Vert[ix_axon+i*total,3]==1.0:
                        if study_number==3:
                            VTA_volumes[i+6]=VTA_volumes[i+6]+np.pi*((((j-5)*0.5)-r_encap)**2)*0.5           # disc volume, starts from 4 mm distance
                        else:
                            VTA_volumes[i+6]=VTA_volumes[i+6]+np.pi*((((j-3)*0.5)-r_encap)**2)*0.5 

                        if create_VTA_mesh==True:                        
                            # r_disc=((j-maximum_r+1.0)*0.5)

                            if study_number==3:
                                r_disc=((j-maximum_r+1.0)*0.5-0.5) 
                            else:
                                r_disc=((j-maximum_r)*0.5) 
                            
                            point_bottom=Point(seeding_point_MRI_coord[0],seeding_point_MRI_coord[1],seeding_point_MRI_coord[2]-maximum_r+z_shift*0.5-0.25)
                            point_top=Point(seeding_point_MRI_coord[0],seeding_point_MRI_coord[1],seeding_point_MRI_coord[2]-maximum_r+z_shift*0.5+0.25)
                            
                            if total_disc==0:
                                total_disc=mshr.Cylinder(point_top,point_bottom,r_disc,r_disc)-mshr.Cylinder(point_top,point_bottom,r_encap,r_encap)                  
                            else:                            
                                total_disc=total_disc+mshr.Cylinder(point_top,point_bottom,r_disc,r_disc)-mshr.Cylinder(point_top,point_bottom,r_encap,r_encap)
                        break


    if create_VTA_mesh==True:
        z_top=seeding_point_MRI_coord[2]+radius_7th
        z_bottom=seeding_point_MRI_coord[2]-radius_7th
        
        point_center=Point(seeding_point_MRI_coord[0],seeding_point_MRI_coord[1],seeding_point_MRI_coord[2])
        point_top=Point(seeding_point_MRI_coord[0],seeding_point_MRI_coord[1],z_top)
        point_bottom=Point(seeding_point_MRI_coord[0],seeding_point_MRI_coord[1],z_bottom)   
        
        #note that this is not a proper approach for study2, where the electrode does not fully pierce the VTA array
        Sphere_VTA=mshr.Sphere(point_center,radius_7th)-mshr.Cylinder(point_top,point_bottom,r_encap,r_encap)        
        total_shapes=total_disc+Sphere_VTA
        
        mesh_VTA=mshr.generate_mesh(total_shapes,40)
        file=File('mesh_VTA.pvd')
        file<<mesh_VTA
        
        mesh_file=File('mesh_VTA.xml.gz')
        mesh_file<<mesh_VTA        
        
    VTA_averaged=np.zeros(6,float)
    for i in range(VTA_averaged.shape[0]):
        #print(VTA_volumes)
        VTA_averaged[i]=(VTA_volumes[i]+VTA_volumes[i+6])/2.0 
    print("Averaged VTA in the vertical planes: ",VTA_averaged)
    
    return True

def check_fiber_intersection(Axons_per_population,Vertices_per_axon):
    
    #these are compratments of the pathway axons
    
    Vertices_get=read_csv('Neuron_model_arrays/Vert_of_Neural_model_NEURON.csv', delimiter=' ', header=None)
    Vertices=Vertices_get.values
    
    Vertices_activ=np.zeros((Vertices.shape[0],4),float)
    Vertices_activ[:,:3]=Vertices
    
    mesh_VTA=Mesh('mesh_VTA.xml.gz')        #from the previous simulation with a VTA array
    #file=File('mesh_VTA_E054.pvd')
    #file<<mesh_VTA    

    VTA_fiber_intersect=np.zeros(len(Axons_per_population),int)
    
    inx_shift=0     #among populations
    for i_population in range(len(Axons_per_population)):  
        fibers_intersected=0
        n_segments_fib_diam_array=Vertices_per_axon[i_population]
        axons_in_population=Axons_per_population[i_population]

        for i_axon in range(axons_in_population):
            for i_segm in range(n_segments_fib_diam_array):
                segm_coords=Point(Vertices[inx_shift+i_axon*n_segments_fib_diam_array+i_segm,0],Vertices[inx_shift+i_axon*n_segments_fib_diam_array+i_segm,1],Vertices[inx_shift+i_axon*n_segments_fib_diam_array+i_segm,2])
                if mesh_VTA.bounding_box_tree().compute_first_entity_collision(segm_coords)<mesh_VTA.num_cells():
                    Vertices_activ[inx_shift+i_axon*n_segments_fib_diam_array:inx_shift+(i_axon+1)*n_segments_fib_diam_array,3]=1.0
                    fibers_intersected+=1
                    break
                
        inx_shift=inx_shift+(i_axon+1)*n_segments_fib_diam_array
        #print(inx_shift)
                
        VTA_fiber_intersect[i_population]=fibers_intersected    
    
    print("Number of intersected per population: ",VTA_fiber_intersect) 
    np.savetxt('Field_solutions/Vertices_activ.csv', Vertices_activ, delimiter=" ")

    return True            


