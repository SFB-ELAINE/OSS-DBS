'''Written by K.Butenko'''
'''The script refines provided mesh in the regions of specified material (CSF)'''


from dolfin import *
from pandas import read_csv
from tissue_dielectrics import DielectricProperties
import numpy as np
import os.path
import pickle

#from Math_module_MPI_only_FEniCS import get_solutions
import subprocess

#from Math_module_new_floating import get_field,get_field_dual,get_field_QS,get_field_on_points
import time as tim

#launch_CSF_refinement is the manager function (called in Launcher)

#mesh,subdomains,boundaries_sol,Field_calc_param,EQS_mode=(0,0,0,0,0)


parameters['linear_algebra_backend']='PETSc'
parameters["refinement_algorithm"] = "plaza_with_parent_facets"
parameters["allow_extrapolation"] = True;





class Field_calc_parameters:
    def __init__(self,default_material,element_order,anisotropy,c_c,CPE,refinement_frequency,EQS_type):
        self.default_material=default_material
        self.element_order=element_order
        self.anisotropy=anisotropy
        self.c_c=c_c
        self.CPE=CPE
        self.frequenc=refinement_frequency      #list
        self.EQS_mode=EQS_type
        
        

#=============================================================================#    
    
def save_mesh_and_kappa_to_h5(mesh_to_h5,subdomains_to_h5,boundaries_to_h5,Field_calc_param):
    print("Number of mesh elements: ",mesh_to_h5.num_cells())

    #due to the glitch with the ghost model and subdomains. Used only for c-c multicontact, so always with floating
    V0_r=FunctionSpace(mesh_to_h5,'DG',0)    
    kappa_r=Function(V0_r)
    default_material=Field_calc_param.default_material
    
    [cond_GM, perm_GM]=DielectricProperties(3).get_dielectrics(Field_calc_param.frequenc)        #3 for grey matter and so on (numeration as in voxel_data)
    [cond_WM, perm_WM]=DielectricProperties(2).get_dielectrics(Field_calc_param.frequenc)
    [cond_CSF, perm_CSF]=DielectricProperties(1).get_dielectrics(Field_calc_param.frequenc)    
    [cond_default,perm_default]=DielectricProperties(default_material).get_dielectrics(Field_calc_param.frequenc)
    
    from GUI_inp_dict import d as d_encap
    [cond_encap, perm_encap]=DielectricProperties(d_encap['encap_tissue_type']).get_dielectrics(Field_calc_param.frequenc) 
    cond_encap=cond_encap*d_encap['encap_scaling_cond']
    perm_encap=perm_encap*d_encap['encap_scaling_perm']
        
    k_val_r=[cond_default*0.001,cond_CSF*0.001,cond_WM*0.001,cond_GM*0.001,cond_encap*0.001,1000.0]
    help = np.asarray(subdomains_to_h5.array(), dtype=np.int32)
    kappa_r.vector()[:] = np.choose(help, k_val_r)
            
    hdf = HDF5File(mesh_to_h5.mpi_comm(), 'Results_adaptive/Mesh_to_solve.h5', 'w')
    hdf.write(mesh_to_h5, "/mesh")
    hdf.write(subdomains_to_h5, "/subdomains")
    hdf.write(boundaries_to_h5, "/boundaries")
    hdf.write(kappa_r, "/kappa_r")
    
    file=File('Results_adaptive/Last_subdomains_map.pvd')
    file<<subdomains_to_h5
    file=File('Results_adaptive/Last_conductivity_map.pvd')
    file<<kappa_r
    
    if Field_calc_param.EQS_mode == 'EQS':
        V0_i=FunctionSpace(mesh_to_h5,'DG',0)    
        kappa_i=Function(V0_i)
        omega_eps0=2*np.pi*130.0*8.854e-12           #2*pi*f*eps0
        k_val_i=[omega_eps0*perm_default*0.001,omega_eps0*perm_CSF*0.001,omega_eps0*perm_WM*0.001,omega_eps0*perm_GM*0.001,1*omega_eps0*perm_encap*0.001,1000000000*omega_eps0]    
        help = np.asarray(subdomains_to_h5.array(), dtype=np.int32)
        kappa_i.vector()[:] = np.choose(help, k_val_i)
        
        hdf.write(kappa_i, "/kappa_i")        
        file=File('Results_adaptive/Last_permittivity_map.pvd')
        file<<kappa_i

    if Field_calc_param.anisotropy == 1:

        #         # Define conductivity expression and matrix
        # c00 = MeshFunction("double", mesh_to_h5, "Tensors/c00.xml.gz")
        # c01 = MeshFunction("double", mesh_to_h5, "Tensors/c01.xml.gz")
        # c02 = MeshFunction("double", mesh_to_h5, "Tensors/c02.xml.gz")
        # c11 = MeshFunction("double", mesh_to_h5, "Tensors/c11.xml.gz")
        # c12 = MeshFunction("double", mesh_to_h5, "Tensors/c12.xml.gz")
        # c22 = MeshFunction("double", mesh_to_h5, "Tensors/c22.xml.gz")

        #we get unscaled tensors and scale them with conductivity here. This approach is needed only for MPI
        c00 = MeshFunction("double", mesh_to_h5, 3, 0.0)
        c01 = MeshFunction("double", mesh_to_h5, 3, 0.0)
        c02 = MeshFunction("double", mesh_to_h5, 3, 0.0)
        c11 = MeshFunction("double", mesh_to_h5, 3, 0.0)
        c12 = MeshFunction("double", mesh_to_h5, 3, 0.0)
        c22 = MeshFunction("double", mesh_to_h5 3, 0.0)

        hdf2 = HDF5File(mesh_to_h5.mpi_comm(), "Results_adaptive/Tensors_to_solve_num_el_"+str(mesh_to_h5.num_cells())+".h5", "r")
        hdf2.read(c00, "/c00")
        hdf2.read(c01, "/c01")
        hdf2.read(c02, "/c02")
        hdf2.read(c11, "/c11")
        hdf2.read(c12, "/c12")
        hdf2.read(c22, "/c22")
        hdf2.close()
        
        cell_Anis = MeshFunction('bool',mesh_to_h5,3)
        cell_Anis.set_all(False)   
        
        for cell in cells(mesh_to_h5):
                
            scale_cond=k_val_r[subdomains_to_h5[cell]]
        
            if c00[cell]!=1.0:
                cell_Anis[cell]=True
        
            c00[cell]=c00[cell]*scale_cond
            c01[cell]=c01[cell]*scale_cond
            c02[cell]=c02[cell]*scale_cond
            c11[cell]=c11[cell]*scale_cond
            c12[cell]=c12[cell]*scale_cond
            c22[cell]=c22[cell]*scale_cond

            
        file=File('Tensors/c00_mapped.pvd')
        file<<c00,mesh_to_h5       
        file=File('Tensors/c01_mapped.pvd')
        file<<c01,mesh_to_h5          
        file=File('Tensors/c02_mapped.pvd')
        file<<c02,mesh_to_h5
        file=File('Tensors/c11_mapped.pvd')
        file<<c11,mesh_to_h5
        file=File('Tensors/c12_mapped.pvd')
        file<<c12,mesh_to_h5
        file=File('Tensors/c22_mapped.pvd')
        file<<c22,mesh_to_h5
        
        file=File('Tensors/Anis_cells.pvd')
        file<<cell_Anis,mesh_to_h5
        
        hdf.write(c00, "/c00")
        hdf.write(c01, "/c01")
        hdf.write(c02, "/c02")
        hdf.write(c11, "/c11")
        hdf.write(c12, "/c12")
        hdf.write(c22, "/c22")

    hdf.close()
                
    return True
    
def Dummy_CSF():
    
    
    mesh = Mesh("Meshes/Mesh_unref.xml")
    boundaries = MeshFunction('size_t',mesh,'Meshes/Mesh_unref_facet_region.xml')
    subdomains_assigned=MeshFunction('size_t',mesh,"Meshes/Mesh_unref_physical_region.xml")
    
    
    mesh_file=File('Results_adaptive/mesh_adapt.xml.gz')
    boundaries_file = File('Results_adaptive/boundaries_adapt.xml')
    subdomains_assigned_file=File('Results_adaptive/subdomains_assigned_adapt.xml')
    
    mesh_file<<mesh
    boundaries_file<<boundaries
    subdomains_assigned_file<<subdomains_assigned    
    
    #Scaling=1.0
    return True

def mesh_refiner(mesh_old,boundaries,subdomains_assigned,cell_markers,Domains,cc_multicontact):
        
    parameters['linear_algebra_backend']='PETSc'
    parameters["refinement_algorithm"] = "plaza_with_parent_facets"
    parameters["allow_extrapolation"] = True;
    
    facets_old = MeshFunction('size_t',mesh_old,2)
    facets_old.set_all(0)
    
    facets_old.array()[boundaries.array() == Domains.Contacts[0]]=1
    
    if cc_multicontact==True and Domains.fi[0]!=0.0:            #because ground contact is always subtracted from the mesh   
        dsSSS=Measure("dS",domain=mesh_old,subdomain_data=facets_old)  
        An_surface_size_old=assemble(1.0*dsSSS(1))
    else:
        dss=Measure("ds",domain=mesh_old,subdomain_data=facets_old)  
        An_surface_size_old=assemble(1.0*dss(1))

    mesh_new = refine(mesh_old, cell_markers)
    subdomains_assigned_new=adapt(subdomains_assigned,mesh_new)
    boundaries_new = adapt(boundaries,mesh_new) # put function space

    facets = MeshFunction('size_t',mesh_new,2)
    facets.set_all(0)
    facets.array()[boundaries_new.array()==Domains.Contacts[0]]=1
    #facets.array()[boundaries_new.array()==Domains.Contacts[1]]=2

    if cc_multicontact==True and Domains.fi[0]!=0.0:
        dsS_new=Measure("dS",domain=mesh_new,subdomain_data=facets)
        An_surface_size_new=assemble(1.0*dsS_new(1))
    else:
        dss_new=Measure("ds",domain=mesh_new,subdomain_data=facets)
        An_surface_size_new=assemble(1.0*dss_new(1))

   
    if (An_surface_size_new-An_surface_size_old)/An_surface_size_new>0.02:
        print((An_surface_size_new-An_surface_size_old)/An_surface_size_new)
        print("Refinement broke the imposed B.C.!")
        exit()

    return (mesh_new,boundaries_new,subdomains_assigned_new)

def index_cell_marker(mesh, index_array ,MRI_param, Scaler):
    
    #cell_ref=MeshFunctionBool(mesh,3)
    cell_ref = MeshFunction('bool',mesh,3)
    #subdomains.set_all(0)
    cell_ref.set_all(False)
    cell_to_ref=0
    cell_processed=0
    c00 = MeshFunction("double", mesh, 3)       #to check, which cells will be refined
    #print("before going over the cells")
    for cell in cells(mesh):
        cell_processed=cell_processed+1
        #print ('Cell N: {}'.format(cell_processed))
        smallest_edge=min([MRI_param.x_vox_size,MRI_param.y_vox_size,MRI_param.z_vox_size])
        
        if np.any(np.isin(index_array,cell.index())) and cell.h()>Scaler*smallest_edge:
            #print "found the cell"
            cell_ref[cell] = True
            cell_to_ref=cell_to_ref+1
            c00[cell]=1.0
        else:
            c00[cell]=0.0
    #print("finished marking cells")

    #if cell_ref.where_equal(True):      #save file only if marked
    #    file=File('CSF_ref/cells_from_indices.pvd')
    #    file<<c00
    
    #print('Number of cells to refine {}'.format(cell_to_ref))
    return cell_ref


#=========================Main================================================#

def Refine_CSF(MRI_param,DTI_param,Scaling,Domains,Field_calc_param,rel_div,CSF_frac_div,CSF_ref_add,cc_multicontact,num_proc,EQS_mode,ref_freq,Best_scaling=0,scaling_old=0):

    #global mesh,subdomains,boundaries_sol
    
    #if cc_multicontact==True:
        #from Math_module_new_floating import compute_field_with_superposition,get_field_on_points,compute_field_with_superposition_QS    
    #else:
        #from Math_module_new import get_field,get_field_QS,get_field_on_points
    
    if cc_multicontact==True:
        from Math_module_floating_MPI import get_solutions,get_field_on_points
    else:
        from Math_module_MPI_only_FEniCS import get_solutions,get_field_on_points 

    
    scaling_old=float(scaling_old)
    Scaling=float(Scaling) 
    Best_scaling=float(Best_scaling)
    
    start_CSF_refinement=tim.clock()    
    
    #from dolfin import *
    if Field_calc_param.anisotropy==1:
        #from Tensor_and_tissue_marking import get_cellmap_tensors
        from Tissue_marking_new import get_cellmap_tensors
        #from Tensor_and_tissue_marking_sundomains import get_cellmap_tensors
        #print("Anisotropy is on")
    else:
        from Tissue_marking_new import get_cellmap
    
    '''load results from the most refined iteration'''
    if Best_scaling!=0:
        Phi_amp_on_neuron_old_get=read_csv('CSF_ref/Field_on_points'+str(Best_scaling)+'.csv', delimiter=' ', header=None)
        Phi_amp_on_neuron_old=Phi_amp_on_neuron_old_get.values
        
        mesh = Mesh('CSF_ref/mesh_adapt_CSF'+str(scaling_old)+'.xml.gz')
        boundaries = MeshFunction('size_t',mesh,'CSF_ref/boundaries_adapt_CSF'+str(scaling_old)+'.xml')
        subdomains_assigned=MeshFunction('size_t',mesh,'CSF_ref/subdomains_assigned_adapt_CSF'+str(scaling_old)+'.xml')  
    else:
        mesh = Mesh("Meshes/Mesh_unref.xml")
        boundaries = MeshFunction('size_t',mesh,'Meshes/Mesh_unref_facet_region.xml')
        subdomains_assigned=MeshFunction('size_t',mesh,"Meshes/Mesh_unref_physical_region.xml")    
            
    Vertices_neur_get=read_csv('Neuron_model_arrays/Vert_of_Neural_model_NEURON.csv', delimiter=' ', header=None)
    Vertices_neur=Vertices_neur_get.values 
    
    '''structure of voxel array'''
    '''rows: x,y,z,material index'''
    '''first voxel has to have coordinates: 0+voxelsize_x,0+voxelsize_y,0+voxelsize_z'''
    '''Indeed, such a structure of the array is not optimal, since x,y,z coordinates can be simply stored in coordinate vectors'''
    '''Any suggestions (incl. minimal working example) will be appreciated'''
    

    #mesh = Mesh("Meshes/Mesh_unref.xml")
    #boundaries = MeshFunction('size_t',mesh,'Meshes/Mesh_unref_facet_region.xml')
    #subdomains_assigned=MeshFunction('size_t',mesh,"Meshes/Mesh_unref_physical_region.xml")
    
    if Best_scaling==0:
        #if os.path.isfile('CSF_ref/Field_on_points'+str(Best_scaling)+'.csv'):
        #    Field_old_get=read_csv('CSF_ref/Field_on_points'+str(Best_scaling)+'.csv', delimiter=' ', header=None)
        #    Field_old=Field_old_get.values
        
       # else:
            
        mesh_file=File('CSF_ref/mesh_adapt_CSF'+str(Best_scaling)+'.xml.gz')
        boundaries_file = File('CSF_ref/boundaries_adapt_CSF'+str(Best_scaling)+'.xml')
        subdomains_assigned_file=File('CSF_ref/subdomains_assigned_adapt_CSF'+str(Best_scaling)+'.xml')

        mesh_file<<mesh
        boundaries_file<<boundaries
        subdomains_assigned_file<<subdomains_assigned              
        
        
        print("Field calculation on the initial mesh")
        if Field_calc_param.anisotropy==1:
            subdomains=get_cellmap_tensors(mesh,subdomains_assigned,Domains,MRI_param,DTI_param,Field_calc_param.default_material)
        else:     
            subdomains=get_cellmap(mesh,subdomains_assigned,Domains,MRI_param,Field_calc_param.default_material)
        
        print("CSF_Subdomains_unref file was created")
        file=File('CSF_ref/CSF_Subdomains_unref.pvd')
        file<<subdomains,mesh

        save_mesh_and_kappa_to_h5(mesh,subdomains,boundaries,Field_calc_param)

        if cc_multicontact==True:
            subprocess.call(["mpirun", "-np", "{}".format(num_proc), "python3", "Math_module_floating_MPI.py"])
        else:
            #run_mpi(8)
            subprocess.call(["mpirun", "-np", "{}".format(num_proc), "python3", "Math_module_MPI_only_FEniCS.py"])

        Phi_r_old,Phi_im_old,Field_r_old,Field_im_old,max_E_old,J_r_old,J_im_old,j_dens_real,j_dens_im=get_solutions(EQS_mode,Field_calc_param.frequenc,Field_calc_param.element_order)
            #raise SystemExit
            #if EQS_mode=='QS':  #QS
            #    Phi_r_old,Phi_im_old,Field_r_old,Field_im_old,max_E_old,J_r_old,J_im_old,j_dens_real,j_dens_im=get_field_QS(mesh,Domains,subdomains,boundaries,Field_calc_param)
            #if EQS_mode=='EQS': #EQS
            #    Phi_r_old,Phi_im_old,Field_r_old,Field_im_old,max_E_old,J_r_old,J_im_old,j_dens_real,j_dens_im=get_field(mesh,Domains,subdomains,boundaries,Field_calc_param)
            
        
           
        file=File('CSF_ref/Field_r_'+str(Best_scaling)+'.pvd')
        file<<Field_r_old
        
        #[Field_old,Av_field_old,Max_field_old]=get_field_on_points(mesh,Field_r_old,Field_im_old,Field_calc_param.c_c,J_r_old,J_im_old)
        Phi_amp_on_neuron_old = get_field_on_points(Phi_r_old,Phi_im_old,Field_calc_param.c_c,J_r_old,J_im_old)
        
        #with open('CSF_ref/Field_on_points'+str(Best_scaling)+'.csv','w') as f_handle:
        #    np.savetxt(f_handle,Field_old)
        np.savetxt('CSF_ref/Field_on_points'+str(Best_scaling)+'.csv', Phi_amp_on_neuron_old, delimiter=" ")
    
    
    #print mesh.num_cells()
    

       
    #=============Math part=======================================================#
    
        
    #mesh_file=File('CSF_ref/mesh_adapt_CSF'+str(Scaling)+'.xml.gz')
    #boundaries_file = File('CSF_ref/boundaries_adapt_CSF'+str(Scaling)+'.xml')
    #subdomains_assigned_file=File('CSF_ref/subdomains_assigned_adapt_CSF'+str(Scaling)+'.xml')
    
    #mesh_file<<mesh
    #boundaries_file<<boundaries
    #subdomains_assigned_file<<subdomains_assigned    
    
    
    if os.path.isfile('MRI_DTI_derived_data/'+MRI_param.name[:-4]+'_voxel_array_CSF_'+str(CSF_ref_add)+'.npy') or os.path.isfile('MRI_DTI_derived_data/'+MRI_param.name[:-7]+'_voxel_array_CSF_'+str(CSF_ref_add)+'.npy'):
    ##if os.path.isfile('MRI_DTI_files/'+'Full_MRI'+'_voxel_array_CSF_'+str(CSF_ref_add)+'.csv'):         #if array was already prepared
        #voxel_array_CSF_get=read_csv('MRI_DTI_files/'+'Full_MRI'+'_voxel_array_CSF_'+str(CSF_ref_add)+'.csv', delimiter=' ', header=None)
        #voxel_array_CSF_get=read_csv('MRI_DTI_derived_data/'+MRI_param.name[:-4]+'_voxel_array_CSF_'+str(CSF_ref_add)+'.csv', delimiter=' ', header=None)
        #voxel_array_CSF=voxel_array_CSF_get.values
        if MRI_param.name[-2:]=='gz':
            voxel_array_CSF=np.load('MRI_DTI_derived_data/'+MRI_param.name[:-7]+'_voxel_array_CSF_'+str(CSF_ref_add)+'.npy') 
        else: 
            voxel_array_CSF=np.load('MRI_DTI_derived_data/'+MRI_param.name[:-4]+'_voxel_array_CSF_'+str(CSF_ref_add)+'.npy') 
        
        print("voxel_array_CSF in ",str(CSF_ref_add)," mm vicinity is loaded")
    else:
        #print "No CSF voxel array("
        start_voxel_array_CSF=tim.clock()
    
        #Tissue_array_get=read_csv('MRI_DTI_derived_data/Tissue_array_MRI.csv', delimiter=' ', header=None)
        #Tissue_array=Tissue_array_get.values
        Tissue_array=np.load('MRI_DTI_derived_data/Tissue_array_MRI.npy')
        
        x_vect=np.genfromtxt('MRI_DTI_derived_data/x_vector_MRI_Box.csv', delimiter=' ')
        y_vect=np.genfromtxt('MRI_DTI_derived_data/y_vector_MRI_Box.csv', delimiter=' ')
        z_vect=np.genfromtxt('MRI_DTI_derived_data/z_vector_MRI_Box.csv', delimiter=' ')
        
        voxel_size_x=abs(round(x_vect[1]-x_vect[0],6))
        voxel_size_y=abs(round(y_vect[1]-y_vect[0],6))
        voxel_size_z=abs(round(z_vect[1]-z_vect[0],6))

        
        voxel_array_CSF=np.zeros((Tissue_array.shape[0],3),float)      #array to store all CSF voxels in the specified ROI
        
        bb = mesh.bounding_box_tree()    
        
        Mx_mri=MRI_param.M_x
        My_mri=MRI_param.M_y
        Mz_mri=MRI_param.M_z
        
        x_neuron_max=max(Vertices_neur[:,0])
        y_neuron_max=max(Vertices_neur[:,1])        
        z_neuron_max=max(Vertices_neur[:,2])
        
        x_neuron_min=min(Vertices_neur[:,0])
        y_neuron_min=min(Vertices_neur[:,1])
        z_neuron_min=min(Vertices_neur[:,2])
        
        for x_coord in x_vect:
            for y_coord in y_vect:
                for z_coord in z_vect:
                    
                    x_pos=x_coord-voxel_size_x/2.0
                    y_pos=y_coord-voxel_size_y/2.0
                    z_pos=z_coord-voxel_size_z/2.0
                    
                    if (x_pos<=x_neuron_max+CSF_ref_add and x_pos>=x_neuron_min-CSF_ref_add and y_pos<=y_neuron_max+CSF_ref_add and y_pos>=y_neuron_min-CSF_ref_add and z_pos<=z_neuron_max+CSF_ref_add and z_pos>=z_neuron_min-CSF_ref_add):
                        
                        xv_mri=int((x_coord)/voxel_size_x-0.000000001)                                  #defines number of steps to get to the voxels containing x[0] coordinate
                        yv_mri=(int((y_coord)/voxel_size_y-0.000000001))*Mx_mri                  #defines number of steps to get to the voxels containing x[0] and x[1] coordinates
                        zv_mri=(int((z_coord)/voxel_size_z-0.000000001))*Mx_mri*My_mri           #defines number of steps to get to the voxels containing x[0], x[1] and x[2] coordinates
                        glob_index=int(xv_mri+yv_mri+zv_mri)
                        #print("glob_index: ",glob_index)
                        
                        
                        pnt=Point(x_pos,y_pos,z_pos)
                
                        if Tissue_array[glob_index]==1 and bb.compute_first_entity_collision(pnt)<mesh.num_cells()*100:
                            voxel_array_CSF[glob_index,0]=x_pos
                            voxel_array_CSF[glob_index,1]=y_pos
                            voxel_array_CSF[glob_index,2]=z_pos
    

        voxel_array_CSF=voxel_array_CSF[~np.all(voxel_array_CSF==0.0,axis=1)]  #deletes all zero enteries

        if MRI_param.name[-2:]=='gz':
            np.save('MRI_DTI_derived_data/'+MRI_param.name[:-7]+'_voxel_array_CSF_'+str(CSF_ref_add), voxel_array_CSF)
            #np.savetxt('MRI_DTI_derived_data/'+MRI_param.name[:-7]+'_voxel_array_CSF_'+str(CSF_ref_add)+'.csv', voxel_array_CSF, delimiter=" ") 
        else:
            np.save('MRI_DTI_derived_data/'+MRI_param.name[:-4]+'_voxel_array_CSF_'+str(CSF_ref_add), voxel_array_CSF)
            #np.savetxt('MRI_DTI_derived_data/'+MRI_param.name[:-4]+'_voxel_array_CSF_'+str(CSF_ref_add)+'.csv', voxel_array_CSF, delimiter=" ")

           
        del Tissue_array
        
        print("----- voxel_array_CSF for ",str(CSF_ref_add)," mm vicinity was prepared in %s seconds -----" % (tim.clock() - start_voxel_array_CSF))             
            
    
    #check collision of voxels with mesh, delete uncollided
    '''Here we pre-refine mesh, where CSF occured'''
    csf_ref=0    
    
    print("refining CSF voxels with scaling ", int(Scaling))
    
    num_cell_old=mesh.num_cells()
    
    if os.path.isfile('CSF_ref/mesh_adapt_CSF'+str(Scaling)+'.xml.gz'):
        print("Mesh was loaded from the refinement on the previous frequency")
        mesh = Mesh('CSF_ref/mesh_adapt_CSF'+str(Scaling)+'.xml.gz')
        boundaries = MeshFunction('size_t',mesh,'CSF_ref/boundaries_adapt_CSF'+str(Scaling)+'.xml')
        subdomains_assigned=MeshFunction('size_t',mesh,'CSF_ref/subdomains_assigned_adapt_CSF'+str(Scaling)+'.xml')            
    else:
        while csf_ref==0:
            
            inx_pnt=0       #to check, how much voxels were processed
        
            cell_index_list=[]
            bb = mesh.bounding_box_tree()
            
           
            for i_csf in range(voxel_array_CSF.shape[0]):
                pnt=Point(voxel_array_CSF[i_csf,0],voxel_array_CSF[i_csf,1],voxel_array_CSF[i_csf,2])
                inx_pnt=inx_pnt+1
                #print(inx_pnt)
                cell_index_list.append(bb.compute_first_entity_collision(pnt))
        
            
            cell_index_array=np.asarray(cell_index_list)
            cell_ref=index_cell_marker(mesh, cell_index_array, MRI_param, Scaling) 
            
            if not(cell_ref.where_equal(True)):     #if any cell was marked for refinement, will return True
                csf_ref=1
                
                mesh_file=File('CSF_ref/mesh_adapt_CSF'+str(Scaling)+'.xml.gz')
                mesh_file<<mesh        
                
                boundaries_file = File('CSF_ref/boundaries_adapt_CSF'+str(Scaling)+'.xml')
                subdomains_assigned_file=File('CSF_ref/subdomains_assigned_adapt_CSF'+str(Scaling)+'.xml')
        
        
                boundaries_file<<boundaries
                subdomains_assigned_file<<subdomains_assigned
                
                print("Number of cells after CSF refinement iteration: ",mesh.num_cells()) 

                if num_cell_old==mesh.num_cells():
                    print("skipping scaling ",Scaling)
                    return 0                 
            else:
                
                mesh_file=File('CSF_ref/mesh_adapt_CSF'+str(Scaling)+'_old.xml.gz')
                mesh_file<<mesh        
                
                boundaries_file = File('CSF_ref/boundaries_adapt_CSF'+str(Scaling)+'_old.xml')
                subdomains_assigned_file=File('CSF_ref/subdomains_assigned_adapt_CSF'+str(Scaling)+'_old.xml')
        
        
                boundaries_file<<boundaries
                subdomains_assigned_file<<subdomains_assigned
                #print("here")
                mesh,boundaries,subdomains_assigned=mesh_refiner(mesh,boundaries,subdomains_assigned,cell_ref,Domains,cc_multicontact)
                #print("there")
                #print("Number of cells: ",mesh.num_cells()) 
                
                
                if mesh.num_cells()>10000000:
                    print("Mesh is too large, will have to check with bigger scaling")
                    csf_refined=-1
                    return csf_refined
            
    #do not need to be redone, actually
    if Field_calc_param.anisotropy==1:
        subdomains=get_cellmap_tensors(mesh,subdomains_assigned,Domains,MRI_param,DTI_param,Field_calc_param.default_material)
    else:     
        subdomains=get_cellmap(mesh,subdomains_assigned,Domains,MRI_param,Field_calc_param.default_material)
    
    print("CSF_Subdomains_refinement file with scaling ",int(Scaling)," was created")
    file=File('CSF_ref/CSF_Subdomains_refinement_'+str(int(Scaling))+'.pvd')
    file<<subdomains,mesh
        
    save_mesh_and_kappa_to_h5(mesh,subdomains,boundaries,Field_calc_param)
    
    if cc_multicontact==True:
        subprocess.call(["mpirun", "-np", "{}".format(num_proc), "python3", "Math_module_floating_MPI.py"])
    else:   
        subprocess.call(["mpirun", "-np", "{}".format(num_proc), "python3", "Math_module_MPI_only_FEniCS.py"])
    
    Phi_r,Phi_im,Field_r,Field_im,max_E,J_r,J_im,j_dens_real,j_dens_im=get_solutions(EQS_mode,Field_calc_param.frequenc,Field_calc_param.element_order)        
        
        #if EQS_mode=='QS':  #QS
        #    Phi_r,Phi_im,Field_r,Field_im,max_E,J_r,J_im,j_dens_real,j_dens_im=get_field_QS(mesh,Domains,subdomains,boundaries,Field_calc_param)
        #if EQS_mode=='EQS': #EQS
        #    Phi_r,Phi_im,Field_r,Field_im,max_E,J_r,J_im,j_dens_real,j_dens_im=get_field(mesh,Domains,subdomains,boundaries,Field_calc_param)
       
    if Scaling==1:
        file=File('CSF_ref/Field_r_'+str(Scaling)+'.pvd')
        file<<Field_r

        print("CSF_Subdomains full refinement was created")
        file=File('CSF_ref/CSF_Subdomains_full_ref.pvd')
        file<<subdomains,mesh        
        

    #[Field,Av_field,Max_field]=get_field_on_points(mesh,Field_r,Field_im,Field_calc_param.c_c,J_r,J_im)
    Phi_amp_on_neuron =get_field_on_points(Phi_r,Phi_im,Field_calc_param.c_c,J_r,J_im)

#    with open('CSF_ref/Field_on_points'+str(Scaling)+'.csv','w') as f_handle:
#        np.savetxt(f_handle,Field)
    np.savetxt('CSF_ref/Field_on_points'+str(Scaling)+'.csv', Phi_amp_on_neuron, delimiter=" ")
   
    csf_refined=1

    max_div=0.0
    #affected_points=np.zeros((Field_old.shape[0],4),float)
    #aff_point_inx=0
        
    if Field_calc_param.c_c == True:
        if EQS_mode=='EQS':
            #not the best approach, but both should be on the Dirichlet BCs
            max_phi_r=max(Phi_r.vector()[:])
            max_phi_im=max(Phi_im.vector()[:])
            
            min_phi_r=min(Phi_r.vector()[:])
            min_phi_im=min(Phi_im.vector()[:])            
            
            phi_error=abs((np.sqrt((max_phi_r-min_phi_r)**2+(max_phi_im-min_phi_im)**2))*CSF_frac_div)
        else:
            phi_error=abs((max(Phi_r.vector()[:])-min(Phi_r.vector()[:]))*CSF_frac_div)   #should be scaled 
    else:    
        Phi_vector=[x for x in Domains.fi if x is not None]
        phi_error=abs((max(Phi_vector)-min(Phi_vector))*CSF_frac_div)      #Absolute potential error defined as a 1% of the maximum potential difference, VC case

        
    for inx in range(Phi_amp_on_neuron_old.shape[0]):
        if Best_scaling==0:             #first iteration
            delim=abs(Phi_amp_on_neuron[inx,3])
        else:
            delim=abs(Phi_amp_on_neuron_old[inx,3])
            
        if max_div<abs(Phi_amp_on_neuron_old[inx,3]-Phi_amp_on_neuron[inx,3]):
            max_div=abs(Phi_amp_on_neuron_old[inx,3]-Phi_amp_on_neuron[inx,3])
            
        if max_div> phi_error:
            print("Deviation threshold: ",phi_error,"V")  
            print("Deviation at least: ", max_div, "V")
            print("At the point: ", Phi_amp_on_neuron_old[inx,0],Phi_amp_on_neuron_old[inx,1],Phi_amp_on_neuron_old[inx,2])
            print("Need further refinement of CSF")
            csf_refined=0
            break

    if csf_refined==1:
        print("Max. deviation: ", max_div, "V")
        print("Deviation threshold: ",phi_error,"V")            
        print("CSF is refined enough")
        if Best_scaling==0:
            "On the initial mesh"
    
    #for inx in xrange(Field_old.shape[0]):
#    if (Field_old[inx,3]-Field[inx,3])/Field_old[inx,3]>0.01:
#        print "Need further refinement of CSF"
        
    del voxel_array_CSF
    #print("---CSF refinement took %s seconds ---" % (time.clock() - start_CSF_refinement))        
       
    minutes=int((tim.clock() - start_CSF_refinement)/60)
    secnds=int(tim.clock() - start_CSF_refinement)-minutes*60
    print("----- CSF refinement iteration took ",minutes," min ",secnds," s -----")
    print(" ")          
       
    return csf_refined
    

def launch_CSF_refinement(d,MRI_param,DTI_param,Domains,anisotrop,cc_multicontact,ref_freqs):


    d["CSF_frac_div"]=2*d["Adaptive_frac_div"]
    if d["CSF_frac_div"]<0.01:
        d["CSF_frac_div"]=0.01
        
    el_order_for_CSF=1      #just for now
    if cc_multicontact==True and el_order_for_CSF==1:
        el_order_for_CSF=2
        "element_order is 1, increasing to 2 for multicontact current-controlled stimulation"
    
    #Field_calc_param=Field_calc_parameters(d["default_material"],el_order_for_CSF,anisotrop,d["current_control"],d["CPE_activ"],d["refinement_frequency"])
    
    print("----- Conducting evaluation of CSF refinement -----")    
    #from CSF_re finement_parallel import Refine_CSF
    #print("CSF module was loaded")
    
    Min_Scaling=d["Min_Scaling"]
    Scaling_results=[]
    
    #global EQS_mode
    #EQS_mode=d["EQS_core"]
    #print("EQS_mode",EQS_mode)
    
    for freq in ref_freqs:
        print("At frequency: ",freq)
               
        #global Field_calc_param
        Field_calc_param=Field_calc_parameters(d["default_material"],el_order_for_CSF,anisotrop,d["current_control"],d["CPE_activ"],freq,d["EQS_core"])
        
        with open('Results_adaptive/Field_calc_param.file', "wb") as f:
            pickle.dump(Field_calc_param, f, pickle.HIGHEST_PROTOCOL)
        
        csf_ref=-1
        '''CSF_ref is 1, when we further refinement of elements with CSF voxels does not significantly change the result'''
        while csf_ref==-1:                  #NOTE: CSF REFINEMENT IS ALWAYS ON THE FIRST ORDER ELEMENTS!
            csf_ref=Refine_CSF(MRI_param,DTI_param,Min_Scaling,Domains,Field_calc_param,d["rel_div_CSF"],d["CSF_frac_div"],d["CSF_ref_reg"],cc_multicontact,d['number_of_processors'],d["EQS_core"],ref_freqs)          #return -1, when solution was not obtained (mesh is too large)
            if csf_ref==-1:                     
                Min_Scaling=Min_Scaling*2           #if the amount of elements is too large, we will have to have bigger tetrs.
    
        #subprocess.call('python Paraview_CSFunref.py', shell=True)
        #subprocess.call('xdg-open "Images/CSF_unref.png"',shell=True)
      
        if csf_ref==1:
            '''resave unref as CSF with scaling 0'''
            Scaling=0
        else:
            iteration_csf=0
            '''I don't like this part'''
            scaling_factors=[32*Min_Scaling,16*Min_Scaling,8*Min_Scaling,4*Min_Scaling,2*Min_Scaling,1*Min_Scaling]
            #Starting_Scaling=float(raw_input('Enter scaling factor for tetrahedron edges, containing CSF voxels (converg. study is recommended)\n'))
            scaling_old=0
            while csf_ref==0:
                print('scaling factor is ',scaling_factors[iteration_csf])
                Scaling=float(scaling_factors[iteration_csf])
                csf_ref=Refine_CSF(MRI_param,DTI_param,Scaling,Domains,Field_calc_param,d["rel_div_CSF"],d["CSF_frac_div"],d["CSF_ref_reg"],cc_multicontact,d['number_of_processors'],d["EQS_core"],ref_freqs,Min_Scaling,scaling_old) 
                iteration_csf=iteration_csf+1           #IF TOO MUCH ITERATIONS, GIVE A WARNING!!!!
                scaling_old=Scaling
    
        Scaling_results.append(Scaling)

    
    Scaling_arr=np.array([min(Scaling_results)])
    np.savetxt('CSF_ref/Scaling_CSF.csv', Scaling_arr, delimiter=" ")
    Scaling=float(Scaling)
    
    return Scaling


