import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

import os
path = os.getcwd()
os.chdir(path)

import importlib.util
spec = importlib.util.spec_from_file_location("d",path+'/GUI_inp_dict.py')
foo = importlib.util.module_from_spec(spec)
spec.loader.exec_module(foo)
d = foo.d

def plot_vta_figs(path_target = path+'/Images/',
                  path_nodes = path+'/Neuron_model_arrays/All_neuron_models.csv',
                  path_pattern = path+'/Neuron_model_arrays/default_pattern.csv',
                  path_activation = path+'/Field_solutions/Activation/Network_status.npy',
                  activate_close_axons=False):
    
    

    args = {'path_target': path_target, 'path_nodes':path_nodes, 'path_pattern': path_pattern, 'path_activation': path_activation, 'activate_close_axons': activate_close_axons}

    if d['VTA_type']=='cylinder':
        plot_cylindircal_vta(**args)
    else:
        plot_grid_vta(**args)


def plot_cylindircal_vta(input_dict = '../GUI_inp_dict.py',
                         path_activation = '../Field_solutions/Activation/Network_status.npy',
                         path_target = '../Images/'):
    """
    Visualizes the VTA cross-section for each surface. 
    """

    # reading the input dict
    #import importlib.util
    #spec = importlib.util.spec_from_file_location("d",input_dict)
    #foo = importlib.util.module_from_spec(spec)
    #spec.loader.exec_module(foo)
    #d = foo.d
    
    Nx = d['Nx']
    Nz_p = d['Nz_p']
    Nz_n = d['Nz_n']
    dx = d['dx']
    dz = d['dz']
    xmin = d['xmin']
    z0 = d['z0']
    n_dirs = d['n_dirs'] 

    # constructing the grid
    xc = np.arange(xmin, xmin+Nx*dx, dx)
    xc = np.sort(np.append(-xc, xc))
    zc = np.arange(z0-Nz_n*dz, z0+(1+Nz_p)*dz, dz)
    x,z = np.meshgrid(xc, zc)
    
    # reading the activation status
    activations = np.load(path_activation)
    if activate_close_axons:
        activations= np.abs(activations)
        
    axon_per_plane = len(xc)*len(zc)
    for n in range(n_dirs):
        vta_plane = activations[n*axon_per_plane: (n+1)*axon_per_plane]
        vta_plane = vta_plane.reshape(len(xc), len(zc))
        
        plt.figure(figsize=(15,4))
        
        vmin = vmin = (activate_close_axons-1)
        ax = plt.pcolormesh(xc, zc, vta_plane.T, shading='nearest', cmap='coolwarm',
                            edgecolors='gray', linewidths=0.25, vmin=vmin, vmax=1)
        plt.colorbar(ax,fraction=0.03,label='activation state')
        ax = plt.gca()
        ax.set_aspect('equal')
        plt.title('VTA profile prependicular to the direction '+ str(n))
        plt.tight_layout()
        plt.savefig(path_target+'Cyl_VTA_dir'+str(n)+'.png', dpi=300)
        plt.close()


def plot_grid_vta(path_target = '../Images/',
                  path_input_dict = '../GUI_input_dict.py',
                  path_nodes = '../Neuron_model_arrays/All_neuron_models.csv',
                  path_pattern = '../Neuron_model_arrays/default_pattern.csv',
                  path_activation = '../Field_solutions/Activation/Network_status.npy',
                  activate_close_axons=False):
    """
    extracts axons directions and the centers of the axon arrays from the 
    `All_neuron_models.csv`. Makes png (and possibly) nifti files to show VTA
    orthogonal to each direction.
    
    If axon centers are seeded in a plane (which will be recognized from the
    input dictionary) then just png files will be generated. Otherwise, for
    each direction, a nifti file is also produced and saved in the same directory.
    
    arguments:
    ----------
        path_target:      loaction to save VTA cross-section images and niftis
        path_input_dict:  path of the input dictionary (.py)
        path_nodes:       path of the full neuron array model (.csv) file 
        path_pattern:     path of the neural model pattern (.csv) file
        path_activation:  path of activation status (.npy) files
        activate_close_axons: True if removed axons are assumed activated (bool)
        
    """
    import nibabel as nib
    from scipy.spatial.transform import Rotation as R

    #import importlib.util
    #spec = importlib.util.spec_from_file_location("d",path_input_dict)
    #foo = importlib.util.module_from_spec(spec)
    #spec.loader.exec_module(foo)
    #d = foo.d
    
    eps = 1e-4
    nodes = np.loadtxt(path_nodes, delimiter=' ')
    pattern = np.loadtxt(path_pattern, delimiter=' ')
    
    ## Preprocessing ## 
    # N is number of nodes in a direction (not num. segments!)
    Nx = d['x_steps']+1 
    Ny = d['y_steps']+1
    Nz = d['z_steps']+1

    # dx is the length of segments
    dx = d['x_step']
    dy = d['y_step']
    dz = d['z_step']

    # electrode location
    elec_loc = np.array([d['Implantation_coordinate_X'],
                         d['Implantation_coordinate_Y'],
                         d['Implantation_coordinate_Z']])

    # eurler angles
    alphas = d['alpha_array_glob']
    betas = d['beta_array_glob']
    gammas = d['gamma_array_glob']

    # some good numbers about neuron array
    n_dir = len(d["alpha_array_glob"]) # number of unique directions
    n_axons_per_dir = Nx*Ny*Nz
    n_axons = n_axons_per_dir*n_dir
    n_nodes_per_axon = len(pattern)#len(nodes)//n_axons # total node per axon (ranvier and non-ranvier) 
    
    # center axons (from positive octant)
    CEN = nodes[:,:3].mean(axis=0)
    nodes[:,:3] -= CEN

    ## Extracting centers and directions ## 
    
    # information about each direction will be stored in a dict
    directions_info = {}
    
    # information about each axons in a dataframe
    centers = []
    direction_ids = []
    for axon_id in range(n_axons):
        vec = np.array(nodes[(axon_id+1)*n_nodes_per_axon-1, :3]-\
                       nodes[axon_id*n_nodes_per_axon, :3]) # axon vector
        cen = vec/2. + nodes[axon_id*n_nodes_per_axon, :3]
        centers.append(cen) # the center location

        # save direction if new
        if axon_id%n_axons_per_dir==0:
            id_ = axon_id//n_axons_per_dir
            info = {'direction': np.round(np.array(vec/np.linalg.norm(vec)),5),
                    'alpha': alphas[id_],
                    'beta': betas[id_],
                    'gamma': gammas[id_]}
            directions_info[id_] = info
            direction_ids += [id_]*n_axons_per_dir

    centers = np.round(np.array(centers), 5)
    activations = np.load(path_activation)
    data = {'center': list(centers), 'active':activations, 'dir_id': direction_ids}
    df = pd.DataFrame(data = data)

    ## Making nii/png files ##
    # centeral node 
    xc = np.arange(-(Nx-1)*dx/2, Nx*dx/2, dx)
    yc = np.arange(-(Ny-1)*dy/2, Ny*dy/2, dy)
    zc = np.arange(-(Nz-1)*dz/2, Nz*dz/2, dz)
    centers = np.meshgrid(xc,yc,zc, indexing='ij')

#     Max = np.stack(df.center.values).max()
#     Min = np.stack(df.center.values).min()
#     M = max(Max, abs(Min))

    for dir_id in sorted(directions_info.keys()):

        VTA = np.zeros_like(centers[0]) #shape = (Nx,Ny,Nz)

        alpha = directions_info[dir_id]['alpha']
        beta = directions_info[dir_id]['beta']
        gamma = directions_info[dir_id]['gamma']

        rx = R.from_euler('x', alpha, degrees=True).as_dcm()#.as_matrix()
        ry = R.from_euler('y', beta, degrees=True).as_dcm()#.as_matrix()
        rz = R.from_euler('z', gamma, degrees=True).as_dcm()#.as_matrix()
        R_tot = np.dot(np.dot(rz,ry), rx)

        # rotate the array
        x,y,z = np.dot(R_tot, np.array(centers).transpose(1,2,0,3))

        # make a ref direction (+y)
#         ref = np.array([[0,0],[0,M*0.8],[0,0]])
#         ref = np.dot(R_tot,ref)

        # finding the activated cells
        
        active_cond = (df.dir_id ==dir_id) & (df.active!=0)
        actives = np.stack(df[active_cond].center)
        status = df[active_cond].active.reset_index(drop=True)
        if activate_close_axons:
            status = np.abs(status)
    
        for idx, row in enumerate(actives):
            x0,y0,z0 = row

            cond_x = abs(x-x0)<eps
            cond_y = abs(y-y0)<eps
            cond_z = abs(z-z0)<eps

            cond =(cond_x * cond_y * cond_z)
            VTA[cond] = status.values[idx]
    
        if Ny>1:
            aff_ = np.zeros((4,4))
            scl = np.diag([dx, dy, dz])
            aff_[:3,:3] = scl
            aff_[:3,3] = np.array([center.min() for center in centers])+ elec_loc 
            nii = nib.nifti1.Nifti1Image(VTA.astype('int8'), aff_,)
            nib.save(nii, path_target+'VTA4nii_dir'+str(dir_id)+'.nii')
        
        
        VTA = VTA[:,0,:] # the section perpendicular to the direction
        
        vmin = (activate_close_axons-1)
        pcm = plt.pcolormesh(xc, zc, VTA, shading='nearest', cmap='coolwarm', 
                            edgecolors='gray', linewidths=1, vmin=vmin, vmax=1,
                           in_layout=False)
        
        ax = plt.gca()
        ax.set_title('VTA profile prependicular to the direction '+ str(dir_id))
        ax.set_aspect('equal')
        
        divider = make_axes_locatable(plt.gca())
        cax = divider.append_axes("right", "5%", pad="3%")
        plt.colorbar(pcm, cax=cax)
        
        plt.tight_layout()
        plt.savefig(path_target+'Grid_VTA_dir'+str(dir_id)+'.png', dpi=300)       
        plt.close()


