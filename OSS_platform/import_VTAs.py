#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'PVD Reader'

for i in range(13):

    mesh_VTA_plane_0pvd = PVDReader(FileName='/home/butenko/OSS-DBS/OSS_platform_het_anis/Mesh_VTA/mesh_VTA_plane_'+str(i)+'.pvd')
    
    # find source
    fiber_activ_plane_Butsoncsv = FindSource('Fiber_activ_plane_Butson.csv')
    
    # find source
    phi_real_unscaled_1300Hzpvd = FindSource('Phi_real_unscaled_130.0Hz.pvd')
    
    # find source
    clip1 = FindSource('Clip1')
    
    # find source
    clip2 = FindSource('Clip2')
    
    # find source
    parallel_Subdomainspvd = FindSource('parallel_Subdomains.pvd')
    
    # find source
    tableToPoints1 = FindSource('TableToPoints1')
    
    # find source
    isoVolume1 = FindSource('IsoVolume1')
    
    # get active view
    renderView1 = GetActiveViewOrCreate('RenderView')
    # uncomment following to set a specific view size
    # renderView1.ViewSize = [1105, 742]
    
    # show data in view
    mesh_VTA_plane_0pvdDisplay = Show(mesh_VTA_plane_0pvd, renderView1)
    # trace defaults for the display properties.
    mesh_VTA_plane_0pvdDisplay.Representation = 'Surface'
    mesh_VTA_plane_0pvdDisplay.ColorArrayName = [None, '']
    mesh_VTA_plane_0pvdDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
    mesh_VTA_plane_0pvdDisplay.SelectOrientationVectors = 'None'
    mesh_VTA_plane_0pvdDisplay.ScaleFactor = 0.6000000000000001
    mesh_VTA_plane_0pvdDisplay.SelectScaleArray = 'None'
    mesh_VTA_plane_0pvdDisplay.GlyphType = 'Arrow'
    mesh_VTA_plane_0pvdDisplay.GlyphTableIndexArray = 'None'
    mesh_VTA_plane_0pvdDisplay.DataAxesGrid = 'GridAxesRepresentation'
    mesh_VTA_plane_0pvdDisplay.PolarAxes = 'PolarAxesRepresentation'
    mesh_VTA_plane_0pvdDisplay.ScalarOpacityUnitDistance = 0.22874548789819044
    
    # update the view to ensure updated data information
    renderView1.Update()
    
    # change solid color
    mesh_VTA_plane_0pvdDisplay.DiffuseColor = [1.0, 0.4549019607843137, 0.09019607843137255]
    
    #### saving camera placements for all active views
    
    # # current camera placement for renderView1
    # renderView1.CameraPosition = [93.99391804484023, 100.1227290553883, 76.49362638859624]
    # renderView1.CameraFocalPoint = [107.6849793823067, 120.29995750964257, 66.12110529851502]
    # renderView1.CameraViewUp = [0.18228033817956385, 0.3488841763942984, 0.9192680293444432]
    # renderView1.CameraParallelScale = 67.55183340389618
    
    #### uncomment the following to render all views
    # RenderAllViews()
    # alternatively, if you want to write images, you can use SaveScreenshot(...).