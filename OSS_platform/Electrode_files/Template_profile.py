# -*- coding: utf-8 -*-

###
### This file is generated automatically by SALOME v8.3.0 with dump python functionality
###### Run with DPS_lead_position_V9.py 

import sys
import salome

salome.salome_init()
theStudy = salome.myStudy

import salome_notebook
notebook = salome_notebook.NoteBook(theStudy)
sys.path.insert( 0, r'/data/trieu/electrode')

###
### GEOM component
###
######################################################################################################
######################################################################################################
########################################### extra code 1 V11 13/1/19##################################
###### This file runs with DBS_lead_position_V10.py
import os
sys.path.insert( 0, r'{}'.format(os.getcwd()))
sys.path.append('/usr/local/lib/python2.7/dist-packages')
from pandas import read_csv
from PrintLog import PrintLog

##### DEFAULT LIST #####

#Lead2nd_Enable = True
#Xt = 0
#Yt = 5
#Zt = 0
#X_2nd = 0
#Y_2nd = 5
#Z_2nd = 0
#OZ_angle = 0
#Xm = 0
#Ym = 0
#Zm = 0
#encap_thickness = 0.1
#ROI_radial = 13
#Vertice_enable = False
#Brain_map = '/home/trieu/electrode_dir/brain_elipse.brep'
#if(Lead2nd_Enable):
#   Xt2 = 0
#   Yt2 = -5
#   Zt2 = 0
#   OX_angle2 = 0
#   OY_angle2 = 0
#   OZ_angle2 = 0

##### VARIABLE LIST #####



########## End of variable list#############

# this is needed to allow rotation by 90 degrees
if Z_2nd == Zt:
    Z_2nd_artif = Zt+1.0 # just to ensure the rotation is possible
else:
    Z_2nd_artif=Z_2nd


if (Vertice_enable):
   Vert_array_get=read_csv('Vert_for_Salome.csv', delimiter=' ', header=None)
   Vert_array=Vert_array_get.values
else:
   Vert_array =[0];
number_vertex = len(Vert_array)
Vert = []
VolumeObject1 = []
ContactObject1 = []
VolumeObject2 = []
ContactObject2 = []
PrintLog ("############################################### DBS340 #######################################################\n")
PrintLog (" DBS_lead's Geometry buid\n")
######################################### end of extra code 1 ########################################
######################################################################################################
import GEOM
from salome.geom import geomBuilder
import math
import SALOMEDS


geompy = geomBuilder.New(theStudy)

O = geompy.MakeVertex(0, 0, 0)
OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)
Sphere_1 = geompy.MakeSphereR(0.05)
Cylinder_1 = geompy.MakeCylinderRH(0.05, 0.2)
tip = geompy.MakeFuseList([Sphere_1, Cylinder_1], True, True)
listFreeFacesIDs = geompy.GetFreeFacesIDs(tip)
[Free_face_1_1,Free_face_1_2,Free_face_1_3] = geompy.SubShapes(tip, [3, 7, 12])
listFreeFacesIDs = geompy.GetFreeFacesIDs(tip)
[Free_face_1_1,Free_face_1_2,Free_face_1_3] = geompy.SubShapes(tip, [3, 7, 12])
Contact_1 = geompy.MakeFuseList([Free_face_1_2, Free_face_1_3], True, True)
geompy.TranslateDXDYDZ(tip, 0, 0, 0.25)
geompy.TranslateDXDYDZ(tip, 0, 0, -0.25)
Cylinder_3 = geompy.MakeCylinderRH(0.165, 0.25)
geompy.TranslateDXDYDZ(Cylinder_3, 0, 0, 0.7)
Cylinder_2 = geompy.MakeCylinderRH(0.07, 1)
geompy.TranslateDXDYDZ(Cylinder_2, 0, 0, 0.2)
Cut_1 = geompy.MakeCutList(Cylinder_3, [Cylinder_2], True)
listFreeFacesIDs = geompy.GetFreeFacesIDs(Cut_1)
[Free_face_1_4,Free_face_1_5,Free_face_1_6,Free_face_1_7] = geompy.SubShapes(Cut_1, [3, 10, 15, 20])
listFreeFacesIDs = geompy.GetFreeFacesIDs(Cut_1)
[Free_face_1_4,Free_face_1_5,Free_face_1_6,Free_face_1_7] = geompy.SubShapes(Cut_1, [3, 10, 15, 20])
Contact_2 = geompy.MakeFuseList([Free_face_1_4, Free_face_1_6], True, True)
Cylinder_4 = geompy.MakeCylinderRH(0.2055, 99)
geompy.TranslateDXDYDZ(Cylinder_4, 0, 0, 0.95)
Sphere_2 = geompy.MakeSphereR(0.05+encap_thickness)
Cylinder_5 = geompy.MakeCylinderRH(0.05+encap_thickness, 0.2)
Cylinder_6 = geompy.MakeCylinderRH(0.07+encap_thickness, 0.5+encap_thickness)
geompy.TranslateDXDYDZ(Cylinder_6, 0, 0, 0.2-encap_thickness)
Cylinder_7 = geompy.MakeCylinderRH(0.165+encap_thickness, 0.25+encap_thickness)
geompy.TranslateDXDYDZ(Cylinder_7, 0, 0, 0.7-encap_thickness)
Cylinder_8 = geompy.MakeCylinderRH(0.2055+encap_thickness, 99+encap_thickness)
geompy.TranslateDXDYDZ(Cylinder_8, 0, 0, 0.95-encap_thickness)
fuse_encap = geompy.MakeFuseList([Sphere_2, Cylinder_5, Cylinder_6, Cylinder_7, Cylinder_8], True, True)

body = geompy.MakeFuseList([Sphere_1, Cylinder_1, Cylinder_3, Cylinder_2, Cylinder_4], True, True)
DBS_lead=geompy.MakeFuseList([Sphere_1, Cylinder_1, Cylinder_3, Cylinder_2, Cylinder_4], True, True)

encap = geompy.MakeCutList(fuse_encap, [body], True)
geompy.TranslateDXDYDZ(Sphere_1, 0, 0, 0.05)
geompy.TranslateDXDYDZ(Cylinder_1, 0, 0, 0.05)
geompy.TranslateDXDYDZ(tip, 0, 0, 0.05)
geompy.TranslateDXDYDZ(Contact_1, 0, 0, 0.05)
geompy.TranslateDXDYDZ(Cylinder_3, 0, 0, 0.05)
geompy.TranslateDXDYDZ(Cylinder_2, 0, 0, 0.05)
geompy.TranslateDXDYDZ(Cut_1, 0, 0, 0.05)
geompy.TranslateDXDYDZ(Contact_2, 0, 0, 0.05)
geompy.TranslateDXDYDZ(Cylinder_4, 0, 0, 0.05)
geompy.TranslateDXDYDZ(Sphere_2, 0, 0, 0.05)
geompy.TranslateDXDYDZ(Cylinder_5, 0, 0, 0.05)
geompy.TranslateDXDYDZ(Cylinder_6, 0, 0, 0.05)
geompy.TranslateDXDYDZ(Cylinder_7, 0, 0, 0.05)
geompy.TranslateDXDYDZ(Cylinder_8, 0, 0, 0.05)
geompy.TranslateDXDYDZ(fuse_encap, 0, 0, 0.05)
geompy.TranslateDXDYDZ(body, 0, 0, 0.05)
geompy.TranslateDXDYDZ(encap, 0, 0, 0.05)
ROI_Sphere_3 = geompy.MakeSphereR(ROI_radial)

encap_outer_ROI = geompy.MakeCutList(encap, [ROI_Sphere_3], True)
encap_inner_ROI = geompy.MakeCutList(encap, [encap_outer_ROI], True)
ROI = geompy.MakeCutList(ROI_Sphere_3, [fuse_encap], True)
Fuse_all=geompy.MakeFuseList([ROI,encap_outer_ROI,encap_inner_ROI,DBS_lead], True, True)
######################################################################################################
########################################### extra code 2 V11 13/1/19#############################################
PrintLog (" Load brain image \n")
if (Brain_map[-4:] == 'brep'):
	brain_solid = geompy.ImportBREP( Brain_map )
elif (Brain_map[-4:] == 'step'):
	brain_solid = geompy.ImportSTEP( Brain_map )
elif (Brain_map[-4:] == 'iges'):
	brain_solid = geompy.ImportIGES( Brain_map )
elif (Brain_map[-4:] == '.stl'):
	brain_solid = geompy.ImportSTL( Brain_map )
else:
	PrintLog (" unknow imported file format")

#################################################### Geometry and extra code interface ###############################################################
Fuse_all_lead_encap_ROI = Fuse_all  
Fuse_all_lead_encap_ROI_no_internal_face = geompy.RemoveInternalFaces(Fuse_all_lead_encap_ROI)
VolumeObject1 = [ encap_outer_ROI,ROI,encap_inner_ROI]         # Declare objects included to partition, encap_outer_ROI always @1st position
Volume_name1  = ['encap_outer_ROI1','ROI1','encap_inner_ROI1'] # Declare name of the group in the partition for volume
ContactObject1 = [Contact_1,Contact_2]   
Contact_name1 = ['Contact1_1','Contact1_2']

if(Lead2nd_Enable): ##################  2nd LEAD ###############################################
  VolumeObject2 = [ROI]*len(VolumeObject1)
  ContactObject2 = [Contact_1]*len(ContactObject1)
  Volume_name2  = [ 'encap_outer_ROI2','ROI2','encap_inner_ROI2']
  Contact_name2 = ['Contact2_1','Contact2_2']
###################################################################################################################################################
  PrintLog ("Position 2nd Fuse all object at [{},{},{}], [{}',{}',{}']\n".format(Xt2,Yt2,Zt2,OX_angle2,OY_angle2,OZ_angle2))
  Fuse_all_lead_encap_ROI_no_internal_face2 = geompy.MakeTranslation(Fuse_all_lead_encap_ROI_no_internal_face,Xt2,Yt2,Zt2)
  OX2 = geompy.MakeTranslation(OX,Xt2,Yt2,Zt2)
  OY2 = geompy.MakeTranslation(OY,Xt2,Yt2,Zt2)
  OZ2 = geompy.MakeTranslation(OZ,Xt2,Yt2,Zt2)
  geompy.Rotate(Fuse_all_lead_encap_ROI_no_internal_face2, OX2,OX_angle2*math.pi/180.0)
  geompy.Rotate(Fuse_all_lead_encap_ROI_no_internal_face2, OY2,OY_angle2*math.pi/180.0)
  geompy.Rotate(Fuse_all_lead_encap_ROI_no_internal_face2, OZ2,OZ_angle2*math.pi/180.0)
  PrintLog (" cut fuse all scaled up with neuron \n")
  for i in range(0,len(VolumeObject1)):
	VolumeObject2[i] = geompy.MakeTranslation(VolumeObject1[i],Xt2,Yt2,Zt2)
	geompy.Rotate(VolumeObject2[i], OX2,OX_angle2*math.pi/180.0)
	geompy.Rotate(VolumeObject2[i], OY2,OY_angle2*math.pi/180.0)
	geompy.Rotate(VolumeObject2[i], OZ2,OZ_angle2*math.pi/180.0)
  for i in range(0,len(ContactObject1)):
	ContactObject2[i] = geompy.MakeTranslation(ContactObject1[i],Xt2,Yt2,Zt2)
	geompy.Rotate(ContactObject2[i], OX2,OX_angle2*math.pi/180.0)
	geompy.Rotate(ContactObject2[i], OY2,OY_angle2*math.pi/180.0)
	geompy.Rotate(ContactObject2[i], OZ2,OZ_angle2*math.pi/180.0)
  PrintLog ("Cut outer ROI2 with brain\n")
  cut_outer_ROI = geompy.MakeCutList(VolumeObject2[0], [brain_solid], True)
  VolumeObject2[0] = geompy.MakeCutList(VolumeObject2[0], [cut_outer_ROI], True)
  PrintLog ("Cut ROI2 with brain\n")
  VolumeObject2[1] = geompy.MakeCommonList([VolumeObject2[1], brain_solid], True) 
  PrintLog ("Group 2nd:volume and area extraction for group ID identification process\n")        
  Volume2_Pro = [geompy.BasicProperties( VolumeObject2[0])]*len(VolumeObject2)
  Contact2_Pro = [geompy.BasicProperties( ContactObject2[0])]*len(ContactObject2)
  for i in range(0,len(VolumeObject2)):
	Volume2_Pro[i] = geompy.BasicProperties( VolumeObject2[i])
  for i in range(0,len(ContactObject2)):
	Contact2_Pro[i] = geompy.BasicProperties( ContactObject2[i])

################## LEAD 1st #############################################################
#print "Position 1st Fuse all object at [{},{},{}], [{}',{}',{}']\n".format(Xt,Yt,Zt,OX_angle,OY_angle,OZ_angle)
geompy.TranslateDXDYDZ(Fuse_all_lead_encap_ROI_no_internal_face,Xt,Yt,Zt)

OX1 = geompy.MakeTranslation(OX,Xt,Yt,Zt)
OY1 = geompy.MakeTranslation(OY,Xt,Yt,Zt)
OZ1 = geompy.MakeTranslation(OZ,Xt,Yt,Zt)

geompy.Rotate(Fuse_all_lead_encap_ROI_no_internal_face, OZ1,OZ_angle*math.pi/180.0)

Vertex_1 = geompy.MakeVertex(X_2nd,Y_2nd,Z_2nd)
Vertex_O = geompy.MakeVertex(Xt,Yt,Zt)
Vertex_3 = geompy.MakeVertex(Xt,Yt,Z_2nd_artif)

if X_2nd!=Xt or Y_2nd!=Yt:
        Fuse_all_lead_encap_ROI_no_internal_face=geompy.MakeRotationThreePoints(Fuse_all_lead_encap_ROI_no_internal_face, Vertex_O, Vertex_3, Vertex_1)


#print "Position 1st Lead at [{},{},{}], [{}',{}',{}']\n".format(Xt,Yt,Zt,OX_angle,OY_angle,OZ_angle)
for i in range(0,len(VolumeObject1)):
    geompy.TranslateDXDYDZ(VolumeObject1[i],Xt,Yt,Zt)
    geompy.Rotate(VolumeObject1[i], OZ1,OZ_angle*math.pi/180.0)
    if X_2nd!=Xt or Y_2nd!=Yt:
        VolumeObject1[i]=geompy.MakeRotationThreePoints(VolumeObject1[i], Vertex_O, Vertex_3, Vertex_1)

for i in range(0,len(ContactObject1)):
    geompy.TranslateDXDYDZ(ContactObject1[i],Xt,Yt,Zt)
    geompy.Rotate(ContactObject1[i], OZ1,OZ_angle*math.pi/180.0)
    if X_2nd!=Xt or Y_2nd!=Yt:
        ContactObject1[i]=geompy.MakeRotationThreePoints(ContactObject1[i], Vertex_O, Vertex_3, Vertex_1)



PrintLog( "Cut outer ROI1 with brain\n")
cut_outer_ROI = geompy.MakeCutList(VolumeObject1[0], [brain_solid], True)
VolumeObject1[0] = geompy.MakeCutList(VolumeObject1[0], [cut_outer_ROI], True)
PrintLog( "Cut ROI1 with brain\n")
VolumeObject1[1] = geompy.MakeCommonList([VolumeObject1[1], brain_solid], True)
PrintLog ("Group 1st:volume and area extraction for group ID identification process\n") 

Volume1_Pro = [geompy.BasicProperties( VolumeObject1[0])]*len(VolumeObject1)
Contact1_Pro = [geompy.BasicProperties( ContactObject1[0])]*len(ContactObject1)
for i in range(0,len(VolumeObject1)):
	Volume1_Pro[i] = geompy.BasicProperties( VolumeObject1[i])
for i in range(0,len(ContactObject1)):
	Contact1_Pro[i] = geompy.BasicProperties( ContactObject1[i])


PrintLog ("Create reference groups for ID identification process\n")

if(Lead2nd_Enable):
  
  Rest = geompy.MakeCutList(brain_solid, [Fuse_all_lead_encap_ROI_no_internal_face,Fuse_all_lead_encap_ROI_no_internal_face2], True)
  Partition_profile = geompy.MakePartition(VolumeObject1+VolumeObject2+ContactObject1+ContactObject2, [], [], [], geompy.ShapeType["SOLID"], 0, [], 0)
  
###reference_volume 
  reference_volume = VolumeObject1 + VolumeObject2
  reference_volume_Pro = Volume1_Pro + Volume2_Pro
  Volume_name = Volume_name1+Volume_name2
### reference_area 
  reference_surface = ContactObject1 + ContactObject2
  reference_surface_Pro = Contact1_Pro + Contact2_Pro
  Contact_name = Contact_name1+Contact_name2
  Group_volume  = [geompy.CreateGroup(Partition_profile, geompy.ShapeType["SOLID"])] * (len(VolumeObject1)+len(VolumeObject2)+1) # +1 is Rest Group
  Group_surface = [geompy.CreateGroup(Partition_profile, geompy.ShapeType["FACE"])] * (len(ContactObject1)+len(ContactObject2))
else:
  
  Rest = geompy.MakeCutList(brain_solid, [Fuse_all_lead_encap_ROI_no_internal_face], True)
  Partition_profile = geompy.MakePartition(VolumeObject1+ContactObject1, [], [], [], geompy.ShapeType["SOLID"], 0, [], 0)
###reference_volume 
  reference_volume = VolumeObject1 
  reference_volume_Pro = Volume1_Pro 
  Volume_name = Volume_name1
### reference_area 
  reference_surface = ContactObject1 
  reference_surface_Pro = Contact1_Pro 
  Contact_name = Contact_name1
  Group_volume  = [geompy.CreateGroup(Partition_profile, geompy.ShapeType["SOLID"])] * (len(VolumeObject1)+1) # +1 is Rest Group
  Group_surface = [geompy.CreateGroup(Partition_profile, geompy.ShapeType["FACE"])] * len(ContactObject1)
### find out subshape and subshape ID 

Group_surface_ListIDs =[]
Group_volume_ListIDs =[]
Group_partition_volume = []
Group_partition_surface = []

### find group volume ID ######################################################################
Partition_volume_IDsList = geompy.SubShapeAllIDs(Partition_profile, geompy.ShapeType["SOLID"]) # list all sub shape volume in Partition
PrintLog ("Partition_volume_IDsList: {} \n".format(Partition_volume_IDsList))

for ref_ind in range (0, len(reference_volume)):
	temp_volume = []
	for sub_ind in range (0, len (Partition_volume_IDsList)):
		subshape = geompy.GetSubShape(Partition_profile, [Partition_volume_IDsList[sub_ind]]) # get subshape
		subshape_Pro = geompy.BasicProperties(subshape)       # extract volume of subshape
		Common_volume = geompy.MakeCommonList([subshape, reference_volume[ref_ind]], True) # check common intersection
		Common_volume_Pro = geompy.BasicProperties(Common_volume) 
		PrintLog ("volume difference {}/{} \n".format(abs(Common_volume_Pro[2]-subshape_Pro[2]),abs(Common_volume_Pro[2]-reference_volume_Pro[ref_ind][2])))
		# if ( common volume = subshape) and (common volume = ref volume) => ref volume = sub shape
		if (abs(Common_volume_Pro[2]-subshape_Pro[2])< 0.00001) and (abs(Common_volume_Pro[2]-reference_volume_Pro[ref_ind][2])<0.00001):
		
			Group_partition_volume.append([Volume_name[ref_ind],Partition_volume_IDsList[sub_ind]])
		# if ( common volume = subshape) and (common volume < ref volume) => sub shape belong to ref volume
		elif (abs(Common_volume_Pro[2]-subshape_Pro[2])< 0.00001) and ((Common_volume_Pro[2] - reference_volume_Pro[ref_ind][2])<-0.00001):
			temp_volume.append( Partition_volume_IDsList[sub_ind] )
	if len(temp_volume) >1 : # the volume is devided
		Group_partition_volume.append([Volume_name[ref_ind],temp_volume ])
		PrintLog (str(Volume_name[ref_ind])+" is devided and has sub IDs:{}\n".format(temp_volume))
if len(reference_volume) != len(Group_partition_volume):
	PrintLog ("Geometry-volume error please check ROI diameter and DBS lead Position {} {}\n".format(len(reference_volume),len(Group_partition_volume)))
PrintLog ('Group_partition_volume: {}\n'.format(Group_partition_volume))

### find group surface ID ######################################################################
Partition_surface_IDsList = geompy.SubShapeAllIDs(Partition_profile, geompy.ShapeType["FACE"]) # list all sub shape face in Partition
PrintLog ('Partition_surface_IDsList: {}\n'.format(Partition_surface_IDsList))
sub_face = [] ## store devided faces
for reff_ind in range (0, len (reference_surface)):
	temp_surface = []
	for subf_ind in range (0, len(Partition_surface_IDsList)):
		subshapef = geompy.GetSubShape(Partition_profile, [Partition_surface_IDsList[subf_ind]]) # get subshape
		Common_face = geompy.MakeCommonList([subshapef, reference_surface[reff_ind]], True) # check common intersection
		Common_face_Pro = geompy.BasicProperties(Common_face) 
		subshapef_Pro = geompy.BasicProperties(subshapef) # extract volume of subshape
		PrintLog ("area difference: {}/{} \n".format(abs(Common_face_Pro[1]-subshapef_Pro[1]),abs(Common_face_Pro[1]-reference_surface_Pro[reff_ind][1])))
		# if ( common face = subface) and (common face = ref face) => ref face = sub face
		if (abs(Common_face_Pro[1]-subshapef_Pro[1])<0.000001 )and (abs(Common_face_Pro[1]-reference_surface_Pro[reff_ind][1])<0.000001): 
			Group_partition_surface.append([ Contact_name[reff_ind],Partition_surface_IDsList[subf_ind] ])
		# if ( common face = subface) and (common face < ref face) => sub face belong to ref face 
		elif (abs(Common_face_Pro[1]-subshapef_Pro[1])<0.000001 ) and ((Common_face_Pro[1] - reference_surface_Pro[reff_ind][1])<-0.000001):
			temp_surface.append(Partition_surface_IDsList[subf_ind])
	if len(temp_surface) >1 : # the face is devided
		Group_partition_surface.append( [Contact_name[reff_ind],temp_surface ])		
		PrintLog (Contact_name[reff_ind] + " is devided and has sub IDs:{}\n".format(temp_surface))
if len(reference_surface) != len(Group_partition_surface): #+len(Group_partition_Multi_surface):
	PrintLog ("Geometry-Surface error please check ROI diameter and DBS lead Position {},{}\n".format(len(reference_surface),len(Group_partition_surface)))
PrintLog ('Group_partition_surface {}\n'.format(Group_partition_surface,'\n'))

if(Lead2nd_Enable):
   Partition_profile = geompy.MakePartition(VolumeObject1+VolumeObject2+ContactObject1+ContactObject2+[Rest], [], [], [], geompy.ShapeType["SOLID"], 0, [], 0)
else:
   Partition_profile = geompy.MakePartition(VolumeObject1+ContactObject1+[Rest], [], [], [], geompy.ShapeType["SOLID"], 0, [], 0)

new_volume_ID= geompy.SubShapeAllIDs(Partition_profile, geompy.ShapeType["SOLID"])
ID= list(set(Partition_volume_IDsList) ^ set (new_volume_ID))
Group_partition_volume.append(['Rest_1',ID[0]])
PrintLog ("REST {}\n".format(ID))
PrintLog ('Group_partition_volume: {}\n'.format(Group_partition_volume))
PrintLog ("Create volume and surface group under partition_profile\n")

for i_solid in range (0,len (Group_partition_volume)):
	Group_volume[i_solid] = geompy.CreateGroup(Partition_profile, geompy.ShapeType["SOLID"])
	if (isinstance (Group_partition_volume[i_solid][1],list) == False):
		geompy.UnionIDs(Group_volume[i_solid], [Group_partition_volume[i_solid][1]])
	if (isinstance (Group_partition_volume[i_solid][1],list) == True):
		geompy.UnionIDs(Group_volume[i_solid], Group_partition_volume[i_solid][1])

#############################################

for i_surface in range (0,len (Group_partition_surface)):
	Group_surface[i_surface] = geompy.CreateGroup(Partition_profile, geompy.ShapeType["FACE"])
	if (isinstance (Group_partition_surface[i_surface][1],list) == False): # not a list
		geompy.UnionIDs(Group_surface[i_surface], [Group_partition_surface[i_surface][1]])
	if (isinstance (Group_partition_surface[i_surface][1],list) == True): #  it is a list
		geompy.UnionIDs(Group_surface[i_surface], Group_partition_surface[i_surface][1])
PrintLog ("Translate whole partition to Xm,Ym,Zm\n")
geompy.TranslateDXDYDZ(Partition_profile, Xm, Ym, Zm)
### add Vertices to geometry
if(Vertice_enable):
   for ver_ind in range (0,number_vertex):
       PrintLog("Add vertices to model\n")
       Vert.append(geompy.MakeVertex(Vert_array[ver_ind][0],Vert_array[ver_ind][1],Vert_array[ver_ind][2]))
       geompy.TranslateDXDYDZ(Vert[ver_ind], Xm, Ym, Zm)       ###Translate vertices to Xm,Ym,Zm
       geompy.addToStudy( Vert[ver_ind], 'Vert_{}'.format(ver_ind))
PrintLog ("add to study\n")
############################################ end of extra code 2 ############################################
#############################################################################################################
geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )
#geompy.addToStudy( Sphere_1, 'Sphere_1' )
#geompy.addToStudy( Cylinder_1, 'Cylinder_1' )
#geompy.addToStudy( tip, 'tip' )
#geompy.addToStudyInFather( tip, Free_face_1_1, 'Free_face_1_1' )
#geompy.addToStudyInFather( tip, Free_face_1_2, 'Free_face_1_2' )
#geompy.addToStudyInFather( tip, Free_face_1_3, 'Free_face_1_3' )
geompy.addToStudy( Contact_1, 'Contact_1' )
#geompy.addToStudy( Cylinder_3, 'Cylinder_3' )
#geompy.addToStudy( Cylinder_2, 'Cylinder_2' )
#geompy.addToStudy( Cut_1, 'Cut_1' )
#geompy.addToStudyInFather( Cut_1, Free_face_1_4, 'Free_face_1_4' )
#geompy.addToStudyInFather( Cut_1, Free_face_1_5, 'Free_face_1_5' )
#geompy.addToStudyInFather( Cut_1, Free_face_1_6, 'Free_face_1_6' )
#geompy.addToStudyInFather( Cut_1, Free_face_1_7, 'Free_face_1_7' )
geompy.addToStudy( Contact_2, 'Contact_2' )
#geompy.addToStudy( Cylinder_4, 'Cylinder_4' )
#geompy.addToStudy( body, 'body' )
#geompy.addToStudy( Sphere_2, 'Sphere_2' )
#geompy.addToStudy( Cylinder_5, 'Cylinder_5' )
#geompy.addToStudy( Cylinder_6, 'Cylinder_6' )
#geompy.addToStudy( Cylinder_7, 'Cylinder_7' )
#geompy.addToStudy( Cylinder_8, 'Cylinder_8' )
#geompy.addToStudy( fuse_encap, 'fuse_encap' )
#geompy.addToStudy( encap, 'encap' )
#geompy.addToStudy( ROI_Sphere_3, 'ROI_Sphere_3' )
geompy.addToStudy( encap_outer_ROI, 'encap_outer_ROI' )
geompy.addToStudy( encap_inner_ROI, 'encap_inner_ROI ' )
geompy.addToStudy( ROI, 'ROI' )
geompy.addToStudy( Fuse_all_lead_encap_ROI, 'Fuse_all_lead_encap_ROI' )
################################################################################################################
####################################### extra code 3 V10  15/12/18##############################################/
#for i in range(0,len(VolumeObject2)):/
#	geompy.addToStudy( VolumeObject2[i], 'VolumeObject2_{}'.format(i) )
#for i in range(0,len(ContactObject2)):
#	geompy.addToStudy( ContactObject2[i], 'ContactObject2_{}'.format(i) )
#for i in range(0,len(VolumeObject1)):
#	geompy.addToStudy( VolumeObject1[i], 'VolumeObject1_{}'.format(i) )
#for i in range(0,len(ContactObject1)):
#	geompy.addToStudy( ContactObject1[i], 'ContactObject1_{}'.format(i) )

geompy.addToStudy( Partition_profile, 'Partition_profile' )
for i_solid1 in range (0,len (Group_partition_volume)):
	geompy.addToStudyInFather( Partition_profile, Group_volume [i_solid1], Group_partition_volume[i_solid1][0])

for i_surface1 in range (0,len (Group_partition_surface)):
	geompy.addToStudyInFather( Partition_profile, Group_surface [i_surface1], Group_partition_surface[i_surface1][0])

##################################### end of extra code 3##########################################
###################################################################################################

#Meshing code

#Uncomment the following lines, if you want to check the  geoetrical surface of the contacts
#anode_surf=Contact1_Pro[0][1]
#cath_surf=Contact1_Pro[1][1]

#Add or delete a contact depending on the number of contacts. For example, if you have 3 contacts, then add Contact1_3=Group_surface[3]
Contact1_1=Group_surface[0]
Contact1_2=Group_surface[1]

encap_inner_ROI1=Group_volume[2]
encap_outer_ROI1=Group_volume[0]
ROI1=Group_volume[1]
Rest_1=Group_volume[3]

if(Lead2nd_Enable):
	#Adjust according to the number of contacts. If you have three contacts on the lead, then Contact2_1=Group_surface[3],Contact2_2=Group_surface[4],Contact2_3=Group_surface[5]
	Contact2_1=Group_surface[2]
	Contact2_2=Group_surface[3]

	encap_inner_ROI2=Group_volume[5]
	encap_outer_ROI2=Group_volume[3]
	ROI2=Group_volume[4]
	Rest_1=Group_volume[6]

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

#define maximum edge length for the submeshess (depends on the geometry of the electrode and the computational domain)
Max_size_rest=5.0
Max_size_ROI=0.75
Max_size_encap_outer=encap_thickness	 #should be not larger than the encap. thickness. Better a half of it.
Max_size_encap_inner=0.1
#define minimum edge length for the submeshes.  (for encap. entities and contacts it is recommended to assign a very small number.)
Min_size_rest=0.01
Min_size_ROI=0.001
Min_size_encap_outer=0.000000001	 #should be not larger than the encap. thickness. Better a half of it.
Min_size_encap_inner=0.00000001
#Another way to control the mesh quality is by changing the Fineness. Remember, this is an initial meshing, its aim is to preserve the shape of the geometries. Physics-based refinement will be conducted later. 


smesh = smeshBuilder.New(theStudy)
Mesh_1 = smesh.Mesh(Partition_profile)
NETGEN_1D_2D_3D = Mesh_1.Tetrahedron(algo=smeshBuilder.NETGEN_1D2D3D)
NETGEN_3D_Parameters_1 = NETGEN_1D_2D_3D.Parameters()
NETGEN_3D_Parameters_1.SetMaxSize( Max_size_rest )
NETGEN_3D_Parameters_1.SetSecondOrder( 0 )
NETGEN_3D_Parameters_1.SetOptimize( 1 )
NETGEN_3D_Parameters_1.SetFineness( 1 )			#increase if the mesh outside of ROI is too coarse
NETGEN_3D_Parameters_1.SetMinSize( Min_size_rest )
NETGEN_3D_Parameters_1.SetUseSurfaceCurvature( 1 )
NETGEN_3D_Parameters_1.SetFuseEdges( 1 )
NETGEN_3D_Parameters_1.SetQuadAllowed( 0 )

#Each domain (ROI,encap_inner_ROI, encap_outer_ROI, every contact and the rest) has its own submesh. Make sure that for each you have a unique "Sub_mesh" entity!
NETGEN_1D_2D = Mesh_1.Triangle(algo=smeshBuilder.NETGEN_1D2D,geom=Contact1_1)
Sub_mesh_1 = NETGEN_1D_2D.GetSubMesh()
NETGEN_2D_Parameters_1 = NETGEN_1D_2D.Parameters()
NETGEN_2D_Parameters_1.SetMaxSize( contact1_max_edge )	#!!! put here the maximum edge length
NETGEN_2D_Parameters_1.SetSecondOrder( 0 )
NETGEN_2D_Parameters_1.SetOptimize( 1 )
NETGEN_2D_Parameters_1.SetFineness( 4 )
NETGEN_2D_Parameters_1.SetMinSize( 1e-08 )				
NETGEN_2D_Parameters_1.SetUseSurfaceCurvature( 1 )
NETGEN_2D_Parameters_1.SetFuseEdges( 1 )
NETGEN_2D_Parameters_1.SetQuadAllowed( 0 )

#If two contacts have different requirements, create an additional submesh
NETGEN_1D_2D_1 = Mesh_1.Triangle(algo=smeshBuilder.NETGEN_1D2D,geom=Contact1_2)
Sub_mesh_2 = NETGEN_1D_2D_1.GetSubMesh()
NETGEN_2D_Parameters_2 = NETGEN_1D_2D_1.Parameters()
NETGEN_2D_Parameters_2.SetMaxSize( contact2_max_edge )  #!!! put here the maximum edge length
NETGEN_2D_Parameters_2.SetSecondOrder( 0 )
NETGEN_2D_Parameters_2.SetOptimize( 1 )
NETGEN_2D_Parameters_2.SetFineness( 4 )
NETGEN_2D_Parameters_2.SetMinSize( 1e-08 )
NETGEN_2D_Parameters_2.SetUseSurfaceCurvature( 1 )
NETGEN_2D_Parameters_2.SetFuseEdges( 1 )
NETGEN_2D_Parameters_2.SetQuadAllowed( 0 )

#if a contact has the same requirements as the first one
NETGEN_1D_2D_2 = Mesh_1.Triangle(algo=smeshBuilder.NETGEN_1D2D,geom=Contact1_3)
Sub_mesh_3 = NETGEN_1D_2D_2.GetSubMesh()
status = Mesh_1.AddHypothesis(NETGEN_2D_Parameters_1,Contact1_4)	#the hypothesis for meshing is the same as for the first contact
###NOTE: if you have only two contacts, then Sub_mesh_3 will be defined for the next domain, i.e. encap_inner_ROI1!!!


NETGEN_1D_2D_3D_1 = Mesh_1.Tetrahedron(algo=smeshBuilder.NETGEN_1D2D3D,geom=encap_inner_ROI1)
Sub_mesh_4 = NETGEN_1D_2D_3D_1.GetSubMesh()
NETGEN_3D_Parameters_2 = NETGEN_1D_2D_3D_1.Parameters()
NETGEN_3D_Parameters_2.SetMaxSize( Max_size_encap_inner )
NETGEN_3D_Parameters_2.SetSecondOrder( 0 )
NETGEN_3D_Parameters_2.SetOptimize( 1 )
NETGEN_3D_Parameters_2.SetFineness( 3 )
NETGEN_3D_Parameters_2.SetMinSize( Min_size_encap_inner )
NETGEN_3D_Parameters_2.SetUseSurfaceCurvature( 1 )
NETGEN_3D_Parameters_2.SetFuseEdges( 1 )
NETGEN_3D_Parameters_2.SetQuadAllowed( 0 )

NETGEN_1D_2D_3D_2 = Mesh_1.Tetrahedron(algo=smeshBuilder.NETGEN_1D2D3D,geom=encap_outer_ROI1)
Sub_mesh_5 = NETGEN_1D_2D_3D_2.GetSubMesh()
NETGEN_3D_Parameters_3 = NETGEN_1D_2D_3D_2.Parameters()
NETGEN_3D_Parameters_3.SetMaxSize( Max_size_encap_outer )
NETGEN_3D_Parameters_3.SetSecondOrder( 0 )
NETGEN_3D_Parameters_3.SetOptimize( 1 )
NETGEN_3D_Parameters_3.SetFineness( 2 )
NETGEN_3D_Parameters_3.SetMinSize( Min_size_encap_outer )
NETGEN_3D_Parameters_3.SetUseSurfaceCurvature( 1 )
NETGEN_3D_Parameters_3.SetFuseEdges( 1 )
NETGEN_3D_Parameters_3.SetQuadAllowed( 0 )

NETGEN_1D_2D_3D_3 = Mesh_1.Tetrahedron(algo=smeshBuilder.NETGEN_1D2D3D,geom=ROI1)
Sub_mesh_6 = NETGEN_1D_2D_3D_3.GetSubMesh()
NETGEN_3D_Parameters_4 = NETGEN_1D_2D_3D_3.Parameters()
NETGEN_3D_Parameters_4.SetMaxSize( Max_size_ROI )
NETGEN_3D_Parameters_4.SetSecondOrder( 0 )
NETGEN_3D_Parameters_4.SetOptimize( 1 )
NETGEN_3D_Parameters_4.SetFineness( 3 )
NETGEN_3D_Parameters_4.SetMinSize( Min_size_ROI )
NETGEN_3D_Parameters_4.SetUseSurfaceCurvature( 1 )
NETGEN_3D_Parameters_4.SetFuseEdges( 1 )
NETGEN_3D_Parameters_4.SetQuadAllowed( 0 )


#Lead2nd_Enable (bipolar stimulation) feature is disabled: such simulation is computationally expensive and pointless as the electrical fields from DBS do not intersect
if(Lead2nd_Enable):	
	NETGEN_1D_2D_1lead2 = Mesh_1.Triangle(algo=smeshBuilder.NETGEN_1D2D,geom=Contact2_1)
	Sub_mesh_1_2 = NETGEN_1D_2D_1lead2.GetSubMesh()
	NETGEN_2D_Parameters_1lead2 = NETGEN_1D_2D_1lead2.Parameters()
	NETGEN_2D_Parameters_1lead2 =NETGEN_2D_Parameters_1

	NETGEN_1D_2D_2lead2 = Mesh_1.Triangle(algo=smeshBuilder.NETGEN_1D2D,geom=Contact2_2)
	Sub_mesh_2_2 = NETGEN_1D_2D_2lead2.GetSubMesh()
	NETGEN_2D_Parameters_2lead2 = NETGEN_1D_2D_2lead2.Parameters()
	NETGEN_2D_Parameters_2lead2 =NETGEN_2D_Parameters_2

	NETGEN_1D_2D_3D_1lead2 = Mesh_1.Tetrahedron(algo=smeshBuilder.NETGEN_1D2D3D,geom=encap_inner_ROI2)
	Sub_mesh_3_2 = NETGEN_1D_2D_3D_1lead2.GetSubMesh()
	NETGEN_3D_Parameters_2lead2 = NETGEN_1D_2D_3D_1lead2.Parameters()
	NETGEN_3D_Parameters_2lead2 =NETGEN_3D_Parameters_2

	NETGEN_1D_2D_3D_2lead2 = Mesh_1.Tetrahedron(algo=smeshBuilder.NETGEN_1D2D3D,geom=encap_outer_ROI2)
	Sub_mesh_4_2 = NETGEN_1D_2D_3D_2lead2.GetSubMesh()
	NETGEN_3D_Parameters_3lead2 = NETGEN_1D_2D_3D_2lead2.Parameters()
	NETGEN_3D_Parameters_3lead2 =NETGEN_3D_Parameters_3

	NETGEN_1D_2D_3D_3lead2 = Mesh_1.Tetrahedron(algo=smeshBuilder.NETGEN_1D2D3D,geom=ROI2)
	Sub_mesh_5_2 = NETGEN_1D_2D_3D_3lead2.GetSubMesh()
	NETGEN_3D_Parameters_4lead2 = NETGEN_1D_2D_3D_3lead2.Parameters()
	NETGEN_3D_Parameters_4lead2 =NETGEN_3D_Parameters_4

#Here the parameters should be as for the Mesh_1
NETGEN_1D_2D_3D_4 = Mesh_1.Tetrahedron(algo=smeshBuilder.NETGEN_1D2D3D,geom=Rest_1)
Sub_mesh_7 = NETGEN_1D_2D_3D_4.GetSubMesh()
NETGEN_3D_Parameters_5 = NETGEN_1D_2D_3D_4.Parameters()
NETGEN_3D_Parameters_5.SetMaxSize( Max_size_rest )
NETGEN_3D_Parameters_5.SetSecondOrder( 0 )
NETGEN_3D_Parameters_5.SetOptimize( 1 )
NETGEN_3D_Parameters_5.SetFineness( 1 )
NETGEN_3D_Parameters_5.SetMinSize( Min_size_rest )
NETGEN_3D_Parameters_5.SetUseSurfaceCurvature( 1 )
NETGEN_3D_Parameters_5.SetFuseEdges( 1 )
NETGEN_3D_Parameters_5.SetQuadAllowed( 0 )


if(Lead2nd_Enable):
	isDone = Mesh_1.SetMeshOrder( [ [ Sub_mesh_1,Sub_mesh_1_2, Sub_mesh_2,Sub_mesh_2_2, Sub_mesh_3,Sub_mesh_3_2, Sub_mesh_4,Sub_mesh_4_2, Sub_mesh_5,Sub_mesh_5_2, Sub_mesh_6 ] ])
else:
#be sure you have all the submeshes here (in ascending order)
	isDone = Mesh_1.SetMeshOrder( [ [ Sub_mesh_1, Sub_mesh_2, Sub_mesh_3, Sub_mesh_4, Sub_mesh_5, Sub_mesh_6 ] ])
isDone = Mesh_1.Compute()

#Modify according to the number of contacts (For example, if there are three, add C1_3 = Mesh_1.GroupOnGeom(Contact1_3,'C1_3',SMESH.FACE))
C1_1 = Mesh_1.GroupOnGeom(Contact1_1,'C1_1',SMESH.FACE)
C1_2 = Mesh_1.GroupOnGeom(Contact1_2,'C1_2',SMESH.FACE)

Encap_rest = Mesh_1.GroupOnGeom(encap_outer_ROI1,'Encap_rest',SMESH.VOLUME)
Encap_contact = Mesh_1.GroupOnGeom(encap_inner_ROI1,'Encap_contact',SMESH.VOLUME)
RegOfInt = Mesh_1.GroupOnGeom(ROI1,'RegOfInt',SMESH.VOLUME)
Rst = Mesh_1.GroupOnGeom(Rest_1,'Rst',SMESH.VOLUME)
if(Lead2nd_Enable):
	C2_1 = Mesh_1.GroupOnGeom(Contact2_1,'C2_1',SMESH.FACE)
	C2_2 = Mesh_1.GroupOnGeom(Contact2_2,'C2_2',SMESH.FACE)
	Encap_rest2 = Mesh_1.GroupOnGeom(encap_outer_ROI2,'Encap_rest2',SMESH.VOLUME)
	Encap_contact2 = Mesh_1.GroupOnGeom(encap_inner_ROI2,'Encap_contact2',SMESH.VOLUME)
	RegOfInt2 = Mesh_1.GroupOnGeom(ROI2,'RegOfInt2',SMESH.VOLUME)

#Uncomment the following lines, if you want to check the surface of the meshed contacts (helpful to check the accuracy of meshing)
#measure = smesh.CreateMeasurements()
#An_meshed_surf=measure.Area(C1_1)
#Cat_meshed_surf=measure.Area(C1_2)
#print "Contact_1 mesh surface: "
#print abs(An_meshed_surf-anode_surf)/anode_surf
#print "Contact_1 mesh surface: "
#print abs(Cat_meshed_surf-cath_surf)/cath_surf


#Uncomment to evaluate the meshed surface against the geometrical for the first two contacts
#while abs(An_meshed_surf-anode_surf)/anode_surf>0.01:
#	contact_mesh_max=contact_mesh_max/2.0
#	NETGEN_2D_Parameters_1.SetMaxSize( anode_mesh_max )
#	isDone = Mesh_1.Compute()

#while abs(Cat_meshed_surf-cath_surf)/cath_surf>0.01:
#	contact_mesh_max=contact_mesh_max/2.0
#	NETGEN_2D_Parameters_2.SetMaxSize( cathode_mesh_max )
#	isDone = Mesh_1.Compute()


## Set names of Mesh objects
smesh.SetName(NETGEN_1D_2D_3D.GetAlgorithm(), 'NETGEN 1D-2D-3D')
smesh.SetName(NETGEN_1D_2D.GetAlgorithm(), 'NETGEN 1D-2D')

#Modify according to the number of contacts and the number of unique hypothesis for their Sub_meshes!
smesh.SetName(C1_1, 'C1_1')
smesh.SetName(C1_2, 'C1_2')
smesh.SetName(NETGEN_2D_Parameters_1, 'NETGEN 2D Parameters_1')
smesh.SetName(NETGEN_2D_Parameters_2, 'NETGEN 2D Parameters_2')
smesh.SetName(Sub_mesh_1, 'Sub-mesh_1')
smesh.SetName(Sub_mesh_2, 'Sub-mesh_2')
#if you have three contacts, add:
#smesh.SetName(C1_3, 'C1_3')
#smesh.SetName(Sub_mesh_3, 'Sub-mesh_3')
#if it has a unique hypothesis, add: 
#smesh.SetName(NETGEN_2D_Parameters_2, 'NETGEN 2D Parameters_3')    

smesh.SetName(NETGEN_3D_Parameters_2, 'NETGEN 3D Parameters_2')
smesh.SetName(NETGEN_3D_Parameters_3, 'NETGEN 3D Parameters_3')
smesh.SetName(NETGEN_3D_Parameters_1, 'NETGEN 3D Parameters_1')

if(Lead2nd_Enable):
    smesh.SetName(NETGEN_2D_Parameters_1lead2, 'NETGEN_2D_Parameters_1lead2')
    smesh.SetName(NETGEN_2D_Parameters_2lead2, 'NETGEN_2D_Parameters_2lead2')
    smesh.SetName(NETGEN_3D_Parameters_2lead2, 'NETGEN_3D_Parameters_2lead2')
    smesh.SetName(NETGEN_3D_Parameters_3lead2, 'NETGEN_3D_Parameters_3lead2')
    smesh.SetName(NETGEN_3D_Parameters_4lead2, 'NETGEN_3D_Parameters_4lead2')

smesh.SetName(NETGEN_3D_Parameters_4, 'NETGEN 3D Parameters_4')
smesh.SetName(NETGEN_3D_Parameters_5, 'NETGEN 3D Parameters_5')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')

smesh.SetName(Rst, 'Rst')
smesh.SetName(RegOfInt, 'RegOfInt')
smesh.SetName(Encap_rest, 'Encap_rest')
smesh.SetName(Encap_contact, 'Encap_contact')


#!!! Check, that the name is set for all previously defined Sub_meshes. For example, you will need to add smesh.SetName(Sub_mesh_7, 'Sub-mesh_7'), if you have 3 contacts on the lead.
smesh.SetName(Sub_mesh_6, 'Sub-mesh_6')
smesh.SetName(Sub_mesh_5, 'Sub-mesh_5')
smesh.SetName(Sub_mesh_3, 'Sub-mesh_3')
smesh.SetName(Sub_mesh_4, 'Sub-mesh_4')

if(Lead2nd_Enable):
	smesh.SetName(C2_1, 'C2_1')
	smesh.SetName(C2_2, 'C2_2')
	smesh.SetName(Sub_mesh_5_2, 'Sub-mesh_5_2')
	smesh.SetName(RegOfInt2, 'RegOfInt2')
	smesh.SetName(Encap_rest2, 'Encap_rest2')
	smesh.SetName(Encap_contact2, 'Encap_contact2')
	smesh.SetName(Sub_mesh_3_2, 'Sub-mesh_3_2')
	smesh.SetName(Sub_mesh_2_2, 'Sub-mesh_2_2')
	smesh.SetName(Sub_mesh_1_2, 'Sub-mesh_1_2')
	smesh.SetName(Sub_mesh_4_2, 'Sub-mesh_4_2')

#to export the mesh
Mesh_1.ExportMED('Meshes/Mesh_unref.med')


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser(True)
