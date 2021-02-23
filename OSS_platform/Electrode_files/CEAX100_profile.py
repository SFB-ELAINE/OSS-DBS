# -*- coding: utf-8 -*-

###
### This file is generated automatically by SALOME v8.3.0 with dump python functionality
###

import sys
import salome
import os

salome.salome_init()
theStudy = salome.myStudy

import salome_notebook
notebook = salome_notebook.NoteBook(theStudy)
sys.path.insert( 0, "r'"+os.getcwd())
sys.path.append('/usr/local/lib/python2.7/dist-packages')

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
#tip_length = 0.0
#contact_shift = 0.0

##### VARIABLE LIST #####


########## End of variable list#############

if Z_2nd == Zt:
    Z_2nd_artif = Zt+1.0 # just to ensure the rotation is possible
else:
    Z_2nd_artif=Z_2nd


####################################################
##       Begin of NoteBook variables section      ##
####################################################
#notebook.set("encap_thickness", 0.2)
#notebook.set("ROI_radial", 2)


### Change the radius (up until 0.05) and the position of the new ground electrode
notebook.set("ground_start", 3.0)              # 2.22 is the lowest (on the border with ROI)
#notebook.set("ground_length", "7.0-ground_start")
notebook.set("ground_length", 2.0)     
         
notebook.set("fillet_radius", 0.01)		# max(R_tip) is 0.01 mm (fixed)
fillet_radius=0.01


notebook.set("encap_thickness", 0.1)
notebook.set("ROI_radial", ROI_radial)
notebook.set("encap_body_radial", "0.2+encap_thickness")
notebook.set("encap_cone1_radial", "0.15+encap_thickness")
notebook.set("encap_cone2_radial", "0.075+encap_thickness")
notebook.set("encap_cone3_radial", "0.01+encap_thickness")
notebook.set("encap_sphere1_radial", "0.01+encap_thickness")

#notebook.set("tip_length", 0.175)		# at least 0.05, at most 0.175
#tip_length=0.175
#contact_shift=0.1
#notebook.set("contact_shift", 0.1)		# at most 0.2
####################################################
##        End of NoteBook variables section       ##
####################################################
###
### GEOM component
###



########## End of variable list#############

#Contact_base_radius=0.15-(0.375*contact_shift)

#R_cut=Contact_base_radius-0.375*tip_length
#encap_cone3_radial=R_cut+(encap_thickness-0.375*encap_thickness)


notebook.set("encap_body_radial", 0.2+encap_thickness)
#notebook.set("encap_cone1_radial", "0.15+encap_thickness")
#notebook.set("encap_body_radial", "0.2+encap_thickness")
notebook.set("encap_cone1_radial", "0.15+encap_thickness")
notebook.set("encap_cone2_radial", "0.075+encap_thickness")
notebook.set("encap_cone3_radial", "0.01+encap_thickness")
notebook.set("encap_sphere1_radial", "0.01+encap_thickness")
#notebook.set("encap_cone1_radial",Contact_base_radius+encap_thickness)
#notebook.set("encap_cone2_radial", 0.075+encap_thickness)
#notebook.set("encap_cone3_radial", 0.06+encap_thickness)
#notebook.set("encap_sphere1_radial", 0.06+encap_thickness)

Vert_array =[0];
number_vertex = len(Vert_array)
Vert = []
VolumeObject1 = []
ContactObject1 = []
VolumeObject2 = []
ContactObject2 = []
print(" DBS_lead's Geometry\n")


import GEOM
from salome.geom import geomBuilder
import math
import SALOMEDS


geompy = geomBuilder.New(theStudy)

O = geompy.MakeVertex(0, 0, 0)
OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)
lead_body = geompy.MakeCylinderRH(0.2, 7.5)
lead_Cone_1 = geompy.MakeConeR1R2H(0.15, 0.2, 0.1333)
geompy.TranslateDXDYDZ(lead_Cone_1, 0, 0, -0.1333)
lead_Cone_2 = geompy.MakeConeR1R2H(0.075, 0.15, 0.2)
geompy.TranslateDXDYDZ(lead_Cone_2, 0, 0, -0.3333)
lead_Cone_3 = geompy.MakeConeR1R2H(0.01, 0.075, 0.17)
geompy.TranslateDXDYDZ(lead_Cone_3, 0, 0, -0.5033)
lead_Sphere_tip = geompy.MakeSphereR(0.01)
encap_body = geompy.MakeCylinderRH("encap_body_radial", 7.5)
encap_Cone_1 = geompy.MakeConeR1R2H("encap_cone1_radial", "encap_body_radial", 0.1333)
geompy.TranslateDXDYDZ(encap_Cone_1, 0, 0, -0.1333)
encap_Cone_2 = geompy.MakeConeR1R2H("encap_cone2_radial", "encap_cone1_radial", 0.2)
geompy.TranslateDXDYDZ(encap_Cone_2, 0, 0, -0.3333)
encap_Cone_3 = geompy.MakeConeR1R2H("encap_cone3_radial", "encap_cone2_radial", 0.17)
geompy.TranslateDXDYDZ(encap_Cone_3, 0, 0, -0.5033)
encap_Sphere_tip = geompy.MakeSphereR("encap_sphere1_radial")
geompy.TranslateDXDYDZ(encap_Sphere_tip, 0, 0, -0.5033)
geompy.TranslateDXDYDZ(lead_Sphere_tip, 0, 0, -0.5033)
encap_cone_object = geompy.MakeFuseList([encap_Cone_1, encap_Cone_2, encap_Cone_3, encap_Sphere_tip], True, True)
lead_cone_all = geompy.MakeFuseList([lead_Cone_1, lead_Cone_2, lead_Cone_3, lead_Sphere_tip], True, True)
encap_cone_layer = geompy.MakeCutList(encap_cone_object, [lead_cone_all], True)

dbs_lead = geompy.MakeFuseList([lead_body, lead_Cone_1, lead_Cone_2, lead_Cone_3, lead_Sphere_tip], True, True)

encap_body_full = geompy.MakeFuseList([encap_body, encap_cone_layer], True, True)
encap_complete_body_layer = geompy.MakeCutList(encap_body_full, [dbs_lead], True)

listFreeFacesIDs = geompy.GetFreeFacesIDs(dbs_lead)
listFreeFacesIDs = geompy.GetFreeFacesIDs(dbs_lead)
[Free_face_1_1,body_surface,Free_face_1_3,Free_face_1_4,Free_face_1_5,Free_face_1_6] = geompy.SubShapes(dbs_lead, [3, 7, 12, 17, 22, 27])
Sphere_ROI = geompy.MakeSpherePntR(O, "ROI_radial")
geompy.TranslateDXDYDZ(lead_body, 0, 0, 0.5133)
geompy.TranslateDXDYDZ(lead_Cone_1, 0, 0, 0.5133)
listFreeFacesIDs = geompy.GetFreeFacesIDs(lead_Cone_1)
listFreeFacesIDs = geompy.GetFreeFacesIDs(lead_Cone_1)
[Contact_2,geomObj_12,geomObj_13] = geompy.SubShapes(lead_Cone_1, [3, 10, 12])
Contact_2 = geompy.MakeFuseList([Contact_2], True, True)
geompy.TranslateDXDYDZ(lead_Cone_2, 0, 0, 0.5133)
geompy.TranslateDXDYDZ(lead_Cone_3, 0, 0, 0.5133)
listFreeFacesIDs = geompy.GetFreeFacesIDs(lead_Cone_3)
listFreeFacesIDs = geompy.GetFreeFacesIDs(lead_Cone_3)
geompy.TranslateDXDYDZ(lead_Sphere_tip, 0, 0, 0.5133)
listFreeFacesIDs = geompy.GetFreeFacesIDs(lead_Sphere_tip)
listFreeFacesIDs = geompy.GetFreeFacesIDs(lead_Sphere_tip)
[sphere_tip_surface] = geompy.SubShapes(lead_Sphere_tip, [3])
geompy.TranslateDXDYDZ(encap_body, 0, 0, 0.5133)
geompy.TranslateDXDYDZ(encap_Cone_1, 0, 0, 0.5133)
geompy.TranslateDXDYDZ(encap_Cone_2, 0, 0, 0.5133)
geompy.TranslateDXDYDZ(encap_Cone_3, 0, 0, 0.5133)
geompy.TranslateDXDYDZ(encap_Sphere_tip, 0, 0, 0.5133)
geompy.TranslateDXDYDZ(encap_cone_object, 0, 0, 0.5133)
geompy.TranslateDXDYDZ(lead_cone_all, 0, 0, 0.5133)
geompy.TranslateDXDYDZ(encap_cone_layer, 0, 0, 0.5133)
geompy.TranslateDXDYDZ(encap_complete_body_layer, 0, 0, 0.5133)
geompy.TranslateDXDYDZ(dbs_lead, 0, 0, 0.5133)
geompy.TranslateDXDYDZ(Sphere_ROI, 0, 0, 0.5133)
Lead_and_encap_fuse = geompy.MakeFuseList([encap_cone_layer, encap_complete_body_layer, dbs_lead], True, True)

ROI = geompy.MakeCutList(Sphere_ROI, [Lead_and_encap_fuse], True)

encap_inner_ROI = geompy.MakeCommonList([encap_complete_body_layer, Sphere_ROI], True)
encap_outer_ROI = geompy.MakeCutList(encap_complete_body_layer, [Sphere_ROI], True)
Contact_1_solid = geompy.MakeFuseList([lead_Cone_3, lead_Sphere_tip], True, True)
listFreeFacesIDs = geompy.GetFreeFacesIDs(Contact_1_solid)
listFreeFacesIDs = geompy.GetFreeFacesIDs(Contact_1_solid)
[geomObj_23,Free_face_1_7,Free_face_1_8] = geompy.SubShapes(Contact_1_solid, [3, 7, 12])
Contact_1 = geompy.MakeFuseList([Free_face_1_7, Free_face_1_8], True, True)

Fuse_all_lead_encap_ROI = geompy.MakeFuseList([lead_body, lead_Cone_1, lead_Cone_2, lead_Cone_3, lead_Sphere_tip, encap_body, encap_Cone_1, encap_Cone_2, encap_Cone_3, encap_Sphere_tip, Sphere_ROI], True, True)
Fuse_all = geompy.MakeFuseList([ROI, encap_inner_ROI, encap_outer_ROI ,dbs_lead ], True, True) #python equivalent

geompy.addToStudy( ROI, 'ROI' )
geompy.addToStudy( encap_inner_ROI, 'encap_inner_ROI' )
geompy.addToStudy( encap_outer_ROI, 'encap_outer_ROI' )
geompy.addToStudy( Contact_2, 'Contact_2_init' )
geompy.addToStudy( Contact_1, 'Contact_1_init' )


######################################################################################################
########################################### extra code 2 V10 15/12/18#############################################
print " Load brain image \n"
if (Brain_map[-4:] == 'brep'):
	brain_solid = geompy.ImportBREP( Brain_map )
elif (Brain_map[-4:] == 'step'):
	brain_solid = geompy.ImportSTEP( Brain_map )
elif (Brain_map[-4:] == 'iges'):
	brain_solid = geompy.ImportIGES( Brain_map )
elif (Brain_map[-4:] == '.stl'):
	brain_solid = geompy.ImportSTL( Brain_map )
else:
	print("unknow imported file format")

Fuse_all_lead_encap_ROI_no_internal_face = geompy.RemoveInternalFaces(Fuse_all_lead_encap_ROI)
#################################################### Geometry and extra code interface ###############################################################
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
  print "Position 2nd Fuse all object at [{},{},{}], [{}',{}',{}']\n".format(Xt2,Yt2,Zt2,OX_angle2,OY_angle2,OZ_angle2)
  Fuse_all_lead_encap_ROI_no_internal_face2 = geompy.MakeTranslation(Fuse_all_lead_encap_ROI_no_internal_face,Xt2,Yt2,Zt2)
  OX2 = geompy.MakeTranslation(OX,Xt2,Yt2,Zt2)
  OY2 = geompy.MakeTranslation(OY,Xt2,Yt2,Zt2)
  OZ2 = geompy.MakeTranslation(OZ,Xt2,Yt2,Zt2)
  geompy.Rotate(Fuse_all_lead_encap_ROI_no_internal_face2, OX2,OX_angle2*math.pi/180.0)
  geompy.Rotate(Fuse_all_lead_encap_ROI_no_internal_face2, OY2,OY_angle2*math.pi/180.0)
  geompy.Rotate(Fuse_all_lead_encap_ROI_no_internal_face2, OZ2,OZ_angle2*math.pi/180.0)
  print "Position 2nd Lead at [{},{},{}], [{}',{}',{}']\n".format(Xt2,Yt2,Zt2,OX_angle2,OY_angle2,OZ_angle2)
  for i in range(0,len(VolumeObject1)):
	VolumeObject2[i] = geompy.MakeTranslation(VolumeObject1[i],Xt2,Yt2,Zt2)
	geompy.Rotate(VolumeObject2[i], OX2,OX_angle2*math.pi/180.0)
	geompy.Rotate(VolumeObject2[i], OY2,OY_angle2*math.pi/180.0)
	geompy.Rotate(VolumeObject2[i], OZ2,OZ_angle2*math.pi/180.0)
  for i in range(0,len(ContactObject1)):
        #geompy.addToStudy( Contact_1, 'Contact_1' )
	ContactObject2[i] = geompy.MakeTranslation(ContactObject1[i],Xt2,Yt2,Zt2)
	geompy.Rotate(ContactObject2[i], OX2,OX_angle2*math.pi/180.0)
	geompy.Rotate(ContactObject2[i], OY2,OY_angle2*math.pi/180.0)
	geompy.Rotate(ContactObject2[i], OZ2,OZ_angle2*math.pi/180.0)
  print "Cut outer ROI2 with brain\n"
  cut_outer_ROI = geompy.MakeCutList(VolumeObject2[0], [brain_solid], True)
  VolumeObject2[0] = geompy.MakeCutList(VolumeObject2[0], [cut_outer_ROI], True)
  print "Cut ROI2 with brain\n"
  VolumeObject2[1] = geompy.MakeCommonList([VolumeObject2[1], brain_solid], True) 
  print "Group 2nd:volume and area extraction for group ID identification process\n"       
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
    geompy.addToStudy(ContactObject1[i], 'Contact_'+str(i+10)) 
    geompy.TranslateDXDYDZ(ContactObject1[i],Xt,Yt,Zt)
    geompy.Rotate(ContactObject1[i], OZ1,OZ_angle*math.pi/180.0)
    if X_2nd!=Xt or Y_2nd!=Yt:
        ContactObject1[i]=geompy.MakeRotationThreePoints(ContactObject1[i], Vertex_O, Vertex_3, Vertex_1)


print "Cut outer ROI1 with brain\n"
cut_outer_ROI = geompy.MakeCutList(VolumeObject1[0], [brain_solid], True)
VolumeObject1[0] = geompy.MakeCutList(VolumeObject1[0], [cut_outer_ROI], True)
print "Cut ROI1 with brain\n"
VolumeObject1[1] = geompy.MakeCommonList([VolumeObject1[1], brain_solid], True)
print "Group 1st:volume and area extraction for group ID identification process\n" 

Volume1_Pro = [geompy.BasicProperties( VolumeObject1[0])]*len(VolumeObject1)
Contact1_Pro = [geompy.BasicProperties( ContactObject1[0])]*len(ContactObject1)
for i in range(0,len(VolumeObject1)):
	Volume1_Pro[i] = geompy.BasicProperties( VolumeObject1[i])
for i in range(0,len(ContactObject1)):
	Contact1_Pro[i] = geompy.BasicProperties( ContactObject1[i])


print "Create reference groups for ID identification process\n"

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
  geompy.addToStudy(Partition_profile, 'Partition_profile_check') 
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
print "Partition_volume_IDsList",Partition_volume_IDsList, '\n'

for ref_ind in range (0, len(reference_volume)):
	temp_volume = []
	for sub_ind in range (0, len (Partition_volume_IDsList)):
		subshape = geompy.GetSubShape(Partition_profile, [Partition_volume_IDsList[sub_ind]]) # get subshape
		subshape_Pro = geompy.BasicProperties(subshape)       # extract volume of subshape
		Common_volume = geompy.MakeCommonList([subshape, reference_volume[ref_ind]], True) # check common intersection
		Common_volume_Pro = geompy.BasicProperties(Common_volume) 
		print "volume difference",abs(Common_volume_Pro[2]-subshape_Pro[2]),"/",abs(Common_volume_Pro[2]-reference_volume_Pro[ref_ind][2])
		# if ( common volume = subshape) and (common volume = ref volume) => ref volume = sub shape
		if (abs(Common_volume_Pro[2]-subshape_Pro[2])< 0.00001) and (abs(Common_volume_Pro[2]-reference_volume_Pro[ref_ind][2])<0.00001):
		
			Group_partition_volume.append([Volume_name[ref_ind],Partition_volume_IDsList[sub_ind]])
		# if ( common volume = subshape) and (common volume < ref volume) => sub shape belong to ref volume
		elif (abs(Common_volume_Pro[2]-subshape_Pro[2])< 0.00001) and ((Common_volume_Pro[2] - reference_volume_Pro[ref_ind][2])<-0.00001):
			temp_volume.append( Partition_volume_IDsList[sub_ind] )
	if len(temp_volume) >1 : # the volume is devided
		Group_partition_volume.append([Volume_name[ref_ind],temp_volume ])
		print Volume_name[ref_ind]," is devided and has sub IDs:{}\n".format(temp_volume)
if len(reference_volume) != len(Group_partition_volume):
	print "Geometry-volume error please check ROI diameter and DBS lead Position ",len(reference_volume),len(Group_partition_volume)
print 'Group_partition_volume',Group_partition_volume,'\n'

### find group surface ID ######################################################################
Partition_surface_IDsList = geompy.SubShapeAllIDs(Partition_profile, geompy.ShapeType["FACE"]) # list all sub shape face in Partition
print 'Partition_surface_IDsList',Partition_surface_IDsList,'\n'
sub_face = [] ## store devided faces
for reff_ind in range (0, len (reference_surface)):
	temp_surface = []
	for subf_ind in range (0, len(Partition_surface_IDsList)):
		subshapef = geompy.GetSubShape(Partition_profile, [Partition_surface_IDsList[subf_ind]]) # get subshape
		Common_face = geompy.MakeCommonList([subshapef, reference_surface[reff_ind]], True) # check common intersection
		Common_face_Pro = geompy.BasicProperties(Common_face) 
		subshapef_Pro = geompy.BasicProperties(subshapef) # extract volume of subshape
		print "area difference",abs(Common_face_Pro[1]-subshapef_Pro[1]),"/",abs(Common_face_Pro[1]-reference_surface_Pro[reff_ind][1])
		# if ( common face = subface) and (common face = ref face) => ref face = sub face
		if (abs(Common_face_Pro[1]-subshapef_Pro[1])<0.000001 )and (abs(Common_face_Pro[1]-reference_surface_Pro[reff_ind][1])<0.000001): 
			Group_partition_surface.append([ Contact_name[reff_ind],Partition_surface_IDsList[subf_ind] ])
		# if ( common face = subface) and (common face < ref face) => sub face belong to ref face 
		elif (abs(Common_face_Pro[1]-subshapef_Pro[1])<0.000001 ) and ((Common_face_Pro[1] - reference_surface_Pro[reff_ind][1])<-0.000001):
			temp_surface.append(Partition_surface_IDsList[subf_ind])
	if len(temp_surface) >1 : # the face is devided
		Group_partition_surface.append( [Contact_name[reff_ind],temp_surface ])		
		print Contact_name[reff_ind]," is devided and has sub IDs:{}\n".format(temp_surface)
if len(reference_surface) != len(Group_partition_surface): #+len(Group_partition_Multi_surface):
	print "Geometry-Surface error please check ROI diameter and DBS lead Position ",len(reference_surface),len(Group_partition_surface),'\n' 

print 'Group_partition_surface',Group_partition_surface,'\n'

if(Lead2nd_Enable):
   Partition_profile = geompy.MakePartition(VolumeObject1+VolumeObject2+ContactObject1+ContactObject2+[Rest], [], [], [], geompy.ShapeType["SOLID"], 0, [], 0)
else:
   Partition_profile = geompy.MakePartition(VolumeObject1+ContactObject1+[Rest], [], [], [], geompy.ShapeType["SOLID"], 0, [], 0)

new_volume_ID= geompy.SubShapeAllIDs(Partition_profile, geompy.ShapeType["SOLID"])
ID= list(set(Partition_volume_IDsList) ^ set (new_volume_ID))
Group_partition_volume.append(['Rest_1',ID[0]])
print "REST ID:",ID
print 'Group_partition_volume',Group_partition_volume,'\n'
print"Create volume and surface group under partition_profile\n"

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
print "Translate whole partition to Xm,Ym,Zm\n"
geompy.TranslateDXDYDZ(Partition_profile, Xm, Ym, Zm)


#additional to put exactly at the imp. point
Tip_point = geompy.CreateGroup(Partition_profile, geompy.ShapeType["VERTEX"])
#geompy.UnionIDs(Tip_point, [85])   #hamster
geompy.UnionIDs(Tip_point, [81])
coords_tip = geompy.PointCoordinates(Tip_point)
#tip_shift=(Yt+Ym)-coords_tip[1]
tip_shift=(Zt+Zm)-coords_tip[2]
#geompy.TranslateDXDYDZ(Partition_profile, 0.0, tip_shift, 0.0)
geompy.TranslateDXDYDZ(Partition_profile, 0.0, 0.0, tip_shift)		

### add Vertices to geometry
if(Vertice_enable):
   for ver_ind in range (0,number_vertex):
       print"Add vertices to model\n"
       Vert.append(geompy.MakeVertex(Vert_array[ver_ind][0],Vert_array[ver_ind][1],Vert_array[ver_ind][2]))
       geompy.TranslateDXDYDZ(Vert[ver_ind], Xm, Ym, Zm)       ###Translate vertices to Xm,Ym,Zm
       geompy.addToStudy( Vert[ver_ind], 'Vert_{}'.format(ver_ind))
print"add to study\n"
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

geompy.addToStudyInFather( Partition_profile, Tip_point, 'Tip_point' )

##################################### end of extra code 3##########################################
###################################################################################################
anode_surf=Contact1_Pro[0][1]
cath_surf=Contact1_Pro[1][1]


Contact1_1=Group_surface[0]
Contact1_2=Group_surface[1]

encap_inner_ROI1=Group_volume[2]
encap_outer_ROI1=Group_volume[0]
ROI1=Group_volume[1]
Rest_1=Group_volume[3]

if(Lead2nd_Enable):
	Contact2_1=Group_surface[2]
	Contact2_2=Group_surface[3]
	encap_inner_ROI2=Group_volume[5]
	encap_outer_ROI2=Group_volume[3]
	ROI2=Group_volume[4]
	Rest_1=Group_volume[6]

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

anode_mesh_max=0.005
cathode_mesh_max=0.1

smesh = smeshBuilder.New(theStudy)
Mesh_1 = smesh.Mesh(Partition_profile)
NETGEN_1D_2D_3D = Mesh_1.Tetrahedron(algo=smeshBuilder.NETGEN_1D2D3D)
NETGEN_3D_Parameters_1 = NETGEN_1D_2D_3D.Parameters()
NETGEN_3D_Parameters_1.SetMaxSize( 2.25 )
NETGEN_3D_Parameters_1.SetSecondOrder( 0 )
NETGEN_3D_Parameters_1.SetOptimize( 1 )
NETGEN_3D_Parameters_1.SetFineness( 0 )
NETGEN_3D_Parameters_1.SetMinSize( 0.1 )
NETGEN_3D_Parameters_1.SetUseSurfaceCurvature( 1 )
NETGEN_3D_Parameters_1.SetFuseEdges( 1 )
NETGEN_3D_Parameters_1.SetQuadAllowed( 0 )
NETGEN_1D_2D = Mesh_1.Triangle(algo=smeshBuilder.NETGEN_1D2D,geom=Contact1_1)
Sub_mesh_1 = NETGEN_1D_2D.GetSubMesh()
NETGEN_2D_Parameters_1 = NETGEN_1D_2D.Parameters()
NETGEN_2D_Parameters_1.SetMaxSize( anode_mesh_max )
NETGEN_2D_Parameters_1.SetSecondOrder( 0 )
NETGEN_2D_Parameters_1.SetOptimize( 1 )
NETGEN_2D_Parameters_1.SetFineness( 4 )
NETGEN_2D_Parameters_1.SetMinSize( 1e-05 )
NETGEN_2D_Parameters_1.SetUseSurfaceCurvature( 1 )
NETGEN_2D_Parameters_1.SetFuseEdges( 1 )
NETGEN_2D_Parameters_1.SetQuadAllowed( 0 )
NETGEN_1D_2D_1 = Mesh_1.Triangle(algo=smeshBuilder.NETGEN_1D2D,geom=Contact1_2)
Sub_mesh_2 = NETGEN_1D_2D_1.GetSubMesh()
NETGEN_2D_Parameters_2 = NETGEN_1D_2D_1.Parameters()
NETGEN_2D_Parameters_2.SetMaxSize( cathode_mesh_max )
NETGEN_2D_Parameters_2.SetSecondOrder( 0 )
NETGEN_2D_Parameters_2.SetOptimize( 1 )
NETGEN_2D_Parameters_2.SetFineness( 4 )
NETGEN_2D_Parameters_2.SetMinSize( 1e-05 )
NETGEN_2D_Parameters_2.SetUseSurfaceCurvature( 1 )
NETGEN_2D_Parameters_2.SetFuseEdges( 1 )
NETGEN_2D_Parameters_2.SetQuadAllowed( 0 )
NETGEN_1D_2D_3D_1 = Mesh_1.Tetrahedron(algo=smeshBuilder.NETGEN_1D2D3D,geom=encap_inner_ROI1)
Sub_mesh_3 = NETGEN_1D_2D_3D_1.GetSubMesh()
NETGEN_3D_Parameters_2 = NETGEN_1D_2D_3D_1.Parameters()
NETGEN_3D_Parameters_2.SetMaxSize( encap_thickness )
NETGEN_3D_Parameters_2.SetSecondOrder( 0 )
NETGEN_3D_Parameters_2.SetOptimize( 1 )
NETGEN_3D_Parameters_2.SetFineness( 2 )
NETGEN_3D_Parameters_2.SetMinSize( 0.0001 )
NETGEN_3D_Parameters_2.SetUseSurfaceCurvature( 1 )
NETGEN_3D_Parameters_2.SetFuseEdges( 1 )
NETGEN_3D_Parameters_2.SetQuadAllowed( 0 )
isDone = Mesh_1.SetMeshOrder( [ [ Sub_mesh_1, Sub_mesh_2, Sub_mesh_3 ] ])
NETGEN_1D_2D_3D_2 = Mesh_1.Tetrahedron(algo=smeshBuilder.NETGEN_1D2D3D,geom=encap_outer_ROI1)
Sub_mesh_4 = NETGEN_1D_2D_3D_2.GetSubMesh()
NETGEN_3D_Parameters_3 = NETGEN_1D_2D_3D_2.Parameters()
NETGEN_3D_Parameters_3.SetMaxSize( encap_thickness )
NETGEN_3D_Parameters_3.SetSecondOrder( 0 )
NETGEN_3D_Parameters_3.SetOptimize( 1 )
NETGEN_3D_Parameters_3.SetFineness( 3 )
NETGEN_3D_Parameters_3.SetMinSize( 0.01 )
NETGEN_3D_Parameters_3.SetUseSurfaceCurvature( 1 )
NETGEN_3D_Parameters_3.SetFuseEdges( 1 )
NETGEN_3D_Parameters_3.SetQuadAllowed( 0 )
isDone = Mesh_1.SetMeshOrder( [ [ Sub_mesh_1, Sub_mesh_2, Sub_mesh_3, Sub_mesh_4 ] ])
NETGEN_1D_2D_3D_3 = Mesh_1.Tetrahedron(algo=smeshBuilder.NETGEN_1D2D3D,geom=ROI1)
Sub_mesh_5 = NETGEN_1D_2D_3D_3.GetSubMesh()
NETGEN_3D_Parameters_4 = NETGEN_1D_2D_3D_3.Parameters()
NETGEN_3D_Parameters_4.SetMaxSize( 0.35 )
NETGEN_3D_Parameters_4.SetSecondOrder( 0 )
NETGEN_3D_Parameters_4.SetOptimize( 1 )
NETGEN_3D_Parameters_4.SetFineness( 2 )
NETGEN_3D_Parameters_4.SetMinSize( 0.00111647 )
NETGEN_3D_Parameters_4.SetUseSurfaceCurvature( 1 )
NETGEN_3D_Parameters_4.SetFuseEdges( 1 )
NETGEN_3D_Parameters_4.SetQuadAllowed( 0 )
isDone = Mesh_1.SetMeshOrder( [ [ Sub_mesh_1, Sub_mesh_2, Sub_mesh_3, Sub_mesh_4, Sub_mesh_5 ] ])


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

NETGEN_1D_2D_3D_4 = Mesh_1.Tetrahedron(algo=smeshBuilder.NETGEN_1D2D3D,geom=Rest_1)
Sub_mesh_6 = NETGEN_1D_2D_3D_4.GetSubMesh()
NETGEN_3D_Parameters_5 = NETGEN_1D_2D_3D_4.Parameters()
NETGEN_3D_Parameters_5.SetMaxSize( 2.25 )
NETGEN_3D_Parameters_5.SetSecondOrder( 0 )
NETGEN_3D_Parameters_5.SetOptimize( 1 )
NETGEN_3D_Parameters_5.SetFineness( 1 )
NETGEN_3D_Parameters_5.SetMinSize( 0.1 )
NETGEN_3D_Parameters_5.SetUseSurfaceCurvature( 1 )
NETGEN_3D_Parameters_5.SetFuseEdges( 1 )
NETGEN_3D_Parameters_5.SetQuadAllowed( 0 )
if(Lead2nd_Enable):
	isDone = Mesh_1.SetMeshOrder( [ [ Sub_mesh_1,Sub_mesh_1_2, Sub_mesh_2,Sub_mesh_2_2, Sub_mesh_3,Sub_mesh_3_2, Sub_mesh_4,Sub_mesh_4_2, Sub_mesh_5,Sub_mesh_5_2, Sub_mesh_6 ] ])
else:
	isDone = Mesh_1.SetMeshOrder( [ [ Sub_mesh_1, Sub_mesh_2, Sub_mesh_3, Sub_mesh_4, Sub_mesh_5, Sub_mesh_6 ] ])
isDone = Mesh_1.Compute()




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

measure = smesh.CreateMeasurements()
An_meshed_surf=measure.Area(C1_1)
Cat_meshed_surf=measure.Area(C1_2)


print "Core_div: "
print abs(An_meshed_surf-anode_surf)/anode_surf

print "Outer_div: "
print abs(Cat_meshed_surf-cath_surf)/cath_surf
#TEST!
#while abs(An_meshed_surf-anode_surf)/anode_surf>0.01:
#	contact_mesh_max=contact_mesh_max/2.0
#	NETGEN_2D_Parameters_1.SetMaxSize( anode_mesh_max )
#	isDone = Mesh_1.Compute()

#while abs(Cat_meshed_surf-cath_surf)/cath_surf>0.01:
#	contact_mesh_max=contact_mesh_max/2.0
#	NETGEN_2D_Parameters_2.SetMaxSize( cathode_mesh_max )
#	isDone = Mesh_1.Compute()

## Set names of Mesh objects
#smesh.SetName(NETGEN_1D_2D_3D.GetAlgorithm(), 'NETGEN 1D-2D-3D')
#smesh.SetName(NETGEN_1D_2D.GetAlgorithm(), 'NETGEN 1D-2D')
#smesh.SetName(NETGEN_2D_Parameters_1, 'NETGEN 2D Parameters_1')
#smesh.SetName(NETGEN_2D_Parameters_2, 'NETGEN 2D Parameters_2')
#smesh.SetName(NETGEN_3D_Parameters_1, 'NETGEN 3D Parameters_1')
#smesh.SetName(NETGEN_3D_Parameters_4, 'NETGEN 3D Parameters_4')
#smesh.SetName(NETGEN_3D_Parameters_5, 'NETGEN 3D Parameters_5')
#smesh.SetName(NETGEN_3D_Parameters_2, 'NETGEN 3D Parameters_2')
#smesh.SetName(NETGEN_3D_Parameters_3, 'NETGEN 3D Parameters_3')
#smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')
#smesh.SetName(Sub_mesh_2, 'Sub-mesh_2')
#smesh.SetName(Sub_mesh_1, 'Sub-mesh_1')
#smesh.SetName(Sub_mesh_5, 'Sub-mesh_5')
#smesh.SetName(Sub_mesh_4, 'Sub-mesh_4')
#smesh.SetName(Sub_mesh_3, 'Sub-mesh_3')
#smesh.SetName(Sub_mesh_6, 'Sub-mesh_6')
## Set names of Mesh objects
smesh.SetName(NETGEN_1D_2D_3D.GetAlgorithm(), 'NETGEN 1D-2D-3D')
smesh.SetName(NETGEN_1D_2D.GetAlgorithm(), 'NETGEN 1D-2D')
smesh.SetName(NETGEN_3D_Parameters_2, 'NETGEN 3D Parameters_2')
smesh.SetName(NETGEN_3D_Parameters_3, 'NETGEN 3D Parameters_3')
smesh.SetName(NETGEN_3D_Parameters_1, 'NETGEN 3D Parameters_1')
smesh.SetName(NETGEN_2D_Parameters_1, 'NETGEN 2D Parameters_1')
smesh.SetName(NETGEN_2D_Parameters_2, 'NETGEN 2D Parameters_2')
if(Lead2nd_Enable):
    smesh.SetName(NETGEN_2D_Parameters_1lead2, 'NETGEN_2D_Parameters_1lead2')
    smesh.SetName(NETGEN_2D_Parameters_2lead2, 'NETGEN_2D_Parameters_2lead2')
    smesh.SetName(NETGEN_3D_Parameters_2lead2, 'NETGEN_3D_Parameters_2lead2')
    smesh.SetName(NETGEN_3D_Parameters_3lead2, 'NETGEN_3D_Parameters_3lead2')
    smesh.SetName(NETGEN_3D_Parameters_4lead2, 'NETGEN_3D_Parameters_4lead2')

smesh.SetName(C1_1, 'C1_1')
smesh.SetName(NETGEN_3D_Parameters_4, 'NETGEN 3D Parameters_4')
smesh.SetName(C1_2, 'C1_2')
smesh.SetName(NETGEN_3D_Parameters_5, 'NETGEN 3D Parameters_5')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')
smesh.SetName(Sub_mesh_6, 'Sub-mesh_6')
smesh.SetName(Sub_mesh_5, 'Sub-mesh_5')
smesh.SetName(Rst, 'Rst')
smesh.SetName(RegOfInt, 'RegOfInt')
smesh.SetName(Encap_rest, 'Encap_rest')
smesh.SetName(Encap_contact, 'Encap_contact')
smesh.SetName(Sub_mesh_3, 'Sub-mesh_3')
smesh.SetName(Sub_mesh_2, 'Sub-mesh_2')
smesh.SetName(Sub_mesh_1, 'Sub-mesh_1')
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


Mesh_1.ExportMED('Meshes/Mesh_unref.med')

#if salome.sg.hasDesktop():
#  salome.sg.updateObjBrowser(True)

import killSalome
killSalome.killAllPorts()
