#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  6 17:42:26 2020

@author: butenko
"""

# This script computes VTA basing on the impendance (Mädler2012)
# See "Explaining Clinical Effects of Deep Brain Stimulation through Simplified Target-Specific Modeling of the Volume of Activated Tissue" by Mädler and Coenen

import numpy as np

def VTA_from_Z(V_across,Z_model):    #only for a 2-contact system

    print(V_across,Z_model)
    
    #using model 10:
    k1=-1.0473
    k3=0.2786
    k4=0.0009856
    
    r_VTA=-1*(k4*Z_model-np.sqrt((k4**2)*(Z_model**2)+2*k1*k4*Z_model+k1**2+4*k3*V_across)+k1)/(2*k3)
        
    print("VTA sphere radius: ",r_VTA)
    print("VTA in the chapter is given after subtracting the electrode and the encapuslation layer from the VTA sphere")
    
    return True