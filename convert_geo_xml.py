#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Fri May 25 17:19:57 20120

@author: monique
"""

import os
from dolfin import *

#Path
dir_name = "/mnt/c/Users/Monique/Documents/Guelph_Dolomite/"
folder = "Malhas/"
path = dir_name + folder
os.chdir(path)

sample = ["YZ","XZ"]
#Sample number
threshold = ["5"]
sl_1 = ["1","2"]
#sl_1 = 150
#s_tot = 10


for am in sample:
    for s in sl_1:
        for t in threshold:
            #File
            file_in = "gd_" + am + "-474-g" + s + "-c" + t + "-0005_00025" 
            file_out = "gd_" + am + "-474-g" + s + "-c" + t 
                
            #Convert
            os.system("gmsh "+ file_in +".geo -2 -o " + file_out + "_5_00025.msh")
            os.system("dolfin-convert " + file_out + "_5_00025.msh " + file_out + "_5_00025.xml")
            
        