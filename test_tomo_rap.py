#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 18:09:26 2017

@author: kostenko
"""
#%%
import tomo_wrap as tw
import os
import numpy
import matplotlib.pyplot as plt

import dxchange

#%% Initialize:
sino = tw.sinogram()


#%% Misc 1:
#sino.what_to_do()

#%% Misc 2:
#print(sino.meta.history)

#%% Read data:

base_data_dir = '/export/scratch2/sarkissi/Data/seashell/Projections'
rel_files_dir = 'Hinf-AMT24-05-SU1/Hinf-AMT24-05-05LS-SU1'
files_dir = os.path.join(base_data_dir, rel_files_dir)
sino.io.read_raw(os.path.join(files_dir,''))   

sino.io.read_meta(os.path.join(files_dir,'logs'),kind='SkyScan') 

#%% Crop:
sino.process.crop([1000, 1330],[1000, 1330])
sino.display.slice(1, dim_num = 1)

# Find the air intensity:
sino.analyse.histogram()

#%%
#roi = numpy.zeros((sino.data.shape(0),sino.data.shape(2),sino.data.shape(2)), dtype=numpy.float32)
#roi[roi.shape[0]//2, :,:] = 1.0
proj_roi = sino.reconstruct.get_proj_ROI(rows=[5,7],cols=[900,1001])

# Display proj_roi
dxchange.writer.write_tiff_stack(proj_roi, 'results/roi', axis=1, overwrite=True)

#%% Prepro
sino.process.log(air_intensity = 41800)

sino.display.slice(1, dim_num = 1)

#%% FDK reconstruction:
sino.process.center_shift(offset=0)
volume = sino.reconstruct.FDK()
volume.display.slice(5, dim_num =0)
#%%

#%%
#sino.process.center_shift(offset = 25)
sino.display.slice(1, dim_num = 1)
#figure()
plt.imshow(numpy.fliplr(sino.data._data[:, 720, :]))   
#%%
prnt = sino.reconstruct

# Initialize ASTRA:            
sz = sino.data.shape()
pixel_size = sino.meta.geometry['det_pixel']
det2obj = sino.meta.geometry['det2obj']
src2obj = sino.meta.geometry['src2obj']
theta = -sino.meta.theta
        
prnt._initialize_astra(sz, pixel_size, det2obj, src2obj, theta)
        
# Run the reconstruction:
epsilon = numpy.pi / 180.0 # 1 degree
short_scan = numpy.abs(theta[-1] - 2*numpy.pi) > epsilon 
vol = prnt._backproject(sino.data._data, algorithm='FDK_CUDA', short_scan=short_scan)
            
plt.figure()
plt.imshow(vol[5, :,:], cmap = 'gray')