#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 18:09:26 2017

@author: kostenko
"""
#%%
import tomo_wrap as tw
import os

#%% Initialize:
sino = tw.sinogram()

#%% Read data:
    
files_dir = '/export/scratch2/ci/data/naturalis/seashell/Projections/Hinf-AMT24-13-SU2-T1/'
sino.io.read_raw(os.path.join(files_dir,'Hinf-AMT24-13-22-SU2-T1'))   

sino.io.read_meta(os.path.join(files_dir,'log'),kind='SkyScan') 

#%% Crop:
sino.process.crop([1000, 1330],[1000, 1330])
sino.display.slice(1, dim_num = 1)

# Find the air intensity:
sino.analyse.histogram()

#%% Prepro
sino.process.log(air_intensity = 41800)

sino.display.slice(1, dim_num = 1)
#%% FDK reconstruction:  

volume = sino.reconstruct.FDK()
volume.display.slice(5, dim_num =0)
#%%

#%%
#sino.process.center_shift(offset = 25)
sino.display.slice(1, dim_num = 1)
figure()
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