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
    
##files_dir = '/export/scratch2/ci/data/naturalis/seashell/Projections/Hinf-AMT24-13-SU2-T1/'
files_dir = '/export/scratch2/ci/data/naturalis/seashell/Projections/Hinf-AMT24-05-SU1/Hinf-AMT24-05-05LS-SU1/'
sino.io.read_raw(os.path.join(files_dir,''))   

#data = sino.data._data.copy()
#sino.data._data = data.copy()

sino.io.read_meta(os.path.join(files_dir,''),kind='SkyScan')

#%% Crop:
sino.process.crop([900, 1000], [900, 1000])
sino.display.slice(1, dim_num = 1)

# Find the air intensity:
sino.analyse.histogram()

#%% Prepro
sino.process.log(air_intensity = 44000)
sino.display.slice(1, dim_num = 1)
#%% Find the center of rotation:
#sino.io.save_backup()
#sino.process.center_shift()
    
#%% Movie through the sinogram:
sino.display.slice_movie(1, 1)
    
#%% FDK reconstruction:  

volume = sino.reconstruct.FDK()
volume.display.slice_movie(0, 1)    
#%% Single slice reconstruction:
plt.figure()
plt.imshow(sino.reconstruct.slice_FDK()[0,:])

#%% Scan reconstruction:

vol = sino.reconstruct.slice_scan(numpy.linspace(-3, 3, 21))
vol.display.slice_movie(0, 1)
