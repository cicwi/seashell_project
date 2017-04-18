#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 18:09:26 2017

@author: kostenko
"""
#%%
import tomo_wrap as tw
import os

import dxchange

#%% Initialize:
sino = tw.sinogram()



base_data_dir = '/export/scratch2/sarkissi/Data/seashell/Projections'
rel_files_dir = 'Hinf-AMT24-05-SU1/Hinf-AMT24-05-05LS-SU1'
files_dir = os.path.join(base_data_dir, rel_files_dir)
sino.io.read_raw(os.path.join(files_dir,''))   

sino.io.read_meta(os.path.join(files_dir,'logs'),kind='SkyScan') 

#%% Crop:
sino.process.crop([0,200],[0,200])
sino.display.slice(1, dim_num = 1)

# Find the air intensity:
#sino.analyse.histogram()

#%%
#roi = numpy.zeros((sino.data.shape(0),sino.data.shape(2),sino.data.shape(2)), dtype=numpy.float32)
#roi[roi.shape[0]//2, :,:] = 1.0
#proj_roi = sino.reconstruct.get_proj_ROI(rows=[5,7],cols=[900,1001])

# Display proj_roi
#dxchange.writer.write_tiff_stack(proj_roi, 'results/roi', axis=1, overwrite=True)

#%% Prepro
sino.process.log(air_intensity = 44000)
sino.display.slice(1, dim_num = 1)


#%% FDK reconstruction:
#sino.process.center_shift(offset=0)

#%% Find the center of rotation:
#sino.io.save_backup()
#sino.process.center_shift()
    
#%% Movie through the sinogram:
#sino.display.slice_movie(1, 1)
    
#%% FDK reconstruction:  

volume = sino.reconstruct.FDK()
volume.display.slice(300,0)
#volume.display.slice_movie(0, 1)


#%% Single slice reconstruction:
plt.figure()
plt.imshow(sino.reconstruct.slice_FDK()[0,:])


#%% Single slice reconstruction:
plt.figure()
plt.imshow(sino.reconstruct.slice_FDK()[0,:])


#%%
#sino.process.center_shift(offset = 25)
sino.display.slice(1, dim_num = 1)
#figure()
plt.imshow(numpy.fliplr(sino.data._data[:, 720, :]))   
#%%
prnt = sino.reconstruct

#%% Scan reconstruction:

vol = sino.reconstruct.slice_scan(numpy.linspace(-3, 3, 21))
vol.display.slice_movie(0, 1)
