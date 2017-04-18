#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 18:09:26 2017

@author: Henri
"""
#%%
import tomo_wrap as tw
import os
import numpy
import gc

#%% Initialize:
sino = tw.sinogram()



base_data_dir = '/export/scratch2/sarkissi/Data/20170329-SpectralLego'

for energy in [26, 23, 20, 17, 14, 11]:
    rel_files_dir = str(energy) + 'kV'
    files_dir = os.path.join(base_data_dir, rel_files_dir)
    sino.io.read_raw(os.path.join(files_dir,''))
    sino.io.read_ref(os.path.join(files_dir,'BH/0um.tif'))
    
    #sino.io.read_meta(os.path.join(files_dir,'logs'),kind='SkyScan') 
    sino.io.read_meta('')
    
    
    #%% Prepro
    sino.process.pixel_calibration()
    sino.process.flat_field()
    sino.process.log()
    sino.process.center_shift()
    sino.process.salt_pepper()    

    #sino.data._data = numpy.ascontiguousarray(sino.data._data, dtype=numpy.float32)
    volume = sino.reconstruct.FDK()
    #volume.display.slice(300,0)
    #volume.display.slice_movie(0, 1)    
    
    #%% Save data
    volume.io.save_tiff(os.path.join('results', str(energy)),'Lego')
    gc.collect()
    
    