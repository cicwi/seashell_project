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


#%% Initialize:
sino = tw.sinogram()

#%% Misc 1:
#sino.what_to_do()

#%% Misc 2:
#print(sino.meta.history)


#%% Preprocess:
#files_dir = '/export/scratch1/home/sarkissi/Development/Data/ASI_20170308'
#files_dir = '/ufs/sarkissi/Development/Data/20170322-Fly/'

#for energy in [8,11,14,17,20,23]:
for offset in range(0,1):
    
    #files_dir = '/export/scratch2/sarkissi/Data/20170329-SpectralLego'
    #files_dir = '/export/scratch2/sarkissi/Data/20170322-Fly/'
    #files_dir_en = os.path.join(files_dir, str(energy) + 'kV')
    #sino.io.read_raw(files_dir_en)
    #sino.io.read_ref(os.path.join(files_dir_en, 'BH/0um.tif'))
    
    #files_dir = '../Data/seashell/Projections/AMT24-25-SU1/AMT24-23-22/'
    files_dir = '../Data/seashell/Projections/Hinf-AMT24-13-SU2-T1/Hinf-AMT24-13-22-SU2-T1'
    sino.io.read_raw(os.path.join(files_dir, 'Cropped'))
    #sino.io.read_ref(os.path.join(files_dir, '0um-100s-1.tif'))
    sino.io.read_meta(path=os.path.join(files_dir,'log'),kind='SkyScan')
    #sino.meta.theta = -sino.meta.theta
    #sino.process.crop_centered((sino.meta.geometry['optical_axis'], sino.data.shape(2)//2),(50,4000))
    
    #
    # For medipix, we acquire flat field 100 times longer than projections, normalization must be performed. 
    #sino.data._ref[sino.data._ref < 10] = 10
    #sino.data._ref = sino.data._ref / 100.0
    
    #sino.data._data[sino.data._data < 0] = 0
    #sino.data._data /= 1000
    #sino.process.flat_field()
    #sino.data._data += 16
    #sino.process.pixel_calibration(10)
    sino.process.log()
    #sino.io.save_tiff('results/sino','before'+str(energy))
    #sino.process.center_shift(offset=4, test_path='results/center', ind = 65, cen_range= [sino.data.shape(2)//2 - 20, sino.data.shape(2)//2+10, 0.2])
    
    #sino.process.center_shift(offset=offset)
    
    #sino.process.salt_pepper()
    #
    ##%% Reconstruct:
    ''' pixel_size = sino.meta.geometry['det_pixel']
    det2obj = sino.meta.geometry['det2obj']
    src2obj = sino.meta.geometry['src2obj']
    theta = sino.meta.theta    
    sino.reconstruct._initialize_astra(sino.data.shape(), pixel_size, det2obj, src2obj, theta)
    vol = sino.reconstruct._backproject(sino.data._data, algorithm='BP3D_CUDA')
    volume = tw.volume(vol)'''
    #sino.process.short_scan_weights(3.0*numpy.pi/180.0)
    #sino.io.save_tiff(os.path.join('results','Sino'),axis=0)
    volume = sino.reconstruct.FDK()
    #volume = sino.reconstruct.SIRT(iterations=100)
    #
    #%% Visualize:
    #    
    #volume.display.slice(100, 0)
    #sino.io.save_tiff(os.path.join(files_dir,'Cropped'),axis=1)
    
    
    volume.io.save_tiff('results','FDK')
    

