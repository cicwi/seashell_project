#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 10 15:39:33 2017

@author: kostenko

 ***********    Pilot for the new tomo_wrap  *************
 
"""


#%% Initialization

import matplotlib.pyplot as plt

#%%

from scipy import misc  # Reading BMPs
import os
import numpy
import re
import pkg_resources
#pkg_resources.require("dxchange==0.1.2")
import dxchange
import time
#from mayavi import mlab
#from tvtk.util.ctf import ColorTransferFunction

   
# **************************************************************
#           Parent class for all sinogram subclasses:
# ************************************************************** 
 
class sino_subclass(object):
    def __init__(self, parent):
        self._parent = parent

# **************************************************************
#           IO class and subroutines
# **************************************************************   
from stat import S_ISREG, ST_CTIME, ST_MODE
import os
import gc

def sort_by_date(files):
    '''
    Sort file entries by date
    '''  
    # get all entries in the directory w/ stats
    entries = [(os.stat(path)[ST_CTIME], path) for path in files]
                  
    return [path for date, path in sorted(entries)]
            
def sort_natural(files):
    '''
    Sort file entries using the natural (human) sorting
    '''  
    # Keys
    keys = [int(re.findall('\d+', f)[-1]) for f in files]
        
    # Sort files using keys:        
    files = [f for (k, f) in sorted(zip(keys, files))]    
            
    return files
    
def read_image_stack(file):
    '''
    Read a stack of some image files (not tiff)
    '''
    # Remove the extention and the last few letters:
    name = os.path.basename(file)[:-8]
    path = os.path.dirname(file)

    # Get the files that are alike and sort:   
    files = [os.path.join(path,x) for x in os.listdir(path) if name in x]    
    #files = sorted(files)         
    files = sort_natural(files)    
    
    # Read the first file:         
    image = misc.imread(files[0], flatten= 0)
    sz = numpy.shape(image)

    data = numpy.zeros((len(files), sz[0], sz[1]), 'float32')

    # Read all files:       
    ii = 0
    for filename in files:
        data[ii, :, :] = misc.imread(filename, flatten= 0)
        ii = ii + 1

    print(ii, 'files were loaded.')
    
    return data
    
def update_path(path, io):
    '''
    Memorize the path if it is provided, otherwise use the one that remember.
    '''
    if path == '':
        path = io.path
    else:
        io.path = path
        
    if path == '':
        io._parent.error('Path to the file was not specified.')
        
    return path    
    
def extract_2d_array(dimension, index, data):
    '''
    Extract a 2d array from 3d.
    '''
    if dimension == 0:
        return data[index, :, :]
    elif dimension == 1:
        return data[:, index, :]
    else:
        return data[:, :, index]

class io(sino_subclass):
    '''
    Static class for loading / saving the data
    '''   

    path = ''
    
    #settings = {'sort_by_date':False}
                  
    def read_raw(self, path = ''):
        '''
        Read projection files automatically.
        This will look for files with numbers in the last 4 letters of their names. 
        '''
        path = update_path(path, self)
         
        # Free memory:
        self._parent.data._data = None
        gc.collect()
        
        # Try to find how many files do we need to read:
        
        # Get the files only:
        files = [x for x in os.listdir(path) if os.path.isfile(os.path.join(path, x))]
        
        # Get the last 4 letters:
        index = [os.path.splitext(x)[0][-4:] for x in files]
        
        # Filter out non-numbers:
        index = [int(re.findall('\d+', x)[0]) for x in index if re.findall('\d+', x)]    
                            
        # Extract a number from the first element of the list:
        first = min(index)
        
        # Extract a number from the first element of the list:
        last = max(index)
        
        print('We found projections in the range from ', first, 'to ', last, flush=True)
        
        # Find the file with the maximum index value:
        filename = [x for x in os.listdir(path) if str(last) in x][0]
                    
        # Find the file with the minimum index:
        filename = sorted([x for x in os.listdir(path) if (filename[:-8] in x)&(filename[-3:] in x)])[0]
        
        # If it's a tiff, use dxchange to read tiff:
        if ((filename[-4:] == 'tiff') | (filename[-3:] == 'tif')):
            
            print('Reading a tiff stack')
            if self._parent:
                self._parent.data._data = dxchange.reader.read_tiff_stack(os.path.join(path,filename), range(first, last + 1), digit=5)
            else:
                return dxchange.reader.read_tiff_stack(os.path.join(path,filename), range(first, last + 1))
        else:
            
            print('Reading a stack of images')
            print('Seed file name is:', filename)
            #srt = self.settings['sort_by_date']AMT24-25-SU1/

            if self._parent:
                self._parent.data._data = (read_image_stack(os.path.join(path,filename)))
            else:
                return read_image_stack(os.path.join(path,filename))
         
        # Transpose to satisfy ASTRA dimensions:
        self._parent.data._data = numpy.transpose(self._parent.data._data, (1, 0, 2)) 
        #self._parent.data._data = numpy.ascontiguousarray(self._parent.data._data, dtype=numpy.float32)
        
        # Cast to float to avoid problems with divisions in the future:
       # self._parent.data._data = numpy.float32(self._parent.data._data, copy=False)
            
        # add record to the history:        
        self._parent.meta.history['io.read_raw'] = path
        
    def read_ref(self, path_file):
        '''
        Read reference flat field.
            
        '''
        if ((path_file[-4:] == 'tiff') | (path_file[-3:] == 'tif')):
            ref = numpy.array(dxchange.reader.read_tiff(path_file))
        else:
            ref  = misc.imread(path_file, flatten= 0)
            
        if self._parent:
            self._parent.data._ref = ref
         
        # Cast to float to avoid problems with divisions in the future:
        self._parent.data._ref = numpy.float32(self._parent.data._ref)
            
        # add record to the history:       
        self._parent.meta.history['io.read_ref'] = path_file

        self._parent.message('Flat field reference image loaded.')
        
    def read_bh(self, path_file):
        '''
        Read reference foil data for signal to equivalent thickness calibration.
            
        ''' 
        if ((path_file[-4:] == 'tiff') | (path_file[-3:] == 'tif')):
            ref = numpy.array(dxchange.reader.read_tiff(path_file))
        else:
            ref  = misc.imread(path_file, flatten= 0)
            
        if self._parent:
            self._parent.data._ref = ref
         
        # Cast to float to avoid problems with divisions in the future:
        self._parent.data._ref = numpy.float32(self._parent.data._ref)    
            
        # add record to the history:        
        self._parent.meta.history['io.read_ref'] = path_file

        self._parent.message('Beam hardening correction reference images loaded.')       
        
    def read_meta(self, path = '', kind=None):
        '''
        Parser for the metadata file that contains information about the acquisition system.
        '''
        path = update_path(path, self)
        geometry = self._parent.meta.geometry 
        
        # TODO: make the actual parser. For now just initialization with default
        #self._parent.meta.geometry['det_pixel'] = 28.0/self._parent.data.shape(0)
        if self._parent.data._data != []:
            geometry['nb_angle'] = self._parent.data.shape(1)
        else:
            self._parent.error('Load the data first. We don`t know how many angles to use')
            
        if kind is None:
            geometry['det_pixel'] = 0.055
            geometry['src2obj'] = 210.0
            geometry['det2obj'] = 50.0
            #geometry['src2det'] = 209.610
            geometry['rot_step'] = 2*numpy.pi / (geometry['nb_angle']-1)
        
        elif (str.lower(kind) == 'skyscan'):
            # Parse the SkyScan log file 
            self._parse_skyscan_meta(path)
        
        print(geometry)
        if (geometry['det2obj'] == 0.0):
            geometry['det2obj'] = geometry['src2det'] - geometry['src2obj']
        
        
        self._parent.meta.theta = numpy.linspace(0, (geometry['nb_angle'] - 1)*geometry['rot_step'], geometry['nb_angle'])
        
        # add record to the history: 
        self._parent.meta.history['io.read_meta'] = path

        self._parent.message('Meta data loaded.')
     
    def save_backup(self):
        '''
        Make a copy of data in memory, just in case.
        '''
        self._parent.data._backup = self._parent.data._data.copy()
        
        # add record to the history: 
        self._parent.meta.history['io.save_backup'] = 'backup saved'

        self._parent.message('Backup saved.')
    
    def load_backup(self):
        '''
        Retrieve a copy of data from the backup.
        '''
        if self._parent.data._backup == []:
            self._parent.error('No backup found in memory!')
            
        self._parent.data._data = self._parent.data._backup.copy()
        
        # Clean memory:
        self._parent.data._backup = None
        gc.collect()
        
        # add record to the history: 
        self._parent.meta.history['io.load_backup'] = 'backup loaded'

        self._parent.message('Backup loaded.')
        
    # **************************************************************
    # Parsers for metadata files
    # **************************************************************
    def _parse_unit(self, string):
            # Look at the inside of trailing parenthesis
            unit = ''
            factor = 1.0
            if string[-1] == ')':
                unit = string[string.rfind('(')+1:-1].lower()
            else:
                return 1.0
                
            # Metric units --> mm
            if unit == 'mm':
                pass 
            elif unit == 'um':
                factor = 0.001
            elif unit == 'cm':
                factor = 10.0
            elif unit == 'm':
                factor = 1000.0
                
            # Angles --> rad
            elif unit == 'rad':
                pass
            elif unit == 'deg':
                factor = numpy.pi / 180.0
                
            # Time --> ms
            elif unit == 'ms':
                pass
            elif unit == 's':
                factor = 1000.0
            elif unit == 'us':
                factor = 0.001
             
            # Energy --> keV
            elif unit == 'kev':
                pass
            elif unit == 'mev':
                factor = 1000.0
            elif unit == 'ev':
                factor = 0.001
                
            # Voltage --> kV
            elif unit == 'kv':
                pass
            elif unit == 'mv':
                factor = 1000.0
            elif unit == 'v':
                factor = 0.001
            
            # Current --> uA
            elif unit == 'ua':
                pass
            elif unit == 'ma':
                factor = 1000.0
            elif unit == 'a':
                factor = 1000000.0

            else:
                self._parent.warning('Unknown unit: ' + unit + '. Skipping.')

            return factor
    
    def _parse_skyscan_meta(self, path = ''):
                
        path = update_path(path,self)
        
        # Try to find the log file in the selected path
        log_file = [x for x in os.listdir(path) if (os.path.isfile(os.path.join(path, x)) and os.path.splitext(os.path.join(path, x))[1] == '.log')]

        if len(log_file) == 0:
            raise FileNotFoundError('Log file not found in path: ' + path)
        if len(log_file) > 1:
            #raise UserWarning('Found several log files. Currently using: ' + log_file[0])
            self._parent.warning('Found several log files. Currently using: ' + log_file[0])
            log_file = os.path.join(path, log_file[0])
        else:
            log_file = os.path.join(path, log_file[0])
            
        #Once the file is found, parse it
        geometry = self._parent.meta.geometry
        physics = self._parent.meta.physics
        
        # Create a dictionary of keywords (skyscan -> our geometry definition):
        geom_dict = {'camera pixel size': 'det_pixel', 'image pixel size': 'img_pixel', 'object to source':'src2obj', 'camera to source':'src2det', 
        'optical axis':'optical_axis', 'rotation step':'rot_step', 'exposure':'exposure', 'source voltage':'voltage', 'source current':'current'}    
        
        with open(log_file, 'r') as logfile:
            for line in logfile:
                name, var = line.partition("=")[::2]
                name = name.strip().lower()
                                
                # If name contains one of the keys (names can contain other stuff like units):
                geom_key = [geom_dict[key] for key in geom_dict.keys() if key in name]

                if geom_key != []:
                    factor = self._parse_unit(name)
                    geometry[geom_key[0]] = float(var)*factor
    
    def save_tiff(self, path = '', fname='data', axis=0):
        '''
        Saves the data to tiff files
        '''
        
        path = update_path(path, self)
        dxchange.writer.write_tiff_stack(self._parent.data.get_data(),fname=os.path.join(path, fname), axis=axis,overwrite=True)
        
# **************************************************************
#           META class and subclasses
# **************************************************************
        
class meta(sino_subclass):
    '''
    This object contains various properties of the imaging system and the history of pre-processing.
    '''         
    geometry = {'src2obj': 0, 'det2obj': 0, 'det_pixel': 0, 'det_offset': 0, 'det_tilt': 0}
    theta =  numpy.linspace(0, 2*numpy.pi, 540)
    
    physics = {'voltage': 0, 'current':0, 'exposure': 0}
    history = {'object initialized': time.time()} 
        
# **************************************************************
#           DISPLAY class and subclasses
# **************************************************************
class display(sino_subclass):   
    '''
    This is a collection of display tools for the raw and reconstructed data
    '''  
    def __init__(self, parent = []):
        self._parent = parent
        self._cmap = 'gray'
    
    def _figure_maker_(self, fig_num):
        '''
        Make a new figure or use old one.
        '''
        if fig_num:
            plt.figure(fig_num)
        else:
            plt.figure()
            
    def slice(self, slice_num, dim_num, fig_num = []):
        '''
        Display a 2D slice of 3D volumel
        '''
        self._figure_maker_(fig_num)
         
        img = extract_2d_array(dim_num, slice_num, self._parent.data.get_data())            
        plt.imshow(img, cmap = self._cmap)
            
    def slice_movie(self, dim_num, fig_num = []):
        '''
        Display a 2D slice of 3D volumel
        '''
        self._figure_maker_(fig_num)
        
        slice_num = 0
        img = extract_2d_array(dim_num, slice_num, self._parent.data.get_data())            
        fig = plt.imshow(img, cmap = self._cmap)
        
        for slice_num in range(1, self._parent.data.shape()[dim_num]):
            img = extract_2d_array(dim_num, slice_num, self._parent.data.get_data())            
            fig.set_data(img)
            plt.show()
            plt.title(slice_num)
            plt.pause(0.0001)
    
    def max_projection(self, dim_num, fig_num = []):
        '''
        Get maximum projection image of the 3d data.
        '''
        self._figure_maker_(fig_num)
        
        img = self._parent.data._data.max(dim_num)
        plt.imshow(img, cmap = self._cmap)
        
    def min_projection(self, dim_num, fig_num = []):
        '''
        Get maximum projection image of the 3d data.
        '''
        self._figure_maker_(fig_num)
        
        img = self._parent.data._data.max(dim_num)
        plt.imshow(img, cmap = self._cmap)
        
    def render(self):
        '''
        Render volume using mayavi routines
        '''
        data = self._parent.data.get_data()
        
        vol = mlab.pipeline.volume(mlab.pipeline.scalar_field(numpy.fliplr(data)), vmin = 0.001, vmax = 0.01)
        mlab.colorbar()

        # Adjust colors:
        ctf = ColorTransferFunction()
            
        for ii in numpy.linspace(0, 1, 10):
            ctf.add_hsv_point(ii * 0.01, 0.99 - ii, 1, 1)
                
        ctf.range= [0, 1]
        vol._volume_property.set_color(ctf)
        vol._ctf = ctf
        vol.update_ctf = True
        
# **************************************************************
#           ANALYSE class and subclasses
# **************************************************************
class analyse(sino_subclass):   
    '''
    This is an anlysis toolbox for the raw and reconstructed data
    '''
    def __init__(self, parent = []):
        self._parent = parent
                
    def mean(self):
        return numpy.mean(self._parent.data.get_data())    
        
    def min(self):
        return numpy.min(self._parent.data.get_data())        
        
    def max(self):
        return numpy.max(self._parent.data.get_data())    
        
    def histogram(self, nbin = 256, plot = True):
        
        mi = self.min()
        ma = self.max()
        
        a, b = numpy.histogram(self._parent.data.get_data(), bins = nbin, range = [mi, ma]) 
        
        # Set bin values to the middle of the bin:
        b = (b[0:-1] + b[1:]) / 2

        if plot:
            plt.figure()
            plt.plot(b, a)

        return a, b   
        
# **************************************************************
#           PROCESS class and subclasses
# **************************************************************
#def interpolate_in_mask(image, mask, kernel):
#    '''
#    Replace masked pixels with interpolated values
#    '''
from scipy import ndimage
from tomopy.recon import rotation
import scipy.ndimage.interpolation as interp

class process(sino_subclass):   
    '''
    Various preprocessing routines
    '''
    def __init__(self, parent = []):
        self._parent = parent
        
    def arbitrary_function(self, func):
        '''
        Apply an arbitrary function:
        '''
        print(func)
        self._parent.data._data = func(self._parent.data._data)

        # add a record to the history: 
        self._parent.meta.history['process.arbitrary_function'] = func.__name__

        self._parent.message('Arbitrary function applied.')        
    
    def pixel_calibration(self, kernel=5):
        '''
        Apply correction to miscalibrated pixels.
        '''
        # Compute mean image of intensity variations that are < 5x5 pixels
        res = self._parent.data._data - ndimage.filters.median_filter(self._parent.data._data, [kernel, 1, kernel])
        res = res.mean(1)    
        self._parent.data._data -= res.reshape((res.shape[0], 1, res.shape[1]))
        self._parent.meta.history['Pixel calibration'] = 1
        self._parent.message('Pixel calibration correction applied.')
        
    def medipix_quadrant_shift(self):
        self._parent.data._data[:,:, 0:self._parent.data.shape(2)//2 - 2] = self._parent.data._data[:,:, 2:self._parent.data.shape(2)//2]
        self._parent.data._data[:,:, self._parent.data.shape(2)//2 + 2:] = self._parent.data._data[:,:, self._parent.data.shape(2)//2:-2]
        
        # Fill in two extra pixels:
        for ii in range(-2,2):
            closest_offset = -3 if (numpy.abs(-3-ii) < numpy.abs(2-ii)) else 2
            self._parent.data._data[:,:, self._parent.data.shape(2)//2 - ii] = self._parent.data._data[:,:, self._parent.data.shape(2)//2 + closest_offset]
        
        
        # Then in columns
        self._parent.data._data[0:self._parent.data.shape(0)//2 - 2,:,:] = self._parent.data._data[2:self.data.shape(0)//2,:,:]
        self._parent.data._data[self.data.shape(0)//2 + 2:, :, :] = self._parent.data._data[self.data.shape(0)//2:-2,:,:]
        
        # Fill in two extra pixels:
        for jj in range(-2,2):
            closest_offset = -3 if (numpy.abs(-3-jj) < numpy.abs(2-jj)) else 2
            self._parent.data._data[self.data.shape(0)//2 - jj,:,:] = self._parent.data._data[self.data.shape(0)//2 + closest_offset,:,:]
        
        self._parent.meta.history['Quadrant shift'] = 1
        self._parent.message('Medipix quadrant shift applied.')

    def flat_field(self):
        '''
        Apply flat field correction.
        '''
        
        if numpy.min(self._parent.data._ref) <= 0:
            self._parent.warning('Flat field reference image contains zero (or negative) values! Will replace those with little tiny numbers.')
            
            tiny = self._parent.data._ref[self._parent.data._ref > 0].min()
            self._parent.data._ref[self._parent.data._ref <= 0] = tiny

        # How many projections:    
        n_proj = self._parent.data.shape(1)

        for ii in range(0, n_proj):
            self._parent.data._data[:, ii, :] = self._parent.data._data[:, ii, :] / self._parent.data._ref                        
        
        # add a record to the history: 
        self._parent.meta.history['process.flat_field'] = 1

        self._parent.message('Flat field correction applied.')
        
    
    def short_scan_weights(self, fan_angle):
        '''
        Apply parker weights correction.
        '''
        def _Parker_window(theta, gamma, fan):
            weight = 0.0
            if (0 <= theta < 2*(gamma+fan)):
                weight = numpy.sin((numpy.pi/4)*(theta/(gamma+fan)))**2
            elif (2*(gamma+fan) <= theta < numpy.pi + 2*gamma):
                weight = 1.0
            elif (numpy.pi + 2*gamma <= theta < numpy.pi + 2*fan):
                weight = numpy.sin((numpy.pi/4)*((numpy.pi + 2*fan - theta)/(gamma+fan)))**2 
            else:
                weight = 0.0
            return weight
        
        weights = numpy.zeros_like(self._parent.data._data, dtype=numpy.float32)
        sdd = self._parent.meta.geometry['src2det']
        for u in range(0,weights.shape[2]):
            weights[:,:,u] = u
            
        weights = weights - weights.shape[2]/2
        weights = self._parent.meta.geometry['det_pixel']*weights
        weights = numpy.arctan(weights/sdd)
        
        theta = self._parent.meta.theta
        for ang in range(0,theta.shape[0]):
            tet = theta[ang]
            for u in range(0, weights.shape[2]):
                weights[:,ang,u] = _Parker_window(theta = tet, gamma = weights[0,ang,u], fan=fan_angle)
        
        self._parent.data._data *= weights
        # add a record to the history: 
        self._parent.meta.history['process.short_scan'] = 1

        self._parent.message('Short scan correction applied.')
        

    def log(self, air_intensity = 1, upper_bound = numpy.log(256)):
        '''
        Apply -log(x) to the sinogram
        '''
        self._parent.data._data = -numpy.log(self._parent.data._data / air_intensity)
        
        # Apply a bound to large values:
        #self._parent.data._data[self._parent.data._data > upper_bound] = upper_bound
        #self._parent.data._data[~numpy.isfinite(self._parent.data._data)] = upper_bound

        self._parent.data._data = numpy.nan_to_num(self._parent.data._data)
        numpy.clip(self._parent.data._data, a_min = -10, a_max = upper_bound, out = self._parent.data._data)
                 
        self._parent.meta.history['process.log(upper_bound)'] = upper_bound                  

    def salt_pepper(self, kernel = 3):
        '''
        Gets rid of nasty speakles
        '''       
        # Make a smooth version of the data and look for outlayers:
        smooth = ndimage.filters.median_filter(self._parent.data._data, [kernel, 1, kernel])
        mask = self._parent.data._data / smooth
        mask = (numpy.abs(mask) > 1.5) | (numpy.abs(mask) < 0.75)

        self._parent.data._data[mask] = smooth[mask]
        
        
        self._parent.meta.history['process.salt_pepper(kernel)'] = kernel
                                  
    def center_shift(self, offset=None, test_path=None, ind=None, cen_range=None):
        '''
        Find the center of the sinogram and apply the shift to corect for the offset
        '''
        sz = self._parent.data.shape(2) // 2
        
        if test_path is not None:
            rotation.write_center(self._parent.data._data, theta=self._parent.meta.theta, dpath=test_path, ind=ind, cen_range=cen_range, sinogram_order=True)
            
        else:
            if offset is None:
                offset = sz-rotation.find_center(self._parent.data._data, self._parent.meta.theta,  sinogram_order=True)[0]
                #offset = sz-rotation.find_center_pc(self._parent.data._data[:,0,:], self._parent.data._data[:,self._parent.data.shape(1)//2,:])p
                
            # Do nothing is the offset is less than 1 pixel
            if (numpy.abs(offset) >= 1):
                self._parent.data._data = interp.shift(self._parent.data._data, (0,0,offset))                            
                self._parent.meta.history['process.center_shift(offset)'] = offset
                              
                self._parent.message('Horizontal offset corrected.')
            else:
                self._parent.warning("Center shift found an offset smaller than 1. Correction won't be applied")
                
                
    def bin_theta(self):
        '''
        Bin angles with a factor of two
        '''
        self._parent.data._data = (self._parent.data._data[:,0:-1:2,:] + self._parent.data._data[:,1::2,:]) / 2
        self._parent.meta.theta = (self._parent.meta.theta[:,0:-1:2,:] + self._parent.meta.theta[:,1::2,:]) / 2                           

    def crop(self, top_left, bottom_right):
        '''
        Crop the sinogram
        '''
        if bottom_right[1] > 0:
            self._parent.data._data = self._parent.data._data[top_left[1]:-bottom_right[1], :, :]
        else:
            self._parent.data._data = self._parent.data._data[top_left[1]:, :, :]                

        if bottom_right[0] > 0:
            self._parent.data._data = self._parent.data._data[:, :, top_left[0]:-bottom_right[0]]
        else:
            self._parent.data._data = self._parent.data._data[:, :, top_left[0]:]

        self._parent.data._data = numpy.ascontiguousarray(self._parent.data._data, dtype=numpy.float32)
        gc.collect()
         
        self._parent.meta.history['process.ccrop(top_left, bottom_right)'] = [top_left, bottom_right]
                                  
        self._parent.message('Sinogram cropped.')
        
    def crop_centered(self, center, dimensions):
        '''
        Crop the sinogram
        '''
        self._parent.data._data = self._parent.data._data[center[0] - dimensions[0]//2:center[0] + dimensions[0]//2, :, center[1] - dimensions[1]//2:center[1] + dimensions[1]//2]
        self._parent.data._data = numpy.ascontiguousarray(self._parent.data._data, dtype=numpy.float32)
        gc.collect()
        
        self._parent.meta.history['process.crop_centered(center, dimensions)'] = [center, dimensions]
                                  
        self._parent.message('Sinogram cropped.')
                                  
# **************************************************************
#           RECONSTRUCTION class and subclasses
# **************************************************************
import astra

class reconstruct(sino_subclass):   
    '''
    Reconstruction algorithms: FDK, SIRT, KL, FISTA etc.
    '''
    # Astra-variables:
    proj_id = []
    rec_id = []
    sinogram_id = []
    
    def __init__(self, parent = []):
        self._parent = parent
    
    def slice_FDK(self, parameter_value = 0, parameter = 'axis_offset'):
        '''
        A quick calculation of a single central slice.
        Returns a numpy array and not a volume object!
        '''
        prnt = self._parent
        
        # Extract 1 pixel thin slice:
        sinogram = prnt.data._data[prnt.data.shape(0)//2, :, :]
        
        # For compatibility purposes make sure that the result is 3D:
        sinogram = numpy.ascontiguousarray(sinogram[None, :])    
        
        # Initialize ASTRA:            
        sz = numpy.array(prnt.data.shape())
        pixel_size = prnt.meta.geometry['img_pixel']
        det2obj = prnt.meta.geometry['det2obj']
        src2obj = prnt.meta.geometry['src2obj']
        theta = prnt.meta.theta
        
        # Temporary change of one of the parameters:      
        if abs(parameter_value) >0:
            
            if parameter == 'axis_offset':
                # Apply shift:
                sinogram = interp.shift(sinogram, (0,0, parameter_value))  
            elif parameter == 'det_pixel':
                pixel_size = parameter
                
            elif parameter == 'det2obj':    
                det2obj = parameter
                
            elif parameter == 'src2obj':   
                src2obj = parameter
                
            else: prnt.error("Can't recognize given parameter.")     
                
        sz[0] = 1
        self._initialize_astra(sz, pixel_size, det2obj, src2obj, theta)
        
        # Run the reconstruction:
        #epsilon = numpy.pi / 180.0 # 1 degree
        
        short_scan = (theta.max() - theta.min()) < (numpy.pi * 1.99)
        vol = self._backproject(sinogram, algorithm='FDK_CUDA', short_scan=short_scan)
            
        return vol
        # No need to make a history record - sinogram is not changed.
        
    def slice_scan(self, scan_range = numpy.linspace(-10, 10, 11), parameter = 'axis_offset'):
        '''
        Create a scan of different rotation axis offsets:
        '''
        sz = self._parent.data.shape()
        
        print('Starting an', parameter, ' scan.')
        
        # Create a volume to put a scan into:
        vol = numpy.zeros([scan_range.shape[0], sz[2], sz[2]])
        
        for ii in scan_range:
            print(ii)
            
            img = self.slice_FDK(ii, parameter)
            vol[ii, :, :] = img

        return volume(vol)
    
    def FDK(self):
        
        prnt = self._parent
        
        # Initialize ASTRA:            
        sz = prnt.data.shape()
        pixel_size = prnt.meta.geometry['det_pixel']
        det2obj = prnt.meta.geometry['det2obj']
        src2obj = prnt.meta.geometry['src2obj']
        theta = prnt.meta.theta
        
        self._initialize_astra(sz, pixel_size, det2obj, src2obj, theta)
        
        # Run the reconstruction:
        #epsilon = numpy.pi / 180.0 # 1 degree - I deleted a part of code here by accident...
        short_scan = (theta.max() - theta.min()) < (numpy.pi * 1.99)
        vol = self._backproject(prnt.data._data, algorithm='FDK_CUDA', short_scan=short_scan)
            
        return volume(vol)
        # No need to make a history record - sinogram is not changed.
        
    def SIRT(self, iterations = 10, min_constraint = None):
        
        prnt = self._parent
        
        # Initialize ASTRA:            
        sz = prnt.data.shape()
        pixel_size = prnt.meta.geometry['det_pixel']
        det2obj = prnt.meta.geometry['det2obj']
        src2obj = prnt.meta.geometry['src2obj']
        theta = prnt.meta.theta
        
        self._initialize_astra(sz, pixel_size, det2obj, src2obj, theta)
        
        # Run the reconstruction:
        vol = self._backproject(prnt.data._data, algorithm = 'SIRT3D_CUDA', iterations = iterations, min_constraint= min_constraint)

        return volume(vol)
        # No need to make a history record - sinogram is not changed.
        
    def SIRT_CPU(self, iterations = 10, min_constraint = None):
        
        prnt = self._parent
        
        # Initialize ASTRA:            
        sz = prnt.data.shape()
        pixel_size = prnt.meta.geometry['det_pixel']
        det2obj = prnt.meta.geometry['det2obj']
        src2obj = prnt.meta.geometry['src2obj']
        theta = prnt.meta.theta
        
        self._initialize_astra(sz, pixel_size, det2obj, src2obj, theta)
        
        # Create a volume containing only ones for forward projection weights
        vol_ones = numpy.ones((sz[0], sz[2], sz[2]), dtype=numpy.float32)
        vol = numpy.zeros_like(vol_ones, dtype=numpy.float32)
        fwd_weights = self._forwardproject(vol_ones)
        fwd_weights = 1.0 / (fwd_weights + (fwd_weights == 0))
        
        bwd_weights = 1.0 / (theta.shape[0])
        
        for ii_iter in range(iterations):
            fwd_proj_vols = self._forwardproject(vol)
            
            residual = prnt.data._data - fwd_proj_vols
            residual *= fwd_weights
            
            vol += bwd_weights * self._backproject(residual, algorithm='BP3D_CUDA')

        
        return volume(vol)
        # No need to make a history record - sinogram is not changed.
        
    def CPLS(self, iterations = 10, min_constraint = None):
        
        prnt = self._parent
        
        # Initialize ASTRA:            
        sz = prnt.data.shape()
        pixel_size = prnt.meta.geometry['det_pixel']
        det2obj = prnt.meta.geometry['det2obj']
        src2obj = prnt.meta.geometry['src2obj']
        theta = prnt.meta.thetaimg
        
        self._initialize_astra(sz, pixel_size, det2obj, src2obj, theta)
        
        # Create a volume containing only ones for forward projection weights
        vol_ones = numpy.ones((sz[0], sz[2], sz[2]), dtype=numpy.float32)
        vol = numpy.zeros_like(vol_ones, dtype=numpy.float32)
        
        sigma = self._forwardproject(vol_ones)
        sigma = 1.0 / (sigma + (sigma == 0))
        sigma_1 = 1.0  / (1.0 + sigma)
        tau = 1.0 / theta.shape[0]
        
        p = numpy.zeros_like(prnt.data._data)
        ehn_sol = vol.copy()

        for ii_iter in range(iterations):
            residual = prnt.data._data - self._forwardproject(ehn_sol)
            p = (p + residual * sigma) * sigma_1

            updates = self._backproject(p, algorithm='BP3D_CUDA')
            old_vol = vol.copy()
            vol += updates * tau
            vol *= (vol > 0)
            
            ehn_sol = vol + (vol - old_vol)

        
        return volume(vol)
        # No need to make a history record - sinogram is not changed.
        
        
        
    def CGLS(self, iterations = 10, min_constraint = None):

        prnt = self._parent
        
        # Initialize ASTRA:            
        sz = prnt.data.shape()
        pixel_size = prnt.meta.geometry['det_pixel']
        det2obj = prnt.meta.geometry['det2obj']
        src2obj = prnt.meta.geometry['src2obj']
        theta = prnt.meta.theta
        
        self._initialize_astra(sz, pixel_size, det2obj, src2obj, theta)
        
        # Run the reconstruction:
        vol = self._backproject(prnt.data._data, algorithm = 'CGLS3D_CUDA', iterations = iterations, min_constraint=min_constraint)
            
        
        return volume(vol)
        # No need to make a history record - sinogram is not changed.
        
        
    def _initialize_astra(self, sz, det_pixel_size, det2obj, src2obj, theta):
        
        # Initialize ASTRA (3D): 
        det_count_x = sz[2]
        det_count_z = sz[0]
        
        # Make volume count x > detector count to include corneres of the object:
        vol_count_x = sz[2]
        vol_count_z = sz[0]
                
        tot_dist = det2obj + src2obj
        
        magnification = tot_dist / src2obj
              
        self.vol_geom = astra.create_vol_geom(vol_count_x, vol_count_x, vol_count_z)
        self.proj_geom = astra.create_proj_geom('cone', magnification, magnification, det_count_z, det_count_x, 
                                           theta, (src2obj*magnification)/det_pixel_size, (det2obj*magnification)/det_pixel_size)

     
    def _backproject(self, y, algorithm = 'FDK_CUDA', iterations=1, min_constraint = None, short_scan=False):
        
      cfg = astra.astra_dict(algorithm)
      cfg['option'] = {'ShortScan' : short_scan}
      if (min_constraint is not None):
          cfg['option']['MinConstraint'] = min_constraint
 
      output = numpy.zeros(astra.functions.geom_size(self.vol_geom), dtype=numpy.float32)
      rec_id = astra.data3d.link('-vol', self.vol_geom, output)
      sinogram_id = astra.data3d.link('-sino', self.proj_geom, y)
          
      cfg['ReconstructionDataId'] = rec_id
      cfg['ProjectionDataId'] = sinogram_id
      #print(cfg['option'])
      alg_id = astra.algorithm.create(cfg)

      #astra.data3d.store(self.sinogram_id, y)
      astra.algorithm.run(alg_id, iterations)
      
      astra.algorithm.delete(alg_id)
      
      astra.data3d.delete([rec_id, sinogram_id])
      
      
      return output #astra.data3d.get(self.rec_id)
    
    def _forwardproject(self, x, algorithm = 'FP3D_CUDA'):
        
      cfg = astra.astra_dict(algorithm)
      
      output = numpy.zeros(astra.functions.geom_size(self.proj_geom), dtype=numpy.float32)
      rec_id = astra.data3d.link('-vol', self.vol_geom, x)
      sinogram_id = astra.data3d.link('-sino', self.proj_geom, output)
      
      cfg['VolumeDataId'] = rec_id
      cfg['ProjectionDataId'] = sinogram_id
        
      alg_id = astra.algorithm.create(cfg)
        
      #astra.data3d.store(self.rec_id, x)
      astra.algorithm.run(alg_id, 1)
      
      astra.data3d.delete([rec_id, sinogram_id])
      astra.algorithm.delete(alg_id)
      
      return output


    def get_proj_ROI(self, rows=[0,512], cols=[0,512], algorithm='FP3D_CUDA'):
        # Computes a mask of minimal projection ROI needed to reconstruct a ROI for FDK
        prnt = self._parent
        
        # Initialize ASTRA:            
        sz = prnt.data.shape()
        pixel_size = prnt.meta.geometry['det_pixel']
        det2obj = prnt.meta.geometry['det2obj']
        src2obj = prnt.meta.geometry['src2obj']
        theta = prnt.meta.theta
        
        
        roi = numpy.zeros((sz[0],sz[2], sz[2]), dtype=numpy.float32)
        roi[rows[0]:rows[1],cols[0]:cols[1],cols[0]:cols[1]] = 1.0
        self._initialize_astra(sz, pixel_size, det2obj, src2obj, theta)
        
        
        mask = self._forwardproject(roi, algorithm=algorithm)
        
        # TODO: Compute the bounds of the minimal non-zero rectangle
        '''
        mask[mask>0]=1.0
        bounds = [[0,0],[0,0]]
        bounds[0][0] = numpy.min(numpy.argmax(numpy.argmax(mask,axis=2),axis=1))
        for row in range(mask.shape[0],-1,-1))
        bounds[0][1] = numpy.argmin(mask,axis=0)
        bounds[1][0] = numpy.argmax(mask,axis=2)
        bounds[1][1] = numpy.argmin(mask,axis=2)
        
        print(bounds)
        '''
        return mask
        
        

# **************************************************************
#           DATA class and subclasses
# **************************************************************
class data(sino_subclass):   
    '''
    Memory allocation, reading and writing the data
    '''
    _data = []  
    _ref = []
    _backup = []

    _isgpu = False

    def __init__(self, parent = []):
        self._parent = parent
        
    def get_data(self):
        '''
        Get sinogram data. Copies data from GPU if needed
        '''
        return self._data
    
    def set_data(self, sino):
        '''
        Set sinogram data. Copies data to GPU if needed
        '''
        self._data = data      
        
        self._parent.meta.history['data.set_data'] = 1

    def shape(self, dim = None):
        if dim is None:
            return self._data.shape
        else:
            return self._data.shape[dim]
            
# **************************************************************
#           VOLUME class
# ************************************************************** 
class volume(object):
    data = []
    io = []
    analyse = []
    display = []

    def __init__(self, vol):
        self.io = io(self) 
        self.display = display(self)
        self.analyse = analyse(self)
        self.data = data(self)    
        
        # Get the data in:
        self.data._data = vol

# **************************************************************
#           SINOGRAM class
# **************************************************************        
import warnings
import logging
import copy

_min_history = ['io.read_raw', 'io.read_ref', 'io.read_meta', 'process.flat_field', 'process.log']

_wisdoms = ['You’d better get busy, though, buddy. The goddam *sands* run out on you \
every time you turn around. I know what I’m talking about. You’re lucky if \
you get time to sneeze in this goddam phenomenal world.', 
'Work done with anxiety about results is far inferior to work done without\
such anxiety, in the calm of self-surrender. Seek refuge in the knowledge\
of Brahman. They who work selfishly for results are miserable.',
'You have the right to work, but for the work`s sake only. You have no right\
to the fruits of work. Desire for the fruits of work must never be your\
motive in working. Never give way to laziness, either.',
'Perform every action with your heart fixed on the Supreme Lord. Renounce\
attachment to the fruits. Be even-tempered [underlined by one of the \
cal-ligraphers] in success and failure; for it is this evenness of temper which is meant by yoga.', 
'God instructs the heart, not by ideas but by pains and contradictions.',
'Sir, we ought to teach the people that they are doing wrong in worshipping\
the images and pictures in the temple.']

class sinogram(object):
    '''
    
    Class that will contain the raw data and links to all operations that we need
    to process and reconstruct it.
    
    '''   
    # Public stuff:
    io = []
    meta = []
    display = []
    analyse = []
    process = []
    reconstruction = []
    data = []  

    def __init__(self):
        self.io = io(self) 
        self.meta = meta(self)
        self.display = display(self)
        self.analyse = analyse(self)
        self.process = process(self)
        self.reconstruct = reconstruct(self)
        self.data = data(self)
            
    def message(self, msg):
        '''
        Send a message to IPython console.
        '''
        log = logging.getLogger()
        log.setLevel(logging.DEBUG)
        log.debug(msg)
        
    def error(self, msg):
        '''
        Throw an error:
        '''
        self.meta.history['error'] = msg
        raise ValueError(msg)
        
    def warning(self, msg):
        '''
        Throw a warning. In their face!
        '''
        self.meta.history['warning'] = msg
        warnings.warn(msg)
        
    def what_to_do(self):
        
        finished = True
        
        for k in _min_history:  
            print((k in meta.history.keys()))
            if ~(k in meta.history.keys()):
                print('You should use ', k, 'as a next step')
                finished = False
                break
                
        if finished:        
            print('All basic processing steps were done. Use "reconstruct.FDK" to compute filtered backprojection.')
                
    def copy(self):
        '''
        Deep copy of the sinogram object:
        '''
        return copy.deepcopy(self)