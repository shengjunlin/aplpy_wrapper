#! /usr/bin/env python
# #########################################################
# Author : Axel Donath (CfA)
# Email : https://gist.github.com/adonath
# Source: https://gist.github.com/adonath/c9a97d2f2d964ae7b9eb
# Description : Provide the ds9 colormaps:
#               'ds9grey'
#               'ds9a'
#               'ds9b'
#               'ds9bb'
#               'ds9he'
#               'ds9i8'
#               'ds9aips0'
#               'ds9heat'
#               'ds9cool'
#               'ds9rainbow'
#               'ds9rainbow_r'
# #########################################################
from matplotlib.cm import register_cmap, cmap_d
import numpy as np

ds9b = {'red': lambda v : 4 * v - 1, 
        'green': lambda v : 4 * v - 2,
        'blue': lambda v : np.select([v < 0.25, v < 0.5, v < 0.75, v <= 1],
                                      [4 * v, -4 * v + 2, 0, 4 * v - 3])}

# Note that this definition slightly differs from ds9cool, but make more sense to me...
ds9cool = {'red': lambda v : 2 * v - 1, 
           'green': lambda v : 2 * v - 0.5,
           'blue': lambda v : 2 * v}

ds9a = {'red': lambda v : np.interp(v, [0, 0.25, 0.5, 1],
                                        [0, 0, 1, 1]),
         'green': lambda v : np.interp(v, [0, 0.25, 0.5, 0.77, 1],
                                          [0, 1, 0, 0, 1]),
         'blue': lambda v : np.interp(v, [0, 0.125, 0.5, 0.64, 0.77, 1],
                                         [0, 0, 1, 0.5, 0, 0])}

ds9i8 = {'red': lambda v : np.where(v < 0.5, 0, 1), 
        'green': lambda v : np.select([v < 1/8., v < 0.25, v < 3/8., v < 0.5,
                                       v < 5/8., v < 0.75, v < 7/8., v <= 1],
                                      [0, 1, 0, 1, 0, 1, 0, 1]),
        'blue': lambda v : np.select([v < 1/8., v < 0.25, v < 3/8., v < 0.5,
                                      v < 5/8., v < 0.75, v < 7/8., v <= 1],
                                      [0, 0, 1, 1, 0, 0, 1, 1])}

ds9aips0 = {'red': lambda v : np.select([v < 1/9., v < 2/9., v < 3/9., v < 4/9., v < 5/9.,
                                        v < 6/9., v < 7/9., v < 8/9., v <= 1],
                                        [0.196, 0.475, 0, 0.373, 0, 0, 1, 1, 1]), 
            'green': lambda v : np.select([v < 1/9., v < 2/9., v < 3/9., v < 4/9., v < 5/9.,
                                        v < 6/9., v < 7/9., v < 8/9., v <= 1],
                                        [0.196, 0, 0, 0.655, 0.596, 0.965, 1, 0.694, 0]),
            'blue': lambda v : np.select([v < 1/9., v < 2/9., v < 3/9., v < 4/9., v < 5/9.,
                                        v < 6/9., v < 7/9., v < 8/9., v <= 1],
                                        [0.196, 0.608, 0.785, 0.925, 0, 0, 0, 0, 0])}

ds9rainbow =   {'red': lambda v : np.interp(v, [0, 0.2, 0.6, 0.8, 1], [1, 0, 0, 1, 1]),
              'green': lambda v : np.interp(v, [0, 0.2, 0.4, 0.8, 1], [0, 0, 1, 1, 0]),
               'blue': lambda v : np.interp(v, [0, 0.4, 0.6, 1], [1, 1, 0, 0])}

ds9rainbow_r = {'red': lambda v : np.interp(v, [0, 0.2, 0.4, 0.8, 1], [1, 1, 0, 0, 1]),
              'green': lambda v : np.interp(v, [0, 0.2, 0.6, 0.8, 1], [0, 1, 1, 0, 0]),
               'blue': lambda v : np.interp(v, [0, 0.4, 0.6, 1], [0, 0, 1, 1])}

# This definition seems a bit strange...
ds9he = {'red': lambda v : np.interp(v, [0, 0.015, 0.25, 0.5, 1],
                                        [0, 0.5, 0.5, 0.75, 1]),
         'green': lambda v : np.interp(v, [0, 0.065, 0.125, 0.25, 0.5, 1],
                                          [0, 0, 0.5, 0.75, 0.81, 1]),
         'blue': lambda v : np.interp(v, [0, 0.015, 0.03, 0.065, 0.25, 1],
                                         [0, 0.125, 0.375, 0.625, 0.25, 1])}

ds9heat = {'red': lambda v : np.interp(v, [0, 0.34, 1], [0, 1, 1]),
           'green': lambda v : np.interp(v, [0, 1], [0, 1]),
           'blue': lambda v : np.interp(v, [0, 0.65, 0.98, 1], [0, 0, 1, 1])}




# Set aliases, where colormap exists in matplotlib
cmap_d['ds9bb'] = cmap_d['afmhot']
cmap_d['ds9grey'] = cmap_d['gray']

# Register all other colormaps
register_cmap('ds9b', data=ds9b)
register_cmap('ds9cool', data=ds9cool)
register_cmap('ds9a', data=ds9a)
register_cmap('ds9i8', data=ds9i8)
register_cmap('ds9aips0', data=ds9aips0)
register_cmap('ds9rainbow', data=ds9rainbow)
register_cmap('ds9rainbow_r', data=ds9rainbow_r)
register_cmap('ds9he', data=ds9he)
register_cmap('ds9heat', data=ds9heat)

