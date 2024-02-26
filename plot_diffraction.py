#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 21 10:09:30 2023

@author: Sam Lambrick

Plots simulated SHeM diffraction data stored in the directory provided as a
command line argument when running the script.
"""

import numpy as np
import SHeMDiffractionAnalysis.shem_spot_profile as ssp
import matplotlib.pyplot as plt
import sys

data_dir = sys.argv[1]
if data_dir[-1] != '/':
    data_dir = data_dir + '/'
if data_dir[0] == '/':
    data_dir = data_dir[1:]

#data_dir = 'simulations/0007_LiF_simple_diffraction_test/'

data_dir = 'simulations/' + data_dir

# New simulated data
sim_data = ssp.SpotProfile.import_ray(data_dir)

print(sim_data.signal)
print(sim_data.DK)
sim_data.shem_raw_plot()
plt.savefig(data_dir + 'diffraction_zs.eps')
plt.savefig(data_dir + 'diffraction_zs.png', dpi=300)



fig, a2, cax = sim_data.shem_diffraction_plot(bar_location='right', 
                                             figsize = [8, 6], 
                                             x_offset=0.08,
                                             logplot=False)
a2.set_xticklabels(['0°', '45°', '$\\alpha=90$°', '135°', '180°',
                          '225°', '270°', '315°'])
a2.set_xlabel(None)
a2.set_rlabel_position(0)
cax.set_ylabel('$\\log_{10} I$')
plt.savefig(data_dir + 'diffraction_Ks.eps')
plt.savefig(data_dir + 'diffraction_Ks.png', dpi=300)
