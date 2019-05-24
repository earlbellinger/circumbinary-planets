#!usr/bin/env python
#### Call planetplanet with quasi-random inputs 
#### Author: Earl Bellinger ( bellinger@phys.au.dk ) 
#### Stellar Astrophysics Centre Aarhus 

import numpy as np
import subprocess
import argparse
import sys 
import random 
import os
import glob
import re
from time import sleep
from tqdm import tqdm

import pandas as pd
import numpy as np

from sys import path
path.append('/home/earl/asteroseismology/scripts')
from sobol_lib import i4_sobol

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

output = parser.add_argument_group('Output')
output.add_argument('--output_dir',      default='simulations', type=str)
output.add_argument('--lightcurves_dir', default='lightcurves', type=str)
output.add_argument('--systems_dir',     default='systems',     type=str)

gen = parser.add_argument_group('General')
gen.add_argument('-N', '--systems', default=1024,    type=int)
gen.add_argument('-s', '--skip',    default=20000,   type=int)
gen.add_argument('--age',           default=[1, 10], type=float)

A = parser.add_argument_group('Primary')
A.add_argument('--A_mass',      default=[0.6, 1.5],   type=float)
A.add_argument('--A_limbdark1', default=[0.6, 0.6],   type=float) # compute from Claret 
A.add_argument('--A_limbdark2', default=[0, 0],       type=float)

B = parser.add_argument_group('Secondary')
B.add_argument('--B_mass',      default=[0.2, 1],    type=float, help='Mass ratio for B component')
B.add_argument('--B_per',       default=[7, 50],      type=float)
B.add_argument('--B_inc',       default=[80, 90],     type=float)
B.add_argument('--B_ecc',       default=[0, 0.7],     type=float)
B.add_argument('--B_w',         default=[0, 360],     type=float)
B.add_argument('--B_Omega',     default=[0, 0],       type=float)
B.add_argument('--B_t0',        default=[0, 1],       type=float) # phase 
A.add_argument('--B_limbdark1', default=[0.6, 0.6],   type=float) # compute from Claret 
A.add_argument('--B_limbdark2', default=[0, 0],       type=float)

b = parser.add_argument_group('Companion')
b.add_argument('--b_exist',     default=[-1,   1],    type=float) # planet is there if b_exist > 0
b.add_argument('--b_radius',    default=[0.09, 1],    type=float) # R_earth to R_jupiter, logarithmic 
b.add_argument('--b_per',       default=[0, 1],       type=float) # 0=B_per, 1=(2 years)
b.add_argument('--b_inc',       default=[1e-4, 10],   type=float) # mutual inclination angle, logarithmic
b.add_argument('--b_ecc',       default=[0, 0.5],     type=float)
b.add_argument('--b_w',         default=[0, 360],     type=float)
b.add_argument('--b_Omega',     default=[1e-7, 1e-2], type=float) # logarithmic
b.add_argument('--b_t0',        default=[0, 1],       type=float) # phase 

logs = ['b_inc', 'b_Omega']

args, unknown = parser.parse_known_args(sys.argv[1:]) 
#args = parser.parse_args(sys.argv[1:])

# collect all the vars with two args so we can pick within the intervals 
flags = []
ranges = []
for flag in vars(args):
    val = vars(args)[flag]
    if isinstance(val, list):
        vals = val
        if flag in logs:
            vals = np.log10(vals)
        ranges += [vals]
        flags += [flag]

shift = np.array(ranges)[:,0]
scale = np.array([(b-a) for a,b in ranges])
nlflag = ' \\\n\t--'


# parse isochrones 
print('Parsing isochrones')
iso_dir = 'MIST_v1.2_feh_p0.00_afe_p0.0_vvcrit0.0_EEPS'
track_files = glob.glob(os.path.join(iso_dir, '*.eep'))
isochrones = {}
for track_file in tqdm(track_files):
    with open(track_file) as f:
        head = [next(f) for x in range(12)]
    header = re.split('\\s+', head[-1].rstrip().strip("#"))[1:]
    DF = pd.read_table(track_file, names=header, sep='\\s+', comment='#')
    mass = int(os.path.basename(track_file).split('M')[0])/100
    isochrones[mass] = DF


# call lightcurve with quasirandomly generated inputs 
for id in range(args.skip, args.systems+args.skip):
    vals = shift+np.array(i4_sobol(len(ranges), id)[0])*scale
    
    flag_dict = {flag: 10**val if flag in logs else val 
        for flag, val in zip(flags, vals)}
    
    # figure out which stellar model to use 
    A_mass = flag_dict['A_mass']
    M_idx = np.argmin([abs(A_mass - M) for M in isochrones.keys()])
    A_mass = list(isochrones.keys())[M_idx]
    track = isochrones[A_mass]
    t_idx = np.argmin(abs(flag_dict['age']*10**9 - track['star_age']))
    flag_dict['age'] = track['star_age'][t_idx] / 10**9
    
    # pick A_radius, A_Teff from the model 
    flag_dict['A_mass'] = A_mass
    flag_dict['A_radius'] = 10**track['log_R'][t_idx]
    flag_dict['A_Teff'] = 10**track['log_Teff'][t_idx]
    
    # compute B_mass from B_mass_ratio 
    ratio = flag_dict['B_mass'] 
    B_mass = flag_dict['A_mass'] * ratio 
    flag_dict['B_mass'] = B_mass
    
    # pick B_radius, B_Teff from a model 
    M_idx = np.argmin([abs(B_mass - M) for M in isochrones.keys()])
    B_mass = list(isochrones.keys())[M_idx]
    track = isochrones[B_mass]
    t_idx = np.argmin(abs(flag_dict['age']*10**9 - track['star_age']))
    
    flag_dict['B_mass'] = B_mass
    flag_dict['B_radius'] = 10**track['log_R'][t_idx]
    flag_dict['B_Teff'] = 10**track['log_Teff'][t_idx]
    
    # compute b_mass from mass-radius relation 
    # Weiss & Marcy 2014
    R_earth = 0.0892141778
    M_earth = 0.00314635457
    b_radius = flag_dict['b_radius']
    if b_radius < 1.5 * R_earth:
        rho = 2.43 + 3.39 * b_radius / R_earth
        flag_dict['b_mass'] = rho * b_radius**3
    else:
        flag_dict['b_mass'] = 2.69 * (b_radius / R_earth)**0.93 * M_earth
    
    if flag_dict['b_exist'] <= 0:
        flag_dict['b_mass'] *= -1 
    
    # compute b_per, where 0=T_c and 1=(2 years) 
    # Holman & Wiegert 1999 
    # T_c is the critical period 
    # a_c is the critical semimajor axis 
    # a_b is the binary semimajor axis 
    # e is the binary eccentricity 
    # mu is the mass ratio 
    e = flag_dict['B_ecc']
    m2 = flag_dict['B_mass']
    m1 = flag_dict['A_mass']
    mu = m2 / (m1 + m2)
    # get a_b from Kepler's third law 
    # T**2 = 4 * pi**2 * a**3 / (G*M)
    # => a = (G*M)^(1/3) * T^(2/3) * (2*pi)^(-2/3)
    B_per = flag_dict['B_per'] 
    a_b = (m1+m2)**(1/3) * B_per**(2/3) * (2*np.pi)**(-2/3)
    a_c = (1.6 + 5.1*e + (-2.22)*e**2 + (4.12)*mu + (-4.27)*e*mu + 
        (-5.09)*mu**2 + 4.61*e**2*mu**2) * a_b
    T_c = np.sqrt(4*np.pi**2 * a_c**3 / (m1+m2))
    two_years = 365*2
    if T_c > two_years:
        print('Critical period too long: skipping') 
        continue 
    flag_dict['b_per'] = T_c + flag_dict['b_per'] * (two_years - T_c)
    
    # compute B_t0 by multiplying by the period 
    factor = flag_dict['B_t0']
    B_t0 = flag_dict['B_per'] * factor 
    #print('Transforming B_t0 from factor', factor, 'to value', B_t0)
    flag_dict['B_t0'] = B_t0
    
    # compute b_t0 by multiplying by the period 
    factor = flag_dict['b_t0']
    b_t0 = flag_dict['b_per'] * factor 
    #print('Transforming b_t0 from factor', factor, 'to value', b_t0)
    flag_dict['b_t0'] = b_t0
    
    # pick a random sign for b_inc and transform from mutual inclination angle 
    flag_dict['b_inc'] *= 1 if random.random() < 0.5 else -1
    flag_dict['b_inc'] += flag_dict['B_inc']
    
    clargs = [nlflag + flag + ' ' + str(flag_dict[flag])
        for flag in flag_dict]
    clargs = ' '.join(clargs)
    
    bash_cmd = "python3 lightcurve.py --id " + str(id) +\
        nlflag + "output_dir "      +      args.output_dir +\
        nlflag + "lightcurves_dir " + args.lightcurves_dir +\
        nlflag + "systems_dir "     +     args.systems_dir +\
        clargs#%\
        #tuple([id] + 
        #      [output_dir, ])
    print(bash_cmd)
    #exit()
    process = subprocess.Popen(bash_cmd.split(), shell=False)
    process.wait()
