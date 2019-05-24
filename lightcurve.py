#!usr/bin/env python
#### Compute the lightcurve of a circumbinary system 
#### Author: Earl Bellinger ( bellinger@phys.au.dk ) 
#### Stellar Astrophysics Centre Aarhus 

from planetplanet import Planet, Star, System
from planetplanet.constants import *
from planetplanet.photo.maps import UniformMap

import numpy as np 
import pandas as pd

import os, sys
import argparse

from astropy.constants import *
M_je = M_jup / M_earth
R_je = R_jup / R_earth

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

output = parser.add_argument_group('Output')
output.add_argument('--output_dir',      default='simulations', type=str)
output.add_argument('--lightcurves_dir', default='lightcurves', type=str)
output.add_argument('--systems_dir',     default='systems',     type=str)

gen = parser.add_argument_group('General')
gen.add_argument('--id',         default=1,                     type=int)
#gen.add_argument('--cadence',    default=MINUTE*(58.84876/60),  type=float)
gen.add_argument('--cadence',    default=MINUTE*29.4244,        type=float)
gen.add_argument('--duration',   default=365*4,                 type=float)
gen.add_argument('--lambda1',    default=0.4,                   type=float)
gen.add_argument('--lambda2',    default=0.865,                 type=float)
gen.add_argument('--integrator', default='ias15',               type=str)

A = parser.add_argument_group('Primary')
A.add_argument('--A_mass',       default=0.6897,   type=float)
A.add_argument('--A_radius',     default=0.6489,   type=float)
A.add_argument('--A_teff',       default=4450,     type=float)
A.add_argument('--A_limbdark1',  default=0.6,      type=float)
A.add_argument('--A_limbdark2',  default=0,        type=float)

B = parser.add_argument_group('Secondary')
B.add_argument('--B_mass',      default=0.20255,   type=float)
B.add_argument('--B_radius',    default=0.22623,   type=float)
B.add_argument('--B_teff',      default=3000,      type=float)
B.add_argument('--B_per',       default=41.079220, type=float)
B.add_argument('--B_inc',       default=90.3401,   type=float)
B.add_argument('--B_ecc',       default=0.15944,   type=float)
B.add_argument('--B_w',         default=263.464,   type=float)
B.add_argument('--B_Omega',     default=0,         type=float)
B.add_argument('--B_t0',        default=31,        type=float)
A.add_argument('--B_limbdark1', default=0.6,       type=float)
A.add_argument('--B_limbdark2', default=0,         type=float)

b = parser.add_argument_group('Companion')
b.add_argument('--b_mass',      default=0.333,     type=float)
b.add_argument('--b_radius',    default=0.7538,    type=float)
b.add_argument('--b_per',       default=228.776,   type=float)
b.add_argument('--b_inc',       default=90.0322,   type=float)
b.add_argument('--b_ecc',       default=0.0069,    type=float)
b.add_argument('--b_w',         default=318,       type=float)
b.add_argument('--b_Omega',     default=0.003,     type=float)
b.add_argument('--b_t0',        default=203,       type=float)

#args = parser.parse_args(sys.argv[1:])
args, unknown = parser.parse_known_args(sys.argv[1:]) 

A = Star('A', m=args.A_mass, r=args.A_radius, teff=args.A_teff, 
    nz=31, color='gold', limbdark=[args.A_limbdark1, args.A_limbdark2])

B = Star('B', m=args.B_mass, r=args.B_radius, teff=args.B_teff,
    per=args.B_per, inc=args.B_inc, ecc=args.B_ecc, w=args.B_w, 
    Omega=args.B_Omega, t0=args.B_t0,
    host=None,
    nz=31, color='darkblue', limbdark=[args.B_limbdark1, args.B_limbdark2])

if args.b_mass > 0:
    planet = Planet('b', m=args.b_mass*M_je, r=args.b_radius*R_je, 
        per=args.b_per, inc=args.b_inc, ecc=args.b_ecc, w=args.b_w, 
        Omega=args.b_Omega, t0=args.b_t0,
        host=None,
        nz=1, color='green', radiancemap=UniformMap())
    
    system = System(A, B, planet, 
        nbody=True, quiet=False,
        integrator=args.integrator,
        timestep=args.cadence)
else:
    system = System(A, B, 
        nbody=True, quiet=False,
        integrator=args.integrator,
        timestep=args.cadence)

time = np.arange(0, args.duration, args.cadence)
system.compute(time, lambda1=args.lambda1, lambda2=args.lambda2)

df = pd.DataFrame(vars(args), index=[0])
flux = system.flux.T[0]/np.max(system.flux.T[0])
nonone = flux != 1
lc = pd.DataFrame(data={'Time': system.time[nonone], 
                        'Flux': flux[nonone]})

for dir in [args.output_dir, 
            os.path.join(args.output_dir, args.lightcurves_dir),
            os.path.join(args.output_dir, args.systems_dir)]:
    if not os.path.isdir(dir):
        os.makedirs(dir)

lc.to_csv(os.path.join(args.output_dir, args.lightcurves_dir, str(args.id)), 
    sep='\t', index=False)

df.to_csv(os.path.join(args.output_dir, args.systems_dir, str(args.id)), 
    sep='\t', index=False)
