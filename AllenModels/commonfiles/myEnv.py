#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 18 10:37:38 2022

@author: fidel
"""

# from bmtk.utils.sim_setup import build_env_bionet

# build_env_bionet(
#     base_dir='sim_ch01',       # Where to save the scripts and config files
#     config_file='config.json', # Where main config will be saved.
#     network_dir='network',     # Location of directory containing network files
#     tstop=2000.0, dt=0.1,      # Run a simulation for 2000 ms at 0.1 ms intervals
#     report_vars=['v'],  # Tells simulator we want to record membrane potential and calcium traces
#     current_clamp={            # Creates a step current from 500.0 ms to 1500.0 ms
#         'amp': 0.120,
#         'delay': 500.0,
#         'duration': 1000.0
#     },
#     overwrite=True,
#     include_examples=False,    # Copies components files for tutorial examples
#     compile_mechanisms=False   # Will try to compile NEURON mechanisms
# )
import sys
temp = sys.argv[1]
Q = sys.argv[2]
iAmp = float(sys.argv[3])
from bmtk.utils.create_environment import create_environment

#iAmp=0.2
create_environment(
    'bionet',
    base_dir='sims/' + Q,                      # Where to save the scripts and config files
    config_file='config'+ temp +'.json',  # 
    network_dir='sims/network',             # Location of directory containing network files
    output_dir='output_iclamp'+ temp + Q ,       # Directory where simulation output will be put into
    tstop=10000.0, dt=0.1,              # Run a simulation for 2000 ms at 0.1 ms intervals
    current_clamp= {
        'amp': iAmp,                  # Current size (pA)
        'delay': 500,                  # Time from start of simulation to onset of step (ms)
        'duration': 9000               # Duration of current step (ms)
    },
    celsius = float(temp),
    report_vars=['v'],                 # Record membrane potential
    compile_mechanisms=False,           # Try to compile NEURON mechanisms
    overwrite=True                  # Overwrite pre-existing config files of same name, default is False
)

# create_environment(
#     'bionet',
#     base_dir='sim_ch34',                      # Where to save the scripts and config files
#     config_file='config.json',  # 
#     network_dir='network',             # Location of directory containing network files
#     output_dir='output_iclamp' ,       # Directory where simulation output will be put into
#     tstop=3000.0, dt=0.1,              # Run a simulation for 2000 ms at 0.1 ms intervals
#     current_clamp= {
#         'amp': iAmp,                  # Current size (pA)
#         'delay': 500,                  # Time from start of simulation to onset of step (ms)
#         'duration': 2000               # Duration of current step (ms)
#     },
#     celsius = 34,
#     report_vars=['v'],                 # Record membrane potential
#     compile_mechanisms=True,           # Try to compile NEURON mechanisms
#     overwrite=True                  # Overwrite pre-existing config files of same name, default is False
# )

# create_environment(
#     'bionet',
#     base_dir='sim_ch40',                      # Where to save the scripts and config files
#     config_file='config.json',  # 
#     network_dir='network',             # Location of directory containing network files
#     output_dir='output_iclamp' ,       # Directory where simulation output will be put into
#     tstop=3000.0, dt=0.1,              # Run a simulation for 2000 ms at 0.1 ms intervals
#     current_clamp= {
#         'amp': iAmp,                  # Current size (pA)
#         'delay': 500,                  # Time from start of simulation to onset of step (ms)
#         'duration': 2000               # Duration of current step (ms)
#     },
#     celsius = 40,
#     report_vars=['v'],                 # Record membrane potential
#     compile_mechanisms=True,           # Try to compile NEURON mechanisms
#     overwrite=True                  # Overwrite pre-existing config files of same name, default is False
# )


# # Q23
# create_environment(
#     'bionet',
#     base_dir='sim_ch21',                      # Where to save the scripts and config files
#     config_file='configQ23.json',  # 
#     network_dir='networkQ23',             # Location of directory containing network files
#     output_dir='output_iclampQ23' ,       # Directory where simulation output will be put into
#     tstop=3000.0, dt=0.1,              # Run a simulation for 2000 ms at 0.1 ms intervals
#     current_clamp= {
#         'amp': iAmp,                  # Current size (pA)
#         'delay': 500,                  # Time from start of simulation to onset of step (ms)
#         'duration': 2000               # Duration of current step (ms)
#     },
#     celsius = 21,
#     report_vars=['v'],                 # Record membrane potential
#     compile_mechanisms=True,           # Try to compile NEURON mechanisms
#     overwrite=True                  # Overwrite pre-existing config files of same name, default is False
# )

# create_environment(
#     'bionet',
#     base_dir='sim_ch34',                      # Where to save the scripts and config files
#     config_file='configQ23.json',  # 
#     network_dir='networkQ23',             # Location of directory containing network files
#     output_dir='output_iclampQ23' ,       # Directory where simulation output will be put into
#     tstop=3000.0, dt=0.1,              # Run a simulation for 2000 ms at 0.1 ms intervals
#     current_clamp= {
#         'amp': iAmp,                  # Current size (pA)
#         'delay': 500,                  # Time from start of simulation to onset of step (ms)
#         'duration': 2000               # Duration of current step (ms)
#     },
#     celsius = 34,
#     report_vars=['v'],                 # Record membrane potential
#     compile_mechanisms=True,           # Try to compile NEURON mechanisms
#     overwrite=True                  # Overwrite pre-existing config files of same name, default is False
# )

# create_environment(
#     'bionet',
#     base_dir='sim_ch40',                      # Where to save the scripts and config files
#     config_file='configQ23.json',  # 
#     network_dir='networkQ23',             # Location of directory containing network files
#     output_dir='output_iclampQ23' ,       # Directory where simulation output will be put into
#     tstop=3000.0, dt=0.1,              # Run a simulation for 2000 ms at 0.1 ms intervals
#     current_clamp= {
#         'amp': iAmp,                  # Current size (pA)
#         'delay': 500,                  # Time from start of simulation to onset of step (ms)
#         'duration': 2000               # Duration of current step (ms)
#     },
#     celsius = 40,
#     report_vars=['v'],                 # Record membrane potential
#     compile_mechanisms=True,           # Try to compile NEURON mechanisms
#     overwrite=True                  # Overwrite pre-existing config files of same name, default is False
# )
