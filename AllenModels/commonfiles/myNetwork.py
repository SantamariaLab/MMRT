#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 18 10:04:52 2022

@author: fidel
"""

from bmtk.builder import NetworkBuilder

# net = NetworkBuilder("network_name")
# net.add_nodes(N=1, model_type='biophysical', ei='exc',
#               model_template='ctdb:Biophys1.hoc',
#               dynamics_parameters='fit_parameters.json',
#               morphology='reconstruction.swc',
#               model_processing='aibs_perisomatic',
#               )

import sys
temp = sys.argv[1]
Q = sys.argv[2]

net = NetworkBuilder('network_name')
net.add_nodes(
    cell_name='myCell',
    potential='exc',
    model_type='biophysical',
    model_template='ctdb:Biophys1.hoc',
    model_processing='aibs_perisomatic',
    dynamics_params='fit_parameters.json',
    morphology='reconstruction.swc'
)

net.build()
net.save_nodes(output_dir='sims/network')




# net = NetworkBuilder('network_name')
# net.add_nodes(
#     cell_name='myCell',
#     potential='exc',
#     model_type='biophysical',
#     model_template='ctdb:Biophys1.hoc',
#     model_processing='aibs_perisomatic',
#     dynamics_params='fit_parameters.json',
#     morphology='reconstruction.swc'
# )

# net.build()
# net.save_nodes(output_dir='sim_ch34/network')

# net = NetworkBuilder('network_name')
# net.add_nodes(
#     cell_name='myCell',
#     potential='exc',
#     model_type='biophysical',
#     model_template='ctdb:Biophys1.hoc',
#     model_processing='aibs_perisomatic',
#     dynamics_params='fit_parameters.json',
#     morphology='reconstruction.swc'
# )

# net.build()
# net.save_nodes(output_dir='sim_ch40/network')

# # Original Q

# net = NetworkBuilder('network_name')
# net.add_nodes(
#     cell_name='myCell',
#     potential='exc',
#     model_type='biophysical',
#     model_template='ctdb:Biophys1.hoc',
#     model_processing='aibs_perisomatic',
#     dynamics_params='fit_parametersQ23.json',
#     morphology='reconstruction.swc'
# )

# net.build()
# net.save_nodes(output_dir='sim_ch21/networkQ23')

# net = NetworkBuilder('network_name')
# net.add_nodes(
#     cell_name='myCell',
#     potential='exc',
#     model_type='biophysical',
#     model_template='ctdb:Biophys1.hoc',
#     model_processing='aibs_perisomatic',
#     dynamics_params='fit_parametersQ23.json',
#     morphology='reconstruction.swc'
# )

# net.build()
# net.save_nodes(output_dir='sim_ch34/networkQ23')

# net = NetworkBuilder('network_name')
# net.add_nodes(
#     cell_name='myCell',
#     potential='exc',
#     model_type='biophysical',
#     model_template='ctdb:Biophys1.hoc',
#     model_processing='aibs_perisomatic',
#     dynamics_params='fit_parametersQ23.json',
#     morphology='reconstruction.swc'
# )

# net.build()
# net.save_nodes(output_dir='sim_ch40/networkQ23')





for node in net.nodes():
    print(node)