#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 18 11:38:39 2022

@author: fidel
"""

from bmtk.analyzer.compartment import plot_traces

_ = plot_traces(config_file='sims/MMRT/config21.json', node_ids=[0], report_name='v_report')
_ = plot_traces(config_file='sims/Q23/config21.json', node_ids=[0], report_name='v_report')

_ = plot_traces(config_file='sims/MMRT/config34.json', node_ids=[0], report_name='v_report')
_ = plot_traces(config_file='sims/Q23/config34.json', node_ids=[0], report_name='v_report')

_ = plot_traces(config_file='sims/MMRT/config38.json', node_ids=[0], report_name='v_report')
_ = plot_traces(config_file='sims/Q23/config38.json', node_ids=[0], report_name='v_report')

_ = plot_traces(config_file='sims/MMRT/config42.json', node_ids=[0], report_name='v_report')
_ = plot_traces(config_file='sims/Q23/config42.json', node_ids=[0], report_name='v_report')

# curve = plot_traces(config_file='sim_ch01/config.json', node_ids=[0], report_name='v_report')