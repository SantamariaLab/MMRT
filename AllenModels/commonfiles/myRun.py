#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 18 10:38:34 2022

@author: fidel
"""
import sys

temp = sys.argv[1]
Q = sys.argv[2]

from bmtk.simulator import bionet


# conf21 = bionet.Config.from_json('sim_ch21/config.json')
# conf21.build_env()
# conf21
# net21 = bionet.BioNetwork.from_config(conf21)
# sim21 = bionet.BioSimulator.from_config(conf21, network=net21)
# sim21.run()

print(['sims/' + Q + '/config' + temp + '.json'])
conf21 = bionet.Config.from_json("sims/" + Q + "/config" + temp + ".json")
conf21.build_env()
conf21
net21 = bionet.BioNetwork.from_config(conf21)
sim21 = bionet.BioSimulator.from_config(conf21, network=net21)
sim21.run()
