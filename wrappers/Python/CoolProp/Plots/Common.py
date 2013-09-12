# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 15:13:20 2013

@author: logan
"""

class BasePlot(object):
    AXIS_LABLES = {'T': ["Temperature", r"[$K$]"],
                   'P': ["Pressure", r"[$kPa$]"],
                   'S': ["Entropy", r"[$kJ/kg K$]"],
                   'H': ["Enthalpy", r"[$kJ/kg$]"],
                   'V': [],
                   'RHO': ["Density", r"[$kg/m^3$]"]}

    COLOR_MAP = {'T': 'Darkred',
                 'P': 'DarkCyan',
                 'H': 'DarkGreen',
                 'D': 'DarkBlue',
                 'S': 'DarkOrange',
                 'Q': 'black',}

    LINE_IDS = {'TS': ['P', 'D'], #'H'],
                'PH': ['S', 'T', 'D'],
                'HS': ['P'], #'T', 'D'],
                'PS': ['H', 'T', 'D'],
                'PD': ['T', 'S', 'H'],
                'TD': ['P'], #'S', 'H'],
                'PT': ['D', 'P', 'S'],}

    def __init__(self, fluid_ref, graph_type, **kwargs):
        if not isinstance(graph_type, str):
            raise TypeError("Invalid graph_type input, expeceted a string")

        graph_type = graph_type.upper()
        if len(graph_type) >= 2 and graph_type[1:len(graph_type)] == 'RHO':
            graph_type = graph_type[0] + graph_type[1:len(graph_type)]

        if graph_type.upper() not in self.LINE_IDS.keys():
            raise ValueError("You have to specify the kind of plot, use " \
                              + str(self.LINE_IDS.keys()))