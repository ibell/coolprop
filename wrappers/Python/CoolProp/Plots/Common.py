# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 15:13:20 2013

@author: logan
"""

import matplotlib


class BasePlot(object):
    #TODO: Simplify / Consolidate dictionary maps
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
                 'Q': 'black'}

    SYMBOL_MAP = {'T' : [r'$T = ','$ K'],
                  'P' : [r'$p = ','$ kPa'],
                  'H' : [r'$h = ','$ kJ/kg'],
                  'D' : [r'$\rho = ','$ kg/m$^3$'],
                  'S' : [r'$s = ','$ kJ/kg-K'],
                  'Q' : [r'$x = ','$']}

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

        self.fluid_ref = fluid_ref
        self.graph_type = graph_type.upper()

        self.axis = kwargs.get('axis', matplotlib.pyplot.gca())


    def __sat_bounds(self, kind, smin=None, smax=None):
        """
        Generates limits for the saturation line in either T or p determined by 'kind'.
        If xmin or xmax are provided, values will be checked against the allowable
        range for the EOS and an error might be generated.

        Returns a tuple containing (xmin,xmax)
        """
        if kind == 'P':
            name = 'pressure'
            minKey = 'ptriple'
        elif kind == 'T':
            name = 'temperature'
            minKey = 'Tmin'

        fluid_crit = CP.Props(self.fluid_ref, str(kind) + 'crit')

        if smin is None:
            smin = CP.Props(self.fluid_ref, str(minKey)) + SMALL
        elif smin > fluid_crit:
            raise ValueError('Minimum ' + name +
                             ' cannot be greater than fluid critical ' +
                             name + '.')

        if smax is None:
            smax = CP.Props(self.fluid_ref, str(kind) + 'crit') - SMALL
        elif smax > fluid_crit:
            raise ValueError('Maximum ' + name +
                             ' cannot be greater than fluid critical ' +
                             name + '.')

        smin = max(smin, CP.Props(self.fluid_ref, minKey) + SMALL)
        smax = min(smax, CP.Props(self.fluid_ref, kind + 'crit') - SMALL)

        return smin, smax

    def _get_fluid_data(self, req_prop, prop1_name,
                        prop2_name, prop1_vals, prop2_vals):
        """
        Calculates lines for constant iName (iVal) over an interval of xName (xVal).
        Returns (x[],y[]) - a tuple of arrays containing the values in x and y dimensions.
        """
        if len(prop1_vals) != len(prop2_vals):
            raise ValueError('We need the same number of x value arrays as iso quantities.')

        y_vals = []
        x_vals = []
        for i, p1_val in enumerate(prop1_vals):
            x_vals.append(prop2_vals[i])
            y_vals.append([CP.Props(req_prop,
                                    prop1_name,
                                    p1_val,
                                    prop2_name,
                                    x,
                                    self.fluid_ref) for x in prop2_vals[i]])
        return [x_vals, y_vals]

    def _get_sat_lines(self, kind='T', smin=None, smax=None, num=500, x=[0., 1.]):
        """
        Calculates bubble and dew line in the quantities for your plot.
        You can specify if you need evenly spaced entries in either
        pressure or temperature by supplying kind='p' and kind='T'
        (default), respectively.
        Limits can be set with kmin (default: minimum from EOS) and
        kmax (default: critical value).
        Returns lines[] - a 2D array of dicts containing 'x' and 'y' coordinates
        for bubble and dew line. Additionally, the dict holds the keys
        'kmax', 'label' and 'opts', those can be used for plotting as well.
        """
        if not kind.upper() in ['T', 'P']:
            raise ValueError("Invalid input for determining the saturation \
                              lines... Expected either 'T' or 'P'")

        smin, smax = self.__sat_bounds(kind, xmin=smin, xmax=smax)
        sat_range = numpy.linspace(smin, smax, num)
        sat_mesh = [sat_range for i in x]

        x_vals = sat_mesh
        y_vals = sat_mesh
        if self.graph_type[1] != kind:
            _, x_vals = self._get_fluid_data(self.graph_type[1], 'Q',
                                             kind, x, sat_mesh)

        if self.graph_type[0] != kind:
            _, y_vals = self._get_fluid_data(self.graph_type[0],'Q',
                                             kind, x, sat_mesh)

        # Merge the two lines, capital Y holds important information.
        # We merge on X values
        # Every entry, eg. Xy, contains two arrays of values.
        sat_lines = []
        for i in range(len(x_vals)): # two dimensions: i = {0,1}
            line = {'x': x_vals[i],
                    'y': y_vals[i],
                    'smax': smax}

            line['label'] = self.SYMBOL_MAP['Q'][0] + str(x[i])
            line['type'] = 'Q'
            line['value'] = x[i]
            line['unit'] = self.SYMBOL_MAP['Q'][1]
            line['opts'] = {'color': self.COLOR_MAP['Q'],
                            'lw': 1.0}

            if x[i] == 0.:
                line['label'] = 'bubble line'
            elif x[i] == 1.:
                line['label'] = 'dew line'
            else:
                line['opts']['lw'] = 0.75
                line['opts']['alpha'] = 0.5

            sat_lines.append(line)

        return sat_lines