# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 18:39:22 2013

@author: logan
"""

from Plots import Graph #TODO: Change to absolute import


def main():
    fluid_ref = 'n-Pentane'
    for graph_type in ['Ts']: #['pt', 'ph', 'ps', 'ts', 'pt', 'prho', 'trho']:
        graph = Graph(fluid_ref, graph_type)
        graph.set_axis_limits([-0.5, 1.5, 300, 530])
        graph.draw_isolines('Q', [0.1, 0.9])
        graph.draw_isolines('P', [100, 2000])
        graph.draw_isolines('D', [2, 600])
        graph.show()

if __name__ == "__main__":
    main()
