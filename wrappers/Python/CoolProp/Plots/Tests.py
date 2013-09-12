# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 18:39:22 2013

@author: logan
"""

from Plots import Graph #TODO: Change to absolute import


def main():
    fluid_ref = 'R245fa'
    for graph_type in ['Ts', 'Ps']: #['pt', 'ph', 'ps', 'ts', 'pt', 'prho', 'trho']:
        graph = Graph(fluid_ref, graph_type)
        graph.draw_isolines('Q', [0.1, 0.9])
        graph.show()

if __name__ == "__main__":
    main()
