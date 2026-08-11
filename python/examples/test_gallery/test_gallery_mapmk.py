#!/usr/bin/env python3
"""Test gallery_mapmk with MVA"""

import sys
sys.path.insert(0, '../gallery')
from gallery_mapmk import gallery_mapmk
from line_solver import MVA

if __name__ == '__main__':
    model = gallery_mapmk()
    solver = MVA(model)
    avg_table = solver.getAvgTable()
    print(f'Model: {model.getName()}')
    print(avg_table)
