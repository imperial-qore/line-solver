#!/usr/bin/env python3
"""Test gallery_hyphyp1_linear with MVA"""

import sys
sys.path.insert(0, '../gallery')
from gallery_hyphyp1_linear import gallery_hyphyp1_linear
from line_solver import MVA

if __name__ == '__main__':
    model = gallery_hyphyp1_linear()
    solver = MVA(model)
    avg_table = solver.getAvgTable()
    print(f'Model: {model.getName()}')
    print(avg_table)
