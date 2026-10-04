#!/usr/bin/env python3
"""Test gallery_gamm1 with MVA"""

import sys
sys.path.insert(0, '../gallery')
from gallery_gamm1 import gallery_gamm1
from line_solver import MVA

if __name__ == '__main__':
    model = gallery_gamm1()
    solver = MVA(model)
    avg_table = solver.getAvgTable()
    print(f'Model: {model.getName()}')
    print(avg_table)
