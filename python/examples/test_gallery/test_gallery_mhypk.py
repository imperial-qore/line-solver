#!/usr/bin/env python3
"""Test gallery_mhypk with MVA"""

import sys
sys.path.insert(0, '../gallery')
from gallery_mhypk import gallery_mhypk
from line_solver import MVA

if __name__ == '__main__':
    model = gallery_mhypk()
    solver = MVA(model)
    avg_table = solver.getAvgTable()
    print(f'Model: {model.getName()}')
    print(avg_table)
