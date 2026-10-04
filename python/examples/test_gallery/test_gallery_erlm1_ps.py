#!/usr/bin/env python3
"""Test gallery_erlm1_ps with MVA"""

import sys
sys.path.insert(0, '../gallery')
from gallery_erlm1_ps import gallery_erlm1_ps
from line_solver import MVA

if __name__ == '__main__':
    model = gallery_erlm1_ps()
    solver = MVA(model)
    avg_table = solver.getAvgTable()
    print(f'Model: {model.getName()}')
    print(avg_table)
