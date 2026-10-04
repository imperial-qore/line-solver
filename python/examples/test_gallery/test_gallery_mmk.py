#!/usr/bin/env python3
"""Test gallery_mmk with MVA"""

import sys
sys.path.insert(0, '../gallery')
from gallery_mmk import gallery_mmk
from line_solver import MVA

if __name__ == '__main__':
    model = gallery_mmk()
    solver = MVA(model)
    avg_table = solver.getAvgTable()
    print(f'Model: {model.getName()}')
    print(avg_table)
