#!/usr/bin/env python3
"""Test gallery_mmapk with MVA"""

import sys
sys.path.insert(0, '../gallery')
from gallery_mmapk import gallery_mmapk
from line_solver import MVA

if __name__ == '__main__':
    model = gallery_mmapk()
    solver = MVA(model)
    avg_table = solver.getAvgTable()
    print(f'Model: {model.getName()}')
    print(avg_table)
