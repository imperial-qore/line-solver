#!/usr/bin/env python3
"""Test gallery_mdk with MVA"""

import sys
sys.path.insert(0, '../gallery')
from gallery_mdk import gallery_mdk
from line_solver import MVA

if __name__ == '__main__':
    model = gallery_mdk()
    solver = MVA(model)
    avg_table = solver.getAvgTable()
    print(f'Model: {model.getName()}')
    print(avg_table)
