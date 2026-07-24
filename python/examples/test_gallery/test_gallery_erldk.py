#!/usr/bin/env python3
"""Test gallery_erldk with MVA"""

import sys
sys.path.insert(0, '../gallery')
from gallery_erldk import gallery_erldk
from line_solver import MVA

if __name__ == '__main__':
    model = gallery_erldk()
    solver = MVA(model)
    avg_table = solver.getAvgTable()
    print(f'Model: {model.getName()}')
    print(avg_table)
