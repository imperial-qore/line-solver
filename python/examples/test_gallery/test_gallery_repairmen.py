#!/usr/bin/env python3
"""Test gallery_repairmen with MVA"""

import sys
sys.path.insert(0, '../gallery')
from gallery_repairmen import gallery_repairmen
from line_solver import MVA

if __name__ == '__main__':
    model = gallery_repairmen()
    solver = MVA(model)
    avg_table = solver.getAvgTable()
    print(f'Model: {model.getName()}')
    print(avg_table)
