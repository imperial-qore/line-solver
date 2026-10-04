#!/usr/bin/env python3
"""Test gallery_cqn with MVA"""

import sys
sys.path.insert(0, '../gallery')
from gallery_cqn import gallery_cqn
from line_solver import MVA

if __name__ == '__main__':
    model = gallery_cqn()
    solver = MVA(model)
    avg_table = solver.getAvgTable()
    print(f'Model: {model.getName()}')
    print(avg_table)
