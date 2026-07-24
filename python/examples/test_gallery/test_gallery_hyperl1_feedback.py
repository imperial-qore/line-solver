#!/usr/bin/env python3
"""Test gallery_hyperl1_feedback with MVA"""

import sys
sys.path.insert(0, '../gallery')
from gallery_hyperl1_feedback import gallery_hyperl1_feedback
from line_solver import MVA

if __name__ == '__main__':
    model = gallery_hyperl1_feedback()
    solver = MVA(model)
    avg_table = solver.getAvgTable()
    print(f'Model: {model.getName()}')
    print(avg_table)
