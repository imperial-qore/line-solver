#!/usr/bin/env python3
"""Test gallery_merl1_tandem with MVA"""

import sys
sys.path.insert(0, '../gallery')
from gallery_merl1_tandem import gallery_merl1_tandem
from line_solver import MVA

if __name__ == '__main__':
    model = gallery_merl1_tandem()
    solver = MVA(model)
    avg_table = solver.getAvgTable()
    print(f'Model: {model.getName()}')
    print(avg_table)
