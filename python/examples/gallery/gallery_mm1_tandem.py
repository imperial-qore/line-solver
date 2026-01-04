#!/usr/bin/env python3
"""Gallery Example: M/M/1-Tandem Queue"""

from line_solver import *
from .gallery_mm1_linear import gallery_mm1_linear

def gallery_mm1_tandem():
    return gallery_mm1_linear(2)

if __name__ == '__main__':
    model = gallery_mm1_tandem()
    print(f"Model: {model.getName()}")
