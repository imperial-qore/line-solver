#!/usr/bin/env python3
"""Gallery Example: M/Hyp/1-Tandem Queue"""

from line_solver import *
from .gallery_mhyp1_linear import gallery_mhyp1_linear

def gallery_mhyp1_tandem():
    return gallery_mhyp1_linear(2)

if __name__ == '__main__':
    model = gallery_mhyp1_tandem()
    print(f"Model: {model.getName()}")
