#!/usr/bin/env python3
"""Gallery Example: Hyp/Hyp/1-Tandem Queue"""

from line_solver import *
from .gallery_hyphyp1_linear import gallery_hyphyp1_linear

def gallery_hyphyp1_tandem():
    return gallery_hyphyp1_linear(2)

if __name__ == '__main__':
    model = gallery_hyphyp1_tandem()
    print(f"Model: {model.getName()}")
