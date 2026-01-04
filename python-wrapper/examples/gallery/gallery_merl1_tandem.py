#!/usr/bin/env python3
"""Gallery Example: M/Erl/1-Tandem Queue"""

from line_solver import *
from .gallery_merl1_linear import gallery_merl1_linear

def gallery_merl1_tandem():
    return gallery_merl1_linear(2)

if __name__ == '__main__':
    model = gallery_merl1_tandem()
    print(f"Model: {model.getName()}")
