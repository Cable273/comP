#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import pandas as pd
import math

def find_index_bisection(s,x):
    "find index of array x which contains s using bisection"
    b_min = 1
    b_max = np.size(x)
    breaker=0

    b=0
    if x[0] == s:
        return 0
    elif x[np.size(x)-1] == s:
        return int(np.size(x)-1)
    else:
        while x[b] != s:
            if (b_max - b_min) % 2 ==0:
                b = int(b_min + (b_max - b_min)/2)
                if s < x[b]:
                    b_max = b-1
                elif s > x[b]:
                    b_min = b+1
            else:
                b = int(b_min + (b_max - b_min-1)/2)
                if s < x[b]:
                    b_max = b-1
                elif s > x[b]:
                    b_min = b+1
        return int(b)

