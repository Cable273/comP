#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import pandas as pd

from comP import op_string_seq,op_string

Hp = op_string_seq()
Hp.update(1,'P-P',4,3)
Hp.update(1,'P+P',4,0)
Hp.update(1,'P+P',4,1)
Hp.update(1,'P+P',4,2)
print("H+")
Hp.print()
print("\n")

Hm = op_string_seq()
Hm.update(1,'P+P',4,3)
Hm.update(1,'P-P',4,0)
Hm.update(1,'P-P',4,1)
Hm.update(1,'P-P',4,2)
print("H-")
Hm.print()
print("\n")

Hz = 1/2 * Hp.comPlin(Hm)
print("Hz")
Hz.print()
print("\n")

print("[Hz,H+]")
error = Hz.comPlin(Hp)
error = error.simplify()
error.print()


