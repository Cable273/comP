#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import pandas as pd

from comP import op_string_seq,op_string

Hp = op_string_seq()
Hp.update(1,'P+P',2,1)
Hp.update(1,'P-P',2,0)
Hp.update(1,'PP+P',2,0,variables=["x"],variable_orders = [1])
Hp.update(1,'P+PP',2,1,variables=["x"],variable_orders = [1])
Hp.update(1,'P-PP',2,0,variables=["x"],variable_orders = [1])
Hp.update(1,'PP-P',2,1,variables=["x"],variable_orders = [1])
print("H+")
Hp.print()
print("\n")

Hm = op_string_seq()
Hm.update(1,'P+P',2,0)
Hm.update(1,'P-P',2,1)
Hm.update(1,'PP-P',2,0,variables=["x"],variable_orders = [1])
Hm.update(1,'P-PP',2,1,variables=["x"],variable_orders = [1])
Hm.update(1,'P+PP',2,0,variables=["x"],variable_orders = [1])
Hm.update(1,'PP+P',2,1,variables=["x"],variable_orders = [1])
print("H-")
Hm.print()
print("\n")

Hz = 1/2 * Hp.comPlin(Hm)
print("Hz")
Hz.print()
print("\n")

print("[Hz,H+]")
error = Hz.comPlin(Hp)
error.print()


