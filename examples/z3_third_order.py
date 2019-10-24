#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import pandas as pd

from comP import op_string_seq,op_string

Hp = op_string_seq()
Hp.update(1,'P-P',3,2)
Hp.update(1,'P+P',3,0)
Hp.update(1,'P+P',3,1)

Hp.update(1,'PP+P',3,2,variables=["x1"],variable_orders = [1])
Hp.update(1,'P+PP',3,1,variables=["x"],variable_orders = [1])
Hp.update(1,'P-PP',3,2,variables=["x1"],variable_orders = [1])
Hp.update(1,'PP-P',3,1,variables=["x1"],variable_orders = [1])

Hp.update(1,'PP-PP',3,1,variables=["x2"],variable_orders = [1])

Hp.update(1,'P-PZP',3,2,variables=["x3"],variable_orders = [1])
Hp.update(1,'PZP-P',3,0,variables=["x3"],variable_orders = [1])

Hp.update(1,'P-PZPP',3,2,variables=["x4"],variable_orders = [1])
Hp.update(1,'PPZP-P',3,2,variables=["x4"],variable_orders = [1])

Hp.update(1,'P+PP',3,0,variables=["x5"],variable_orders = [1])
Hp.update(1,'PP+P',3,0,variables=["x5"],variable_orders = [1])

Hp.update(1,'PP+PP',3,2,variables=["x6"],variable_orders = [1])
Hp.update(1,'PP+PP',3,0,variables=["x6"],variable_orders = [1])

Hp.update(1,'P+PPP',3,0,variables=["x7"],variable_orders = [1])
Hp.update(1,'PPP+P',3,2,variables=["x7"],variable_orders = [1])

Hp.update(1,'PZP+P',3,1,variables=["x8"],variable_orders = [1])
Hp.update(1,'P+PZP',3,1,variables=["x8"],variable_orders = [1])

Hp.update(1,'P+PZP',3,0,variables=["x9"],variable_orders = [1])
Hp.update(1,'PZP+P',3,2,variables=["x9"],variable_orders = [1])

Hp.update(1,'P+PQP',3,0,variables=["x10"],variable_orders = [1])
Hp.update(1,'PQP+P',3,2,variables=["x10"],variable_orders = [1])

Hp.update(1,'PP+PPP',3,2,variables=["x11"],variable_orders = [1])
Hp.update(1,'PPP+PP',3,2,variables=["x11"],variable_orders = [1])

Hp.update(1,'PP+PZP',3,2,variables=["x12"],variable_orders = [1])
Hp.update(1,'PZP+PP',3,2,variables=["x12"],variable_orders = [1])

Hp.update(1,'PP+PQP',3,2,variables=["x13"],variable_orders = [1])
Hp.update(1,'PQP+PP',3,2,variables=["x13"],variable_orders = [1])

Hp.update(1,'PP+PZPP',3,2,variables=["x14"],variable_orders = [1])
Hp.update(1,'PPZP+PP',3,1,variables=["x14"],variable_orders = [1])

Hp.update(1,'PP+PZPP',3,2,variables=["x15"],variable_orders = [1])
Hp.update(1,'PPZP+PP',3,1,variables=["x15"],variable_orders = [1])

Hp.update(1,'PPP+P',3,1,variables=["x16"],variable_orders = [1])
Hp.update(1,'P+PPP',3,1,variables=["x16"],variable_orders = [1])


print("H+")
Hp.print()
print("\n")

Hm = op_string_seq()
Hm.update(1,'P+P',3,2)
Hm.update(1,'P-P',3,0)
Hm.update(1,'P-P',3,1)
Hm.update(1,'PP-P',3,2,variables=["x"],variable_orders = [1])

Hm.update(1,'P-PP',3,1,variables=["x"],variable_orders = [1])
Hm.update(1,'P+PP',3,2,variables=["x"],variable_orders = [1])
Hm.update(1,'PP+P',3,1,variables=["x"],variable_orders = [1])

Hp.update(1,'PP+PP',3,1,variables=["x2"],variable_orders = [1])

Hp.update(1,'P+PZP',3,2,variables=["x3"],variable_orders = [1])
Hp.update(1,'PZP+P',3,0,variables=["x3"],variable_orders = [1])

Hp.update(1,'P+PZPP',3,2,variables=["x4"],variable_orders = [1])
Hp.update(1,'PPZP+P',3,2,variables=["x4"],variable_orders = [1])

Hp.update(1,'P-PP',3,0,variables=["x5"],variable_orders = [1])
Hp.update(1,'PP-P',3,0,variables=["x5"],variable_orders = [1])

Hp.update(1,'PP-PP',3,2,variables=["x6"],variable_orders = [1])
Hp.update(1,'PP-PP',3,0,variables=["x6"],variable_orders = [1])

Hp.update(1,'P-PPP',3,0,variables=["x7"],variable_orders = [1])
Hp.update(1,'PPP-P',3,2,variables=["x7"],variable_orders = [1])

Hp.update(1,'PZP-P',3,1,variables=["x8"],variable_orders = [1])
Hp.update(1,'P-PZP',3,1,variables=["x8"],variable_orders = [1])

Hp.update(1,'P-PZP',3,0,variables=["x9"],variable_orders = [1])
Hp.update(1,'PZP-P',3,2,variables=["x9"],variable_orders = [1])

Hp.update(1,'P-PQP',3,0,variables=["x10"],variable_orders = [1])
Hp.update(1,'PQP-P',3,2,variables=["x10"],variable_orders = [1])

Hp.update(1,'PP-PPP',3,2,variables=["x11"],variable_orders = [1])
Hp.update(1,'PPP-PP',3,2,variables=["x11"],variable_orders = [1])

Hp.update(1,'PP-PZP',3,2,variables=["x12"],variable_orders = [1])
Hp.update(1,'PZP-PP',3,2,variables=["x12"],variable_orders = [1])

Hp.update(1,'PP-PQP',3,2,variables=["x13"],variable_orders = [1])
Hp.update(1,'PQP-PP',3,2,variables=["x13"],variable_orders = [1])

Hp.update(1,'PP-PZPP',3,2,variables=["x14"],variable_orders = [1])
Hp.update(1,'PPZP-PP',3,1,variables=["x14"],variable_orders = [1])

Hp.update(1,'PP-PZPP',3,2,variables=["x15"],variable_orders = [1])
Hp.update(1,'PPZP-PP',3,1,variables=["x15"],variable_orders = [1])

Hp.update(1,'PPP-P',3,1,variables=["x16"],variable_orders = [1])
Hp.update(1,'P-PPP',3,1,variables=["x16"],variable_orders = [1])

print("H-")
Hm.print()
print("\n")

Hz =  1/2 * Hp.comPlin(Hm)
print("Hz")
Hz.print()
print("\n")

print("[Hz,H+]")
error = Hz.comPlin(Hp)
error.print()
