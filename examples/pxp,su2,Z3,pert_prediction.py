#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import pandas as pd

from comP import op_string_seq,op_string
from relations import su2_2d_relations,su3_8d_relations

# relation_object = su3_8d_relations()
relation_object = su2_2d_relations()

Hp = op_string_seq()
Hp.update(1,'P-P',3,2)
Hp.update(1,'P+P',3,0)
Hp.update(1,'P+P',3,1)

print("H+")
Hp.print()
print("\n")

Hm = op_string_seq()
Hm.update(1,'P+P',3,2)
Hm.update(1,'P-P',3,0)
Hm.update(1,'P-P',3,1)

print("H-")
Hm.print()

Hz =  1/2 * Hp.comPlin(Hm,relation_object)
print("\nHz")
print("\n")
Hz.print()
print("\n")

error = Hz.comPlin(Hp,relation_object)

print("Error Terms:")
from comP import collect_terms
collect_terms(error)

