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
Hp.print()
print("\n")

Hm = op_string_seq()
Hm.update(1,'P+P',2,0)
Hm.update(1,'P-P',2,1)
Hm.print()
print("\n")

Hz = 1/2* Hp.comPlin(Hm)
Hz.print()
print("\n")

error = Hz.comPlin(Hp)
error.print()


