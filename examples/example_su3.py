#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import pandas as pd

from comP import op_string_seq,op_string
from relations import su2_2d_relations,su3_8d_relations

relation_object = su3_8d_relations()

#broken su3 cartan weyl basis
Ip = op_string_seq()
Ip.update(1,'PIP',2,1)
Ip.update(-1,'PiP',2,0)

Kp = op_string_seq()
Kp.update(1,'PKP',2,1)
Kp.update(-1,'PkP',2,0)

Lp = op_string_seq()
Lp.update(1,'PLP',2,1)
Lp.update(-1,'PlP',2,0)

Im = op_string_seq()
Im.update(1,'PiP',2,1)
Im.update(-1,'PIP',2,0)

Km = op_string_seq()
Km.update(1,'PkP',2,1)
Km.update(-1,'PKP',2,0)

Lm = op_string_seq()
Lm.update(1,'PlP',2,1)
Lm.update(-1,'PLP',2,0)

#"diagonal" broken su(3) generators
I3 = 1/2 * Ip.comPlin(Im,relation_object)
g8 = 1/(2*np.power(3,0.5))*(Kp.comPlin(Km,relation_object) + Lp.comPlin(Lm,relation_object))

#broken lie algebra
Ip_root1 = I3.comPlin(Ip,relation_object)
Ip_root2 = g8.comPlin(Ip,relation_object)
Im_root1 = I3.comPlin(Im,relation_object)
Im_root2 = g8.comPlin(Im,relation_object)

Kp_root1 = I3.comPlin(Kp,relation_object)
Kp_root2 = g8.comPlin(Kp,relation_object)
Km_root1 = I3.comPlin(Km,relation_object)
Km_root2 = g8.comPlin(Km,relation_object)

Lp_root1 = I3.comPlin(Lp,relation_object)
Lp_root2 = g8.comPlin(Lp,relation_object)
Lm_root1 = I3.comPlin(Lm,relation_object)
Lm_root2 = g8.comPlin(Lm,relation_object)

from comP import collect_terms
print("\n")
print("[I3,Ip]")
collect_terms(Ip_root1)
print("\n")
print("[g8,Ip]")
collect_terms(Ip_root2)
print("\n")
print("[I3,Kp]")
collect_terms(Kp_root1)
print("\n")
print("[g8,Kp]")
collect_terms(Kp_root2)
print("\n")
print("[I3,Lp]")
collect_terms(Lp_root1)
print("\n")
print("[g8,Lp]")
collect_terms(Lp_root2)
