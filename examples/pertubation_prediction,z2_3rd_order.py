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

Hp.update(1,'PP+P',2,0,variables=["x"],variable_orders=[1])
Hp.update(1,'P+PP',2,1,variables=["x"],variable_orders=[1])
Hp.update(1,'P-PP',2,0,variables=["x"],variable_orders=[1])
Hp.update(1,'PP-P',2,1,variables=["x"],variable_orders=[1])

Hp.update(1,'PZP+P',2,1,variables=["y1"],variable_orders=[1])
Hp.update(1,'P+PZP',2,1,variables=["y1"],variable_orders=[1])
Hp.update(1,'PZP-P',2,0,variables=["y1"],variable_orders=[1])
Hp.update(1,'P-PZP',2,0,variables=["y1"],variable_orders=[1])

Hp.update(1,'P+PPP',2,1,variables=["y2"],variable_orders=[1])
Hp.update(1,'PPP+P',2,1,variables=["y2"],variable_orders=[1])
Hp.update(1,'P-PPP',2,0,variables=["y2"],variable_orders=[1])
Hp.update(1,'PPP-P',2,0,variables=["y2"],variable_orders=[1])

Hp.update(1,'PP+PP',2,0,variables=["y3"],variable_orders=[1])
Hp.update(1,'PP-PP',2,1,variables=["y3"],variable_orders=[1])

Hp.update(1,'PP+PZP',2,0,variables=["y4"],variable_orders=[1])
Hp.update(1,'PZP+PP',2,1,variables=["y4"],variable_orders=[1])
Hp.update(1,'PP-PZP',2,1,variables=["y4"],variable_orders=[1])
Hp.update(1,'PZP-PP',2,0,variables=["y4"],variable_orders=[1])

Hp.update(1,'PPP+PP',2,1,variables=["y5"],variable_orders=[1])
Hp.update(1,'PP+PPP',2,0,variables=["y5"],variable_orders=[1])
Hp.update(1,'PPP-PP',2,0,variables=["y5"],variable_orders=[1])
Hp.update(1,'PP-PPP',2,1,variables=["y5"],variable_orders=[1])

Hp.update(1,'P+PZPP',2,1,variables=["y6"],variable_orders=[1])
Hp.update(1,'PPZP+P',2,0,variables=["y6"],variable_orders=[1])
Hp.update(1,'PPZP-P',2,1,variables=["y6"],variable_orders=[1])
Hp.update(1,'P-PZPP',2,0,variables=["y6"],variable_orders=[1])

Hp.update(1,'PPPP+P',2,0,variables=["y7"],variable_orders=[1])
Hp.update(1,'P+PPPP',2,1,variables=["y7"],variable_orders=[1])
Hp.update(1,'PPPP-P',2,1,variables=["y7"],variable_orders=[1])
Hp.update(1,'P-PPPP',2,0,variables=["y7"],variable_orders=[1])

Hp.update(1,'PP+PZPP',2,0,variables=["y8"],variable_orders=[1])
Hp.update(1,'PPZP+PP',2,0,variables=["y8"],variable_orders=[1])
Hp.update(1,'PP-PZPP',2,1,variables=["y8"],variable_orders=[1])
Hp.update(1,'PPZP-PP',2,1,variables=["y8"],variable_orders=[1])

print("H+")
Hp.print()
print("\n")

Hm = op_string_seq()
Hm.update(1,'P-P',2,1)
Hm.update(1,'P+P',2,0)

Hm.update(1,'PP-P',2,0,variables=["x"],variable_orders=[1])
Hm.update(1,'P-PP',2,1,variables=["x"],variable_orders=[1])
Hm.update(1,'P+PP',2,0,variables=["x"],variable_orders=[1])
Hm.update(1,'PP+P',2,1,variables=["x"],variable_orders=[1])

Hm.update(1,'PZP-P',2,1,variables=["y1"],variable_orders=[1])
Hm.update(1,'P-PZP',2,1,variables=["y1"],variable_orders=[1])
Hm.update(1,'PZP+P',2,0,variables=["y1"],variable_orders=[1])
Hm.update(1,'P+PZP',2,0,variables=["y1"],variable_orders=[1])

Hm.update(1,'P-PPP',2,1,variables=["y2"],variable_orders=[1])
Hm.update(1,'PPP-P',2,1,variables=["y2"],variable_orders=[1])
Hm.update(1,'P+PPP',2,0,variables=["y2"],variable_orders=[1])
Hm.update(1,'PPP+P',2,0,variables=["y2"],variable_orders=[1])

Hm.update(1,'PP-PP',2,0,variables=["y3"],variable_orders=[1])
Hm.update(1,'PP+PP',2,1,variables=["y3"],variable_orders=[1])

Hm.update(1,'PP-PZP',2,0,variables=["y4"],variable_orders=[1])
Hm.update(1,'PZP-PP',2,1,variables=["y4"],variable_orders=[1])
Hm.update(1,'PP+PZP',2,1,variables=["y4"],variable_orders=[1])
Hm.update(1,'PZP+PP',2,0,variables=["y4"],variable_orders=[1])

Hm.update(1,'PPP-PP',2,1,variables=["y5"],variable_orders=[1])
Hm.update(1,'PP-PPP',2,0,variables=["y5"],variable_orders=[1])
Hm.update(1,'PPP+PP',2,0,variables=["y5"],variable_orders=[1])
Hm.update(1,'PP+PPP',2,1,variables=["y5"],variable_orders=[1])

Hm.update(1,'P-PZPP',2,1,variables=["y6"],variable_orders=[1])
Hm.update(1,'PPZP-P',2,0,variables=["y6"],variable_orders=[1])
Hm.update(1,'PPZP+P',2,1,variables=["y6"],variable_orders=[1])
Hm.update(1,'P+PZPP',2,0,variables=["y6"],variable_orders=[1])

Hm.update(1,'PPPP-P',2,0,variables=["y7"],variable_orders=[1])
Hm.update(1,'P-PPPP',2,1,variables=["y7"],variable_orders=[1])
Hm.update(1,'PPPP+P',2,1,variables=["y7"],variable_orders=[1])
Hm.update(1,'P+PPPP',2,0,variables=["y7"],variable_orders=[1])

Hm.update(1,'PP-PZPP',2,0,variables=["y8"],variable_orders=[1])
Hm.update(1,'PPZP-PP',2,0,variables=["y8"],variable_orders=[1])
Hm.update(1,'PP+PZPP',2,1,variables=["y8"],variable_orders=[1])
Hm.update(1,'PPZP+PP',2,1,variables=["y8"],variable_orders=[1])

print("H-")
Hm.print()

Hz =  1/2 * Hp.comPlin(Hm)
print("Hz calculated algebra calculated")
print(Hz.length)
error = Hz.comPlin(Hp)
print("Lie algebra calculated")
print(error.length)

from term_structures import string_term,term_same_coef
            
# collect same terms, ie same string, period and pos
same_terms = dict()
for n in range(0,len(error.string_seq)):
    hash_key = error.string_seq[n].string
    hash_key += str(error.string_seq[n].period)
    hash_key += str(error.string_seq[n].loc)
    term_hash = hash(hash_key)

    if term_hash not in list(same_terms.keys()):
        same_terms[term_hash] = string_term()
    same_terms[term_hash].add_term(error.string_seq[n])

#group differint string terms with same coef
keys = list(same_terms.keys())
found = []
same_coef = dict()
c=0
for n in range(0,np.size(keys,axis=0)):
    if n not in found:
        same_coef[c] = term_same_coef(same_terms[keys[n]])
        found = np.append(found,n)
        for m in range(n,np.size(keys,axis=0)):
            if m not in found:
                init_length = same_coef[c].length
                same_coef[c].add_term(same_terms[keys[m]])
                new_length = same_coef[c].length
                if init_length != new_length:
                    found = np.append(found,m)
        c += 1

#filter terms which have single spin flips and print perts
print("\n")
print("Possible pertubations (restricted to single spin flips)")
for n in range(0,len(same_coef)):
    string = list(same_coef[n].terms[0].string)
    spin_flips = 0
    for u in range(0,np.size(string,axis=0)):
        if string[u] == "+" or string[u] == "-":
            spin_flips += 1
    if spin_flips <=1:
        print(same_coef[n].coef)
        for m in range(0,same_coef[n].length):
            print(same_coef[n].terms[m].string,same_coef[n].terms[m].period,same_coef[n].terms[m].loc)
        paint("\n")