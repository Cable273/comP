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
Hp.update(1,'PP-P',2,1,variables=["x"],variable_orders = [1])
Hp.update(1,'P-PP',2,0,variables=["x"],variable_orders = [1])

a = -1+np.power(2,0.5)
print("H+")
Hp.print()
print("\n")

Hm = op_string_seq()
Hm.update(1,'P+P',2,0)
Hm.update(1,'P-P',2,1)
Hm.update(1,'PP-P',2,0,variables=["x"],variable_orders = [1])
Hm.update(1,'P-PP',2,1,variables=["x"],variable_orders = [1])
Hm.update(1,'PP+P',2,1,variables=["x"],variable_orders = [1])
Hm.update(1,'P+PP',2,0,variables=["x"],variable_orders = [1])

print("H-")
Hm.print()
print("\n")

Hz =  1/2 * Hp.comPlin(Hm)
# print("Hz")
# Hz.print()
# print("\n")

# print("[Hz,H+]")
error = Hz.comPlin(Hp)
# error.print()

class string_term:
    def __init__(self):
        self.terms = dict()
        self.length = 0

    def add_term(self,op_string):
        self.terms[self.length] = op_string
        self.length += 1

        self.string = op_string.string
        self.period = op_string.period
        self.loc = op_string.loc

    def coef_eval(self):
        term_coef = ""
        for n in range(0,self.length):
            inserted_term = 0
            if np.size(self.terms[n].variables) == 0:
                term_coef += str(self.terms[n].coef)
                inserted_term = 1
            else:
                for m in range(0,np.size(self.terms[n].variables,axis=0)):
                    term_coef += str(self.terms[n].coef)
                    term_coef += self.terms[n].variables[m]
                    term_coef += "^"
                    term_coef += str(self.terms[n].variable_orders[m])
                    inserted_term = 1
            if inserted_term != 0:
                if n!= self.length-1:
                    term_coef += " + "
        return term_coef

            
# collect same terms
same_terms = dict()
for n in range(0,len(error.string_seq)):
    hash_key = error.string_seq[n].string
    hash_key += str(error.string_seq[n].period)
    hash_key += str(error.string_seq[n].loc)
    term_hash = hash(hash_key)

    if term_hash not in list(same_terms.keys()):
        same_terms[term_hash] = string_term()
    same_terms[term_hash].add_term(error.string_seq[n])

class term_same_coef:
    def __init__(self,init_term):
        self.coef = init_term.coef_eval()
        self.terms = dict()
        self.terms[0] = init_term
        self.length = 1

    def add_term(self,term):
        #check same coef first
        coef_compare = term.coef_eval()

        coef_compare_list = list(coef_compare)
        self_compare_list = list(self.coef)
        coef_compare_list.sort()
        self_compare_list.sort()

        if coef_compare_list == self_compare_list:
            self.terms[self.length] = term
            self.length += 1
        # else:
            # print("ERROR: Not same coef")

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
        print("\n")
