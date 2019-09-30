#!/usr/bin/env python
# -*- coding: utf-8 -*-

import copy
import math
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import pandas as pd
from Search_functions import find_index_bisection
from relations import relations

class op_string:
    def __init__(self,coef,uc_string,period,loc):
        self.coef = coef
        self.string = uc_string
        self.length = len(self.string)
        self.period = period
        self.loc = loc

#dictionary of op strings. Def method for commuting sum of terms linearly
class op_string_seq:
    def __init__(self,string_seq=None):
        if string_seq is not None:
            self.string_seq = string_seq
            self.length = len(self.string_seq)
        else:
            self.string_seq = dict()
            self.length = 0

    #new turm to sum of op strings 
    def update(self,coef,uc_string,period,loc):
        self.string_seq[self.length] = op_string(coef,uc_string,period,loc)
        self.length += 1

    #scalar multiplication
    def __mul__(self,a):
        new_string_seq = copy.deepcopy(self.string_seq)
        for n in range(0,len(self.string_seq)):
            new_string_seq[n].coef = new_string_seq[n].coef*a
        return op_string_seq(new_string_seq)
    def __rmul__(self,a):
        new_string_seq = copy.deepcopy(self.string_seq)
        for n in range(0,len(self.string_seq)):
            new_string_seq[n].coef = new_string_seq[n].coef*a
        return op_string_seq(new_string_seq)

    def print(self):
        for n in range(0,len(self.string_seq)):
            print(self.string_seq[n].coef,self.string_seq[n].string,self.string_seq[n].period,self.string_seq[n].loc)

    #[sum_n A_n , sum_m B_m] = sum_{n,m} [A_n,B_m]
    def comPlin(self,B):
        string_seq = dict()
        c = 0
        for n in range(0,len(self.string_seq)):
            for m in range(0,len(B.string_seq)):
                term_strings = comP(self.string_seq[n],B.string_seq[m])
                for u in range(0,len(term_strings.string_seq)):
                    string_seq[c] = term_strings.string_seq[u]
                    c += 1
        return op_string_seq(string_seq)

#commute two op strings
def comP(A,B):
    #array with unit cell of A/B indices (eg 2n+1,2n+2,2n+3->1,2,3)
    A_loc = np.arange(A.loc,A.loc+A.length)
    B_loc = np.arange(B.loc,B.loc+B.length)

    #RHS indices = B.period*n+B.loc
    #find n_vals=[...] such that end of B string touch A string
    max_found = 0
    counter=0
    while max_found == 0:
        #B.period*n+B.loc = (last A site) - counter
        n_max = (A_loc[np.size(A_loc)-1]-B.loc-counter)/B.period
        if n_max.is_integer() == True:
            max_found = 1
        else:
            counter += 1
    min_found = 0
    counter=0
    while min_found == 0:
        #B.period*n+(last B site) = (first A site) + counter
        n_min = (A.loc-B_loc[np.size(B_loc)-1]+counter)/B.period
        if n_min.is_integer() == True:
            min_found = 1
        else:
            counter += 1
    n_vals = np.arange(n_min,n_max+1)

    #find all overlapping RH terms, store loc in 2d array
    rhs_loc = np.zeros((np.size(n_vals),B.length))
    for n in range(0,np.size(n_vals,axis=0)):
        rhs_loc[n] = np.arange(B.period*n_vals[n]+B.loc,B.period*n_vals[n]+B.loc+B.length)
            
    #for each rhs loc form "sandwich" of overlapping sites,store each term in nested dictionary
    #eg strings like "P(+P)(-P)P = P+-P"
    rhs_commutes = dict()
    new_coef = np.zeros(np.size(n_vals))
    for n in range(0,np.size(rhs_loc,axis=0)):
        min_loc = int(np.min(np.append(A_loc,rhs_loc[n])))
        max_loc = int(np.max(np.append(A_loc,rhs_loc[n])))

        rhs_commutes[n] = dict()
        for m in range(min_loc,max_loc+1):
            temp = ''
            if m in rhs_loc[n] and m in A_loc:
                A_loc_index = find_index_bisection(m,A_loc)
                rhs_index = find_index_bisection(m,rhs_loc[n])
                temp += A.string[A_loc_index]
                temp += B.string[rhs_index]
                rhs_commutes[n][m] =  temp
            elif m in rhs_loc[n]:
                rhs_index = find_index_bisection(m,rhs_loc[n])
                temp = B.string[rhs_index]
                rhs_commutes[n][m] =  temp
            else:
                A_loc_index = find_index_bisection(m,A_loc)
                temp = A.string[A_loc_index]
                rhs_commutes[n][m] =  temp
        new_coef[n] = A.coef * B.coef

    #replace known products from relation class
    for n in range(0,np.size(rhs_loc,axis=0)):
        #identify all pairs of producted operators
        keys = list(rhs_commutes[n])
        paired_keys = []
        for m in range(0,np.size(keys,axis=0)):
            if len(rhs_commutes[n][keys[m]]) > 1:
                paired_keys = np.append(paired_keys,keys[m])
        #replace any PP pairs with P
        for m in range(0,np.size(paired_keys,axis=0)):
            if rhs_commutes[n][paired_keys[m]] == "PP":
                rhs_commutes[n][paired_keys[m]] = "P"
        #identify remaining pairs
        paired_keys = []
        for m in range(0,np.size(keys,axis=0)):
            if len(rhs_commutes[n][keys[m]]) > 1:
                paired_keys = np.append(paired_keys,keys[m])
            
        #if no. paired keys == 1 can do regular commutator
        if np.size(paired_keys) == 1:
            simplified_commutator = relations.commutator(rhs_commutes[n][paired_keys[0]])
            new_coef[n] = new_coef[n] * simplified_commutator.coef
            rhs_commutes[n][paired_keys[0]] = simplified_commutator.string
        else:
            #AB term
            AB_coef = np.copy(new_coef[n])
            simplified_products_AB = dict()
            for m in range(0,np.size(paired_keys,axis=0)):
                simplified_products_AB[m] = relations.product(rhs_commutes[n][paired_keys[m]])
            for m in range(0,len(simplified_products_AB)):
                AB_coef = AB_coef * simplified_products_AB[m].coef

            #BA term (flip string)
            BA_coef = np.copy(new_coef[n])
            simplified_products_BA = dict()
            for m in range(0,np.size(paired_keys,axis=0)):
                simplified_products_BA[m] = relations.product(rhs_commutes[n][paired_keys[m]][::-1])
            for m in range(0,len(simplified_products_BA)):
                BA_coef = BA_coef * simplified_products_BA[m].coef

            if AB_coef != 0 and BA_coef != 0:
                #change first term to be AB, append BA term with negative coefficient
                for m in range(0,len(simplified_products_AB)):
                    rhs_commutes[n][paired_keys[m]] = simplified_products_AB[m].string
                new_coef[n] = AB_coef
                new_index = len(rhs_commutes)
                rhs_commutes[new_index] = copy.copy(rhs_commutes[n])
                for m in range(0,len(simplified_products_BA)):
                    rhs_commutes[new_index][paired_keys[m]] = simplified_products_BA[m].string
                new_coef = np.append(new_coef,-BA_coef)

            elif AB_coef != 0:
                for m in range(0,len(simplified_products_AB)):
                    rhs_commutes[n][paired_keys[m]] = simplified_products_AB[m].string
                new_coef[n] = AB_coef
            elif BA_coef != 0:
                for m in range(0,len(simplified_products_BA)):
                    rhs_commutes[n][paired_keys[m]] = simplified_products_BA[m].string
                new_coef[n] = -BA_coef
            else:
                new_coef[n] = 0

    #append simplified site ops to form strings, only considering those with coef>0
    rhs_strings = dict()
    new_loc = []
    new_coef_trimmed = []
    c=0
    for n in range(0,len(rhs_commutes)):
        if np.abs(new_coef[n]) > 1e-5: #!=0
            keys = list(rhs_commutes[n].keys())
            rhs_strings[c] = rhs_commutes[n][keys[0]]
            for m in range(1,np.size(keys)):
                if len(rhs_commutes[n][keys[m]])>1:
                    temp = "("
                    temp += rhs_commutes[c][keys[m]]
                    temp += ")"
                    rhs_strings[c] += temp
                else:
                    rhs_strings[c] += rhs_commutes[n][keys[m]]
            loc = np.min(list(rhs_commutes[n].keys()))
            new_loc = np.append(new_loc,loc)
            new_coef_trimmed = np.append(new_coef_trimmed,new_coef[n])
            c += 1
    #simplify final terms - cancel any same terms/replace commutators (eg P(+-)P - P(-+)P = 2*PZP)
    #dictionary storing string + coef using hash of string as key (so coef sum can be updated, unique key)
    new_terms = dict()
    new_term_coef = dict()
    keys = list(rhs_strings.keys())
    #init
    for n in range(0,np.size(keys,axis=0)):
        new_term_coef[hash(rhs_strings[keys[n]])]  = 0
        new_terms[hash(rhs_strings[keys[n]])] = rhs_strings[keys[n]]
    #update coef
    for n in range(0,np.size(keys,axis=0)):
        new_term_coef[hash(rhs_strings[keys[n]])]  += new_coef_trimmed[n]

    #form new op_string_seq
    keys = list(new_terms.keys())
    non_zero_strings = dict()
    non_zero_coef = dict()
    non_zero_loc = []
    c=0
    string_seq = dict()
    for n in range(0,np.size(keys,axis=0)):
        if np.abs(new_term_coef[keys[n]])>1e-5:
            string_seq[c] = op_string(new_term_coef[keys[n]],new_terms[keys[n]],A.period,new_loc[n]%A.period)
            c += 1
    return op_string_seq(string_seq)
