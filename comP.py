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
from progressbar import ProgressBar
from copy import deepcopy

class op_string:
    def __init__(self,coef,uc_string,period,loc,variables=[],variable_orders=[]):
        self.coef = coef
        self.string = uc_string
        self.length = len(self.string)
        self.period = period
        self.loc = loc
        self.variables = variables
        self.variable_orders = variable_orders
        self.update_hash()

    #form unique hash to identify op_string type (to sum terms in an op_string_seq)
    #hash(string+period_loc+variables+order)
    def update_hash(self):
        string_to_hash = self.string
        string_to_hash += str(self.period)
        string_to_hash += str(self.loc)
        string_to_hash += str(self.variables)
        string_to_hash += str(self.variable_orders)
        key = hash(string_to_hash)
        self.key = key

    #multiply variables and keep track of the order (power) of each variable
    def multiply_variables(self,B):
        if np.size(self.variables)>0 and np.size(B.variables)>0:
            #find any shared variables
            shared_A_indices = []
            shared_B_indices = []
            for n in range(0,np.size(self.variables,)):
                for m in range(0,np.size(B.variables)):
                    if self.variables[n] == B.variables[m]:
                        shared_A_indices = np.append(shared_A_indices,n)
                        shared_B_indices = np.append(shared_B_indices,m)
            #sum there orders
            new_orders = np.zeros(np.size(shared_A_indices))
            for n in range(0,np.size(new_orders,axis=0)):
                new_orders[n] = self.variable_orders[int(shared_A_indices[n])] + B.variable_orders[int(shared_B_indices[n])]
            #form new variable array and new order array
            variables = []
            orders = []
            #append shared variables and there order
            for n in range(0,np.size(shared_A_indices,axis=0)):
                variables = np.append(variables,self.variables[int(shared_A_indices[n])])
                orders = np.append(orders,new_orders[n])
            #loop through self and append any variables not shared
            for n in range(0,np.size(self.variables,axis=0)):
                if n not in shared_A_indices:
                    variables = np.append(variables,self.variables[n])
                    orders = np.append(orders,self.variable_orders[n])
            #loop through B and append any variables not shared
            for n in range(0,np.size(B.variables,axis=0)):
                if n not in shared_B_indices:
                    variables = np.append(variables,B.variables[n])
                    orders = np.append(orders,B.variable_orders[n])

            #need a unique way of sorting variables so same terms have same ordering
            #hash each variable and then order the hashes
            var_hash = np.zeros(np.size(variables))
            for n in range(0,np.size(variables,axis=0)):
                var_hash[n] = hash(variables[n])
            indices = np.arange(0,np.size(var_hash))
            var_hash_sorted,indices_sorted = (list(t) for t in zip(*sorted(zip(var_hash,indices))))
            final_variables = []
            final_orders = []
            for n in range(0,np.size(indices_sorted,axis=0)):
                final_variables = np.append(final_variables,variables[indices_sorted[n]])
                final_orders = np.append(final_orders,orders[indices_sorted[n]])
            final_orders= final_orders.astype(int)
            return final_variables,final_orders
        #just return A (multiply by 1)
        elif np.size(self.variables)>0: 
            return self.variables,self.variable_orders
        #just return B (multiply by 1)
        elif np.size(B.variables)>0: 
            return B.variables,B.variable_orders
        else:
            return [],[]

#dictionary of op strings. Def method for commuting sum of terms linearly
class op_string_seq:
    def __init__(self,string_seq=None):
        if string_seq is not None:
            self.string_seq = string_seq
            self.length = len(self.string_seq)
        else:
            self.string_seq = dict()
            self.length = 0

    #new term to sum of op strings 
    def update(self,coef,uc_string,period,loc,variables=[],variable_orders=[]):
        self.string_seq[self.length] = op_string(coef,uc_string,period,loc,variables = variables, variable_orders = variable_orders)
        self.length += 1

    def comPlin(self,B,relation_object):
        string_seq = dict()
        c = 0
        pbar = ProgressBar(maxval=self.length*B.length)
        pbar.start()
        pbar_counter = 0
        for n in range(0,len(self.string_seq)):
            for m in range(0,len(B.string_seq)):
                term_strings = comP(self.string_seq[n],B.string_seq[m],relation_object)
                pbar_counter += 1
                pbar.update(pbar_counter)
                for u in range(0,len(term_strings.string_seq)):
                    string_seq[c] = term_strings.string_seq[u]
                    c += 1
        new_string_seq = op_string_seq(string_seq)
        new_string_seq = new_string_seq.simplify()
        new_string_seq = new_string_seq.reorder()
        return new_string_seq

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
            if np.size(self.string_seq[n].variables)==0:
                print(self.string_seq[n].coef,self.string_seq[n].string,self.string_seq[n].period,self.string_seq[n].loc)
            else:
                var_string = ''
                for m in range(0,np.size(self.string_seq[n].variables,axis=0)):
                    var_string += self.string_seq[n].variables[m]
                    var_string += '^'
                    var_string += str(self.string_seq[n].variable_orders[m])
                    var_string += " "
                print(self.string_seq[n].coef,self.string_seq[n].string,self.string_seq[n].period,self.string_seq[n].loc,var_string)

    # sum any repeated terms and form new op string seq with reduced term count
    def simplify(self):
        # find terms which share a hash (same term)
        shared_term_loc = dict()
        pbar=ProgressBar()
        for n in pbar(range(0,len(self.string_seq))):
            for m in range(n,len(self.string_seq)):
                if self.string_seq[n].key == self.string_seq[m].key:
                    key = self.string_seq[n].key
                    #init dict entry
                    if key not in list(shared_term_loc.keys()):
                        shared_term_loc[key] = np.array([n,m])
                    #append dict entry
                    else:
                        shared_term_loc[key] = np.append(shared_term_loc[key],np.array([n,m]))
        keys = list(shared_term_loc.keys())
        for n in range(0,np.size(keys,axis=0)):
            shared_term_loc[keys[n]] = np.unique(np.sort(shared_term_loc[keys[n]]))

        #sum shared terms
        string_seq = dict()
        c=0
        for n in range(0,np.size(keys,axis=0)):
            shared_indices = shared_term_loc[keys[n]]
            coef = 0
            for m in range(0,np.size(shared_indices,axis=0)):
                coef = coef + self.string_seq[shared_indices[m]].coef
            if np.abs(coef)>1e-5:
                string_seq[c] = copy.deepcopy(self.string_seq[shared_indices[0]])
                string_seq[c].coef = coef
                c+=1
        # terms not shared
        for n in range(0,len(self.string_seq)):
            if self.string_seq[n].key not in keys:
                string_seq[c] = self.string_seq[n]
                c+=1
        return op_string_seq(string_seq)
    
    #reorder first by string length, then by similar strings
    #reorder string seq so all similar strings (same string,period,loc, different variables) are printed in order when using print
    def reorder(self):
        #reorder by string length
        string_lengths = np.zeros(len(self.string_seq))
        for n in range(0,len(self.string_seq)):
            string_lengths[n]  = self.string_seq[n].length
        indices = np.arange(0,len(self.string_seq))
        ordered_lengths,reordered_indices = list((t) for t in zip(*sorted(zip(string_lengths,indices))))

        new_string_seq = dict()
        for n in range(0,np.size(reordered_indices,axis=0)):
            new_string_seq[n] = self.string_seq[reordered_indices[n]]
            
        reordered_indices = []
        for n in range(0,len(new_string_seq)):
            if n not in reordered_indices:
                reordered_indices = np.append(reordered_indices,n)
                for m in range(n,len(new_string_seq)):
                    if new_string_seq[n].string == new_string_seq[m].string:
                            if m not in reordered_indices:
                                reordered_indices = np.append(reordered_indices,m)
        final_string_seq = dict()
        for n in range(0,np.size(reordered_indices,axis=0)):
            final_string_seq[n] = new_string_seq[int(reordered_indices[n])]
        return op_string_seq(final_string_seq)

class commuted_string_term:
    def __init__(self,string,coef,loc):
        self.string = string
        self.coef = coef
        self.loc = loc
class commuted_string_seq:
    def __init__(self):
        self.terms = dict()
        self.no_terms = len(self.terms)
    def update(self,string_dict,new_coef):
        # if np.abs(new_coef) > 1e-5: #!=0
        keys = list(string_dict.keys())
        appended_string = string_dict[keys[0]]
        for m in range(1,np.size(keys)):
            if len(string_dict[keys[m]])>1:
                temp = "("
                temp += string_dict[keys[m]]
                temp += ")"
                appended_string += temp
            else:
                appended_string += string_dict[keys[m]]
        loc = np.min(list(string_dict.keys()))
        self.terms[self.no_terms] = commuted_string_term(appended_string,new_coef,loc)
        self.no_terms += 1
     
#commute two op strings
def comP(A,B,relation_object):
    #array with unit cell of A/B indices (eg 2n+1,2n+2,2n+3->1,2,3)
    A_loc = np.arange(A.loc,A.loc+A.length)
    B_loc = np.arange(B.loc,B.loc+B.length)

    #find n_vals=[...] such that end of B string touch A string
    #RHS indices = B.period*n+B.loc
    max_found = 0
    counter=0
    while max_found == 0:
        #B.period*n+B.loc = (last A site) - counter
        n_max = (A_loc[np.size(A_loc)-1]-B.loc-counter)/B.period
        if n_max.is_integer() == True:
            max_found = 1
            break
        else:
            counter += 1
    min_found = 0
    counter=0
    while min_found == 0:
        #B.period*n+(last B site) = (first A site) + counter
        n_min = (A.loc-B_loc[np.size(B_loc)-1]+counter)/B.period
        if n_min.is_integer() == True:
            min_found = 1
            break
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
    simplified_commutes = commuted_string_seq()
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
        #replace any QQ pairs with Q
        for m in range(0,np.size(paired_keys,axis=0)):
            if rhs_commutes[n][paired_keys[m]] == "QQ":
                rhs_commutes[n][paired_keys[m]] = "Q"

        #identify remaining pairs
        paired_keys = []
        for m in range(0,np.size(keys,axis=0)):
            if len(rhs_commutes[n][keys[m]]) > 1:
                paired_keys = np.append(paired_keys,keys[m])
            
        #if no. paired keys == 1 can do regular commutator
        if np.size(paired_keys) == 1:
            #desired syntax
            # simplified_commutators= relations.commutator(rhs_commutes[n][paired_keys[0]])
            simplified_commutators= relation_object.commutator(rhs_commutes[n][paired_keys[0]])
            for u in range(0,simplified_commutators.length):
                temp_string_dict = deepcopy(rhs_commutes[n])
                temp_string_dict[paired_keys[0]] = simplified_commutators.entry[u].string

                temp_coef = new_coef[n] * simplified_commutators.entry[u].coef
                simplified_commutes.update(temp_string_dict,temp_coef)
        else:
            AB_coef_orig = np.copy(new_coef[n])
            BA_coef_orig = np.copy(new_coef[n])

            simplified_products_AB = dict()
            simplified_products_BA = dict()
            for m in range(0,np.size(paired_keys,axis=0)):
                simplified_products_AB[m] = relation_object.product(rhs_commutes[n][paired_keys[m]])
                simplified_products_BA[m] = relation_object.product(rhs_commutes[n][paired_keys[m]][::-1])

            index_array_AB = []
            for m in range(0,len(simplified_products_AB)):
                temp = list(simplified_products_AB[m].entry.keys())
                index_array_AB.append(temp)
            from itertools import product
            indices_AB = np.array(list(product(*index_array_AB)))

            index_array_BA = []
            for m in range(0,len(simplified_products_BA)):
                temp = list(simplified_products_BA[m].entry.keys())
                index_array_BA.append(temp)
            from itertools import product
            indices_BA = np.array(list(product(*index_array_BA)))

            if np.size(indices_AB)>0:
                for u in range(0,np.size(indices_AB,axis=0)):
                    temp_string_dict = deepcopy(rhs_commutes[n])
                    paired_key_index = 0
                    temp_coef = AB_coef_orig
                    for v in range(0,np.size(indices_AB[u],axis=0)):
                        temp_string_dict[paired_keys[paired_key_index]] = simplified_products_AB[v].entry[indices_AB[u][v]].string
                        temp_coef = temp_coef * simplified_products_AB[v].entry[indices_AB[u][v]].coef
                        paired_key_index += 1
                    simplified_commutes.update(temp_string_dict,temp_coef)

            if np.size(indices_BA)>0:
                for u in range(0,np.size(indices_BA,axis=0)):
                    temp_string_dict = deepcopy(rhs_commutes[n])
                    paired_key_index = 0
                    temp_coef = BA_coef_orig
                    for v in range(0,np.size(indices_BA[u],axis=0)):
                        temp_string_dict[paired_keys[paired_key_index]] = simplified_products_BA[v].entry[indices_BA[u][v]].string
                        temp_coef = temp_coef * simplified_products_BA[v].entry[indices_BA[u][v]].coef
                        paired_key_index += 1
                    simplified_commutes.update(temp_string_dict,-temp_coef)

    rhs_strings = dict()
    new_loc = []
    new_coef_trimmed = []
    c=0
    for n in range(0,simplified_commutes.no_terms):
        if np.abs(simplified_commutes.terms[n].coef)>1e-5:
            rhs_strings[c] = simplified_commutes.terms[n].string
            new_loc = np.append(new_loc,simplified_commutes.terms[n].loc)
            new_coef_trimmed = np.append(new_coef_trimmed,simplified_commutes.terms[n].coef)
            c+=1

    #keep non zero terms
    keys = list(rhs_strings.keys())
    non_zero_strings = dict()
    non_zero_coef = dict()
    non_zero_loc = []
    c=0
    string_seq = dict()
    for n in range(0,np.size(keys,axis=0)):
        if np.abs(new_coef_trimmed[keys[n]])>1e-5:
            string_seq[c] = op_string(new_coef_trimmed[keys[n]],rhs_strings[keys[n]],A.period,new_loc[n]%A.period)
            c += 1

    #multiply out variables and orders (all terms will have same variables)
    new_variables,new_orders = A.multiply_variables(B)
    for n in range(0,len(string_seq)):
        string_seq[n].variables = new_variables
        string_seq[n].variable_orders = new_orders
        string_seq[n].update_hash()
    return op_string_seq(string_seq)
