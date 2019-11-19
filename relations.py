#!/usr/bin/env python
# -*- coding: utf-8 -*-

import copy
import math
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import pandas as pd

def com(a,b):
    return np.dot(a,b)-np.dot(b,a)
def frob_product(A,B):
    return np.trace(np.dot(A,np.conj(np.transpose(B))))
def project_component(element,gen):
    return frob_product(gen,element)/frob_product(gen,gen)

#for sending output of commutators back to comP in a form it can unwrap
class simplified_string:
    def __init__(self,string,coef):
        self.string = string
        self.coef = coef
class simplified_strings:
    def __init__(self):
        self.entry = dict()
        self.length = 0
    def update(self,string,coef):
        self.entry[self.length] = simplified_string(string,coef)
        self.length += 1

#take products/commutators and expand them as a linear superposition of generators of a given lie algebra
class relations:
    def product(self,string_pair):
        site_ops1 = string_pair[0]
        site_ops2 = string_pair[1]

        product = np.dot(self.ops[site_ops1],self.ops[site_ops2] )

        new_string = simplified_strings()
        #check if product not a projector first
        is_a_projector = 0
        for n in range(0,np.size(self.proj_keys,axis=0)):
            overlap = frob_product(self.projectors[self.proj_keys[n]],product)/np.power(frob_product(product,product),0.5)
            if overlap == 1 or overlap == -1:
                coef = project_component(product,self.projectors[self.proj_keys[n]])
                new_string.update(self.proj_keys[n],coef)
                is_a_projector = 1
                break

        if is_a_projector == 0:
            #find expansion of product in terms of generators
            for n in range(0,np.size(self.keys,axis=0)):
                coef = project_component(product,self.generators[self.keys[n]])
                if np.abs(coef)>1e-5:
                    new_string.update(self.keys[n],coef)
        return new_string

    def commutator(self,string_pair):
        site_ops1 = string_pair[0]
        site_ops2 = string_pair[1]

        commutator = com(self.ops[site_ops1],self.ops[site_ops2])

        new_string = simplified_strings()
        #check if com not a projector first
        is_a_projector = 0
        for n in range(0,np.size(self.proj_keys,axis=0)):
            overlap = frob_product(self.projectors[self.proj_keys[n]],commutator)/np.power(frob_product(commutator,commutator),0.5)
            if overlap == 1 or overlap == -1:
                coef = project_component(commutator,self.projectors[self.proj_keys[n]])
                new_string.update(self.proj_keys[n],coef)
                is_a_projector = 1
                break

        if is_a_projector == 0:
            #find expansion of product in terms of generators
            for n in range(0,np.size(self.keys,axis=0)):
                coef = project_component(commutator,self.generators[self.keys[n]])
                if np.abs(coef)>1e-5:
                    new_string.update(self.keys[n],coef)
        return new_string

class su2_2d_relations(relations):
    def __init__(self):
        #generators defined sep to expand products/commutators as lin supp of the generators
        self.generators = dict()
        self.generators['+'] = np.array([[0,0],[1,0]])
        self.generators['-'] = np.array([[0,1],[0,0]])
        self.generators['Z'] = np.array([[-1/2,0],[0,1/2]])
        self.generators['I'] = np.array([[1,0],[0,1]])
        self.keys = list(self.generators.keys())

        #all poss ops (including projectors P etc)
        self.ops = dict()
        self.ops['+'] = np.array([[0,0],[1,0]])
        self.ops['-'] = np.array([[0,1],[0,0]])
        self.ops['Z'] = np.array([[-1/2,0],[0,1/2]])
        self.ops['I'] = np.array([[1,0],[0,1]])

        self.ops['P'] = np.array([[1,0],[0,0]])
        self.ops['Q'] = np.array([[0,0],[0,1]])

        #projector ops - store these terms sep to check by hand if [a,b], ab == projector 
        #cant include in linear superposition of group generators as not orthonormal wrt frob norm
        self.projectors = dict()
        self.projectors['P'] = np.array([[1,0],[0,0]])
        self.projectors['Q'] = np.array([[0,0],[0,1]])
        self.proj_keys = list(self.projectors.keys())

class su3_8d_relations(relations):
    def __init__(self):
        root3 = np.power(3,0.5)
        self.generators = dict()
        self.generators['I'] = np.array([[0,1,0],[0,0,0],[0,0,0]])
        self.generators['i'] = np.array([[0,0,0],[1,0,0],[0,0,0]])
        self.generators['K'] = np.array([[0,0,1],[0,0,0],[0,0,0]])
        self.generators['k'] = np.array([[0,0,0],[0,0,0],[1,0,0]])
        self.generators['L'] = np.array([[0,0,0],[0,0,1],[0,0,0]])
        self.generators['l'] = np.array([[0,0,0],[0,0,0],[0,1,0]])
        self.generators['M'] = np.array([[1/2,0,0],[0,-1/2,0],[0,0,0]])
        self.generators['D'] = np.array([[1/(2*root3),0,0],[0,1/(2*root3),0],[0,0,-2/(2*root3)]])
        self.generators['e'] = np.array([[1,0,0],[0,1,0],[0,0,1]])
        self.keys = list(self.generators.keys())

        #all poss ops (including projectors P etc)
        self.ops = dict()
        self.ops['I'] = np.array([[0,1,0],[0,0,0],[0,0,0]])
        self.ops['i'] = np.array([[0,0,0],[1,0,0],[0,0,0]])
        self.ops['K'] = np.array([[0,0,1],[0,0,0],[0,0,0]])
        self.ops['k'] = np.array([[0,0,0],[0,0,0],[1,0,0]])
        self.ops['L'] = np.array([[0,0,0],[0,0,1],[0,0,0]])
        self.ops['l'] = np.array([[0,0,0],[0,0,0],[0,1,0]])
        self.ops['M'] = np.array([[1/2,0,0],[0,-1/2,0],[0,0,0]])
        self.ops['D'] = np.array([[1/(2*root3),0,0],[0,1/(2*root3),0],[0,0,-2/(2*root3)]])
        self.ops['e'] = np.array([[1,0,0],[0,1,0],[0,0,1]])

        self.ops['P'] = np.array([[1,0,0],[0,0,0],[0,0,0]])
        self.ops['Q'] = np.array([[0,0,0],[0,1,0],[0,0,0]])
        self.ops['T'] = np.array([[0,0,0],[0,0,0],[0,1,0]])

        #projector ops - store these terms sep to check by hand if [a,b], ab == projector 
        #cant include in linear superposition of group generators as not orthonormal wrt frob norm
        self.projectors = dict()
        self.projectors['P'] = np.array([[1,0,0],[0,0,0],[0,0,0]])
        self.projectors['Q'] = np.array([[0,0,0],[0,1,0],[0,0,0]])
        self.projectors['T'] = np.array([[0,0,0],[0,0,0],[0,1,0]])
        self.proj_keys = list(self.projectors.keys())
