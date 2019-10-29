#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import pandas as pd

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

