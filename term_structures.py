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

    def gen_coef_list(self):
        coef_list = []
        for n in range(0,self.length):
            term_coef = ""
            if np.size(self.terms[n].variables) == 0:
                term_coef += str(self.terms[n].coef)
                inserted_term = 1
            else:
                for m in range(0,np.size(self.terms[n].variables,axis=0)):
                    term_coef += str(self.terms[n].coef)
                    term_coef += self.terms[n].variables[m]
                    term_coef += "^"
                    term_coef += str(self.terms[n].variable_orders[m])
            coef_list.append(term_coef)
        return coef_list
            

class term_same_coef:
    def __init__(self,init_term):
        self.coef = init_term.coef_eval()
        self.terms = dict()
        self.terms[0] = init_term
        self.length = 1
        self.init_term = init_term
        self.coef_list = self.terms[0].gen_coef_list()

    def check_same_coef(self,term):
        if len(term.terms) == len(self.init_term.terms):
            all_terms_identical = 1
            #scan through coef of new term
            coef_matched = []
            for n in range(0,len(term.terms)):
                coef = term.terms[n].coef
                variables = term.terms[n].variables
                variable_orders = term.terms[n].variable_orders
                #see if this coef is in initial term coef
                term_found = 0
                for m in range(0,len(self.init_term.terms)):
                    if variables == self.init_term.terms[m].variables and variable_orders == self.init_term.terms[m].variable_orders and coef == self.init_term.terms[m].coef:
                        term_found = 1
                        break

                if term_found == 0:
                    all_terms_identical = 0
                    break
            if all_terms_identical == 1:
                return True
            else:
                return False
        else:
            return False
                    

    def add_term(self,term):
        #check same coef first
        # coef_compare = term.coef_eval()
        # coef_compare_list = list(coef_compare)
        # self_compare_list = list(self.coef)
        # coef_compare_list.sort()
        # self_compare_list.sort()

        # compare_coef_list = term.gen_coef_list()

        # if coef_compare_list == self_compare_list:
        if self.check_same_coef(term):
            self.terms[self.length] = term
            self.length += 1

