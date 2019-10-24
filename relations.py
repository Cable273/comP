#!/usr/bin/env python
# -*- coding: utf-8 -*-

import copy
import math
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import pandas as pd

class simplified_string:
    def __init__(self,string,coef):
        self.string = string
        self.coef = coef

class relations:
    def product(string):
        if string == "+P":
            return simplified_string("+",1)
        elif string == "P+":
            return simplified_string("I",0)

        elif string == "P-":
            return simplified_string("-",1)
        elif string == "-P":
            return simplified_string("I",0)

        elif string == "PP":
            return simplified_string("P",1)

        elif string == "-+":
            return simplified_string("P",1)
        elif string == "+-":
            return simplified_string("Q",1)

        elif string == "PZ":
            return simplified_string("P",-1/2)
        elif string == "ZP":
            return simplified_string("P",-1/2)

        elif string == "Z+":
            return simplified_string("+",1/2)
        elif string == "+Z":
            return simplified_string("+",-1/2)
        elif string == "Z-":
            return simplified_string("-",-1/2)
        elif string == "-Z":
            return simplified_string("-",1/2)

        elif string == "++":
            return simplified_string("I",0)
        elif string == "--":
            return simplified_string("I",0)
        elif string == "ZZ":
            return simplified_string("Z",1/2)
        elif string == "QP":
            return simplified_string("I",0)
        elif string == "PQ":
            return simplified_string("I",0)

        elif string == "Q+":
            return simplified_string("+",1)
        elif string == "+Q":
            return simplified_string("+",-1)
        elif string == "-Q":
            return simplified_string("-",1)
        elif string == "Q-":
            return simplified_string("-",-1)

        elif string == "QZ":
            return simplified_string("Q",1/2)
        elif string == "ZQ":
            return simplified_string("Q",1/2)

        elif string == "QQ":
            return simplified_string("Q",1)
        else:
            print(string)
            print("Product not known, add to product relations \n")

        # else:
            # return simplified_string(string,1)

    def commutator(string):
        if string == "+-":
            return simplified_string("Z",2)
        elif string == "-+":
            return simplified_string("Z",-2)
        elif string == "Z+":
            return simplified_string("+",1)
        elif string == "+Z":
            return simplified_string("+",-1)
        elif string == "Z-":
            return simplified_string("-",-1)
        elif string == "-Z":
            return simplified_string("-",1)
        elif string == "PP":
            return simplified_string("I",0)
        elif string == "PP":
            return simplified_string("I",0)
        elif string == "+P":
            return simplified_string("+",1)
        elif string == "P+":
            return simplified_string("+",-1)
        elif string == "-P":
            return simplified_string("-",-1)
        elif string == "P-":
            return simplified_string("-",1)
        elif string == "ZP":
            return simplified_string("I",0)
        elif string == "++":
            return simplified_string("I",0)
        elif string == "--":
            return simplified_string("I",0)

        elif string == "Q+":
            return simplified_string("+",1)
        elif string == "+Q":
            return simplified_string("+",-1)
        elif string == "Q-":
            return simplified_string("-",-1)
        elif string == "-Q":
            return simplified_string("-",1)
        elif string == "PQ":
            return simplified_string("I",0)
        elif string == "QP":
            return simplified_string("I",0)
        elif string == "ZZ":
            return simplified_string("I",0)
        elif string == "PZ":
            return simplified_string("I",0)
        elif string == "ZP":
            return simplified_string("I",0)
        elif string == "ZQ":
            return simplified_string("I",0)
        elif string == "QZ":
            return simplified_string("I",0)

        elif string == "QQ":
            return simplified_string("I",0)
        else:
            print(string)
            print("Commutator not known, add to commutator relations \n")

