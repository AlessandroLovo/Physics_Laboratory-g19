# -*- coding: utf-8 -*-
"""
Created on Fri May 11 14:48:22 2018

@author: athenis
"""

import shelve
import dumbdbm

def dumbdbm_shelve(filename,flag="c"):
    return shelve.Shelf(dumbdbm.open(filename,flag))

out_shelf=dumbdbm_shelve(".anlang.dumbdbm")
in_shelf=shelve.open(".anlang")

key_list=in_shelf.keys()
for key in key_list:
    out_shelf[key]=in_shelf[key]

out_shelf.close()
in_shelf.close()