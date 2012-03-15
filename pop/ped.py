#!/usr/bin/env python

import sys

class pedformat:
    
    def __init__(self,field):
        
        self.fields = field[:]
        self.family = field[0]
        self.name = field[1]
        self.paternal = field[2]
        self.maternal = field[3]
        self.sex = field[4]
        self.phenotype = field[5]
        self.ethnicity = field[6] 
    
    def __str__(self):
        return ",".join([self.family, self.name, self.paternal, self.maternal, self.sex, self.phenotype, self.ethnicity])    