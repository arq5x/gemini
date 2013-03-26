#!/usr/bin/env python

import sys


class pedformat:

    def __init__(self, fields):

        self.fields = fields[:]
        self.family = fields[0]
        self.name = fields[1]
        self.paternal = fields[2]
        self.maternal = fields[3]
        self.sex = fields[4]
        self.phenotype = fields[5]
        self.ethnicity = fields[6]

    def __str__(self):
        return ",".join([self.family, self.name, self.paternal, self.maternal, self.sex, self.phenotype, self.ethnicity])
