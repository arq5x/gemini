#!/usr/bin/env python

import sys


class pedformat:

    def __init__(self, fields):

        self.fields = fields[:]
        self.family = self._validate_field(fields[0])
        self.name = self._validate_field(fields[1])
        self.paternal = self._validate_field(fields[2])
        self.maternal = self._validate_field(fields[3])
        self.sex = self._validate_field(fields[4])
        self.phenotype = self._validate_field(fields[5])
        
        if len(fields) > 6:
            self.ethnicity = self._validate_field(fields[6])
        else:
            self.ethnicity = None

    def _validate_field(self, field):
        if len(field) > 0:
            return field
        else:
            return None

    def __str__(self):
        return ",".join([self.family, self.name, self.paternal, self.maternal, self.sex, self.phenotype, self.ethnicity])
