#!/usr/bin/python
# -*- coding: utf-8 -*-
# author: Neo
import re


class atom_mass(object):
    sci_num = "\d\.\d+[eE][-\+]\d+"
    mua = 1.660539040 * 1e-27  # 1mua=1.66...kg

    def __init__(self, abundance="abundance", mass="mass_rmd.mas95"):
        class sub(object):
            pass

        self.amu = sub()
        self.kg = sub()
        f = open(abundance)
        self._abundance_file = f.readlines()
        f.close()
        f = open(mass)
        self._mass_file = f.readlines()
        f.close()
        self.__mass_attribute()
        self.__abund_attribute()

    def __str2num(self, string):
        try:
            return int(string)
        except ValueError:
            pass
        try:
            return float(string)
        except ValueError:
            return string.strip()

    def _parse_abund(self):
        sci_num = self.sci_num
        pattern = r'([\s\d]{3})\s([A-Z]+)\s+(\d+)\s+(%s)\s+(%s)\s+(%s)\s+([FM])\s+(stable|nostable)' % (
        sci_num, sci_num, sci_num)
        self.abund = []
        for line in self._abundance_file:
            match = re.findall(pattern, line)
            if match:
                atom = {i: self.__str2num(j) for i, j in \
                        zip(("A", "sym", "Z", "abund_noyau", 'abund_mass', 'mass', 'state', 'stablility'), match[0])}
                self.abund.append(atom)

    def _parse_mass(self):
        num = r"\d+\.?\d*"
        pattern = r'0?\s+-?\d+\s+(\d+)\s+\d+\s+(\d+)\s+([a-zA-Z]{1,2})\s.*?%s#?\s+%s#?\s*-?%s#?\s+%s#?\s*B[-\+]\s.{20}\s([\s\d]{3})\s*(%s)#?\s*%s#?\s+' \
                  % (num, num, num, num, num, num)
        self._mass = []
        for line in self._mass_file:
            match = re.findall(pattern, line)
            if match:
                A, Z, sym, mass_int, mass_d = match[0]
                mass = float(mass_int) + float(mass_d) * 1e-6
                atom = {i: self.__str2num(j) for i, j in \
                        zip(("Z", "A", "sym"), (A, Z, sym))}
                atom['mass'] = mass
                self._mass.append(atom)

    def __abund_attribute(self):
        if not hasattr(self, "abund"):
            self._parse_abund()
        if not hasattr(self.amu, "H1"):
            self.__mass_attribute()
        element = set()
        for isotope in self.abund:
            element.add(isotope["sym"])
        for atom in element:
            atomname = atom + "_NAT"
            atomname = atomname.upper()
            mass = 0.0
            for isotope in self.abund:
                if isotope['sym'] == atom:
                    mass += isotope["mass"] * isotope["abund_noyau"]
            self.amu.__dict__[atomname] = mass
            self.kg.__dict__[atomname] = mass * self.mua

    def __mass_attribute(self):
        if not hasattr(self, "_mass"):
            self._parse_mass()
        for atom in self._mass:
            atomname = atom['sym'] + str(atom['A'])
            atomname = atomname.upper()
            mass = atom['mass']
            self.amu.__dict__[atomname] = mass
            self.kg.__dict__[atomname] = mass * self.mua

    def __getitem__(self, key):
        key = key.replace('-', '_')
        key = key.upper()
        return self.amu.__dict__[key]

    def getKg(self, key):
        key = key.replace('-', '_')
        return self.kg.__dict__[key]
