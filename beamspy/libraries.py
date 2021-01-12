#!/usr/bin/env python
# -*- coding: utf-8 -*-

from collections import OrderedDict
from beamspy.auxiliary import order_composition_by_hill


class Adducts:
    def __init__(self, ion_mode=None, e=0.0005486):

        self.e = e
        if ion_mode == "pos":
            self.lib = OrderedDict()
        elif ion_mode == "neg":
            self.lib = OrderedDict()
        elif ion_mode is None:
            self.lib = OrderedDict()

    def add(self, name, mass, charge):
        self.lib[name] = OrderedDict([("mass", float(mass)), ("charge", int(charge))])
        self.lib = OrderedDict(sorted(self.lib.items(), key=lambda x: x[1]['mass']))

    def remove(self, name):
        if name == "*":
            self.lib = OrderedDict()
        else:
            if name in self.lib:
                self.lib.remove(name)
            else:
                raise IOError("Entry not in library: {}".format(name))

    def __str__(self):
        out = "Adducts in library\n"
        out += "-----------------\n"
        out += "name\texact_mass\n"
        for key in self.lib:
            out += "%s\t%s\n" % (key, self.lib[key])
        return out


class Isotopes:

    def __init__(self, ion_mode=None):

        self.ion_mode = ion_mode
        self.lib = [OrderedDict([("C", {"abundance": 100.0}), ("(13C)", {"abundance": 1.1}),
                                 ("mass_difference", 1.003355),
                                 ("charge", 1)]),
                    OrderedDict([("S", {"abundance": 100.0}), ("(34S)", {"abundance": 4.21}),
                                 ("mass_difference", 1.995796),
                                 ("charge", 1)])]

        if self.ion_mode == "pos":
            self.lib.append(OrderedDict([("K", {"abundance": 100.0}), ("(41K)", {"abundance": 6.73}),
                                         ("mass_difference", 1.998117), ("charge", 1)]))
            #self.lib.append(OrderedDict([("(6Li)", {"abundance": 7.42}), ("Li", {"abundance": 1.0}), ("mass_difference", 1.000882)]))

        elif self.ion_mode == "neg":
            self.lib.append(OrderedDict([("Cl", {"abundance": 100.0}), ("(37Cl)", {"abundance": 24.23}),
                                         ("mass_difference", 1.997050), ("charge", 1)]))
        self.lib = sorted(self.lib, key=lambda k: k['mass_difference'])

    def add(self, label_x, label_y, mx_abundance, my_abundance, mass_difference, charge):
        self.lib.append(OrderedDict([(label_x, {"abundance": float(mx_abundance)}),
                                     (label_y, {"abundance": float(my_abundance)}),
                                     ("mass_difference", float(mass_difference)),
                                     ("charge", int(charge))]))
        self.lib = sorted(self.lib, key=lambda k: k['mass_difference'])

    def remove(self, label_x="*", label_y="*"):
        if label_x == "*" and label_y == "*":
            self.lib = []
        else:
            for item in self.lib:
                if label_x in item or label_y in item:
                    self.lib.remove(item)
            else:
                print("Entry not in library")

    def __str__(self):
        out = "Isotopes in library:\n"
        out += "--------------------------------------------\n"
        out += "label_x\tlabel_y\tmass_difference\tcharge\tabundance_x\tabundance_y\n"
        for item in self.lib:
            label_x = list(item.items())[0][0]
            label_y = list(item.items())[1][0]
            out += "{}\t{}\t{}\t{}\t{}\t{}\n".format(label_x, label_y,
                                                 item["mass_difference"],
                                                 item["charge"],
                                                 item[label_x]["abundance"], item[label_y]["abundance"])
        return out


class NeutralLosses:

    def __init__(self):
        self.lib = []

    def add(self, label, mass_difference):
        self.lib.append(OrderedDict([("label", label), ("mass_difference", mass_difference)]))
        self.lib = sorted(self.lib, key=lambda k: k['mass_difference'])

    def remove(self, label="*"):
        if label == "*":
            self.lib = []
        else:
            for item in self.lib:
                if label in self.lib:
                    self.lib.remove(item)
            else:
                print("Entry not in library")

    def __str__(self):
        out = "Neutral losses in library:\n"
        out += "--------------------------------------------\n"
        out += "label\tmass_difference\n"
        for d in self.lib:
            out += "{}\t{}\n".format(d["label"], d["mass_difference"])
        return out


class MassDifferences:

    def __init__(self, ion_mode=None):

        self.ion_mode = ion_mode
        self.lib = []

    def add(self, label_x, label_y, mass_difference, charge_x=1, charge_y=1):
        self.lib.append(OrderedDict([(label_x, {"charge": float(charge_x)}),
                                     (label_y, {"charge": float(charge_y)}),
                                     ("mass_difference", mass_difference)]))
        self.lib = sorted(self.lib, key=lambda k: k['mass_difference'])

    def remove(self, label_x="*", label_y="*"):
        if label_x == "*" and label_y == "*":
            self.lib = []
        else:
            for item in self.lib:
                if label_x in item or label_y in item:
                    self.lib.remove(item)
            else:
                print("Entry not in library")

    def __str__(self):
        out = "Mass differences in library:\n"
        out += "--------------------------------------------\n"
        out += "label_x\tlabel_y\tmass_difference\tcharge_x\tcharge_y\n"
        for item in self.lib:
            label_x, label_y = list(item.items())[0][0], list(item.items())[1][0]
            out += "{}\t{}\t{}\t{}\t{}\n".format(label_x, label_y, item["mass_difference"],
                                                 item[label_x]["charge"], item[label_y]["charge"])
        return out
