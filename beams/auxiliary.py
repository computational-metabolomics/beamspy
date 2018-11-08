#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import copy
from collections import OrderedDict
from pandas import read_csv
import pyteomics
from pyteomics.mass import nist_mass


def order_composition_by_hill(composition):
    symbols = set(composition)
    if 'C' in symbols:
        symbols.remove('C')
        yield 'C'
        if 'H' in symbols:
            symbols.remove('H')
            yield 'H'
    for symbol in sorted(symbols):
        yield symbol


def composition_to_string(composition):
    molecular_formula = ""
    for atom in order_composition_by_hill(composition):
        if composition[atom] > 1:
            molecular_formula += atom + str(composition[atom])
        elif composition[atom] == 1:
            molecular_formula += atom
    return molecular_formula


def double_bond_equivalents(composition):
    c = {}
    X = sum([composition[h] for h in ["F", "Cl", "Br", "I", "At"] if h in composition])
    for e in ["C", "H", "N"]:
        if e in composition:
            c[e] = composition[e]
        else:
            c[e] = 0

    return float(c["C"]) - (float(c["H"])/2) - (float(X)/2) + (float(c["N"])/2) + 1


def HC_HNOPS_rules(molecular_formula):

    composition = pyteomics.mass.mass.Composition(molecular_formula)

    rules = {"HC": 0, "NOPSC": 0}

    if "C" not in composition or "H" not in composition:
        rules["HC"] = 0
    elif "C" not in composition and "H" not in composition:
        rules["HC"] = 0
    elif "C" in composition and "H" in composition:
        if float(composition['H']) / float((composition['C'])) > 0 and float(composition['H'] / (composition['C'])) < 6:
            rules["HC"] = 1
        if float(composition['H']) / float((composition['C'])) >= 6:
            rules["HC"] = 0

    NOPS_check = []
    for element in ['N', 'O', 'P', 'S']:
        if element in composition and "C" in composition:
            NOPS_check.append(float(float(composition[element])) / float((composition['C'])))
        else:
            NOPS_check.append(float(0))

    if NOPS_check[0] >= float(0) and \
       NOPS_check[0] <= float(4) and \
       NOPS_check[1] >= float(0) and \
       NOPS_check[1] <= float(3) and \
       NOPS_check[2] >= float(0) and \
       NOPS_check[2] <= float(2) and \
       NOPS_check[3] >= float(0) and \
       NOPS_check[3] <= float(3):
        rules["NOPSC"] = 1

    if NOPS_check[0] > float(4) or NOPS_check[1] > float(3) or NOPS_check[2] > float(2) or NOPS_check[3] > float(3):
        rules["NOPSC"] = 0
    return rules


def lewis_senior_rules(molecular_formula):

    valence = {'C': 4, 'H': 1, 'N': 3, 'O': 2, 'P': 3, 'S': 2}

    composition = pyteomics.mass.mass.Composition(molecular_formula)

    rules = {"lewis": 0, "senior": 0}

    lewis_sum = 0
    for element in valence:
        if element in composition:
            lewis_sum += valence[element] * composition[element]

    if lewis_sum % 2 == 0:
        rules["lewis"] = 1
    if lewis_sum % 2 != 0:
        rules["lewis"] = 0
    if lewis_sum >= ((sum(composition.values()) - 1) * 2):
        rules["senior"] = 1
    if lewis_sum < ((sum(composition.values()) - 1) * 2):
        rules["senior"] = 0

    return rules


def update_and_sort_nist_mass(path="", digits=6):

    if os.path.isfile(path):
        df = read_csv(path, sep="\t")
    else:
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'elements.txt')
        df = read_csv(path, sep="\t")

    nist_mass_copy = copy.deepcopy(nist_mass)
    es = list(order_composition_by_hill(nist_mass_copy.keys()))
    nist_mass_updated = OrderedDict((k, nist_mass_copy[k]) for k in es)

    for index, row in df.iterrows():
        exact_mass = nist_mass_copy[row["name"]][int(round(row["exact_mass"], 0))][0]
        for isotope, data in nist_mass_updated[row["name"]].items():
            abundance = nist_mass_updated[row["name"]][isotope][1]
            if data[0] == exact_mass:
                nist_mass_updated[row["name"]][isotope] = (round(float(row["exact_mass"]), digits), abundance)
            else:
                nist_mass_updated[row["name"]][isotope] = (round(nist_mass_updated[row["name"]][isotope][0], digits), abundance)
    return nist_mass_updated
