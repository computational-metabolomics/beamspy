#!/usr/bin/env python
# -*- coding: utf-8 -*-

from collections import OrderedDict
import gzip
import sqlite3
import pandas as pd
from pyteomics import mass as pyteomics_mass
from beamspy.db_parsers import parse_nist_database


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

    composition = pyteomics_mass.Composition(molecular_formula)

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

    composition = pyteomics_mass.Composition(molecular_formula)

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


def nist_database_to_pyteomics(fn, skip_lines=10):

    """
    :param fn: text file (NISTs Linearized ASCII Output)
    :param skip_lines: the number of lines of the data file to skip before beginning to read data.
    :return: Ordered dictionary containing NIST records compatible with 'Pyteomics'
    """

    def add_record(r, nm):
        if r["Atomic Symbol"] not in nm:
            nm[r["Atomic Symbol"]] = OrderedDict([(0, (0.0, 0.0))]) # update after all records have been added
            nm[r["Atomic Symbol"]][r["Mass Number"]] = (r["Relative Atomic Mass"][0], r["Isotopic Composition"][0])
        else:
            nm[r["Atomic Symbol"]][r["Mass Number"]] = (r["Relative Atomic Mass"][0], r["Isotopic Composition"][0])
        return nm

    lib = OrderedDict()
    for record in parse_nist_database(fn, skip_lines=skip_lines):
        if record["Atomic Symbol"] in ["D", "T"]:
            lib = add_record(record, lib)
            record["Atomic Symbol"] = "H"
            lib = add_record(record, lib)
        else:
            lib = add_record(record, lib)

    for element in list(lib.keys()):
        lib_sorted = sorted(lib[element].items(), key=lambda e: e[1][1], reverse=True)
        if lib_sorted[0][1][0] > 0.0:
            lib[element][0] = (lib_sorted[0][1][0], 1.0)
        elif len(lib_sorted) == 2:
            lib[element][0] = (lib_sorted[1][1][0], 1.0)
        else:
            del lib[element]

    es = list(order_composition_by_hill(lib.keys()))
    return OrderedDict((k, lib[k]) for k in es)


def convert_sql_to_text(path_sql, table_name, path_out, separator="\t"):

    if path_sql.endswith(".gz"):
        db_dump = gzip.GzipFile(path_sql, mode='rb')
    else:
        db_dump = open(path_sql, mode='rb')

    conn = sqlite3.connect(":memory:")
    cursor = conn.cursor()
    cursor.executescript(db_dump.read().decode('utf-8'))
    conn.commit()

    df = pd.read_sql_query("select * from " + table_name, conn)
    df.to_csv(path_out, sep=separator)

    conn.close()
    db_dump.close()
