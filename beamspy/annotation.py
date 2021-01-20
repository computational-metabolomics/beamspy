#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import itertools
import gzip
import sqlite3
from collections import OrderedDict
from urllib.parse import urlparse
import requests
import difflib
import pandas as pd
import numpy as np
import networkx as nx
from pyteomics import mass as pyteomics_mass
from beamspy.in_out import read_molecular_formulae
from beamspy.in_out import read_compounds
from beamspy.auxiliary import nist_database_to_pyteomics
from beamspy.auxiliary import composition_to_string


_sql_create_table_isotopes_ = """
                   CREATE TABLE isotopes (
                   peak_id_a TEXT DEFAULT NULL,
                   peak_id_b TEXT DEFAULT NULL,
                   label_a TEXT DEFAULT NULL,
                   label_b TEXT DEFAULT NULL,
                   atoms REAL DEFAULT NULL,
                   exact_mass_diff REAL DEFAULT NULL,
                   ppm_error REAL DEFAULT NULL,
                   charge INT DEFAULT 1,
                   PRIMARY KEY (peak_id_a, peak_id_b, label_a, label_b));
                   """

_sql_create_table_neutral_losses_ = """
                   CREATE TABLE neutral_losses (
                   peak_id_a TEXT DEFAULT NULL,
                   peak_id_b TEXT DEFAULT NULL,
                   label TEXT DEFAULT NULL,
                   exact_mass_diff REAL DEFAULT NULL,
                   ppm_error REAL DEFAULT NULL,
                   PRIMARY KEY (peak_id_a, peak_id_b, label));
                   """


def calculate_mz_tolerance(mass, ppm):
    min_tol = mass - (mass * 0.000001 * ppm)
    max_tol = mass + (mass * 0.000001 * ppm)
    return min_tol, max_tol


def calculate_rt_tolerance(rt, rt_tol):
    return rt - rt_tol, rt + rt_tol


def calculate_ppm_error(mass, theo_mass):
    return float(mass - theo_mass) / (theo_mass * 0.000001)


def _remove_elements_from_compositions(records, keep):

    path_nist_database = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'nist_database.txt')
    nist_database = nist_database_to_pyteomics(path_nist_database)

    elements = [e for e in nist_database if e not in keep]
    for record in records:
        for e in elements:
            if "composition" in record:
                record["composition"].pop(e, None)
            else:
                record.pop(e, None)
    return records


def _flatten_composition(records):
    for record in records:
        record.update(record["composition"])
        del record["composition"]
    return records


def _prep_lib(lib):
    lib_pairs = []
    if isinstance(lib, OrderedDict):
        combs = list(itertools.combinations(lib, 2))
        for pair in combs:
            if lib[pair[0]]["charge"] == lib[pair[1]]["charge"]:
                if lib[pair[0]]["mass"] < lib[pair[1]]["mass"]:
                    # print("yes", lib[pair[0]]["mass"], lib[pair[1]]["mass"], lib[pair[0]]["mass"] < lib[pair[1]]["mass"])
                    lib_pairs.append(OrderedDict([(pair[0], {"mass": lib[pair[0]]["mass"], "charge": lib[pair[0]]["charge"]}),
                                                  (pair[1], {"mass": lib[pair[1]]["mass"], "charge": lib[pair[1]]["charge"]})]))
                else:
                    lib_pairs.append(OrderedDict([(pair[1], {"mass": lib[pair[1]]["mass"], "charge": lib[pair[1]]["charge"]}),
                                                  (pair[0], {"mass": lib[pair[0]]["mass"], "charge": lib[pair[0]]["charge"]})]))
            elif lib[pair[0]]["charge"] > lib[pair[1]]["charge"]:
                lib_pairs.append(OrderedDict([(pair[0], {"mass": lib[pair[0]]["mass"], "charge": lib[pair[0]]["charge"]}),
                                              (pair[1], {"mass": lib[pair[1]]["mass"], "charge": lib[pair[1]]["charge"]})]))
            else:
                lib_pairs.append(OrderedDict([(pair[1], {"mass": lib[pair[1]]["mass"], "charge": lib[pair[1]]["charge"]}),
                                              (pair[0], {"mass": lib[pair[0]]["mass"], "charge": lib[pair[0]]["charge"]})]))

        #lib_pairs = sorted(lib_pairs, key=lambda pair: (list(pair.items())[0][1]["mass"] - list(pair.items())[1][1]["mass"]), reverse=True)
        return lib_pairs

    elif isinstance(lib, list) and isinstance(lib[0], OrderedDict):
        if "mass_difference" in lib[0]:
            return sorted(lib, key=lambda d: d["mass_difference"], reverse=True)
        else:
            raise ValueError("Format library incorrect")
        #else:
        #    return sorted(lib_pairs, key=lambda pair: (list(pair.items())[0][1]["mass"] - list(pair.items())[1][1]["mass"]), reverse=True)
    else:
        raise ValueError("Incorrect format for library: {}".format(type(lib)))


def _annotate_artifacts(peaklist, diff=0.02):
    n = peaklist.iloc[:, 1]
    for i in range(n):
        for j in range(i + 1, n):
            mz_diff = peaklist.iloc[i,1] - peaklist.iloc[j,1]
            ppm_error = calculate_ppm_error(peaklist.iloc[i,1], peaklist.iloc[j,1])
            if abs(mz_diff) < diff:
                yield i, j, mz_diff, ppm_error


def _check_tolerance(mz_x, mz_y, lib_pair, ppm, charge):
    min_tol_a, max_tol_a = calculate_mz_tolerance(mz_x, ppm)
    min_tol_b, max_tol_b = calculate_mz_tolerance(mz_y, ppm)
    if "mass_difference" in lib_pair.keys():
        min_tol_b = min_tol_b - (lib_pair["mass_difference"])
        max_tol_b = max_tol_b - (lib_pair["mass_difference"])
    elif "mass" in list(lib_pair.items())[0][1]:
        charge_a = list(lib_pair.items())[0][1]["charge"]
        charge_b = list(lib_pair.items())[1][1]["charge"]
        mass_a = list(lib_pair.items())[0][1]["mass"]
        mass_b = list(lib_pair.items())[1][1]["mass"]

        min_tol_a -= mass_a
        max_tol_a -= mass_a
        min_tol_b -= mass_b
        max_tol_b -= mass_b

        if charge_a > 1:
            min_tol_a += (charge_a - 1) * (mz_x - mass_a)
            max_tol_a += (charge_a - 1) * (mz_x - mass_a)
        if charge_b > 1:
            min_tol_b += (charge_b - 1) * (mz_y - mass_b)
            max_tol_b += (charge_b - 1) * (mz_y - mass_b)

    else:
        raise ValueError("Incorrect format: {}".format(lib_pair))
    #if min_tol_b > min_tol_a and min_tol_b > max_tol_a:
    #    return -1

    # x1 <= mass <= x2
    # y1 <= mass <= y2
    # x1 <= y2 AND y1 <= x2

    if min_tol_a < max_tol_b and min_tol_b < max_tol_a:
        return 1

    return 0


def _annotate_pairs_from_graph(G, ppm, lib_pairs, charge):

    for e in G.edges(data=True):

        mz_x = G.nodes[e[0]]["mz"]
        mz_y = G.nodes[e[1]]["mz"]

        for lib_pair in lib_pairs:

            ct = _check_tolerance(mz_x, mz_y, lib_pair, ppm, charge)

            if ct:

                if "charge" in list(lib_pair.items())[0][1]:
                    charge_a = list(lib_pair.items())[0][1]["charge"]
                    charge_b = list(lib_pair.items())[1][1]["charge"]
                elif charge:
                    charge_a = charge
                    charge_b = charge
                else:
                    charge_a = 1
                    charge_b = 1

                if "mass_difference" in lib_pair:
                    charge = lib_pair["charge"]
                    ppm_error = calculate_ppm_error(
                        mz_x,
                        mz_y - (lib_pair["mass_difference"]))
                    exact_mass_diff = lib_pair["mass_difference"]
                else:
                    ppm_error = calculate_ppm_error(
                        (mz_x - list(lib_pair.items())[0][1]["mass"]) * charge_a,
                        (mz_y - list(lib_pair.items())[1][1]["mass"]) * charge_b)
                    exact_mass_diff = list(lib_pair.items())[1][1]["mass"] - list(lib_pair.items())[0][1]["mass"]

                yield OrderedDict([("peak_id_a", e[0]), ("peak_id_b", e[1]),
                                   ("label_a", list(lib_pair.keys())[0]),
                                   ("label_b", list(lib_pair.keys())[1]),
                                   ('charge_a', charge_a),
                                   ('charge_b', charge_b),
                                   ('exact_mass_diff', exact_mass_diff),
                                   ('ppm_error', round(ppm_error, 2))])


def _annotate_pairs_from_peaklist(peaklist, ppm, lib_pairs, charge):
    n = len(peaklist.iloc[:, 1])
    for i in range(n):
        for j in range(i + 1, n):

            for lib_pair in lib_pairs:

                ct = _check_tolerance(peaklist.iloc[i, 1], peaklist.iloc[j, 1], lib_pair, ppm, charge)

                if ct:
                    if "charge" in list(lib_pair.items())[0][1]:
                        charge_a = list(lib_pair.items())[0][1]["charge"]
                        charge_b = list(lib_pair.items())[1][1]["charge"]
                    elif charge:
                        charge_a = charge
                        charge_b = charge
                    else:
                        charge_a = 1
                        charge_b = 1

                    if "mass_difference" in lib_pair:
                        ppm_error = calculate_ppm_error(
                            peaklist.iloc[i, 1],
                            peaklist.iloc[j, 1] - lib_pair["mass_difference"])
                        exact_mass_diff = lib_pair["mass_difference"]
                        charge_a = lib_pair["charge"]
                        charge_b = lib_pair["charge"]
                    else:
                        ppm_error = calculate_ppm_error(
                            (peaklist.iloc[i, 1] - list(lib_pair.items())[0][1]["mass"]) * charge_a,
                            (peaklist.iloc[j, 1] - list(lib_pair.items())[1][1]["mass"]) * charge_b)
                        exact_mass_diff = list(lib_pair.items())[1][1]["mass"] - list(lib_pair.items())[0][1]["mass"]

                    yield OrderedDict([("peak_id_a", peaklist.iloc[i,0]), ("peak_id_b", peaklist.iloc[j,0]),
                                       ("label_a", list(lib_pair.keys())[0]),
                                       ("label_b", list(lib_pair.keys())[1]),
                                       ('charge_a', charge_a),
                                       ('charge_b', charge_b),
                                       ('exact_mass_diff', exact_mass_diff),
                                       ('ppm_error', round(ppm_error, 2))])


class DbCompoundsMemory:

    def __init__(self, filename, lib_adducts=[]):

        self.filename = filename
        self.lib_adducts = lib_adducts

        records = read_compounds(self.filename, lib_adducts=self.lib_adducts)

        if "retention_time" in list(records[0].keys()):
            rt_column = "retention_time REAL DEFAULT NULL,"
        else:
            rt_column = ""

        self.conn = sqlite3.connect(":memory:")
        self.cursor = self.conn.cursor()

        self.cursor.execute("""CREATE TABLE COMPOUNDS(
                            compound_id TEXT PRIMARY KEY NOT NULL,
                            compound_name TEXT,
                            exact_mass REAL,
                            {}
                            C INTEGER DEFAULT 0,
                            H INTEGER DEFAULT 0,
                            N INTEGER DEFAULT 0,
                            O INTEGER DEFAULT 0,
                            P INTEGER DEFAULT 0,
                            S INTEGER DEFAULT 0,
                            CHNOPS INTEGER DEFAULT NULL,
                            molecular_formula TEXT DEFAULT NULL,
                            adduct TEXT DEFAULT NULL
                            );""".format(rt_column))

        records = _remove_elements_from_compositions(records, keep=["C", "H", "N", "O", "P", "S"])
        records = _flatten_composition(records)
        for record in records:
            columns = ",".join(map(str, list(record.keys())))
            qms = ', '.join(['?'] * len(record.values()))
            query = """insert into COMPOUNDS ({}) values ({})""".format(columns, qms)
            self.cursor.execute(query, list(record.values()))

        self.cursor.execute("""CREATE INDEX IDX_EXACT_MASS ON COMPOUNDS (exact_mass);""")
        self.conn.commit()

        self.cursor.execute("select * from COMPOUNDS")

    def select_compounds(self, min_tol, max_tol, min_rt=None, max_rt=None):
        col_names = ["compound_id", "compound_name", "exact_mass", "C", "H", "N", "O", "P", "S", "CHNOPS", "molecular_formula", "adduct"]
        if min_rt:
            col_names.insert(3, "retention_time")
            sql_rt = " and retention_time >= {} and retention_time <= {}".format(min_rt, max_rt)
        else:
            sql_rt = ""
        sql_str = """SELECT {} FROM COMPOUNDS WHERE exact_mass >= {} and exact_mass <= {}{}""".format(",".join(map(str, col_names)), min_tol, max_tol, sql_rt)
        self.cursor.execute(sql_str)
        return [OrderedDict(zip(col_names, list(record))) for record in self.cursor.fetchall()]

    def close(self):
        self.conn.close()


class DbMolecularFormulaeMemory:

    def __init__(self, filename):

        self.filename = filename
        self.conn = sqlite3.connect(":memory:")
        self.cursor = self.conn.cursor()
        self.cursor.execute("""CREATE TABLE MF(
                            exact_mass REAL,
                            C INTEGER DEFAULT 0,
                            H INTEGER DEFAULT 0,
                            N INTEGER DEFAULT 0,
                            O INTEGER DEFAULT 0,
                            P INTEGER DEFAULT 0,
                            S INTEGER DEFAULT 0,
                            CHNOPS INTEGER DEFAULT NULL,
                            HC INTEGER DEFAULT NULL,
                            NOPSC INTEGER DEFAULT NULL,
                            lewis INTEGER DEFAULT NULL,
                            senior INTEGER DEFAULT NULL,
                            double_bond_equivalents REAL,
                            primary key (C,H,N,O,P,S,exact_mass)
                            );""")

        records = read_molecular_formulae(self.filename)
        records = _remove_elements_from_compositions(records, keep=["C", "H", "N", "O", "P", "S"])
        records = _flatten_composition(records)
        for record in records:
            columns = ",".join(map(str, list(record.keys())))
            qms = ', '.join(['?'] * len(record.values()))
            query = """insert into mf ({}) values ({})""".format(columns, qms)
            self.cursor.execute(query, list(record.values()))

        self.cursor.execute("""CREATE INDEX IDX_EXACT_MASS ON MF (exact_mass);""")
        self.cursor.execute("""CREATE INDEX IDX_EXACT_MASS_RULES ON MF (exact_mass, HC, NOPSC, LEWIS, SENIOR);""")
        self.conn.commit()

    def select_mf(self, min_tol, max_tol, rules):

        if rules:
            sql_filters = " and lewis = 1 and senior = 1 and HC = 1 and NOPSC = 1"
        else:
            sql_filters = ""

        col_names = ["exact_mass", "C", "H", "N", "O", "P", "S",
                     "double_bond_equivalents", "LEWIS", "SENIOR", "HC", "NOPSC"]

        self.cursor.execute("""SELECT exact_mass, C, H, N, O, P, S,
                            double_bond_equivalents, LEWIS, SENIOR, HC, NOPSC
                            from mf where exact_mass >= {} and exact_mass <= {}{}
                            """.format(min_tol, max_tol, sql_filters))

        return [OrderedDict(zip(col_names, list(record))) for record in self.cursor.fetchall()]

    def close(self):
        self.conn.close()


def annotate_adducts(source, db_out, ppm, lib, add=False):

    conn = sqlite3.connect(db_out)
    cursor = conn.cursor()

    if not add:
        cursor.execute("DROP TABLE IF EXISTS adduct_pairs")

        cursor.execute("""CREATE TABLE adduct_pairs (
                       peak_id_a TEXT DEFAULT NULL,
                       peak_id_b TEXT DEFAULT NULL,
                       label_a TEXT DEFAULT NULL,
                       label_b TEXT DEFAULT NULL,
                       charge_a INTEGER DEFAULT NULL,
                       charge_b INTEGER DEFAULT NULL,
                       exact_mass_diff REAL DEFAULT NULL,
                       ppm_error REAL DEFAULT NULL,
                       PRIMARY KEY (peak_id_a, peak_id_b, label_a, label_b));""")

    lib_pairs = _prep_lib(lib.lib)

    if isinstance(source, nx.classes.digraph.DiGraph):
        source = list(source.subgraph(c) for c in nx.weakly_connected_components(source))

    if isinstance(source, list) and len(source) > 0 and isinstance(source[0], nx.classes.digraph.DiGraph):
        for i, graph in enumerate(source):
            for assignment in _annotate_pairs_from_graph(graph, lib_pairs=lib_pairs, ppm=ppm, charge=None):
                cursor.execute("""INSERT OR REPLACE into adduct_pairs (peak_id_a, peak_id_b, 
                               label_a, label_b, 
                               charge_a, charge_b,
                               exact_mass_diff, ppm_error)
                               values (?,?,?,?,?,?,?,?)""", (str(assignment["peak_id_a"]), str(assignment["peak_id_b"]),
                                                         assignment["label_a"], assignment["label_b"],
                                                         assignment["charge_a"], assignment["charge_b"],
                                                         float(assignment["exact_mass_diff"]),
                                                         float(assignment["ppm_error"])))

    elif isinstance(source, pd.core.frame.DataFrame):
        for assignment in _annotate_pairs_from_peaklist(source, lib_pairs=lib_pairs, ppm=ppm, charge=None):
            cursor.execute("""INSERT OR REPLACE into adduct_pairs (peak_id_a, peak_id_b, 
                           label_a, label_b,
                           charge_a, charge_b,
                           exact_mass_diff, ppm_error)
                           values (?,?,?,?,?,?,?,?)""", (assignment["peak_id_a"], assignment["peak_id_b"],
                                                     assignment["label_a"], assignment["label_b"],
                                                     assignment["charge_a"], assignment["charge_b"],
                                                     float(assignment["exact_mass_diff"]),
                                                     float(assignment["ppm_error"])))

    cursor.execute("""CREATE INDEX IDX_peak_id_a ON adduct_pairs (peak_id_a);""")
    cursor.execute("""CREATE INDEX IDX_peak_id_b ON adduct_pairs (peak_id_b);""")

    conn.commit()
    conn.close()
    return


def annotate_isotopes(source, db_out, ppm, lib):

    conn = sqlite3.connect(db_out)
    cursor = conn.cursor()

    cursor.execute("DROP TABLE IF EXISTS isotopes")

    cursor.execute(_sql_create_table_isotopes_)

    lib_pairs = _prep_lib(lib.lib)

    abundances = {}
    for pair in lib.lib:
        abundances[list(pair.items())[0][0]] = list(pair.items())[0][1]
        abundances[list(pair.items())[1][0]] = list(pair.items())[1][1]

    if isinstance(source, nx.classes.digraph.DiGraph):
        source = list(source.subgraph(c) for c in nx.weakly_connected_components(source))

    if isinstance(source, list) and len(source) > 0 and isinstance(source[0], nx.classes.digraph.DiGraph):

        for graph in source:

            peaklist = graph.nodes(data=True)

            for assignment in _annotate_pairs_from_graph(G=graph, lib_pairs=lib_pairs, ppm=ppm, charge=None):

                if abundances[assignment["label_a"]]["abundance"] < abundances[assignment["label_b"]]["abundance"]:
                    # Lithium
                    y = abundances[assignment["label_a"]]['abundance'] * peaklist[assignment["peak_id_b"]]["intensity"]
                    x = 100.0 * peaklist[assignment["peak_id_a"]]["intensity"]
                else:
                    y = 100.0 * peaklist[assignment["peak_id_b"]]["intensity"]
                    x = abundances[assignment["label_b"]]['abundance'] * peaklist[assignment["peak_id_a"]]["intensity"]

                if x == 0.0 or y == 0.0:
                    atoms = None
                elif abundances[assignment["label_a"]]["abundance"] < abundances[assignment["label_b"]]["abundance"]:
                    atoms = x/y
                else:
                    atoms = y/x

                cursor.execute("""insert into isotopes (peak_id_a, peak_id_b, label_a, label_b, 
                               atoms, exact_mass_diff, ppm_error, charge)
                               values (?,?,?,?,?,?,?,?)""", (str(assignment["peak_id_a"]), str(assignment["peak_id_b"]),
                               assignment["label_a"], assignment["label_b"], float(atoms),
                               float(assignment["exact_mass_diff"]), float(assignment["ppm_error"]), assignment["charge_a"]))

    elif isinstance(source, pd.core.frame.DataFrame):
        for assignment in _annotate_pairs_from_peaklist(peaklist=source, lib_pairs=lib_pairs, ppm=ppm, charge=None):

            if abundances[assignment["label_a"]]["abundance"] < abundances[assignment["label_b"]]["abundance"]:
                # Lithium
                y = abundances[assignment["label_a"]]["abundance"] * source.loc[source['name'] == assignment["peak_id_b"]]["intensity"].iloc[0]
                x = 100.0 * source.loc[source['name'] == assignment["peak_id_a"]]["intensity"].iloc[0]
            else:
                y = 100.0 * source.loc[source['name'] == assignment["peak_id_b"]]["intensity"].iloc[0]
                x = abundances[assignment["label_b"]]["abundance"] * source.loc[source['name'] == assignment["peak_id_a"]]["intensity"].iloc[0]

            if x == 0.0 or y == 0.0:
                atoms = None
            elif abundances[assignment["label_a"]]["abundance"] < abundances[assignment["label_b"]]["abundance"]:
                atoms = x/y
            else:
                atoms = y/x

            cursor.execute("""insert into isotopes (peak_id_a, peak_id_b, label_a, label_b, 
                           atoms, exact_mass_diff, ppm_error, charge)
                           values (?,?,?,?,?,?,?,?)""", (assignment["peak_id_a"], assignment["peak_id_b"],
                           assignment["label_a"], assignment["label_b"], atoms,
                           float(assignment["exact_mass_diff"]), assignment["ppm_error"], assignment["charge_a"]))

    cursor.execute("""CREATE INDEX IDX_isotopes_peak_id_a ON isotopes (peak_id_a);""")
    cursor.execute("""CREATE INDEX IDX_isotopes_peak_id_b ON isotopes (peak_id_b);""")

    conn.commit()
    conn.close()
    return


def annotate_oligomers(source, db_out, ppm, lib, maximum=3):

    conn = sqlite3.connect(db_out)
    cursor = conn.cursor()

    cursor.execute("DROP TABLE IF EXISTS oligomers")

    cursor.execute("""CREATE TABLE oligomers (
                   peak_id_a TEXT DEFAULT NULL,
                   peak_id_b TEXT DEFAULT NULL,
                   mz_a REAL DEFAULT NULL,
                   mz_b REAL DEFAULT NULL,
                   label_a TEXT DEFAULT NULL,
                   label_b TEXT DEFAULT NULL,
                   charge_a INTEGER DEFAULT NULL,
                   charge_b INTEGER DEFAULT NULL,
                   mz_ratio REAL DEFAULT NULL,
                   ppm_error REAL DEFAULT NULL,
                   PRIMARY KEY (peak_id_a, peak_id_b, label_a, label_b));""")

    if isinstance(source, nx.classes.digraph.DiGraph):
        source = list(source.subgraph(c) for c in nx.weakly_connected_components(source))

    if isinstance(source, list) and len(source) > 0 and isinstance(source[0], nx.classes.digraph.DiGraph):

        for graph in source:

            for n in graph.nodes():

                neighbors = list(graph.neighbors(n))

                for d in range(1, maximum):

                    for nn in neighbors:

                        mz_x = graph.nodes[n]["mz"]
                        mz_y = graph.nodes[nn]["mz"]

                        if mz_x < mz_y:

                            for adduct in lib.lib.keys():

                                if lib.lib[adduct]["charge"] > 1:
                                    continue

                                min_tol_a, max_tol_a = calculate_mz_tolerance(mz_x + ((mz_x - lib.lib[adduct]["mass"]) * d), ppm)
                                min_tol_b, max_tol_b = calculate_mz_tolerance(mz_y, ppm)

                                if (min_tol_b > max_tol_a and max_tol_b > max_tol_a):# or (min_tol_a < min_tol_b and max_tol_a < min_tol_b):
                                    #print(source.iloc[i][1], source.iloc[j][1], adduct)
                                    break

                                min_tol_a -= lib.lib[adduct]["mass"]
                                max_tol_a -= lib.lib[adduct]["mass"]

                                min_tol_b -= lib.lib[adduct]["mass"]
                                max_tol_b -= lib.lib[adduct]["mass"]

                                if min_tol_a < max_tol_b and min_tol_b < max_tol_a:

                                    a = (mz_x - lib.lib[adduct]["mass"]) + (mz_x - lib.lib[adduct]["mass"]) * d
                                    b = mz_y - lib.lib[adduct]["mass"]

                                    ratio = (mz_y - lib.lib[adduct]["mass"]) / (mz_x - lib.lib[adduct]["mass"])
                                    ppm_error = calculate_ppm_error(a, b)

                                    if "M" in adduct:
                                        adduct_oligo = adduct.replace("M", "{}M".format(int(round(ratio))))
                                    else:
                                        adduct_oligo = "{}{}".format(int(round(ratio)), adduct)

                                    cursor.execute("""insert into oligomers (peak_id_a, peak_id_b, mz_a, mz_b, label_a, label_b, charge_a, charge_b, mz_ratio, ppm_error)
                                                   values (?,?,?,?,?,?,?,?,?,?)""", (n, nn, mz_x, mz_y, adduct, adduct_oligo,
                                                                                     lib.lib[adduct]["charge"], lib.lib[adduct]["charge"], round(ratio, 2), round(ppm_error, 2)))
    elif isinstance(source, pd.core.frame.DataFrame):

        n = len(source.iloc[:,0])
        for adduct in lib.lib.keys():

            if lib.lib[adduct]["charge"] > 1:
                continue

            for i in range(n):
                for d in range(1, maximum):
                    for j in range(i + 1, n):

                        min_tol_a, max_tol_a = calculate_mz_tolerance(source.iloc[i][1] + ((source.iloc[i][1] - lib.lib[adduct]["mass"]) * d), ppm)
                        min_tol_b, max_tol_b = calculate_mz_tolerance(source.iloc[j][1], ppm)

                        if (min_tol_b > max_tol_a and max_tol_b > max_tol_a):# or (min_tol_a < min_tol_b and max_tol_a < min_tol_b):
                            #print(source.iloc[i][1], source.iloc[j][1], adduct)
                            break

                        min_tol_a -= lib.lib[adduct]["mass"]
                        max_tol_a -= lib.lib[adduct]["mass"]

                        min_tol_b -= lib.lib[adduct]["mass"]
                        max_tol_b -= lib.lib[adduct]["mass"]

                        if min_tol_a < max_tol_b and min_tol_b < max_tol_a:

                            a = (source.iloc[i][1] - lib.lib[adduct]["mass"]) + (source.iloc[i][1] - lib.lib[adduct]["mass"]) * d
                            b = source.iloc[j][1] - lib.lib[adduct]["mass"]

                            ratio = (source.iloc[j][1] - lib.lib[adduct]["mass"]) / (source.iloc[i][1] - lib.lib[adduct]["mass"])
                            ppm_error = calculate_ppm_error(a, b)

                            if "M" in adduct:
                                adduct_oligo = adduct.replace("M", "{}M".format(int(round(ratio))))
                            else:
                                adduct_oligo = "{}{}".format(int(round(ratio)), adduct)

                            cursor.execute("""insert into oligomers (peak_id_a, peak_id_b, mz_a, mz_b, label_a, label_b, charge_a, charge_b, mz_ratio, ppm_error)
                                           values (?,?,?,?,?,?,?,?,?,?)""", (source.iloc[i][0], source.iloc[j][0], source.iloc[i][1], source.iloc[j][1], adduct, adduct_oligo,
                                                                             lib.lib[adduct]["charge"], lib.lib[adduct]["charge"], round(ratio, 2), round(ppm_error, 2)))

    cursor.execute("""CREATE INDEX IDX_oligomers_peak_id_a ON oligomers (peak_id_a);""")
    cursor.execute("""CREATE INDEX IDX_oligomers_peak_id_b ON oligomers (peak_id_b);""")

    conn.commit()
    conn.close()
    return


def annotate_neutral_losses(source, db_out, ppm, lib):

    conn = sqlite3.connect(db_out)
    cursor = conn.cursor()

    cursor.execute("DROP TABLE IF EXISTS neutral_losses")

    cursor.execute(_sql_create_table_neutral_losses_)

    if isinstance(source, nx.classes.digraph.DiGraph):
        source = list(source.subgraph(c) for c in nx.weakly_connected_components(source))

    if isinstance(source, list) and len(source) > 0 and isinstance(source[0], nx.classes.digraph.DiGraph):

        for graph in source:

            for e in graph.edges(data=True):

                mz_x = graph.nodes[e[0]]["mz"]
                mz_y = graph.nodes[e[1]]["mz"]

                for nl in lib.lib:

                    ct = _check_tolerance(mz_x, mz_y, nl, ppm, charge=1)
                    if ct:
                        ppm_error = calculate_ppm_error(
                            mz_x, mz_y - nl["mass_difference"])

                        cursor.execute("""insert into neutral_losses (peak_id_a, peak_id_b, label,
                                       exact_mass_diff, ppm_error)
                                       values (?,?,?,?,?)""", (e[0], e[1],
                                       nl["label"], nl["mass_difference"], ppm_error))

    elif isinstance(source, pd.core.frame.DataFrame):

        n = len(source.iloc[:, 1])
        for i in range(n):
            for j in range(i + 1, n):
                for nl in lib.lib:
                    ct = _check_tolerance(source.iloc[i, 1], source.iloc[j, 1], nl, ppm, charge=1)
                    if ct:
                        ppm_error = calculate_ppm_error(
                            source.iloc[i, 1],
                            source.iloc[j, 1] - nl["mass_difference"])
                        cursor.execute("""insert into neutral_losses (peak_id_a, peak_id_b, label,
                                       exact_mass_diff, ppm_error)
                                       values (?,?,?,?,?)""", (source.iloc[i, 0], source.iloc[j, 0],
                                       nl["label"], nl["mass_difference"], ppm_error))

    cursor.execute("""CREATE INDEX IDX_nls_peak_id_a ON neutral_losses (peak_id_a);""")
    cursor.execute("""CREATE INDEX IDX_nls_peak_id_b ON neutral_losses (peak_id_b);""")

    conn.commit()
    conn.close()
    return


def annotate_artifacts(source, db_out, diff):

    conn = sqlite3.connect(db_out)
    cursor = conn.cursor()

    cursor.execute("DROP TABLE IF EXISTS artifacts")

    cursor.execute("""CREATE TABLE artifacts (
                   peak_id_a TEXT DEFAULT NULL,
                   peak_id_b TEXT DEFAULT NULL,
                   mz_diff REAL DEFAULT NULL,
                   ppm_error REAL DEFAULT NULL,
                   PRIMARY KEY (peak_id_a, peak_id_b));""")

    if isinstance(source, nx.classes.digraph.DiGraph):
        source = list(source.subgraph(c) for c in nx.weakly_connected_components(source))

    if (isinstance(source, list) or isinstance(source, np.ndarray)) and isinstance(source[0], nx.classes.graph.Graph):
        for graph in source:
            peaklist = graph.nodes(data=True)
            for assignment in _annotate_artifacts(peaklist, diff=diff):
                cursor.execute("""insert into artifacts (peak_id_a, peak_id_b, mz_diff, ppm_error)
                               values (?,?,?,?)""", (assignment["peak_id_a"], assignment["peak_id_b"], assignment["label_a"], assignment["label_b"]))

    elif isinstance(source, pd.core.frame.DataFrame):
        for assignment in _annotate_artifacts(source, diff=diff):
            cursor.execute("""insert into artifacts (peak_id_a, peak_id_b, mz_diff, ppm_error)
                           values (?,?,?,?)""", (assignment["peak_id_a"], assignment["peak_id_b"], assignment["label_a"], assignment["label_b"]))

    conn.commit()
    return


def _select_unions_peak_patterns(cursor):

    cursor.execute("SELECT name FROM sqlite_master WHERE type='table'")
    tables = cursor.fetchall()

    # Prepare (empty) sqlite tables if not exist to make union queries more consistent
    if ("isotopes",) not in tables:
        cursor.execute(_sql_create_table_isotopes_.replace("CREATE", "CREATE TEMP"))
    if ("neutral_losses",) not in tables:
        cursor.execute(_sql_create_table_neutral_losses_.replace("CREATE", "CREATE TEMP"))

    records = []
    sql_unions = []
    excl_ids_pp = []
    for t in ["adduct_pairs", "oligomers"]:
        if (t,) in tables:

            if t == "oligomers":
                col_mz_ratio = "ap.mz_ratio AS mz_ratio"
            else:
                col_mz_ratio = "1 AS mz_ratio"

            sql_str = """
                      SELECT ap.peak_id_a      AS peak_id_a,
                             iso.peak_id_b     AS peak_id_b,
                             nl_ap.peak_id_a   AS peak_id_a_nl,
                             nl_ap.peak_id_b   AS peak_id_aa_nl,
                             nl_iso.peak_id_a  AS peak_id_b_nl,
                             nl_iso.peak_id_b  AS peak_id_bb_nl,
                             ap.label_a        AS adduct_label,
                             IFNULL(ap.peak_id_b = iso.peak_id_b, 0) AS flag,
                             iso.label_a       AS iso_label_a,
                             iso.label_b       AS iso_label_b,
                             ap.charge_a       AS ap_charge,
                             iso.charge        AS iso_charge,
                             1                 AS mz_ratio,
                             IFNULL(iso.exact_mass_diff, 0.0)    AS iso_exact_mass_diff,
                             IFNULL(nl_ap.exact_mass_diff, 0.0)  AS nl_exact_mass_diff,
                             nl_ap.label       AS nl_label
                      FROM   {}                AS ap
                      LEFT JOIN isotopes       AS iso
                      ON     ap.peak_id_a = iso.peak_id_a AND ap.charge_a = iso.charge
                      LEFT JOIN neutral_losses AS nl_ap
                      ON     (ap.peak_id_a = nl_ap.peak_id_a OR ap.peak_id_a = nl_ap.peak_id_b)
                      LEFT JOIN neutral_losses AS nl_iso
                      ON     (iso.peak_id_b = nl_iso.peak_id_a OR iso.peak_id_b = nl_iso.peak_id_b)
                      UNION
                      SELECT ap.peak_id_b      AS peak_id_a,
                             iso.peak_id_b     AS peak_id_b,
                             nl_ap.peak_id_a   AS peak_id_a_nl,
                             nl_ap.peak_id_b   AS peak_id_aa_nl,
                             nl_iso.peak_id_a  AS peak_id_b_nl,
                             nl_iso.peak_id_b  AS peak_id_bb_nl,
                             ap.label_b        AS adduct_label,
                             0                 AS flag,
                             iso.label_a       AS iso_label_a,
                             iso.label_b       AS iso_label_b,
                             ap.charge_b       AS ap_charge,
                             iso.charge        AS iso_charge,
                             {},
                             IFNULL(iso.exact_mass_diff, 0)    AS iso_exact_mass_diff,
                             IFNULL(nl_ap.exact_mass_diff, 0)  AS nl_exact_mass_diff,
                             nl_ap.label       AS nl_label
                      FROM   {}                AS ap
                      LEFT JOIN isotopes       AS iso
                      ON     ap.peak_id_b = iso.peak_id_a AND ap.charge_b = iso.charge
                      LEFT JOIN neutral_losses AS nl_ap
                      ON     (ap.peak_id_b = nl_ap.peak_id_a OR ap.peak_id_b = nl_ap.peak_id_b)
                      LEFT JOIN neutral_losses AS nl_iso
                      ON     (iso.peak_id_b = nl_iso.peak_id_a OR iso.peak_id_b = nl_iso.peak_id_b)
                      """.format(t, col_mz_ratio, t)

            sql_unions.append(sql_str)
            excl_ids_pp.append("""SELECT peak_id_a FROM {} UNION SELECT peak_id_b FROM {}""".format(t, t))

    cursor.execute(" UNION ".join(map(str, sql_unions)))
    records.extend([dict(zip([c[0] for c in cursor.description], record)) for record in cursor.fetchall()])

    if len(excl_ids_pp) > 0:
        excl_ids_adducts_oligomers = " union ".join(map(str, excl_ids_pp))
        sql_excl = """AND iso.peak_id_a NOT IN ({})
                      AND iso.peak_id_b NOT IN ({})""".format(excl_ids_adducts_oligomers, excl_ids_adducts_oligomers)
    else:
        sql_excl = ""

              # WHERE  (nl_ap.label = nl_iso.label
              #        OR IFNULL(nl_ap.label, nl_iso.label) is NOT NULL
              #        OR IFNULL(nl_ap.label, nl_iso.label) is NULL)

    sql_str = """
                      SELECT iso.peak_id_a     AS peak_id_a,
                             iso.peak_id_b     AS peak_id_b,
                             nl_ap.peak_id_a   AS peak_id_a_nl,
                             nl_ap.peak_id_b   AS peak_id_aa_nl,
                             nl_iso.peak_id_a  AS peak_id_b_nl,
                             nl_iso.peak_id_b  AS peak_id_bb_nl,
                             NULL              AS adduct_label,
                             0                 AS flag,
                             iso.label_a       AS iso_label_a,
                             iso.label_b       AS iso_label_b,
                             NULL              AS ap_charge,
                             iso.charge        AS iso_charge,
                             1                 AS mz_ratio,
                             IFNULL(iso.exact_mass_diff, 0.0)    AS iso_exact_mass_diff,
                             IFNULL(nl_ap.exact_mass_diff, 0.0)  AS nl_exact_mass_diff,
                             nl_ap.label       AS nl_label
                      FROM   isotopes          AS iso
                      LEFT JOIN neutral_losses AS nl_ap
                      ON     (iso.peak_id_a = nl_ap.peak_id_a OR iso.peak_id_a = nl_ap.peak_id_b)
                      LEFT JOIN neutral_losses AS nl_iso
                      ON     (iso.peak_id_b = nl_iso.peak_id_a OR iso.peak_id_b = nl_iso.peak_id_b)

              {}""".format(sql_excl)

    cursor.execute(sql_str)
    records.extend([dict(zip([c[0] for c in cursor.description], record)) for record in cursor.fetchall()])

    excl_ids_pp.append("""SELECT peak_id_a FROM isotopes UNION SELECT peak_id_b FROM isotopes""")
    excl_ids = " UNION ".join(map(str, excl_ids_pp))
    sql_excl = """WHERE peak_id_a NOT IN ({})
                  AND peak_id_a_nl NOT IN ({})""".format(excl_ids, excl_ids)

    sql_str = """
              SELECT                    peak_id_a,
                     NULL            AS peak_id_b,
                     peak_id_b       AS peak_id_a_nl,
                     NULL            AS peak_id_b_nl,
                     NULL            AS adduct_label,
                     0               AS flag,
                     NULL            AS iso_label_a,
                     NULL            AS iso_label_b,
                     1               AS ap_charge,
                     1               AS iso_charge,
                     1               AS mz_ratio,
                     0.0             AS iso_exact_mass_diff,
                     exact_mass_diff AS nl_exact_mass_diff,
                     label           AS nl_label
              FROM   neutral_losses
              {}
              UNION
              SELECT peak_id_b       AS peak_id_a,
                     NULL            AS peak_id_b,
                     peak_id_a       AS peak_id_a_nl,
                     NULL            AS peak_id_b_nl,      
                     NULL            AS adduct_label,
                     0               AS flag,
                     NULL            AS iso_label_a,
                     NULL            AS iso_label_b,
                     1               AS ap_charge,
                     1               AS iso_charge,
                     1               AS mz_ratio,
                     0.0             AS iso_exact_mass_diff,
                     exact_mass_diff AS nl_exact_mass_diff,
                     label           AS nl_label
              FROM   neutral_losses
              {}""".format(sql_excl, sql_excl)
    cursor.execute(sql_str)
    records.extend([dict(zip([c[0] for c in cursor.description], record)) for record in cursor.fetchall()])

    return records


def annotate_molecular_formulae(peaklist, lib_adducts, ppm, db_out, db_in="https://mfdb.bham.ac.uk", patterns=True, rules=True, max_mz=None):

    conn = sqlite3.connect(db_out)
    cursor = conn.cursor()

    cursor.execute("DROP TABLE IF EXISTS molecular_formulae")

    cursor.execute("""CREATE TABLE molecular_formulae (
                    id TEXT DEFAULT NULL,
                    mz REAL DEFAULT NULL,
                    exact_mass REAL DEFAULT NULL,
                    ppm_error REAL DEFAULT NULL,
                    adduct TEXT DEFAULT NULL,
                    isotope TEXT DEFAULT '',
                    neutral_loss TEXT DEFAULT '',
                    C INTEGER DEFAULT 0,
                    H INTEGER DEFAULT 0,
                    N INTEGER DEFAULT 0,
                    O INTEGER DEFAULT 0,
                    P INTEGER DEFAULT 0,
                    S INTEGER DEFAULT 0,
                    CHNOPS INTEGER DEFAULT NULL,
                    molecular_formula TEXT DEFAULT NULL,
                    HC INTEGER DEFAULT NULL,
                    NOPSC INTEGER DEFAULT NULL,
                    lewis INTEGER DEFAULT NULL,
                    senior INTEGER DEFAULT NULL,
                    double_bond_equivalents REAL DEFAULT NULL,
                    primary key (id, mz, molecular_formula, adduct, isotope, neutral_loss)
                    );""")

    if os.path.isfile(db_in):
        conn_mem = DbMolecularFormulaeMemory(db_in)
        source = "sqlite"
        max_mz = None
    else:
        source = "api"
        url = '{}/api/formula/mass_range'.format(db_in)
        url_test = '{}/api/formula/mass?mass=180.06339&tol=0.0&tol_unit=ppm&rules=1'.format(db_in)
        o = urlparse(url)
        if o.scheme != "http" and o.netloc != "mfdb.bham.ac.uk":
            raise ValueError("No database or local db available")
        else:
            r = requests.get(url_test)
            r.raise_for_status()

    path_nist_database = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'nist_database.txt')
    nist_database = nist_database_to_pyteomics(path_nist_database)

    def _select_mfs(source, peak_id, mz, ppm, adducts, isotope, neutral_loss, n_oligo, exact_mass_diff=0.0, rules=True):

        min_mz, max_mz = calculate_mz_tolerance(mz, ppm)

        mf_records = []
        for adduct in adducts:

            if mz - lib_adducts.lib[adduct]["mass"] > 0.5:
                if source == "api":
                    params = {"lower": (min_mz - lib_adducts.lib[adduct]["mass"] - exact_mass_diff) *
                                        lib_adducts.lib[adduct]["charge"] / n_oligo,
                              "upper": (max_mz - lib_adducts.lib[adduct]["mass"] - exact_mass_diff) *
                                        lib_adducts.lib[adduct]["charge"] / n_oligo,
                              "rules": int(rules)}
                    response = requests.get(url, params=params)
                    records = response.json()["records"]
                else:
                    records = conn_mem.select_mf((min_mz - lib_adducts.lib[adduct]["mass"] - exact_mass_diff) *
                                                  lib_adducts.lib[adduct]["charge"] / n_oligo,
                                                 (max_mz - lib_adducts.lib[adduct]["mass"] - exact_mass_diff) *
                                                 lib_adducts.lib[adduct]["charge"] / n_oligo, rules=rules)
                for record in records:
                    record["id"] = peak_id
                    if "CHNOPS" not in record:  # MFdb API specific
                        record["CHNOPS"] = True  # MFdb API specific
                    if "rules" in record:
                        record.update(record["rules"])
                        del record["rules"]
                    if "atoms" in record:
                        record.update(record["atoms"])
                        del record["atoms"]
                    record["exact_mass"] = (record["exact_mass"] / lib_adducts.lib[adduct]["charge"] * n_oligo) + \
                                           (float(lib_adducts.lib[adduct]["mass"]) + exact_mass_diff)
                    record["adduct"] = adduct
                    record["mz"] = mz
                    record["ppm_error"] = calculate_ppm_error(mz, record["exact_mass"])
                    comp = OrderedDict([(item, record[item]) for item in record if item in nist_database.keys()])
                    record["molecular_formula"] = composition_to_string(comp)
                    record["adduct"] = adduct
                    if isotope:
                        record["isotope"] = isotope
                    else:
                        record["isotope"] = ""
                    if neutral_loss:
                        record["neutral_loss"] = neutral_loss
                    else:
                        record["neutral_loss"] = ""
                records = _remove_elements_from_compositions(records, keep=["C", "H", "N", "O", "P", "S"])
                mf_records.extend(records)

        return mf_records

    rows = _select_unions_peak_patterns(cursor)
    names_to_skip = []

    for row in rows:

        if row["peak_id_a"] in names_to_skip and row["peak_id_b"] in names_to_skip and not row["adduct_label"]:
            continue

        records = []
        match = None

        if row["adduct_label"]:
            if row["mz_ratio"] > 1: # oligomers
                match = difflib.get_close_matches(row["adduct_label"], lib_adducts.lib.keys(), n=1)
                adducts = [match[0]]
            else:
                adducts = [str(row["adduct_label"])]
        else:
            if row["iso_charge"]:
                adducts = []
                for a in lib_adducts.lib.keys():
                    if lib_adducts.lib[a]["charge"] == row["iso_charge"]:
                        adducts.append(a)
            else:
                adducts = lib_adducts.lib.keys()

        index_name = peaklist["name"].tolist().index(str(row["peak_id_a"]))
        mz = peaklist["mz"].iloc[index_name]

        if max_mz is not None and mz > max_mz:
            continue

        # (source, peak_id, mz, ppm, adducts, isotope, neutral_loss, n_oligo, exact_mass_diff=0.0, rules=True):
        records_a = _select_mfs(source, row["peak_id_a"], mz, ppm, adducts, row["iso_label_a"], row["nl_label"], row["mz_ratio"], 0.0, rules=rules)
        if row["nl_label"] is not None: # Neutral Loss

            if row["peak_id_a_nl"] is not None:

                if row["peak_id_a_nl"] == row["peak_id_a"]:
                    peak_id_nl = row["peak_id_aa_nl"]
                    nl_exact_mass_diff = -row["nl_exact_mass_diff"]
                    a_nl_exact_mass_diff = row["nl_exact_mass_diff"]
                else:
                    peak_id_nl = row["peak_id_a_nl"]
                    nl_exact_mass_diff = row["nl_exact_mass_diff"]
                    a_nl_exact_mass_diff = -row["nl_exact_mass_diff"]

                records_nl = _select_mfs(source, row["peak_id_a"], mz, ppm, adducts, row["iso_label_a"], row["nl_label"], row["mz_ratio"], nl_exact_mass_diff, rules=rules)
                records_a.extend(records_nl)  # Neutral Loss

                index_name = peaklist["name"].tolist().index(str(peak_id_nl))
                mz = peaklist["mz"].iloc[index_name]
                records_a_nl = _select_mfs(source, peak_id_nl, mz, ppm, adducts, row['iso_label_a'], row["nl_label"], row["mz_ratio"], 0.0, rules=rules) # Neutral Loss
                records_nl = _select_mfs(source, peak_id_nl, mz, ppm, adducts, row['iso_label_a'], row["nl_label"], row["mz_ratio"], a_nl_exact_mass_diff, rules=rules)# Neutral Loss
                records_a_nl.extend(records_nl) # Neutral Loss

                for record_a in reversed(records_a):  # list changes during iteration
                    for record_a_nl in records_a_nl:
                        if record_a["molecular_formula"] == record_a_nl["molecular_formula"]:
                            records_a.append(record_a_nl)

                names_to_skip.append(peak_id_nl)

        if row["peak_id_b"]:

            index_name = peaklist["name"].tolist().index(str(row["peak_id_b"]))
            mz = peaklist["mz"].iloc[index_name]

            # if row["flag"]:  # different adducts - label_a and label_b?
            #     exact_mass_diff = 0.0  # adduct == isotope e.g. K / (41K) and [M+K]+ / [M+(41K)]+
            # else:
            exact_mass_diff = float(row["iso_exact_mass_diff"])

            records_b = _select_mfs(source, row["peak_id_b"], mz, ppm, adducts, row["iso_label_b"], row["nl_label"], row["mz_ratio"], exact_mass_diff, rules=rules)

            if row["nl_label"] is not None:  # Neutral Loss

                if row["peak_id_b_nl"] is not None:

                    if row["peak_id_b_nl"] == row["peak_id_b"]:
                        peak_id_nl = row["peak_id_bb_nl"]
                        nl_exact_mass_diff = -(row["nl_exact_mass_diff"] - exact_mass_diff)
                        b_nl_exact_mass_diff = exact_mass_diff + row["nl_exact_mass_diff"]
                    else:
                        peak_id_nl = row["peak_id_b_nl"]
                        nl_exact_mass_diff = (row["nl_exact_mass_diff"] + exact_mass_diff)
                        b_nl_exact_mass_diff = -row["nl_exact_mass_diff"] + exact_mass_diff

                    records_nl = _select_mfs(source, row["peak_id_b"], mz, ppm, adducts, row["iso_label_b"], row["nl_label"], row["mz_ratio"], nl_exact_mass_diff, rules=rules)
                    records_b.extend(records_nl)  # Neutral Loss

                    index_name = peaklist["name"].tolist().index(str(peak_id_nl))
                    mz = peaklist["mz"].iloc[index_name]
                    records_b_nl = _select_mfs(source, peak_id_nl, mz, ppm, adducts, row['iso_label_b'], row["nl_label"], row["mz_ratio"], exact_mass_diff, rules=rules)  # Neutral Loss
                    records_nl = _select_mfs(source, peak_id_nl, mz, ppm, adducts, row['iso_label_b'], row["nl_label"], row["mz_ratio"], b_nl_exact_mass_diff, rules=rules)  # Neutral Loss
                    records_b_nl.extend(records_nl)  # Neutral Loss

                    for record_b in reversed(records_b):
                        for record_b_nl in records_b_nl:
                            if record_b["molecular_formula"] == record_b_nl["molecular_formula"]:
                                records_b.append(record_b_nl)

                    names_to_skip.append(peak_id_nl)

            for record_a in records_a:
                for record_b in records_b:
                    if record_a["molecular_formula"] == record_b["molecular_formula"]:
                        if record_a not in records:
                            records.append(record_a)
                        if record_b not in records:
                            records.append(record_b)

            names_to_skip.append(row["peak_id_a"])
            names_to_skip.append(row["peak_id_b"])
        else:
            names_to_skip.append(row["peak_id_a"])
            records.extend(records_a)

        if len(records) > 0:

            if match:
                for record in records:
                    record["adduct"] = row["adduct_label"]

            sql_str = """INSERT OR IGNORE INTO molecular_formulae ({}) VALUES (:{})
                      """.format(",".join(map(str, records[0].keys())), ", :".join(map(str, records[0].keys())))
            cursor.executemany(sql_str, records)
            conn.commit()

    cursor.execute("select id, molecular_formula from molecular_formulae")
    mfs_subset = cursor.fetchall()

    for i in range(len(peaklist.iloc[:, 0])):

        mz = float(peaklist["mz"].iloc[i])
        name = str(peaklist["name"].iloc[i])

        if name in names_to_skip and patterns:
            continue

        if max_mz is not None and mz > max_mz:
            continue

        records = _select_mfs(source, name, mz, ppm, lib_adducts.lib.keys(), None, None, 1.0, 0.0, rules=rules)

        records_filt = [record for record in records if (record["id"], record["molecular_formula"]) not in mfs_subset]
        if len(records_filt) > 0:
            sql_str = """INSERT INTO molecular_formulae ({}) VALUES (:{})
                      """.format(",".join(map(str, records_filt[0].keys())), ", :".join(map(str, records_filt[0].keys())))
            cursor.executemany(sql_str, records_filt)

    if source != "api":
        conn_mem.close()

    conn.commit()
    conn.close()
    return


def annotate_compounds(peaklist, lib_adducts, ppm, db_out, db_name, patterns=True, db_in="", rt_tol=None):

    if db_in is None or db_in == "":
        conn_cpds = None
        path_dbs = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'databases')
        for db_local in os.listdir(path_dbs):
            if db_name == db_local.replace(".sql.gz", ""):

                with gzip.GzipFile(os.path.join(path_dbs, db_local), mode='rb') as db_dump:

                    conn_cpds = sqlite3.connect(":memory:")
                    cursor_cpds = conn_cpds.cursor()
                    cursor_cpds.executescript(db_dump.read().decode('utf-8'))
                    conn_cpds.commit()

                    cursor_cpds.execute("CREATE INDEX idx_exact_mass ON {} (exact_mass)".format(db_name.replace(".sql.gz", "")))

                    cursor_cpds.execute("SELECT name FROM sqlite_master WHERE type='table'")
                    if (db_name.replace(".sql.gz", ""), ) not in cursor_cpds.fetchall():
                        raise ValueError("Database {} not available".format(db_name))
                    break

        if conn_cpds is None:
            raise ValueError("Database {} not available".format(db_name))

    elif os.path.isfile(db_in):
        with open(db_in, 'rb') as fd:
            if fd.read(100)[:16].decode() == 'SQLite format 3\x00':
                conn_cpds = sqlite3.connect(db_in)
                cursor_cpds = conn_cpds.cursor()
                cursor_cpds.execute("SELECT name FROM sqlite_master WHERE type='table'")
                if not (db_name, ) in cursor_cpds.fetchall():
                    raise ValueError("Database {} not available".format(db_name))
            else:
                cursor_cpds = DbCompoundsMemory(db_in, lib_adducts=lib_adducts)

    else:
        raise IOError("[Errno 2] No such file or directory: {}".format(db_in))

    conn = sqlite3.connect(db_out)
    cursor = conn.cursor()

    cursor.execute("DROP TABLE IF EXISTS compounds_{}".format(db_name))
    cursor.execute("""CREATE TABLE compounds_{} (
                   id TEXT NOT NULL,
                   mz REAL DEFAULT NULL,
                   exact_mass REAL NOT NULL,
                   ppm_error REAL DEFAULT NULL,
                   rt_diff REAL DEFAULT NULL,
                   adduct TEXT NOT NULL,
                   isotope TEXT DEFAULT '',
                   neutral_loss TEXT DEFAULT '',
                   C INTEGER DEFAULT 0,
                   H INTEGER DEFAULT 0,
                   N INTEGER DEFAULT 0,
                   O INTEGER DEFAULT 0,
                   P INTEGER DEFAULT 0,
                   S INTEGER DEFAULT 0,
                   CHNOPS INTEGER DEFAULT NULL,
                   molecular_formula TEXT DEFAULT NULL,
                   compound_id TEXT NOT NULL,
                   compound_name TEXT DEFAULT NULL,
                   primary key (id, compound_id, adduct, isotope, neutral_loss) 
                   );""".format(db_name))
    conn.commit()

    def _select_compounds(db_cursor, peak_id, mz, ppm, adducts, isotope, neutral_loss, min_rt=None, max_rt=None, n_oligo=1, exact_mass_isotope=0.0):

        cpd_records = []

        min_mz, max_mz = calculate_mz_tolerance(mz, ppm)

        if min_rt and max_rt:
            cpd_records = db_cursor.select_compounds(min_mz, max_mz, min_rt=min_rt, max_rt=max_rt)
        else:
            for adduct in adducts:

                if min_mz - lib_adducts.lib[adduct]["mass"] < 0.5:
                    continue

                if isinstance(db_cursor, sqlite3.Cursor):
                    col_names = ["compound_id", "C", "H", "N", "O", "P", "S", "CHNOPS",
                                 "molecular_formula", "compound_name", "exact_mass"]
                    db_cursor.execute("""SELECT id, C, H, N, O, P, S, CHNOPS,
                                         molecular_formula, name, exact_mass
                                         from {} where exact_mass >= {} and exact_mass <= {}
                                         """.format(db_name,
                                                    (min_mz - lib_adducts.lib[adduct]["mass"] - exact_mass_isotope) *
                                                    lib_adducts.lib[adduct]["charge"] / n_oligo,
                                                    (max_mz - lib_adducts.lib[adduct]["mass"] - exact_mass_isotope) *
                                                    lib_adducts.lib[adduct]["charge"] / n_oligo))
                    cpd_records_subset = [OrderedDict(zip(col_names, list(record))) for record in db_cursor.fetchall()]

                else:
                    if min_rt and max_rt:
                        cpd_records_subset = db_cursor.select_compounds(min_mz - lib_adducts.lib[adduct]["mass"] - exact_mass_isotope,
                                                                        max_mz - lib_adducts.lib[adduct]["mass"] - exact_mass_isotope,
                                                                        min_rt=min_rt, max_rt=max_rt)
                    else:
                        cpd_records_subset = db_cursor.select_compounds((min_mz - lib_adducts.lib[adduct]["mass"] - exact_mass_isotope) *
                                                                        lib_adducts.lib[adduct]["charge"] / n_oligo,
                                                                        (max_mz - lib_adducts.lib[adduct]["mass"] - exact_mass_isotope) *
                                                                        lib_adducts.lib[adduct]["charge"] / n_oligo)
                for record in cpd_records_subset:
                    record["exact_mass"] = (record["exact_mass"] / lib_adducts.lib[adduct]["charge"] * n_oligo) + \
                                           (float(lib_adducts.lib[adduct]["mass"]) + exact_mass_isotope)
                    record["adduct"] = adduct

                cpd_records.extend(cpd_records_subset)

        for record in cpd_records:
            record["id"] = peak_id
            record["mz"] = mz
            record["ppm_error"] = calculate_ppm_error(mz, record["exact_mass"])
            if isotope:
                record["isotope"] = isotope
            else:
                record["isotope"] = "" # NULL/None is not used here as isotope is part of the primary key
            if neutral_loss:
                record["neutral_loss"] = neutral_loss
            else:
                record["neutral_loss"] = "" # NULL/None is not used here as isotope is part of the primary key
            if "retention_time" in record:
                record["rt_diff"] = rt - float(record["retention_time"])
                del record["retention_time"]

        return cpd_records

    rows = _select_unions_peak_patterns(cursor)

    names_to_skip = []
    for row in rows:

        if row["peak_id_a"] in names_to_skip and row["peak_id_b"] in names_to_skip and not row["adduct_label"]:
            continue

        records = []
        match = None

        if row["adduct_label"]:
            if row["mz_ratio"] > 1:  # oligomers
                match = difflib.get_close_matches(row["adduct_label"], lib_adducts.lib.keys(), n=1)
                adducts = [match[0]]
            else:
                adducts = [str(row["adduct_label"])]
        else:
            if row["iso_charge"]:
                adducts = []
                for a in lib_adducts.lib.keys():
                    if int(lib_adducts.lib[a]["charge"]) == int(row["iso_charge"]):
                        adducts.append(a)
            else:
                adducts = lib_adducts.lib.keys()

        index_name = peaklist["name"].tolist().index(str(row["peak_id_a"]))
        mz = peaklist["mz"].iloc[index_name]

        records_a = _select_compounds(cursor_cpds, row["peak_id_a"], mz, ppm, adducts, row['iso_label_a'], row["nl_label"], None, None, row["mz_ratio"], 0.0)
        if row["nl_label"] is not None: # Neutral Loss

            if row["peak_id_a_nl"] is not None:

                if row["peak_id_a_nl"] == row["peak_id_a"]:
                    peak_id_nl = row["peak_id_aa_nl"]
                    nl_exact_mass_diff = -row["nl_exact_mass_diff"]
                    a_nl_exact_mass_diff = row["nl_exact_mass_diff"]
                else:
                    peak_id_nl = row["peak_id_a_nl"]
                    nl_exact_mass_diff = row["nl_exact_mass_diff"]
                    a_nl_exact_mass_diff = -row["nl_exact_mass_diff"]

                records_nl = _select_compounds(cursor_cpds, row["peak_id_a"], mz, ppm, adducts, row['iso_label_a'], row["nl_label"], None, None, row["mz_ratio"], nl_exact_mass_diff) # Neutral Loss
                records_a.extend(records_nl) # Neutral Loss

                index_name = peaklist["name"].tolist().index(str(peak_id_nl))
                mz = peaklist["mz"].iloc[index_name]
                records_a_nl = _select_compounds(cursor_cpds, peak_id_nl, mz, ppm, adducts, row['iso_label_a'], row["nl_label"], None, None, row["mz_ratio"], 0.0) # Neutral Loss
                records_nl = _select_compounds(cursor_cpds, peak_id_nl, mz, ppm, adducts, row['iso_label_a'], row["nl_label"], None, None, row["mz_ratio"], a_nl_exact_mass_diff) # Neutral Loss
                records_a_nl.extend(records_nl) # Neutral Loss

                for record_a in reversed(records_a):  # list changes during iteration
                    for record_a_nl in records_a_nl:
                        if record_a["compound_id"] == record_a_nl["compound_id"]:
                            records_a.append(record_a_nl)

                names_to_skip.append(peak_id_nl)

        if row["peak_id_b"]:

            index_name = peaklist["name"].tolist().index(str(row["peak_id_b"]))
            mz = peaklist["mz"].iloc[index_name]

            #if row["flag"]:  # different adducts - label_a and label_b?
            #    exact_mass_diff = 0.0  # adduct == isotope e.g. K / (41K) and [M+K]+ / [M+(41K)]+
            #else:
            exact_mass_diff = float(row["iso_exact_mass_diff"])

            records_b = _select_compounds(cursor_cpds, row["peak_id_b"], mz, ppm, adducts, row['iso_label_b'], row["nl_label"], None, None, row["mz_ratio"], exact_mass_diff)
            if row["nl_label"] is not None:  # Neutral Loss

                if row["peak_id_b_nl"] is not None:

                    if row["peak_id_b_nl"] == row["peak_id_b"]:
                        peak_id_nl = row["peak_id_bb_nl"]
                        nl_exact_mass_diff = -(row["nl_exact_mass_diff"] - exact_mass_diff)
                        b_nl_exact_mass_diff = exact_mass_diff + row["nl_exact_mass_diff"]
                    else:
                        peak_id_nl = row["peak_id_b_nl"]
                        nl_exact_mass_diff = (row["nl_exact_mass_diff"] + exact_mass_diff)
                        b_nl_exact_mass_diff = -row["nl_exact_mass_diff"] + exact_mass_diff

                    records_nl = _select_compounds(cursor_cpds, row["peak_id_b"], mz, ppm, adducts, row['iso_label_b'], row["nl_label"], None, None, row["mz_ratio"], nl_exact_mass_diff)  # Neutral Loss
                    records_b.extend(records_nl) # Neutral Loss

                    index_name = peaklist["name"].tolist().index(str(peak_id_nl))
                    mz = peaklist["mz"].iloc[index_name]

                    records_b_nl = _select_compounds(cursor_cpds, peak_id_nl, mz, ppm, adducts, row['iso_label_b'], row["nl_label"], None, None, row["mz_ratio"], exact_mass_diff)  # Neutral Loss
                    records_nl = _select_compounds(cursor_cpds, peak_id_nl, mz, ppm, adducts, row['iso_label_b'], row["nl_label"], None, None, row["mz_ratio"], b_nl_exact_mass_diff)  # Neutral Loss
                    records_b_nl.extend(records_nl)  # Neutral Loss
                    for record_b in reversed(records_b):  # list changes during iteration
                        for record_b_nl in records_b_nl:
                            if record_b["compound_id"] == record_b_nl["compound_id"]:
                                records_b.append(record_b_nl)

                names_to_skip.append(peak_id_nl)

            for record_a in records_a:
                for record_b in records_b:
                    if record_a["compound_id"] == record_b["compound_id"]:
                        if record_a not in records:
                            records.append(record_a)
                        if record_b not in records:
                            records.append(record_b)

            names_to_skip.append(row["peak_id_b"])
        else:
            names_to_skip.append(row["peak_id_a"])
            records.extend(records_a)

        if len(records) > 0:
            if match:
                for record in records:
                    record["adduct"] = row["adduct_label"]

            sql_str = """INSERT OR IGNORE INTO compounds_{} ({}) VALUES (:{})
                      """.format(db_name, ",".join(map(str, records[0].keys())), ", :".join(map(str, records[0].keys())))
            cursor.executemany(sql_str, records)
            conn.commit()

    cursor.execute("select id, compound_id from compounds_{}".format(db_name))
    cpds_subset = cursor.fetchall()

    for i in range(len(peaklist.iloc[:, 0])):

        mz = float(peaklist["mz"].iloc[i])
        name = str(peaklist["name"].iloc[i])
        rt = float(peaklist["rt"].iloc[i])

        records = []
        if rt_tol:
            min_rt, max_rt = calculate_rt_tolerance(rt, rt_tol)
        else:
            min_rt, max_rt = None, None

        if min_rt and max_rt:
            records = _select_compounds(cursor_cpds, name, mz, ppm, None, None, None, min_rt, max_rt, 1.0, 0.0)
        else:
            if name in names_to_skip and patterns:
                continue
            records = _select_compounds(cursor_cpds, name, mz, ppm, lib_adducts.lib.keys(), None, None, None, None, 1.0, 0.0)

        records_filt = [record for record in records if (record["id"], record["compound_id"]) not in cpds_subset]
        if len(records_filt) > 0:
            sql_str = """INSERT INTO compounds_{} ({}) VALUES (:{})
                      """.format(db_name, ",".join(map(str, records_filt[0].keys())), ", :".join(map(str, records_filt[0].keys())))
            cursor.executemany(sql_str, records_filt)

    if isinstance(cursor_cpds, sqlite3.Cursor):
        conn_cpds.close()
    else:
        cursor_cpds.close()

    conn.commit()
    conn.close()
    return


def predict_drug_products(smiles, phase1_cycles, phase2_cycles):

    try:
        from rdkit import Chem
        import sygma
    except ImportError:
        raise ImportError('Install RDKit and/or SyGMa')

    # sygma/rules/phase1.txt
    # sygma/rules/phase2.txt

    # Each step in a scenario lists the ruleset and the number of reaction cycles to be applied
    scenario = sygma.Scenario([
        [sygma.ruleset['phase1'], phase1_cycles],
        [sygma.ruleset['phase2'], phase2_cycles]])

    # An rdkit molecule, optionally with 2D coordinates, is required as parent molecule
    parent = Chem.MolFromSmiles(smiles)

    metabolic_tree = scenario.run(parent)
    metabolic_tree.calc_scores()
    return metabolic_tree


class DbDrugCompoundsMemory:

    def __init__(self):
        self.conn = sqlite3.connect(":memory:")
        self.cursor = self.conn.cursor()
        self.cursor.execute("""CREATE TABLE predicted_drug_products (
                            compound_id TEXT PRIMARY KEY  NOT NULL,
                            compound_name TEXT,
                            smiles TEXT,
                            exact_mass decimal(15,7),
                            C INTEGER DEFAULT 0,
                            H INTEGER DEFAULT 0,
                            N INTEGER DEFAULT 0,
                            O INTEGER DEFAULT 0,
                            P INTEGER DEFAULT 0,
                            S INTEGER DEFAULT 0,
                            CHNOPS INTEGER DEFAULT NULL,
                            molecular_formula TEXT DEFAULT NULL,
                            sygma_score decimal(15,7),
                            sygma_pathway TEXT,
                            parent TEXT
                            );""")

        self.cursor.execute("""CREATE INDEX IDX_EXACT_MASS ON predicted_drug_products (exact_mass);""")
        self.conn.commit()

    def insert(self, records):
        for record in records:
            columns = ",".join(map(str, list(record.keys())))
            qms = ', '.join(['?'] * len(record.values()))
            query = """insert into predicted_drug_products ({}) values ({})""".format(columns, qms)
            self.cursor.execute(query, list(record.values()))
        self.conn.commit()

    def select(self, min_tol, max_tol):
        col_names = ["compound_id", "compound_name", "smiles", "sygma_score", "sygma_pathway", "parent", "exact_mass", "C", "H", "N", "O", "P", "S", "CHNOPS", "molecular_formula"]
        self.cursor.execute("""SELECT {} FROM predicted_drug_products WHERE 
                            exact_mass >= {} and exact_mass <= {}
                            """.format(",".join(map(str, col_names)), min_tol, max_tol))
        return [OrderedDict(zip(col_names, list(record))) for record in self.cursor.fetchall()]

    def close(self):
        self.conn.close()


def annotate_drug_products(peaklist, db_out, list_smiles, lib_adducts, ppm, phase1_cycles, phase2_cycles):

    try:
        from rdkit import Chem
        import sygma
    except ImportError:
        raise ImportError('Install RDKit and/or SyGMa')

    conn = sqlite3.connect(db_out)
    cursor = conn.cursor()

    cursor.execute("DROP TABLE IF EXISTS drug_products")
    cursor.execute("""CREATE TABLE drug_products (
                   id TEXT DEFAULT NULL,
                   mz REAL DEFAULT NULL,
                   exact_mass REAL DEFAULT NULL,
                   ppm_error REAL DEFAULT NULL,
                   adduct TEXT DEFAULT NULL,
                   C INTEGER DEFAULT 0,
                   H INTEGER DEFAULT 0,
                   N INTEGER DEFAULT 0,
                   O INTEGER DEFAULT 0,
                   P INTEGER DEFAULT 0,
                   S INTEGER DEFAULT 0,
                   CHNOPS INTEGER DEFAULT NULL,
                   molecular_formula TEXT DEFAULT NULL,
                   compound_id TEXT DEFAULT NULL,
                   compound_name TEXT DEFAULT NULL,
                   smiles TEXT,
                   sygma_score REAL DEFAULT 0.0,
                   sygma_pathway TEXT,
                   parent TEXT,
                   primary key (id, adduct, compound_id)
                   );""")

    path_nist_database = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'nist_database.txt')
    nist_database = nist_database_to_pyteomics(path_nist_database)

    records = []
    for smiles_parent in list_smiles:
        metabolic_tree = predict_drug_products(smiles_parent, phase1_cycles, phase2_cycles)
        for entry in metabolic_tree.to_list():
            smiles_product = Chem.MolToSmiles(entry['SyGMa_metabolite'])
            record = OrderedDict()
            record["compound_id"] = smiles_product
            record["compound_name"] = smiles_product
            record["sygma_pathway"] = entry["SyGMa_pathway"]
            record["parent"] = Chem.MolToSmiles(entry["parent"])
            mf = Chem.rdMolDescriptors.CalcMolFormula(Chem.MolFromSmiles(smiles_product))
            record["smiles"] = smiles_product
            record["sygma_score"] = entry['SyGMa_score']
            comp = pyteomics_mass.Composition(mf)
            record.update(comp)
            record["molecular_formula"] = composition_to_string(comp)
            record["exact_mass"] = round(pyteomics_mass.calculate_mass(formula=str(mf), mass_data=nist_database), 6)
            record["CHNOPS"] = sum([comp[e] for e in comp if e in ["C", "H", "N", "O", "P", "S"]]) == sum(list(comp.values()))
            records.append(record)

    conn_mem = DbDrugCompoundsMemory()
    records = _remove_elements_from_compositions(records, keep=["C", "H", "N", "O", "P", "S"])
    conn_mem.insert(records)

    for i in range(len(peaklist.iloc[:, 0])):
        mz = float(peaklist["mz"].iloc[i])
        name = str(peaklist["name"].iloc[i])
        min_tol, max_tol = calculate_mz_tolerance(mz, ppm)
        for adduct in lib_adducts.lib:

            if mz - lib_adducts.lib[adduct]["mass"] > 0.5:

                records = conn_mem.select(min_tol - lib_adducts.lib[adduct]["mass"], max_tol - lib_adducts.lib[adduct]["mass"])

                for record in records:
                    record["id"] = name
                    record["exact_mass"] = record["exact_mass"] + float(lib_adducts.lib[adduct]["mass"])
                    record["mz"] = mz
                    record["ppm_error"] = calculate_ppm_error(mz, record["exact_mass"])
                    record["adduct"] = adduct
                    cursor.execute("""insert into drug_products ({}) values (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
                                   """.format(",".join(map(str, list(record.keys())))), list(record.values()))
    conn.commit()
    conn.close()
    return


def summary(df, db, single_row=False, single_column=False, convert_rt=None, ndigits_mz=None):

    conn = sqlite3.connect(db)
    cursor = conn.cursor()

    cursor.execute("DROP TABLE IF EXISTS peaklist")
    df[["name", "mz", "rt", "intensity"]].to_sql("peaklist", conn, index=False)
    cursor.execute("CREATE INDEX idx_name on peaklist (name)")

    cursor.execute("SELECT name FROM sqlite_master WHERE type='table'")
    tables = cursor.fetchall()

    tables_pp = ["adduct_pairs", "oligomers", "isotopes", "neutral_losses"]  # TODO: make more efficient
    tables_to_union = []
    for tn in tables:
        if tn[0] in tables_pp:  # TODO - add addtional tables (neutral losses)
            tables_to_union.append(str(tn[0]))

    flag_groups = ("groups",) in tables

    if len(tables_to_union) > 0 and flag_groups:

        if len(tables_to_union) > 1:
            query = "select peak_id_a, peak_id_b from "
            query += " union select peak_id_a, peak_id_b from ".join(map(str, tables_to_union))
        elif len(tables_to_union) == 1:
            query = "select peak_id_a, peak_id_b from {}".format(tables_to_union[0])
        else:
            query = ""
        cursor.execute(query)

        records = [(str(record[0]), str(record[1])) for record in cursor.fetchall()]

        G = nx.OrderedDiGraph()
        G.add_edges_from(records)

        graphs = list(G.subgraph(c) for c in nx.weakly_connected_components(G))

        to_add = []
        for i, g in enumerate(graphs):
            for n in g.nodes():
                to_add.append([i+1, n, g.degree(n), g.number_of_nodes(), g.number_of_edges()])

        cursor.execute("""CREATE TEMP TABLE sub_groups (
                           sub_group_id INTEGER DEFAULT NULL,
                           peak_id INTEGER DEFAULT NULL,
                           degree INTEGER DEFAULT NULL,             
                           n_nodes INTEGER DEFAULT NULL,
                           n_edges INTEGER DEFAULT NULL,
                           PRIMARY KEY (sub_group_id, peak_id));""")

        cursor.executemany("""INSERT INTO sub_groups (sub_group_id, peak_id, degree, n_nodes, n_edges) 
                           VALUES (?,?,?,?,?)""", to_add)

        columns_groupings = ["peak_id", "group_id", "degree_cor", "sub_group_id", "degree", "n_nodes", "n_edges"]

        query_groupings = """SELECT DISTINCT gr.peak_id AS peak_id, gr.group_id AS group_id, degree_cor,
                          sub_groups.sub_group_id AS sub_group_id, sub_groups.degree AS degree,
                          sub_groups.n_nodes AS n_nodes, sub_groups.n_edges AS n_edges
                          FROM (SELECT group_id, peak_id_a AS peak_id, degree_a AS degree_cor FROM groups
                          UNION
                          SELECT group_id, peak_id_b AS peak_id, degree_b AS degree_cor FROM groups) AS gr
                          LEFT JOIN sub_groups
                          ON gr.peak_id = sub_groups.peak_id"""
    else:
        query_groupings = ""
        columns_groupings = []

    columns_adducts_oligo, columns_isotopes, columns_nls = [], [], []
    query_adducts_oligo, query_isotopes, query_nls = "", "", ""

    # flag_pp = len([tl for tl in tables_pp if (tl,) in tables]) > 0
    flag_adducts_oligo = ("adduct_pairs",) in tables or ("oligomers",) in tables
    flag_isotopes = ("isotopes",) in tables
    flag_neutral_losses = ("neutral_losses",) in tables

    if flag_adducts_oligo:
        sub_queries = []
        for tl in tables_pp:
            if (tl,) in tables:
                if tl == "adduct_pairs":
                    sub_queries.append("""select peak_id_a AS peak_id_pp, label_a AS label, charge_a AS charge, 1 AS oligomer from adduct_pairs
                    union
                    select peak_id_b AS peak_id_pp, label_b AS label, charge_b AS charge, 1 AS oligomer from adduct_pairs""")
                elif tl == "oligomers":
                    sub_queries.append("""select peak_id_a AS peak_id_pp, label_a AS label, charge_a AS charge, 1 AS oligomer from oligomers
                    union
                    select peak_id_b AS peak_id_pp, label_b AS label, charge_b AS charge, cast(round(mz_ratio) AS integer) AS oligomer from oligomers""")
        columns_adducts_oligo = ["label", "charge", "oligomer"]
        query_adducts_oligo = " union ".join(map(str, sub_queries))

    if flag_isotopes:
        columns_isotopes = ["isotope_labels_a", "isotope_ids", "isotope_labels_b", "isotope_charges", "atoms"]
        query_isotopes = """SELECT peak_id_a as isotope_peak_id_a, group_concat(label_a) AS isotope_labels_a,
                            group_concat(peak_id_b, ",") AS isotope_ids,
                            group_concat(label_b) AS isotope_labels_b, group_concat(charge) AS isotope_charges,
                            group_concat(round(atoms,1), ",") AS atoms
                            from (select peak_id_a, label_a, peak_id_b, label_b, charge, atoms, ppm_error from isotopes
                            union
                            select peak_id_b AS peak_id_a, label_b AS label_a,
                            peak_id_a AS peak_id_b, label_a AS label_b, charge, atoms, ppm_error
                            from isotopes                    
                            ) group by isotope_peak_id_a"""

    if flag_neutral_losses:
        columns_nls = ["nl_labels", "nl_ids"]
        query_nls = """SELECT peak_id_a as nl_peak_id_a, group_concat(label) AS nl_labels,
                            group_concat(peak_id_b, ",") AS nl_ids
                            from (select peak_id_a, label, peak_id_b, ppm_error 
                            from neutral_losses
                            union
                            select peak_id_b as peak_id_a, label, peak_id_a AS peak_id_b, ppm_error
                            from neutral_losses                            
                            ) group by nl_peak_id_a"""

    cursor.execute("DROP TABLE IF EXISTS peak_labels")
    sql_str_index = "CREATE INDEX idx_peak_id ON peak_labels (peak_id)"

    if flag_adducts_oligo or flag_isotopes or flag_neutral_losses:
        columns = ", ".join(map(str, columns_groupings + columns_adducts_oligo + columns_isotopes + columns_nls))
        query = "CREATE TABLE peak_labels AS "
        if query_groupings:
            query += "SELECT {} FROM ({})".format(columns, query_groupings)
        else:
            query += "SELECT peaklist.name AS peak_id, {} FROM peaklist ".format(columns)

        if flag_adducts_oligo:
            query += " LEFT JOIN ({}) ON peak_id = peak_id_pp ".format(query_adducts_oligo)
        if flag_isotopes:
            query += " LEFT JOIN ({}) ON peak_id = isotope_peak_id_a ".format(query_isotopes)
        if flag_neutral_losses:
            query += " LEFT JOIN ({}) ON peak_id = nl_peak_id_a ".format(query_nls)

        cursor.execute(query)
        cursor.execute(sql_str_index)

    if flag_adducts_oligo:

        # Add dummy row for features to make joint statement possible where label and adduct do not match.
        cursor.execute('PRAGMA table_info("peak_labels")')
        columns = cursor.fetchall()

        set_to_NULL = ["label", "charge", "oligomer"]
        columns_to_select = []
        for cn in columns:
            if cn[1] in set_to_NULL:
                columns_to_select.append("NULL")
            else:
                columns_to_select.append(cn[1])

        query = """INSERT INTO peak_labels
                   SELECT {} FROM peak_labels where label is not NULL""".format(", ".join(map(str, columns_to_select)))
        cursor.execute(query)
        conn.commit()

    cpd_tables = [tn[0] for tn in tables if "compound" in tn[0]]

    flag_mf = ("molecular_formulae",) in tables
    flag_cpd = len(cpd_tables) > 0

    columns = ["exact_mass", "ppm_error", "rt_diff", "adduct", "C", "H", "N", "O", "P", "S", "molecular_formula"]

    cpd_t_flags = []
    for i in range(0, len(cpd_tables)):
        b = [0] * len(cpd_tables)
        b[i] = 1
        cpd_t_flags.append(["{} AS {}".format(flag[0], flag[1]) for flag in list(zip(b, cpd_tables))])

    if len(cpd_tables) > 1:
        unions_cpd_sub_query = "LEFT JOIN (select *, {} from {}".format(", ".join(map(str, cpd_t_flags[0])), cpd_tables[0])
        for i, cpd_t in enumerate(cpd_tables[1:], start=1):
            unions_cpd_sub_query += " union select *, {} from {}".format(", ".join(map(str, cpd_t_flags[i])), cpd_t)
        unions_cpd_sub_query += ") AS ct "
    elif len(cpd_tables) == 1:
        unions_cpd_sub_query = "LEFT JOIN (select *, 1 AS {} from {}) AS ct".format(cpd_tables[0], cpd_tables[0])
    else:
        unions_cpd_sub_query = ""

    if flag_mf and flag_cpd:

        # unions_cpd_query = "CREATE TEMP TABLE compounds AS select * from "
        # unions_cpd_query += " union select * from ".join(map(str, cpd_tables))

        unions_cpd_query = "CREATE TEMP TABLE compounds AS select *, {} from {}".format(", ".join(map(str, cpd_t_flags[0])), cpd_tables[0])
        for i, cpd_t in enumerate(cpd_tables[1:], start=1):
            unions_cpd_query += " union select *, {} from {}".format(", ".join(map(str, cpd_t_flags[i])), cpd_t)

        cursor.execute(unions_cpd_query)
        unions_cpd_sub_query = ""

        query = """CREATE TEMP TABLE mf_cd as
                   SELECT mf.id, mf.exact_mass, mf.ppm_error, cpds.rt_diff, mf.adduct, 
                   mf.C, mf.H, mf.N, mf.O, mf.P, mf.S,
                   mf.molecular_formula, cpds.compound_name, cpds.compound_id, NULL AS compound_count, {}
                   FROM molecular_formulae AS mf
                   LEFT JOIN compounds AS cpds
                   ON mf.molecular_formula = cpds.molecular_formula AND mf.adduct = cpds.adduct
                   UNION
                   SELECT cpds.id, cpds.exact_mass, cpds.ppm_error, cpds.rt_diff, cpds.adduct, cpds.C, 
                   cpds.H, cpds.N, cpds.O, cpds.P, cpds.S,
                   cpds.molecular_formula, cpds.compound_name, cpds.compound_id, NULL AS compound_count, {}
                   FROM compounds AS cpds
                   LEFT JOIN molecular_formulae AS mf
                   ON mf.molecular_formula = cpds.molecular_formula AND mf.adduct = cpds.adduct
                   WHERE mf.molecular_formula IS NULL
                """.format(", ".join(map(str, ["cpds." + ct for ct in cpd_tables])),
                           ", ".join(map(str, ["cpds." + ct for ct in cpd_tables])))

        cursor.execute(query)
        mf_cpc_columns = "".join(map(str, [", mf_cd.{} AS {}".format(c, c) for c in columns]))
        mf_cpc_columns += ", mf_cd.compound_name AS compound_name, mf_cd.compound_id AS compound_id, NULL AS compound_count, "
        mf_cpc_columns += ", ".join(map(str, cpd_tables))
        if flag_adducts_oligo:
            union_mf_sub_query = "LEFT JOIN mf_cd ON (peaklist.name = mf_cd.id and peak_labels.label = mf_cd.adduct)"
            union_mf_sub_query += " OR (peaklist.name = mf_cd.id AND peak_labels.label is NULL and not exists (select 1 from peak_labels where peak_id = mf_cd.id and label = mf_cd.adduct))"
        else:
            union_mf_sub_query = "LEFT JOIN mf_cd ON peaklist.name = mf_cd.id"

    elif not flag_mf and flag_cpd:
        mf_cpc_columns = "".join(map(str,[", ct.{} AS {}".format(c, c) for c in columns]))
        mf_cpc_columns += ", compound_name AS compound_name, compound_id AS compound_id, NULL AS compound_count, "
        mf_cpc_columns += ", ".join(map(str, cpd_tables))
        if flag_adducts_oligo:
            unions_cpd_sub_query += " ON (peaklist.name = ct.id AND peak_labels.label = adduct)"
            unions_cpd_sub_query += " OR (peaklist.name = ct.id AND peak_labels.label is NULL and not exists (select 1 from peak_labels where peak_id = ct.id and label = ct.adduct))"
        else:
            unions_cpd_sub_query += " ON peaklist.name = ct.id"
        union_mf_sub_query = ""

    elif flag_mf and not flag_cpd:
        mf_cpc_columns = "".join(map(str, [", mf.{} AS {}".format(c, c) for c in columns if c != "rt_diff"]))
        if flag_adducts_oligo:
            union_mf_sub_query = "LEFT JOIN molecular_formulae AS mf"
            union_mf_sub_query += " ON (peaklist.name = mf.id AND peak_labels.label = mf.adduct)"
            union_mf_sub_query += " OR (peaklist.name = mf.id AND peak_labels.label is NULL and not exists (select 1 from peak_labels where peak_id = mf.id and label = mf.adduct))"
        else:
            union_mf_sub_query = "LEFT JOIN molecular_formulae AS mf"
            union_mf_sub_query += " ON peaklist.name = mf.id"
    else:
        mf_cpc_columns = ""
        union_mf_sub_query = ""

    cursor.execute('PRAGMA table_info("peak_labels")')
    columns_peak_labels = cursor.fetchall()

    if len(columns_peak_labels) == 0 and mf_cpc_columns == "":
        raise ValueError("No annotation results available to create summary from")

    exclude_cns = ["peak_id"]

    if len(columns_peak_labels) > 0:
        pl_columns = ", " + ", ".join(map(str, ["peak_labels.{}".format(cn[1]) for cn in columns_peak_labels if cn[1] not in exclude_cns]))
        join_peak_labels = """
                           LEFT JOIN
                           peak_labels
                           ON peaklist.name = peak_labels.peak_id
                           """
    else:
        pl_columns = ""
        join_peak_labels = ""

    sql_str_order = "ORDER BY peaklist.rowid"
    if ".label," in pl_columns:
        sql_str_order += ", label is NULL, label"
    if "isotope" in pl_columns:
        sql_str_order += ", isotope_labels_a is NULL, isotope_labels_a"
    if "ppm_error" in mf_cpc_columns:
        sql_str_order += ", abs(ppm_error) is NULL, abs(ppm_error)"
    if "compound_name" in mf_cpc_columns:
        sql_str_order += ", compound_name is NULL, compound_name"

    query = """
            CREATE TABLE summary AS SELECT distinct
            peaklist.name, peaklist.mz, peaklist.rt, peaklist.intensity{}{}
            FROM peaklist
            {}
            {}
            {}
            {}
            """.format(pl_columns, mf_cpc_columns, join_peak_labels, union_mf_sub_query, unions_cpd_sub_query, sql_str_order)

    cursor.execute("DROP TABLE IF EXISTS summary")
    # print(query)
    cursor.execute(query)
    conn.commit()

    # build where statement to remove dummy rows from the summary table # TODO: refactor code block
    cursor.execute('PRAGMA table_info("summary")')
    columns_summary = cursor.fetchall()

    where_str = ""
    for cn in columns_summary:
        if cn[1] in ["label", "isotope_labels_a", "neutral_loss_labels_a", "adduct"]:
            where_str += " AND {} is NULL".format(cn[1])

    query_d = """
              SELECT name
              FROM summary
              GROUP BY name
              HAVING COUNT(name) > 1
              """
    cursor.execute(query_d)
    r = cursor.fetchall()

    query_d = """
              SELECT name
              FROM summary
              WHERE name IS NOT NULL{}
              """.format(where_str)
    cursor.execute(query_d)
    rr = cursor.fetchall()

    query_d = """
              DELETE FROM summary
              WHERE name IN ("{}"){}
              """.format('","'.join(map(str, [name[0] for name in set(r) & set(rr)])), where_str)
    cursor.execute(query_d)
    conn.commit()

    if columns_groupings and flag_cpd:
        if "sub_group_id" in columns_groupings:
            grt = "sub_group_id"
        else:
            grt = "group_id"

        query = """            
                UPDATE summary
                SET compound_count = 
                (SELECT scs.c FROM
                    (SELECT {}, compound_id, COUNT(DISTINCT name) AS c FROM summary AS s
                    GROUP BY {}, compound_id) AS scs
                    WHERE scs.compound_id = summary.compound_id 
                    AND (scs.{} = summary.{} and scs.{} IS NOT NULL AND summary.{} IS NOT NULL)
                )
                """.format(grt, grt, grt, grt, grt, grt)

        cursor.execute(query)
        conn.commit()
        query = """            
                UPDATE summary
                SET compound_count = 1 
                WHERE compound_id is NOT NULL AND {} IS NULL
                """.format(grt)

        cursor.execute(query)
        conn.commit()
    elif flag_cpd:
        query = """  
                UPDATE summary
                SET compound_count = 
                (SELECT COUNT(DISTINCT name) FROM summary AS s
                    WHERE s.compound_id = summary.compound_id 
                    AND summary.compound_id IS NOT NULL
                    ) where summary.compound_id IS NOT NULL
                """
        cursor.execute(query)
    conn.commit()

    if single_row:

        columns_to_select = []
        if ("groups",) in tables:
            columns_to_select.extend(["group_id", "degree_cor", "sub_group_id", "degree", "n_nodes", "n_edges"])
        if ("adduct_pairs",) in tables or ("oligomers",) in tables:
            columns_to_select.append("""(select group_concat(label || '::' || charge || '::' || oligomer, '||')
            from (select distinct label, charge, oligomer from summary AS s where summary.name = s.name)
            ) AS label_charge_oligomer""")
        if ("isotopes",) in tables:
            columns_to_select.extend(["isotope_labels_a", "isotope_ids", "isotope_labels_b", "isotope_charges", "atoms"])
        if ("neutral_losses",) in tables:
            columns_to_select.extend(["nl_labels", "nl_ids"])

        if flag_cpd:
            if single_column:
                for cpd_t in cpd_tables:

                    cursor.execute("SELECT COUNT(*) FROM {} WHERE rt_diff is not NULL".format(cpd_t))
                    if int(cursor.fetchone()[0]) > 0:
                        rt_col = """|| '::' || ifnull(round(rt_diff, 2), "None")"""
                    else:
                        rt_col = ""

                    columns_to_select.append("""
                        group_concat(
                            CASE WHEN {} = 1 THEN 
                                molecular_formula || '::' || 
                                adduct || '::' || 
                                compound_name || '::' || 
                                compound_id  || '::' || 
                                compound_count || '::' || 
                                exact_mass || '::' || 
                                round(ppm_error, 2)
                                {}
                            ELSE NULL END 
                            , '||' 
                        ) AS {} 
                        """.format(cpd_t, rt_col, cpd_t))
            else:
                cursor.execute("SELECT COUNT(*) FROM summary WHERE rt_diff is not NULL")
                if int(cursor.fetchone()[0]) > 0:
                    rt_col = """, group_concat(ifnull(round(rt_diff, 2), "None"), '||') AS rt_diff"""
                else:
                    rt_col = ""

                columns_to_select.append("""
                    group_concat(molecular_formula, '||') AS molecular_formula,
                    group_concat(adduct, '||') AS adduct, 
                    group_concat(compound_name, '||') AS compound_name, 
                    group_concat(compound_id, '||') AS compound_id,
                    group_concat(compound_count, '||') AS compound_count,
                    group_concat(exact_mass, '||') AS exact_mass,
                    group_concat(round(ppm_error, 2), '||') AS ppm_error
                    {}
                    """.format(rt_col))
        elif flag_mf:
            if single_column:
                columns_to_select.append("""
                    group_concat(
                        molecular_formula || '::' || 
                        adduct || '::' || 
                        exact_mass || '::' || 
                        round(ppm_error, 2) ,
                        '||'
                    ) AS annotation
                    """)
            else:
                columns_to_select.append("""
                    group_concat(molecular_formula, '||') AS molecular_formula, 
                    group_concat(adduct, '||') AS adduct,
                    group_concat(exact_mass, '||') AS exact_mass,
                    group_concat(round(ppm_error, 2), '||') AS ppm_error
                    """)

        query = """
                SELECT DISTINCT name, mz, rt, intensity, {}
                from summary
                GROUP BY NAME
                ORDER BY rowid
                """.format(", ".join(map(str, columns_to_select)))

        df_out = pd.read_sql(query, conn)
        df_out.columns = [name.replace("peaklist.", "").replace("peak_labels.", "") for name in list(df_out.columns.values)]

        if flag_cpd:
            if not single_column:
                df_out["compound_id"] = df_out["compound_id"].replace({"None": ""})
                df_out["compound_name"] = df_out["compound_name"].replace({"None": ""})
            else:
                for cpd_t in cpd_tables:
                    df_out[cpd_t] = df_out[cpd_t].replace({"None": ""})
    else:
        df_out = pd.read_sql("select * from summary", conn)
        df_out.columns = [name.replace("peaklist.", "").replace("peak_labels.", "") for name in list(df_out.columns.values)]

    if convert_rt == "min" and "rt" in df_out.columns.values:
        rt_min = df_out["rt"] / 60.0
        df_out.insert(loc=df_out.columns.get_loc("rt")+1, column='rt_min', value=rt_min.round(2))
    elif convert_rt == "sec" and "rt" in df_out.columns.values:
        rt_sec = df_out["rt"] * 60.0
        df_out.insert(loc=df_out.columns.get_loc("rt")+1, column='rt_sec', value=rt_sec.round(1))
    elif convert_rt is not None:
        raise ValueError("Provide min, sec or None for convert_rt")

    if isinstance(ndigits_mz, int):
        df_out["mz"] = df_out["mz"].round(ndigits_mz)
    elif ndigits_mz is not None:
        raise ValueError("Provide integer or None for ndigits_mz")

    # Workaround for Pandas casting INT fo Float when Nan is present
    for c in df_out.columns:
        columns = ["charge", "oligomer", "group_id", "degree_cor", "sub_group_id", "degree", "n_nodes", "n_edges", "C", "H", "N", "O", "P", "S"]
        columns.extend(cpd_tables) # include compound tables

        if not single_row and not single_column:
            if c in columns:
                df_out[c] = df_out[c].astype('Int64')

    conn.close()
    return df_out
