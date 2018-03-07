#!/usr/bin/python
# -*- coding: utf-8 -*-

import itertools
import os
import networkx as nx
import sqlite3
from collections import OrderedDict
from in_out import read_molecular_formulae
from in_out import read_compounds
import numpy as np
from urlparse import urlparse
import requests
from in_out import mf_dict_to_str
import pandas as pd


def calculate_mz_tolerance(mass, ppm):
    min_tol = mass - (mass * 0.000001 * ppm)
    max_tol = mass + (mass * 0.000001 * ppm)
    return min_tol, max_tol


def calculate_ppm_error(mass, theo_mass):
    return float(theo_mass - mass) / (theo_mass * 0.000001)


def _prep_lib(lib):
    lib_pairs = []
    if isinstance(lib, OrderedDict):
        combs = list(itertools.combinations(lib, 2))
        for pair in combs:
            if isinstance(lib[pair[0]], float):
                lib_pairs.append(OrderedDict([(pair[0], {"mass": lib[pair[0]], "charge": 1}),
                                              (pair[1], {"mass": lib[pair[1]], "charge": 1})]))
            else:
                lib_pairs.append(OrderedDict([(pair[0], {"mass": lib[pair[0]]["mass"], "charge": lib[pair[0]]["charge"]}),
                                              (pair[1], {"mass": lib[pair[1]]["mass"], "charge": lib[pair[1]]["charge"]})]))
        lib_pairs = sorted(lib_pairs, key=lambda pair: (pair.items()[0][1]["mass"] - pair.items()[1][1]["mass"]), reverse=True)
        return lib_pairs
    elif isinstance(lib, list) and isinstance(lib[0], OrderedDict):
        if "mass_difference" in lib[0]:
            return sorted(lib, key=lambda d: d["mass_difference"], reverse=True)
        else:
            raise ValueError("Format library incorrect")
        #else:
        #    return sorted(lib_pairs, key=lambda pair: (pair.items()[0][1]["mass"] - pair.items()[1][1]["mass"]), reverse=True)
    else:
        raise ValueError("Incorrect format for library: {}".format(type(lib)))


def _annotate_artifacts(peaklist, diff=0.02):
    n = peaklist.iloc[:,1]
    for i in range(n):
        for j in range(i + 1, n):
            mz_diff = peaklist.iloc[i,1] - peaklist.iloc[j,1]
            ppm_error = calculate_ppm_error(peaklist.iloc[i,1], peaklist.iloc[j,1])
            if abs(mz_diff) < diff:
                yield i, j, mz_diff, ppm_error


def _check_tolerance(mz_x, mz_y, lib_pair, ppm):

    min_tol_a, max_tol_a = calculate_mz_tolerance(mz_x, ppm)
    min_tol_b, max_tol_b = calculate_mz_tolerance(mz_y, ppm)
    if "mass_difference" in lib_pair.keys():
        # Need to fix the order, charge is one
        min_tol_b = (min_tol_b - lib_pair["mass_difference"])
        max_tol_b = (max_tol_b - lib_pair["mass_difference"])
    elif "mass" in lib_pair.items()[0][1]:
        # Need to fix the order
        min_tol_a = (min_tol_a - lib_pair.items()[0][1]["mass"]) * lib_pair.items()[0][1]["charge"]
        max_tol_a = (max_tol_a - lib_pair.items()[0][1]["mass"]) * lib_pair.items()[0][1]["charge"]

        min_tol_b = (min_tol_b - lib_pair.items()[1][1]["mass"]) * lib_pair.items()[1][1]["charge"]
        max_tol_b = (max_tol_b - lib_pair.items()[1][1]["mass"]) * lib_pair.items()[1][1]["charge"]
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


def _annotate_pairs_from_graph(G, ppm, lib_pairs):

    for e in G.edges(data=True):
        #if G.node[e[0]]["mz"] < G.node[e[1]]["mz"]:
        #    mz_x = G.node[e[0]]["mz"]
        #    mz_y = G.node[e[1]]["mz"]
        #else:
        mz_x = G.node[e[0]]["mz"]
        mz_y = G.node[e[1]]["mz"]

        for lib_pair in lib_pairs:
            ct = _check_tolerance(mz_x, mz_y, lib_pair, ppm)
            if ct == 1 or ct == True:

                if "charge" in lib_pair.items()[0][1]:
                    charge_a = lib_pair.items()[0][1]["charge"]
                    charge_b = lib_pair.items()[1][1]["charge"]
                else:
                    charge_a = 1
                    charge_b = 1

                if "mass_difference" in lib_pair:
                    ppm_error = calculate_ppm_error(
                        mz_x,
                        mz_y - lib_pair["mass_difference"])
                else:
                    ppm_error = calculate_ppm_error(
                        (mz_x - lib_pair.items()[0][1]["mass"]) * charge_a,
                        (mz_y - lib_pair.items()[1][1]["mass"]) * charge_b)

                yield OrderedDict([("peak_id_a", e[0]), ("peak_id_b", e[1]),
                                   ("label_a", lib_pair.keys()[0]),
                                   ("label_b", lib_pair.keys()[1]),
                                   ('charge_a', charge_a),
                                   ('charge_b', charge_b),
                                   ('ppm_error', round(ppm_error, 2))])


def _annotate_pairs_from_peaklist(peaklist, ppm, lib_pairs):
    ct = 0
    n = len(peaklist.iloc[:,1])
    for i in range(n):
        for j in range(i + 1, n):

            for lib_pair in lib_pairs:
                ct = _check_tolerance(peaklist.iloc[i,1], peaklist.iloc[j,1], lib_pair, ppm)

                if ct == 1:

                    if "charge" in lib_pair.items()[0][1]:
                        charge_a = lib_pair.items()[0][1]["charge"]
                        charge_b = lib_pair.items()[1][1]["charge"]
                    else:
                        charge_a = 1
                        charge_b = 1

                    if "mass_difference" in lib_pair:
                        ppm_error = calculate_ppm_error(
                            peaklist.iloc[i,1],
                            peaklist.iloc[j,1] - lib_pair["mass_difference"])

                    else:
                        ppm_error = calculate_ppm_error(
                            (peaklist.iloc[i,1] - lib_pair.items()[0][1]["mass"]) * lib_pair.items()[0][1]["charge"],
                            (peaklist.iloc[j,1] - lib_pair.items()[1][1]["mass"]) * lib_pair.items()[1][1]["charge"])

                    yield OrderedDict([("peak_id_a", i), ("peak_id_b", j),
                                       ("label_a", lib_pair.keys()[0]),
                                       ("label_b", lib_pair.keys()[1]),
                                       ('charge_a', charge_a),
                                       ('charge_b', charge_b),
                                       ('ppm_error', round(ppm_error,2))])


class DbCompoundsMemory:

    def __init__(self, filename):

        self.filename = filename
        self.conn = sqlite3.connect(":memory:")
        self.cursor = self.conn.cursor()
        self.cursor.execute("""CREATE TABLE COMPOUNDS(
                            compound_id TEXT PRIMARY KEY  NOT NULL,
                            compound_name TEXT,
                            exact_mass decimal(15,7),
                            C int(3),
                            H int(3),
                            N int(3),
                            O int(3),
                            P int(3),
                            S int(3),
                            molecular_formula text DEFAULT NULL
                            );""")

        records = read_compounds(self.filename)
        for record in records:
            self.cursor.execute("""insert into COMPOUNDS ({}) values (?,?,?,?,?,?,?,?,?,?)""".format(
                           ",".join(map(str, record.keys()))), record.values())
        self.cursor.execute("""CREATE INDEX IDX_EXACT_MASS ON COMPOUNDS (exact_mass);""")
        self.conn.commit()

    def select_compounds(self, min_tol, max_tol):
        col_names = ["compound_id", "compound_name", "exact_mass", "C", "H", "N", "O", "P", "S", "molecular_formula"]
        self.cursor.execute("""SELECT {} FROM COMPOUNDS WHERE 
                            exact_mass >= {} and exact_mass <= {}
                            """.format(",".join(map(str, col_names)), min_tol, max_tol))
        return [OrderedDict(zip(col_names, list(record))) for record in self.cursor.fetchall()]

    def close(self):
        self.conn.close()


class DbMolecularFormulaeMemory:

    def __init__(self, filename):

        self.filename = filename
        self.conn = sqlite3.connect(":memory:")
        self.cursor = self.conn.cursor()
        self.cursor.execute("""CREATE TABLE MF(
                            ExactMass decimal(15,7),
                            C int(3),
                            H int(3),
                            N int(3),
                            O int(3),
                            P int(3),
                            S int(3),
                            HC INTEGER DEFAULT NULL,
                            NOPSC INTEGER DEFAULT NULL,
                            lewis INTEGER DEFAULT NULL,
                            senior INTEGER DEFAULT NULL,
                            DoubleBondEquivalents int(3),
                            primary key (C,H,N,O,P,S,ExactMass)
                            );""")

        records = read_molecular_formulae(self.filename)
        for record in records:
            self.cursor.execute("""insert into mf ({}) values (?,?,?,?,?,?,?,?,?,?,?)""".format(
                           ",".join(map(str, record.keys()))), record.values())

        self.cursor.execute("""CREATE INDEX IDX_EXACT_MASS ON MF (exactmass);""")
        self.cursor.execute("""CREATE INDEX IDX_EXACT_MASS_RULES ON MF (exactmass, HC, NOPSC, LEWIS, SENIOR);""")
        self.conn.commit()

    def select_mf(self, min_tol, max_tol, rules):

        if rules:
            sql_filters = " and lewis = 1 and senior = 1 and HC = 1 and NOPSC = 1"
        else:
            sql_filters = ""

        col_names = ["ExactMass", "C", "H", "N", "O", "P", "S",
                     "DoubleBondEquivalents", "LEWIS", "SENIOR", "HC", "NOPSC"]

        self.cursor.execute("""SELECT ExactMass, C, H, N, O, P, S,
                            DoubleBondEquivalents, LEWIS, SENIOR, HC, NOPSC
                            from mf where ExactMass >= {} and ExactMass <= {}{}
                            """.format(min_tol, max_tol, sql_filters))

        return [OrderedDict(zip(col_names, list(record))) for record in self.cursor.fetchall()]


def annotate_adducts(source, db_out, ppm, lib, add=False):

    conn = sqlite3.connect(db_out)
    cursor = conn.cursor()

    if not add:
        cursor.execute("DROP TABLE IF EXISTS adduct_pairs")

        cursor.execute("""CREATE TABLE adduct_pairs (
                       peak_id_a int(11) DEFAULT NULL,
                       peak_id_b int(11) DEFAULT NULL,
                       label_a char(15) DEFAULT NULL,
                       label_b char(15) DEFAULT NULL,
                       ppm_error float DEFAULT NULL,
                       PRIMARY KEY (peak_id_a, peak_id_b, label_a, label_b));""")

    lib_pairs = _prep_lib(lib.lib)

    if isinstance(source, nx.classes.digraph.DiGraph):
        source = list(nx.weakly_connected_component_subgraphs(source))

    if isinstance(source, list) and isinstance(source[0], nx.classes.digraph.DiGraph):
        for i, graph in enumerate(source):
            for assignment in _annotate_pairs_from_graph(graph, lib_pairs=lib_pairs, ppm=ppm):
                cursor.execute("""INSERT OR REPLACE into adduct_pairs (peak_id_a, peak_id_b, label_a, label_b, ppm_error)
                               values (?,?,?,?,?)""", (str(assignment["peak_id_a"]), str(assignment["peak_id_b"]),
                                                       assignment["label_a"], assignment["label_b"], float(assignment["ppm_error"])))

    elif isinstance(source, pd.core.frame.DataFrame):
        print lib_pairs
        raw_input()
        for assignment in _annotate_pairs_from_peaklist(source, lib_pairs=lib_pairs, ppm=ppm):
            cursor.execute("""INSERT OR REPLACE into adduct_pairs (peak_id_a, peak_id_b, label_a, label_b, ppm_error)
                           values (?,?,?,?,?)""", (source.iloc[assignment["peak_id_a"]][0], source.iloc[assignment["peak_id_b"]][0],
                                                   assignment["label_a"], assignment["label_b"], assignment["ppm_error"]))
    conn.commit()
    conn.close()
    return


def annotate_isotopes(source, db_out, ppm, lib):

    conn = sqlite3.connect(db_out)
    cursor = conn.cursor()

    cursor.execute("DROP TABLE IF EXISTS isotopes")

    cursor.execute("""CREATE TABLE isotopes (
                   peak_id_a int(11) DEFAULT NULL,
                   peak_id_b int(11) DEFAULT NULL,
                   label_a char(15) DEFAULT NULL,
                   label_b char(15) DEFAULT NULL,
                   atoms float DEFAULT NULL,
                   ppm_error float DEFAULT NULL,
                   PRIMARY KEY (peak_id_a, peak_id_b, label_a, label_b));""")

    lib_pairs = _prep_lib(lib.lib)

    abundances = {}
    for pair in lib.lib:
        abundances[pair.items()[0][0]] = pair.items()[0][1]
        abundances[pair.items()[1][0]] = pair.items()[1][1]

    if isinstance(source, nx.classes.digraph.DiGraph):
        source = list(nx.weakly_connected_component_subgraphs(source))

    if isinstance(source, list) and isinstance(source[0], nx.classes.digraph.DiGraph):

        for graph in source:

            peaklist = graph.nodes(data=True)

            for assignment in _annotate_pairs_from_graph(graph, lib_pairs=lib_pairs, ppm=ppm):

                y = abundances[assignment["label_a"]]['abundance'] * peaklist[assignment["peak_id_b"]]["intensity"]
                x = abundances[assignment["label_b"]]['abundance'] * peaklist[assignment["peak_id_a"]]["intensity"]
                atoms = (y / x)

                cursor.execute("""insert into isotopes (peak_id_a, peak_id_b, label_a, label_b, atoms, ppm_error)
                               values (?,?,?,?,?,?)""", (str(assignment["peak_id_a"]), str(assignment["peak_id_b"]),
                               assignment["label_a"], assignment["label_b"], float(atoms), float(assignment["ppm_error"])))

    elif isinstance(source, pd.core.frame.DataFrame):

        for assignment in _annotate_pairs_from_peaklist(source, lib_pairs=lib_pairs, ppm=ppm):

            y = abundances[assignment["label_a"]]["abundance"] * source.iloc[assignment["peak_id_b"]][3]
            x = abundances[assignment["label_b"]]["abundance"] * source.iloc[assignment["peak_id_a"]][3]
            atoms = y/x

            cursor.execute("""insert into isotopes (peak_id_a, peak_id_b, label_a, label_b, atoms, ppm_error)
                           values (?,?,?,?,?,?)""", (source.iloc[assignment["peak_id_a"]][0], source.iloc[assignment["peak_id_b"]][0],
                           assignment["label_a"], assignment["label_b"], atoms, assignment["ppm_error"]))
            conn.commit()

    conn.commit()
    conn.close()
    return


def annotate_oligomers(peaklist, db_out, ppm, lib, maximum=2):

    conn = sqlite3.connect(db_out)
    cursor = conn.cursor()

    cursor.execute("DROP TABLE IF EXISTS oligomers")

    cursor.execute("""CREATE TABLE oligomers (
                   peak_id_a int(11) DEFAULT NULL,
                   peak_id_b int(11) DEFAULT NULL,
                   mz_a float DEFAULT NULL,
                   mz_b float DEFAULT NULL,
                   label_a char(15) DEFAULT NULL,
                   label_b char(15) DEFAULT NULL,
                   mz_ratio float DEFAULT NULL,
                   ppm_error float DEFAULT NULL,
                   PRIMARY KEY (peak_id_a, peak_id_b));""")

    n = len(peaklist.iloc[:,0])
    for adduct in lib.lib.keys():

        for i in range(n):

            for d in range(1, maximum):
                for j in range(i + 1, n):

                    min_tol_a, max_tol_a = calculate_mz_tolerance(peaklist.iloc[i][1] + ((peaklist.iloc[i][1] - lib.lib[adduct]) * d), ppm)
                    min_tol_b, max_tol_b = calculate_mz_tolerance(peaklist.iloc[j][1], ppm)

                    if (min_tol_b > max_tol_a and max_tol_b > max_tol_a):# or (min_tol_a < min_tol_b and max_tol_a < min_tol_b):
                        #print peaklist.iloc[i][1], peaklist.iloc[j][1], adduct
                        break

                    min_tol_a = min_tol_a - lib.lib[adduct]
                    max_tol_a = max_tol_a - lib.lib[adduct]

                    min_tol_b = min_tol_b - lib.lib[adduct]
                    max_tol_b = max_tol_b - lib.lib[adduct]

                    if min_tol_a < max_tol_b and min_tol_b < max_tol_a:

                        a = (peaklist.iloc[i][1] - lib.lib[adduct]) + (peaklist.iloc[i][1] - lib.lib[adduct]) * d
                        b = peaklist.iloc[j][1] - lib.lib[adduct]

                        ratio = (peaklist.iloc[j][1] - lib.lib[adduct]) / (peaklist.iloc[i][1] - lib.lib[adduct])
                        ppm_error = calculate_ppm_error(a, b)

                        if "M" in adduct:
                            adduct_oligo = adduct.replace("M", "{}M".format(int(round(ratio))))
                        else:
                            adduct_oligo = "{}{}".format(int(round(ratio)), adduct)

                        cursor.execute("""insert into oligomers (peak_id_a, peak_id_b, mz_a, mz_b, label_a, label_b, mz_ratio, ppm_error)
                                       values (?,?,?,?,?,?,?,?)""", (i, j, peaklist.iloc[i][1], peaklist.iloc[j][1], adduct, adduct_oligo, ratio, ppm_error))
    conn.commit()
    conn.close()
    return


def annotate_artifacts(source, db_out, diff):

    conn = sqlite3.connect(db_out)
    cursor = conn.cursor()

    cursor.execute("DROP TABLE IF EXISTS artifacts")

    cursor.execute("""CREATE TABLE artifacts (
                   peak_id_a int(11) DEFAULT NULL,
                   peak_id_b int(11) DEFAULT NULL,
                   mz_diff float DEFAULT NULL,
                   ppm_error float DEFAULT NULL,
                   PRIMARY KEY (peak_id_a, peak_id_b));""")

    if isinstance(source, nx.classes.digraph.DiGraph):
        source = list(nx.weakly_connected_component_subgraphs(source))

    if (isinstance(source, list) or isinstance(source, np.ndarray)) and isinstance(source[0], nx.classes.graph.Graph):
        for graph in source:
            peaklist = graph.nodes(data=True)
            for assignment in _annotate_artifacts(peaklist, diff=diff):
                cursor.execute("""insert into artifacts (peak_id_a, peak_id_b, mz_diff, ppm_error)
                               values (?,?,?,?)""", (source[assignment["peak_id_a"]][0], source[assignment["peak_id_b"]][0], assignment["label_a"], assignment["label_b"]))

    elif isinstance(source, pd.core.frame.DataFrame):
        for assignment in _annotate_artifacts(source, diff=diff):
            cursor.execute("""insert into artifacts (peak_id_a, peak_id_b, mz_diff, ppm_error)
                           values (?,?,?,?)""", (source[assignment["peak_id_a"]][0], source[assignment["peak_id_b"]][0], assignment["label_a"], assignment["label_b"]))

    conn.commit()
    return


def annotate_multiple_charged_ions(source, db_out, ppm, lib, add=False):

    conn = sqlite3.connect(db_out)
    cursor = conn.cursor()

    if not add:
        cursor.execute("DROP TABLE IF EXISTS multiple_charged_ions")

        cursor.execute("""CREATE TABLE multiple_charged_ions (
                       peak_id_a int(11) DEFAULT NULL,
                       peak_id_b int(11) DEFAULT NULL,
                       label_a char(15) DEFAULT NULL,
                       label_b char(15) DEFAULT NULL,
                       charge_a int(11) DEFAULT NULL,
                       charge_b int(11) DEFAULT NULL,
                       ppm_error float DEFAULT NULL,
                       PRIMARY KEY (peak_id_a, peak_id_b, label_a, label_b, charge_a, charge_b));""")

    lib_pairs = _prep_lib(lib.lib)

    if isinstance(source, nx.classes.digraph.DiGraph):
        source = list(nx.weakly_connected_component_subgraphs(source))

    if (isinstance(source, list) or isinstance(source, np.ndarray)) and isinstance(source[0], nx.classes.graph.Graph):
        for graph in source:
            for assignment in _annotate_pairs_from_graph(graph, lib_pairs=lib_pairs, ppm=ppm):
                cursor.execute("""INSERT OR REPLACE into multiple_charged_ions (peak_id_a, peak_id_b, label_a, label_b, charge_a, charge_b, ppm_error)
                               values (?,?,?,?,?,?,?)""", (assignment["peak_id_a"], assignment["peak_id_b"], assignment["label_a"], assignment["label_b"],
                                                           assignment["charge_a"], assignment["charge_b"], assignment["ppm_error"]))

    elif isinstance(source, pd.core.frame.DataFrame):
        for assignment in _annotate_pairs_from_peaklist(source, lib_pairs=lib_pairs, ppm=ppm):
            cursor.execute("""INSERT OR REPLACE into multiple_charged_ions (peak_id_a, peak_id_b, label_a, label_b, charge_a, charge_b, ppm_error)
                           values (?,?,?,?,?,?,?)""", (source[assignment["peak_id_a"]][0], source[assignment["peak_id_b"]][0],
                                                       assignment["label_a"], assignment["label_b"], assignment["charge_a"], assignment["charge_b"], assignment["ppm_error"]))
    conn.commit()
    conn.close()
    return


def annotate_molecular_formulae(peaklist, lib_adducts, ppm, db_out, db_in="http://multiomics-int.cs.bham.ac.uk", rules=True, max_mz=None):

    conn = sqlite3.connect(db_out)
    cursor = conn.cursor()

    cursor.execute("DROP TABLE IF EXISTS molecular_formulae")

    cursor.execute("""CREATE TABLE molecular_formulae (
                    id int(11) DEFAULT NULL,
                    mz decimal(12,7) DEFAULT NULL,
                    exact_mass decimal(12,7) DEFAULT NULL,
                    ppm_error decimal(12,2) DEFAULT NULL,
                    adduct text DEFAULT NULL,
                    C INTEGER DEFAULT NULL,
                    H INTEGER DEFAULT NULL,
                    N INTEGER DEFAULT NULL,
                    O INTEGER DEFAULT NULL,
                    P INTEGER DEFAULT NULL,
                    S INTEGER DEFAULT NULL,
                    molecular_formula text DEFAULT NULL,
                    HC INTEGER DEFAULT NULL,
                    NOPSC INTEGER DEFAULT NULL,
                    lewis INTEGER DEFAULT NULL,
                    senior INTEGER DEFAULT NULL,
                    double_bond_equivalents INTEGER DEFAULT NULL,
                    primary key  (mz, adduct, C, H, N, O, P, S)
                    );""")

    if os.path.isfile(db_in):
        conn_mem = DbMolecularFormulaeMemory(db_in)
    else:
        url = '{}/mass_range'.format(db_in)
        url_test = '{}//mass?mass=180.06339&tol=0.0&unit=ppm&rules=1'.format(db_in)
        o = urlparse(url)
        if o.scheme != "http" and o.netloc != "multiomics-int.cs.bham.ac.uk":
            raise ValueError("No database or local db available")
        else:
            r = requests.get(url_test)
            r.raise_for_status()

    for i in range(len(peaklist.iloc[:, 0])):
        mz = float(peaklist["mz"].iloc[i])
        name = str(peaklist["name"].iloc[i])

        min_tol, max_tol = calculate_mz_tolerance(mz, ppm)

        if max_mz is not None and mz > max_mz:  # TODO
            continue

        for adduct in lib_adducts.lib:

            if mz - lib_adducts.lib[adduct] > 0.5:

                if "conn_mem" in locals():
                    records = conn_mem.select_mf(min_tol - lib_adducts.lib[adduct], max_tol - lib_adducts.lib[adduct], rules)
                else:
                    params = {"lower": min_tol - lib_adducts.lib[adduct], "upper": max_tol - lib_adducts.lib[adduct], "rules": int(rules)}
                    response = requests.get(url, params=params)
                    records = response.json()["records"]

                for record in records:
                    record["id"] = name
                    record["exact_mass"] = record["ExactMass"] + lib_adducts.lib[adduct]
                    record["double_bond_equivalents"] = record["DoubleBondEquivalents"]
                    del record["DoubleBondEquivalents"]
                    record["mz"] = mz
                    record["ppm_error"] = calculate_ppm_error(mz, record["exact_mass"])
                    del record["ExactMass"]
                    record["molecular_formula"] = mf_dict_to_str(record)
                    record["adduct"] = adduct
                    cursor.execute("""insert into molecular_formulae ({}) values (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)""".format(
                                   ",".join(map(str, record.keys()))), record.values())
    conn.commit()
    conn.close()
    return


def annotate_compounds(peaklist, lib_adducts, ppm, db_out, db_in, db_name):

    conn = sqlite3.connect(db_out)
    cursor = conn.cursor()

    cursor.execute("DROP TABLE IF EXISTS compounds_{}".format(db_name))
    cursor.execute("""CREATE TABLE compounds_{} (
                   id int(11) DEFAULT NULL,
                   mz decimal(12,7) DEFAULT NULL,
                   exact_mass decimal(12,7) DEFAULT NULL,
                   ppm_error decimal(12,2) DEFAULT NULL,
                   adduct text DEFAULT NULL,
                   C INTEGER DEFAULT NULL,
                   H INTEGER DEFAULT NULL,
                   N INTEGER DEFAULT NULL,
                   O INTEGER DEFAULT NULL,
                   P INTEGER DEFAULT NULL,
                   S INTEGER DEFAULT NULL,
                   molecular_formula text DEFAULT NULL,
                   compound_id char(6) DEFAULT NULL,
                   compound_name text DEFAULT NULL,
                   primary key (id, adduct, compound_id)
                   );""".format(db_name))

    if os.path.isfile(db_in):
        with open(db_in, 'rb') as fd:
            if fd.read(100)[:16] == 'SQLite format 3\x00':
                conn_local = sqlite3.connect(db_in)
                cursor_local = conn_local.cursor()
            else:
                conn_mem = DbCompoundsMemory(db_in)
    else:
        raise IOError("[Errno 2] No such file or directory: {}".format(db_in))

    for i in range(len(peaklist.iloc[:, 0])):
        mz = float(peaklist["mz"].iloc[i])
        name = str(peaklist["name"].iloc[i])
        min_tol, max_tol = calculate_mz_tolerance(mz, ppm)

        for adduct in lib_adducts.lib:

            if mz - lib_adducts.lib[adduct] > 0.5:

                if "conn_mem" in locals():
                    records = conn_mem.select_compounds(min_tol - lib_adducts.lib[adduct], max_tol - lib_adducts.lib[adduct])
                elif "conn_local" in locals():
                    col_names = ["compound_id", "C", "H", "N", "O", "P", "S", "molecular_formula", "compound_name", "exact_mass"]
                    cursor_local.execute("""SELECT ID, C, H, N, O, P, S,
                                            molecular_formula, name, exact_mass
                                            from {} where exact_mass >= {} and exact_mass <= {}
                                            """.format(db_name, min_tol - lib_adducts.lib[adduct], max_tol - lib_adducts.lib[adduct]))
                    records = [OrderedDict(zip(col_names, list(record))) for record in cursor_local.fetchall()]

                for record in records:
                    record["id"] = name
                    record["exact_mass"] = record["exact_mass"] +float(lib_adducts.lib[adduct])
                    record["mz"] = mz
                    record["ppm_error"] = calculate_ppm_error(mz, record["exact_mass"])
                    record["adduct"] = adduct
                    cursor.execute("""insert into compounds_{} ({}) values (?,?,?,?,?,?,?,?,?,?,?,?,?,?)
                                   """.format(db_name, ",".join(map(str, record.keys()))), record.values())
    conn.commit()
    conn.close()
    return


def summary(df, db):

    conn = sqlite3.connect(db)
    cursor = conn.cursor()

    cursor.execute("DROP TABLE IF EXISTS peaklist")
    df[['name', 'mz', 'rt', "intensity"]].sort_values(by=["rt", "mz"]).to_sql('peaklist', conn, index=False)

    cursor.execute("SELECT name FROM sqlite_master WHERE type='table'")
    tables = cursor.fetchall()

    tables_amo = ["adduct_pairs", "multiple_charged_ions", "oligomers", "isotopes"]
    tables_to_union = []
    for tn in tables:
        if tn[0] in tables_amo:
            tables_to_union.append(str(tn[0]))

    if len(tables_to_union) > 0 and ("groups",) in tables:

        if len(tables_to_union) > 1:
            query = "select peak_id_a, peak_id_b from "
            query += " union select peak_id_a, peak_id_b from ".join(map(str, tables_to_union))
        elif len(tables_to_union) == 1:
            query = "select peak_id_a, peak_id_b from {}".format(tables_to_union[0])
        cursor.execute(query)

        records = [(str(record[0]), str(record[1])) for record in cursor.fetchall()]

        G = nx.OrderedDiGraph()
        G.add_edges_from(records)

        graphs = list(nx.weakly_connected_component_subgraphs(G))

        to_add = []
        for i, g in enumerate(graphs):
            for n in g.nodes():
                to_add.append([i+1, n, g.degree(n), g.number_of_nodes(), g.number_of_edges()])

        cursor.execute("""CREATE TEMPORARY TABLE sub_groups (
                           sub_group_id int(11) DEFAULT NULL,
                           peak_id int(11) DEFAULT NULL,
                           degree int(11) DEFAULT NULL,             
                           n_nodes int(11) DEFAULT NULL,
                           n_edges int(11) DEFAULT NULL,
                           PRIMARY KEY (sub_group_id, peak_id));""")

        cursor.executemany("""insert into sub_groups (sub_group_id, peak_id, degree, n_nodes, n_edges) 
                           values (?,?,?,?,?)""", to_add)

        columns_groupings = """peak_id, group_id, degree_cor, sub_group_id, degree, n_nodes, n_edges"""

        query_groupings = """select distinct gr.peak_id as peak_id, gr.group_id as group_id, degree_cor,
                          sub_groups.sub_group_id as sub_group_id, sub_groups.degree as degree,
                          sub_groups.n_nodes as n_nodes, sub_groups.n_edges as n_edges
                          from (select group_id, peak_id_a as peak_id, degree_a as degree_cor from groups
                          union
                          select group_id, peak_id_b as peak_id, degree_b as degree_cor from groups) AS gr
                          LEFT JOIN sub_groups
                          ON gr.peak_id = sub_groups.peak_id"""
    else:
        query_groupings = ""
        columns_groupings = ""

    flag_amo = len([tl for tl in tables_amo if (tl,) in tables]) > 0
    flag_isotopes = ("isotopes",) in tables

    if flag_amo:
        sub_queries = []
        for tl in tables_amo:
            if (tl,) in tables:
                if tl == "adduct_pairs":
                    sub_queries.append("""select peak_id_a as peak_id_amo, label_a as label, 1 as charge, 1 as oligomer from adduct_pairs
                    union
                    select peak_id_b as peak_id_amo, label_b as label, 1 as charge, 1 as oligomer from adduct_pairs""")
                elif tl == "multiple_charged_ions":
                    sub_queries.append("""select peak_id_a as peak_id_amo, label_a as label, charge_a as charge, 1 as oligomer from multiple_charged_ions
                    union
                    select peak_id_b as peak_id_amo, label_b as label, charge_b as charge, 1 as oligomer from multiple_charged_ions""")
                elif tl == "oligomers":
                    sub_queries.append("""select peak_id_a as peak_id_amo, label_a as label, 1 as charge, 1 as oligomer from oligomers
                    union
                    select peak_id_b as peak_id_amo, label_b as label, 1 as charge, round(mz_ratio) as oligomer from oligomers""")
        columns_amo = ", label, charge, oligomer"
        query_amo = " union ".join(map(str, sub_queries))

    if flag_isotopes:
        columns_isotopes = ", isotope_labels_a, isotope_ids, isotope_labels_b, atoms"
        query_isotopes = """SELECT peak_id_a, group_concat(label_a) as isotope_labels_a,
                            group_concat(peak_id_b, ",") as isotope_ids,
                            group_concat(label_b) as isotope_labels_b, group_concat(round(atoms,1), ",") as atoms
                            from (select peak_id_a, label_a, peak_id_b, label_b, atoms, ppm_error from isotopes
                            union
                            select peak_id_b as peak_id_a, label_b as label_a,
                            peak_id_a as peak_id_b, label_a as label_b, atoms, ppm_error
                            from isotopes
                            ) group by peak_id_a"""

    cursor.execute("DROP TABLE IF EXISTS peak_labels")
    if flag_amo and flag_isotopes:
        query = "CREATE TABLE peak_labels as "
        if query_groupings != "":
            query += "SELECT {}{}{} from """.format(columns_groupings, columns_amo, columns_isotopes)
            query += "({}) LEFT JOIN ({}) ON peak_id = peak_id_amo LEFT JOIN ({}) ON peak_id = peak_id_a".format(query_groupings, query_amo, query_isotopes)
        else:
            query += "SELECT peaklist.name as peak_id{}{} from ".format(columns_amo, columns_isotopes)
            query += "peaklist LEFT JOIN ({}) ON peaklist.name = peak_id LEFT JOIN ({}) ON peaklist.name = peak_id_a".format(query_amo, query_isotopes)
            query = query.replace("peak_id_amo", "peak_id")
        cursor.execute(query)
    elif flag_isotopes and not flag_amo:
        query = "CREATE TABLE peak_labels as "
        if query_groupings != "":
            query += "select {}{} from ".format(columns_groupings, columns_isotopes)
            query += "({}) LEFT JOIN ({}) ON peak_id = peak_id_a".format(query_groupings, query_isotopes)
        else:
            query += query_isotopes
        cursor.execute(query)
    elif not flag_isotopes and flag_amo:
        query = "CREATE TABLE peak_labels as "
        if query_groupings != "":
            query += """select {}{} from """.format(columns_groupings, columns_amo)
            query += """({}) LEFT JOIN ({}) ON peak_id = peak_id_amo""".format(query_groupings, query_amo)
        else:
            query += query_amo.replace("peak_id_amo", "peak_id")
        cursor.execute(query)

    if flag_amo:
        cursor.execute('PRAGMA table_info("peak_labels")')
        columns = cursor.fetchall()

        set_to_NULL = ["label", "charge", "oligomer"]
        columns_to_select = []
        for cn in columns:
            if cn[1] in set_to_NULL:
                columns_to_select.append("NULL")
            else:
                columns_to_select.append(cn[1])

        query = "INSERT INTO peak_labels"
        query += " SELECT {} FROM peak_labels where label is not NULL".format(", ".join(map(str, columns_to_select)))

        cursor.execute(query)
        conn.commit()

    cpd_tables = [tn[0] for tn in tables if "compound" in tn[0]]

    flag_mf = ("molecular_formulae",) in tables
    flag_cpd = len(cpd_tables) > 0

    columns = ["exact_mass", "ppm_error", "adduct", "C", "H", "N", "O", "P", "S", "molecular_formula"]

    if len(cpd_tables) > 1:
        unions_cpd_sub_query = "LEFT JOIN (select * from "
        unions_cpd_sub_query += " union select * from ".join(map(str, cpd_tables))
        unions_cpd_sub_query += ") as ct "
    elif len(cpd_tables) == 1:
        unions_cpd_sub_query = "LEFT JOIN (select * from {}) as ct".format(cpd_tables[0])
    else:
        unions_cpd_sub_query = ""

    if flag_mf and flag_cpd:
        mf_cpc_columns = "".join(map(str, [", mf.{} as {}".format(c, c) for c in columns]))
        mf_cpc_columns += ", ct.compound_name as compound_name, ct.compound_id as compound_id"
        unions_cpd_sub_query += " ON mf.molecular_formula = ct.molecular_formula AND mf.adduct = ct.adduct"
        if flag_amo:
            union_mf_sub_query = "LEFT JOIN molecular_formulae AS mf ON (peaklist.name = mf.id and peak_labels.label = mf.adduct)"
            union_mf_sub_query += " OR (peaklist.name = mf.id AND peak_labels.label is NULL)"
        else:
            union_mf_sub_query = "LEFT JOIN molecular_formulae AS mf ON peaklist.name = mf.id"

    elif not flag_mf and flag_cpd:
        mf_cpc_columns = "".join(map(str,[", ct.{} as {}".format(c, c) for c in columns]))
        mf_cpc_columns += ", compound_name as compound_name, compound_id as compound_id"
        if flag_amo:
            unions_cpd_sub_query += " ON (peaklist.name = ct.id AND peak_labels.label = adduct)"
            unions_cpd_sub_query += " OR (peaklist.name = ct.id AND peak_labels.label is NULL)"
        else:
            unions_cpd_sub_query += " ON peaklist.name = ct.id"
        union_mf_sub_query = ""

    elif flag_mf and not flag_cpd:
        mf_cpc_columns = "".join(map(str, [", mf.{} as {}".format(c, c) for c in columns]))
        if flag_amo:
            union_mf_sub_query = "LEFT JOIN molecular_formulae AS mf"
            union_mf_sub_query += " ON (peaklist.name = mf.id AND peak_labels.label = mf.adduct)"
            union_mf_sub_query += " OR (peaklist.name = mf.id AND peak_labels.label is NULL)"
        else:
            union_mf_sub_query = "LEFT JOIN molecular_formulae AS mf"
            union_mf_sub_query += " ON peaklist.name = mf.id"
    else:
        mf_cpc_columns = ""
        union_mf_sub_query = ""

    cursor.execute('PRAGMA table_info("peak_labels")')
    columns = cursor.fetchall()

    if len(columns) == 0 and mf_cpc_columns == "":
        raise ValueError("No annotation results available to create summary from")

    exclude_cns = ["peak_id"]
    pl_columns = ", " + ", ".join(map(str, ["peak_labels.{}".format(cn[1]) for cn in columns if cn[1] not in exclude_cns]))
    query = """CREATE TABLE summary AS SELECT
               peaklist.name, peaklist.mz, peaklist.rt{}{}
               FROM peaklist
               LEFT JOIN
               peak_labels
               ON peaklist.name = peak_labels.peak_id
               {}
               {}
               ORDER BY peaklist.rt""".format(pl_columns, mf_cpc_columns, union_mf_sub_query, unions_cpd_sub_query)
    cursor.execute("DROP TABLE IF EXISTS summary")
    cursor.execute(query)

    conn.commit()

    df_out = pd.read_sql("select * from summary", conn)
    df_out.columns = [name.replace("peaklist.", "").replace("peak_labels.", "") for name in list(df_out.columns.values)]


    conn.close()
    return df_out

