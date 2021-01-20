#!/usr/bin/env python
# -*- coding: utf-8 -*-

import copy
import os
import collections
import numpy as np
from pandas import read_csv
import pandas as pd
from pyteomics import mass as pyteomics_mass
from beamspy import libraries
from beamspy.auxiliary import nist_database_to_pyteomics
from beamspy.auxiliary import order_composition_by_hill
from beamspy.auxiliary import composition_to_string
from beamspy.auxiliary import double_bond_equivalents
from beamspy.auxiliary import HC_HNOPS_rules
from beamspy.auxiliary import lewis_senior_rules


def read_adducts(filename, ion_mode, separator="\t"):
    df = read_csv(filename, sep=separator, float_precision="round_trip")
    adducts = libraries.Adducts()
    adducts.remove("*")
    for index, row in df.iterrows():
        if "ion_mode" not in row:
            adducts.add(row["label"], row["exact_mass"], row["charge"])
        elif (row["ion_mode"] == "pos" or row["ion_mode"] == "both") and ion_mode == "pos":
            adducts.add(row["label"], row["exact_mass"], row["charge"])
        elif (row["ion_mode"] == "neg" or row["ion_mode"] == "both") and ion_mode == "neg":
            adducts.add(row["label"], row["exact_mass"], row["charge"])
    return adducts


def read_isotopes(filename, ion_mode, separator="\t"):
    df = read_csv(filename, sep=separator, float_precision="round_trip")
    isotopes = libraries.Isotopes()
    isotopes.remove("*")
    for index, row in df.iterrows():
        if "ion_mode" not in row:
            isotopes.add(row["label_x"], row["label_y"], row["abundance_x"], row["abundance_y"],
                         row["mass_difference"], row["charge"])
        elif (row["ion_mode"] == "pos" or row["ion_mode"] == "both") and ion_mode == "pos":
            isotopes.add(row["label_x"], row["label_y"], row["abundance_x"], row["abundance_y"],
                         row["mass_difference"], row["charge"])
        elif (row["ion_mode"] == "neg" or row["ion_mode"] == "both") and ion_mode == "neg":
            isotopes.add(row["label_x"], row["label_y"], row["abundance_x"], row["abundance_y"],
                         row["mass_difference"], row["charge"])
    return isotopes


def read_molecular_formulae(filename, separator="\t", calculate=True, filename_atoms=""):

    if calculate:
        path_nist_database = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'nist_database.txt')
        nist_database = nist_database_to_pyteomics(path_nist_database)

    df = read_csv(filename, sep=separator, float_precision="round_trip")
    records = []
    for index, row in df.iterrows():
        record = collections.OrderedDict()
        comp = pyteomics_mass.Composition(str(row.molecular_formula))
        if comp:
            record["composition"] = collections.OrderedDict((k, comp[k]) for k in order_composition_by_hill(comp.keys()))
            sum_CHNOPS = sum([comp[e] for e in comp if e in ["C", "H", "N", "O", "P", "S"]])
            record["CHNOPS"] = sum_CHNOPS == sum(list(comp.values()))
            if calculate:
                record["exact_mass"] = round(pyteomics_mass.mass.calculate_mass(formula=str(row.molecular_formula), mass_data=nist_database), 6)
            else:
                record["exact_mass"] = float(row.exact_mass)
            record.update(HC_HNOPS_rules(str(row.molecular_formula)))
            record.update(lewis_senior_rules(str(row.molecular_formula)))
            record["double_bond_equivalents"] = double_bond_equivalents(record["composition"])
            records.append(record)
        else:
            Warning("{} Skipped".format(row))

    return records


def read_compounds(filename, separator="\t", calculate=True, lib_adducts=[], filename_atoms=""):

    if calculate:
        path_nist_database = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'nist_database.txt')
        nist_database = nist_database_to_pyteomics(path_nist_database)

    df = read_csv(filename, sep=separator, float_precision="round_trip")
    records = []
    for index, row in df.iterrows():
        record = collections.OrderedDict()
        comp = pyteomics_mass.Composition(str(row.molecular_formula))
        if comp:
            record["composition"] = collections.OrderedDict((k, comp[k]) for k in order_composition_by_hill(comp.keys()))
            sum_CHNOPS = sum([comp[e] for e in comp if e in ["C", "H", "N", "O", "P", "S"]])
            record["CHNOPS"] = sum_CHNOPS == sum(list(comp.values()))
            if calculate:
                record["exact_mass"] = round(pyteomics_mass.calculate_mass(formula=str(str(row.molecular_formula)), mass_data=nist_database),6)
            else:
                record["exact_mass"] = float(row.exact_mass)

            record["compound_id"] = row.compound_id
            record["compound_name"] = row.compound_name
            comp = pyteomics_mass.Composition(str(row.molecular_formula))
            record["molecular_formula"] = composition_to_string(comp)

            if "retention_time" in df.columns:
                record["retention_time"] = row.retention_time
            elif "rt" in df.columns:
                record["retention_time"] = row.rt
            if "adduct" in df.columns:
                record["adduct"] = row.adduct
                if lib_adducts and calculate:
                    record["exact_mass"] += lib_adducts.lib[row.adduct]["mass"]

            records.append(record)
        else:
            Warning("{} Skipped".format(row))

    return records


def read_mass_differences(filename, ion_mode, separator="\t"):
    df = read_csv(filename, sep=separator, float_precision="round_trip")
    mass_differences = libraries.MassDifferences()
    for index, row in df.iterrows():
        if "charge_x" in row:
            charge_x = row["charge_x"]
            charge_y = row["charge_y"]
        else:
            charge_x = 1
            charge_y = 1
        if "ion_mode" not in row:
            mass_differences.add(row["label_x"], row["label_y"], row["mass_difference"], charge_x, charge_y)
        elif (row["ion_mode"] == "pos" or row["ion_mode"] == "both") and ion_mode == "pos":
            mass_differences.add(row["label_x"], row["label_y"], row["mass_difference"], charge_x, charge_y)
        elif (row["ion_mode"] == "neg" or row["ion_mode"] == "both") and ion_mode == "neg":
            mass_differences.add(row["label_x"], row["label_y"], row["mass_difference"], charge_x, charge_y)
    return mass_differences


def read_neutral_losses(filename, separator="\t"):
    df = read_csv(filename, sep=separator, float_precision="round_trip")
    nls = libraries.NeutralLosses()
    for index, row in df.iterrows():
        nls.add(row["label"], row["mass_difference"])
    return nls


def read_xset_matrix(fn_matrix, first_sample, separator="\t", mapping={"mz": "mz", "rt": "rt", "name": "name"}, samples_in_columns=True):
    if "mz" not in mapping and "rt" not in mapping and "name" not in mapping:
        raise ValueError("Incorrect column mapping: provide column names for mz, and name")

    df = pd.read_csv(fn_matrix, header=0, sep=separator, dtype={"name": str}, float_precision="round_trip")
    df.replace(0, np.nan, inplace=True)

    if not samples_in_columns:
        df = df.T

    df_peaklist = df[[mapping["name"], mapping["mz"], mapping["rt"]]]
    df_matrix = df.iloc[:, df.columns.get_loc(first_sample):]
    df_peaklist = df_peaklist.assign(intensity=pd.Series(df_matrix.median(axis=1, skipna=True).values))
    df_peaklist.columns = ["name", "mz", "rt", "intensity"]
    return pd.concat([df_peaklist, df_matrix], axis=1)


def combine_peaklist_matrix(fn_peaklist, fn_matrix, separator="\t", median_intensity=True,
                            mapping={"name": "name", "mz": "mz", "rt": "rt", "intensity": "intensity"},
                            merge_on="name", samples_in_columns=True):
    if "mz" not in mapping and "rt" not in mapping and "name" not in mapping:
        raise ValueError("Incorrect column mapping: provide column names for mz, and name")

    df_peaklist = pd.read_csv(fn_peaklist, header=0, sep=separator, dtype={"name": str}, float_precision="round_trip")
    df_matrix = pd.read_csv(fn_matrix, header=0, sep=separator, dtype={"name": str}, float_precision="round_trip")
    df_matrix.replace(0, np.nan, inplace=True)

    if not samples_in_columns:
        df_matrix = df_matrix.T

    if mapping["mz"] in df_peaklist.columns and mapping["name"] not in df_peaklist.columns and mapping["mz"] in df_matrix.columns:
        df_peaklist = read_peaklist(fn_peaklist, separator=separator)
        df_peaklist = df_peaklist[[mapping["name"], mapping["mz"], mapping["rt"], "intensity"]]

        df_matrix = df_matrix.rename(columns={"mz": 'name'})
        df_matrix["name"] = [str(x).replace(".", "_") for x in df_matrix["name"]]
    else:
        df_peaklist = df_peaklist[[mapping["name"], mapping["mz"], mapping["rt"]]]
        df_peaklist.columns = ["name", "mz", "rt"]

        df_matrix = df_matrix.rename(columns={mapping["name"]: 'name'})

    if mapping["intensity"] not in df_peaklist.columns:
        if median_intensity:
            df_peaklist["intensity"] = pd.Series(df_matrix.median(axis=1, skipna=True), index=df_matrix.index)
        else:
            df_peaklist["intensity"] = pd.Series(df_matrix.mean(axis=1, skipna=True), index=df_matrix.index)

    if len(df_peaklist[mapping["name"]].unique()) != len(df_peaklist[mapping["name"]]):
        raise ValueError("Peaklist: Values column '{}' are not unique".format(mapping["name"]))
    if len(df_matrix[mapping["name"]].unique()) != len(df_matrix[mapping["name"]]):
        raise ValueError("Matrix: Values column '{}' are not unique".format(mapping["name"]))

    return pd.merge(df_peaklist, df_matrix, how='left', left_on=merge_on, right_on=merge_on)


def read_peaklist(fn_peaklist, separator="\t",
                  mapping={"name": "name", "mz": "mz", "rt": "rt", "intensity": "intensity"}):

    df_peaklist = pd.read_csv(fn_peaklist, header=0, sep=separator, dtype={"name": str}, float_precision="round_trip")
    if mapping["mz"] not in df_peaklist.columns.values or mapping["intensity"] not in df_peaklist.columns.values:
        raise ValueError("Incorrect mapping of columns: {}".format(str(mapping)))

    if ("rt" in mapping and mapping["rt"] not in df_peaklist.columns.values) or "rt" not in mapping:
        if mapping["name"] not in df_peaklist.columns.values:
            df_peaklist = pd.read_csv(fn_peaklist, header=0, sep=separator, dtype={"mz": str})
            df_peaklist = df_peaklist[[mapping["mz"], mapping["intensity"]]]
            df_peaklist.columns = ["mz", "intensity"]
            df_peaklist.insert(0, "name", [str(x).replace(".","_") for x in df_peaklist[mapping["mz"]]])
            df_peaklist["mz"] = df_peaklist["mz"].astype(float)
            df_peaklist["intensity"] = df_peaklist["intensity"].astype(float)
        else:
            df_peaklist = df_peaklist[[mapping["name"], mapping["mz"], mapping["intensity"]]]
            df_peaklist.columns = ["name", "mz", "intensity"]
            df_peaklist["mz"] = df_peaklist["mz"].astype(float)
            df_peaklist["intensity"] = df_peaklist["intensity"].astype(float)
        df_peaklist.insert(2, "rt", 0.0)
    elif "rt" in mapping:
        if mapping["name"] in df_peaklist.columns.values:
            df_peaklist = df_peaklist[[mapping["name"], mapping["mz"], mapping["rt"], mapping["intensity"]]]
            df_peaklist.columns = ["name", "mz", "rt", "intensity"]
        else:
            df_peaklist = df_peaklist[[mapping["mz"], mapping["rt"], mapping["intensity"]]]
            df_peaklist.columns = ["mz", "rt", "intensity"]

            uids = df_peaklist["mz"].round().astype(int).astype(str).str.cat(df_peaklist["rt"].round().astype(int).astype(str), sep="T")
            ms = pd.Series(['M'] * len(uids))
            names = ms.str.cat(uids, sep='')

            for n in names.copy():
                idxs = names.index[names == n].tolist()
                if len(idxs) > 1:
                    for i, idx_t in enumerate(idxs):
                        names[idx_t] = names[idx_t] + "_" + str(i + 1)
            df_peaklist.insert(0, "name", names)
    else:
        df_peaklist = df_peaklist[[mapping["name"], mapping["mz"], mapping["rt"], mapping["intensity"]]]
        df_peaklist.columns = ["name", "mz", "rt", "intensity"]

    df_peaklist["mz"] = df_peaklist["mz"].astype(float)
    df_peaklist["rt"] = df_peaklist["rt"].astype(float)
    df_peaklist["intensity"] = df_peaklist["intensity"].astype(float)

    return df_peaklist
