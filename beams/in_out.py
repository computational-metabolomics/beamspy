from pandas import read_csv
import libraries
import collections
import pandas as pd


def mf_dict_to_str(d):
    mf = ""
    for atom in ["C", "H", "N", "O", "P", "S"]:
        if d[atom] > 0:
            mf += atom + str(d[atom])
        elif d[atom] == 1:
            mf += str(d[atom])
    return mf


def read_atoms(filename, separator="\t"):
    df = read_csv(filename, sep=separator)
    atoms = libraries.Atoms()
    for index, row in df.iterrows():
        atoms.add(row["name"],  row["mass"])
    return atoms


def read_adducts(filename, ion_mode, separator="\t"):
    df = read_csv(filename, sep=separator)
    adducts = libraries.Adducts()
    adducts.delete("*")
    for index, row in df.iterrows():
        if "ion_mode" not in row:
            adducts.add(row["label"], row["mass"])
        elif (row["ion_mode"] == "pos" or row["ion_mode"] == "both") and ion_mode == "pos":
            adducts.add(row["label"], row["mass"])
        elif (row["ion_mode"] == "neg" or row["ion_mode"] == "both") and ion_mode == "neg":
            adducts.add(row["label"], row["mass"])
    return adducts


def read_isotopes(filename, ion_mode, separator="\t"):
    df = read_csv(filename, sep=separator)
    isotopes = libraries.Isotopes()
    isotopes.delete("*")
    for index, row in df.iterrows():
        if "ion_mode" not in row:
            isotopes.add(row["label_x"], row["label_y"], row["abundance_x"], row["abundance_y"], row["mass_difference"])
        elif (row["ion_mode"] == "pos" or row["ion_mode"] == "both") and ion_mode == "pos":
            isotopes.add(row["label_x"], row["label_y"], row["abundance_x"], row["abundance_y"], row["mass_difference"])
        elif (row["ion_mode"] == "neg" or row["ion_mode"] == "both") and ion_mode == "neg":
            isotopes.add(row["label_x"], row["label_y"], row["abundance_x"], row["abundance_y"], row["mass_difference"])
    return isotopes


def read_salts(filename, ion_mode, separator="\t"):
    df = read_csv(filename, sep=separator)
    return df


def read_molecular_formulae(filename, separator="\t"):

    df = read_csv(filename, sep=separator)
    atoms = libraries.Atoms()
    ions = libraries.Adducts()
    isotopes = libraries.Isotopes()
    records = []
    for index, row in df.iterrows():
        f = libraries.Formula(str(row.mf), atoms, ions, isotopes)
        record = collections.OrderedDict()
        atom_counts = f.count()
        if atom_counts:
            record.update(atom_counts[0])
            record["ExactMass"] = f.calc_mass()[0]
            record.update(f.HNOPS())
            record.update(f.lewis_senior())
            records.append(record)
        else:
            Warning("{} Skipped".format(row))

    return records


def read_compounds(filename, separator="\t"):
    df = read_csv(filename, sep=separator)
    atoms = libraries.Atoms()
    adducts = libraries.Adducts()
    isotopes = libraries.Isotopes()
    records = []
    for index, row in df.iterrows():
        f = libraries.Formula(str(row.mf), atoms, adducts, isotopes)
        record = collections.OrderedDict()
        atom_counts = f.count()
        if atom_counts:
            record.update(atom_counts[0])
            record["exact_mass"] = f.calc_mass()[0]
            record["compound_id"] = row.compound_id
            record["compound_name"] = row.name
            record["molecular_formula"] = mf_dict_to_str(record)
            records.append(record)
        else:
            Warning("{} Skipped".format(row))

    return records


def read_multiple_charged_ions(filename, ion_mode, separator="\t"):
    df = read_csv(filename, sep=separator)
    multiple_charges = libraries.MultipleChargedIons()
    multiple_charges.delete("*")
    for index, row in df.iterrows():
        if "ion_mode" not in row:
            multiple_charges.add(row["label"], row["mass"], row["charge"])
        elif (row["ion_mode"] == "pos" or row["ion_mode"] == "both") and ion_mode == "pos":
            multiple_charges.add(row["label"], row["mass"], row["charge"])
        elif (row["ion_mode"] == "neg" or row["ion_mode"] == "both") and ion_mode == "neg":
            multiple_charges.add(row["label"], row["mass"], row["charge"])
    return multiple_charges


def read_mass_differences(filename, ion_mode, separator="\t"):
    df = read_csv(filename, sep=separator)
    mass_differences = libraries.MassDifferences()
    for index, row in df.iterrows():
        if "charge_x" in row:
            charge_x = row["charge_x"]
            charge_y = row["charge_y"]
        else:
            charge_x = 1
            charge_y = 2
        if "ion_mode" not in row:
            mass_differences.add(row["label_x"], row["label_y"], row["mass_difference"], charge_x, charge_y)
        elif (row["ion_mode"] == "pos" or row["ion_mode"] == "both") and ion_mode == "pos":
            mass_differences.add(row["label_x"], row["label_y"], row["mass_difference"], charge_x, charge_y)
        elif (row["ion_mode"] == "neg" or row["ion_mode"] == "both") and ion_mode == "neg":
            mass_differences.add(row["label_x"], row["label_y"], row["mass_difference"], charge_x, charge_y)
    return mass_differences


def read_xset_matrix(fn_matrix, first_sample, separator="\t", mapping={"mz": "mz", "rt": "rt", "name": "name"}, samples_in_columns=True):
    if "mz" not in mapping and "rt" not in mapping and "name" not in mapping:
        raise ValueError("Incorrect column mapping: provide column names for mz, and name")

    df = pd.read_csv(fn_matrix, header=0, sep=separator)

    if not samples_in_columns:
        df = df.T

    df_peaklist = df[[mapping["name"], mapping["mz"], mapping["rt"]]]
    df_peaklist["name"] = df_peaklist["name"].astype(str)
    df_matrix = df.iloc[:, df.columns.get_loc(first_sample):]
    df_peaklist = df_peaklist.assign(intensity=pd.Series(df_matrix.median(axis=1, skipna=True).values))
    return pd.concat([df_peaklist, df_matrix], axis=1)


def combine_peaklist_matrix(fn_peaklist, fn_matrix, separator="\t", mapping={"name": "name", "mz": "mz", "rt": "rt"}, merge_on="name", samples_in_columns=True):
    if "mz" not in mapping and "rt" not in mapping and "name" not in mapping:
        raise ValueError("Incorrect column mapping: provide column names for mz, and name")

    df_peaklist = pd.read_csv(fn_peaklist, header=0, sep=separator)
    df_matrix = pd.read_csv(fn_matrix, header=0, sep=separator)

    if not samples_in_columns:
        df_matrix = df_matrix.T

    df_peaklist = df_peaklist[[mapping["name"], mapping["mz"], mapping["rt"]]]
    df_peaklist.columns = ["name", "mz", "rt"]
    df_peaklist["name"] = df_peaklist["name"].astype(str)

    df_matrix = df_matrix.rename(columns={mapping["name"]: 'name'})
    df_matrix["name"] = df_matrix["name"].astype(str)

    df_peaklist["intensity"] = pd.Series(df_matrix.median(axis=1, skipna=True), index=df_matrix.index)
    # pd.merge(df_peaklist, df_matrix, how='left', left_on=merge_on, right_on=merge_on).to_csv("test_out_df.txt", sep="\t")
    return pd.merge(df_peaklist, df_matrix, how='left', left_on=merge_on, right_on=merge_on)


def read_peaklist(fn_peaklist, separator="\t", mapping={"name": "name", "mz": "mz", "rt": "rt", "intensity": "intensity"}):

    df_peaklist = pd.read_csv(fn_peaklist, header=0, sep=separator)
    if "rt" in mapping and mapping["rt"] is not None:
        if mapping["name"] in df_peaklist.columns.values:
            df_peaklist = df_peaklist[[mapping["name"], mapping["mz"], mapping["intensity"]]]
            df_peaklist.columns = ["name", "mz", "intensity"]
            df_peaklist["name"] = df_peaklist["name"].astype(str)
        else:
            df_peaklist = df_peaklist[[mapping["mz"], mapping["intensity"]]]
            df_peaklist.columns = ["mz", "intensity"]
            df_peaklist.insert(0, "name", df_peaklist[mapping["mz"]].astype(str))
        df_peaklist.insert(2, "rt", 0.0)
    else:
        df_peaklist = df_peaklist[[mapping["name"], mapping["mz"], mapping["rt"], mapping["intensity"]]]
        df_peaklist["name"] = df_peaklist["name"].astype(str)
    return df_peaklist
