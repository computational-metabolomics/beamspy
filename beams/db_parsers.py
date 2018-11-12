#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
from collections import OrderedDict
import io
import xml.etree.ElementTree as etree
import csv
from Bio.KEGG import Compound
import re


def parse_delimited(source, delimiter):
    with open(source, 'r') as inp:
        reader = csv.DictReader(inp, delimiter=delimiter)
        for row in reader:
            yield row


def parse_kegg_compound(source, sdf=False):

    with open(source, "r") as inp:
        ref_line = inp.readline()

    with open(source, "r") as inp:
        if "ENTRY" in ref_line and "Reaction" not in ref_line:
            for record in Compound.parse(inp):
                record_out = OrderedDict()
                if "C" in record.entry or "D" in record.entry:
                    for attribute in dir(record):
                        if "_" not in attribute:
                            record_out[attribute.upper()] = ""
                            record_out[attribute.upper()] = getattr(record, attribute.lower())

                    if sdf:
                        record_out["SDF"] = REST.GetMol(record_out["ENTRY"])

                    yield record_out

        elif "ENTRY" in ref_line and "Reaction" in ref_line:
            for record in Reaction.parse(inp):
                record_out = OrderedDict()
                for attribute in dir(record):
                    if "_" not in attribute:
                        record_out[attribute.upper()] = getattr(record, attribute.lower())
                yield record_out


def parse_xml(source, encoding="utf8"):

    with io.open(source, "r", encoding=encoding) as inp:
        record_out = OrderedDict()

        xmldec = inp.readline()
        xmldec2 = inp.readline()

        xml_record = ""
        path = []

        for line in inp:
            xml_record += line
            if line == "</metabolite>\n" or line == "</drug>\n":

                if sys.version_info[0] == 3:
                    inp = io.StringIO(xml_record)
                else:
                    inp = io.BytesIO(xml_record.encode('utf-8').strip())

                for event, elem in etree.iterparse(inp, events=("start", "end")):
                    if event == 'end':
                        path.pop()

                    if event == 'start':
                        path.append(elem.tag)
                        if elem.text != None:
                            if elem.text.replace(" ", "") != "\n":

                                path_elem = ".".join(map(str, path[1:]))
                                if path_elem in record_out:
                                    if type(record_out[path_elem]) != list:
                                        record_out[path_elem] = [record_out[path_elem]]
                                    record_out[path_elem].append(elem.text)
                                else:
                                    record_out[path_elem] = elem.text

                xml_record = ""
                yield record_out
                record_out = OrderedDict()


def parse_sdf(source):

    with open(source, "r") as inp:
        record_out = OrderedDict()
        c = 0
        temp = ""
        for line in inp:
            line = line.replace("'", "").replace('"', "")
            if "$$$$" in line:
                temp = temp.split("> <")
                for attribute_value in temp[1:]:
                    attribute_value = attribute_value.split(">\n")
                    if len(attribute_value) == 1:
                        record_out[list(record_out.keys())[-1]] += attribute_value[0].strip("\n")
                    else:
                        record_out[attribute_value[0]] = attribute_value[1].strip("\n")
                c += 1
                record_out["SDF"] = temp[0]
                yield record_out
                temp = ""
                record_out = OrderedDict()
            else:
                temp += line


def parse_biocyc(source):

    with open(source, "r") as inp:
        record_out = OrderedDict()
        temp = ""

        for line in inp:
            line = line.replace("'", "").replace('"', "")
            if "//" in line:

                temp_attribute = ""
                extra_line = 0
                extra_line_added = 1
                for line_temp in temp.split("\n")[0:-1]:
                    extra_line += 1  # did it pass the attribute?
                    attribute_value = line_temp.split(" - ", 1)

                    if attribute_value[0][0] != "/":

                        if attribute_value[0] not in record_out:

                            try:
                                record_out[attribute_value[0]] = float(attribute_value[1])
                            except:
                                record_out[attribute_value[0]] = attribute_value[1].replace('"', "'")
                            temp_attribute = attribute_value[0]
                            extra_line = 1

                        elif attribute_value[0] in record_out:
                            if type(record_out[temp_attribute]) != list:
                                record_out[temp_attribute] = [record_out[temp_attribute]]
                            record_out[temp_attribute].append(attribute_value[1].replace('"', "'"))
                            temp_attribute = attribute_value[0]
                            extra_line = 1

                    elif attribute_value[0][0] == "/" and len(attribute_value[0]) > 1:
                        if temp_attribute in record_out and (extra_line == 2 or extra_line_added >= 1):
                            if type(record_out[temp_attribute]) != list:
                                record_out[temp_attribute] = [record_out[temp_attribute]]
                            index_to_add = len(record_out[temp_attribute]) - 1
                            record_out[temp_attribute][index_to_add] = record_out[temp_attribute][index_to_add] + attribute_value[0].replace("/", "", 1).replace('"', "'")
                            extra_line_added += 1

                # PRINT FORMULA BIOCYC IN CORRECT FORMAT ######
                if "CHEMICAL-FORMULA" in record_out:
                    formula = ""
                    for atom in record_out["CHEMICAL-FORMULA"]:
                        if " 1)" in atom:
                            atom = atom.replace(" 1)", "")
                            formula += atom.replace("(", "")
                        else:
                            formula += atom.replace(" ", "")[1:-1]
                    record_out["CHEMICAL-FORMULA"] = formula
                # PRINT FORMULA BIOCYC IN CORRECT FORMAT ######

                yield record_out
                temp = ""
                record_out = OrderedDict()

            elif line[0] != "#":
                temp += line


def parse_nist_database(fn, skip_lines=10):

    """
    :param fn: text file (NISTs Linearized ASCII Output)
    :param skip_lines: the number of lines of the data file to skip before beginning to read data.
    :return: Ordered dictionary containing the parsed records
    """

    with open(fn, "r") as inp:
        for i in range(skip_lines):
            inp.readline()
        for e in inp.read().split("\n\n"):
            record = OrderedDict()
            for line in e.strip().split("\n"):
                kv = line.split(" =")
                if kv[0] == "Relative Atomic Mass":
                    record[kv[0]] = re.findall(r'\d+(?:\.\d+)?', kv[1])
                    record[kv[0]][0] = float(record[kv[0]][0])
                    record[kv[0]][1] = int(record[kv[0]][1])

                elif kv[0] == "Isotopic Composition":
                    matches = re.findall(r'\d+(?:\.\d+)?', kv[1])
                    if len(matches) > 0:
                        record[kv[0]] = matches
                        if len(matches) > 1:
                            record[kv[0]][0] = float(record[kv[0]][0])
                            record[kv[0]][1] = int(record[kv[0]][1])
                        else:
                            record[kv[0]] = [float(record[kv[0]][0]), None]
                    else:
                        record[kv[0]] = [0.0, None]
                elif kv[0] == "Atomic Number" or kv[0] == "Mass Number":
                    record[kv[0]] = int(kv[1])
                elif kv[0] == "Standard Atomic Weight":
                    matches = re.findall(r'\d+(?:\.\d+)?', kv[1])
                    record[kv[0]] = matches
                else:
                    record[kv[0]] = kv[1].strip()
            yield record
