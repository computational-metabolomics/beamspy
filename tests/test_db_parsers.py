#!/usr/bin/env python
#  -*- coding: utf-8 -*-

import os
import unittest
from collections import OrderedDict
from beams.db_parsers import *
from tests.utils import to_test_data
import numpy as np


class DbParsersTestCase(unittest.TestCase):

    def setUp(self):
        self.path, f = os.path.split(os.path.dirname(os.path.abspath(__file__)))

    def test_parse_biocyc(self):
        records = list(parse_biocyc(to_test_data("biocyc_record.txt")))
        self.assertEqual(records[0]['UNIQUE-ID'], 'PANTOTHENATE')

    def test_parse_sdf(self):
        records = list(parse_sdf(to_test_data("sdf_record.sdf")))
        self.assertEqual(records[0]['ChEBI ID'], 'CHEBI:90')
        self.assertEqual(records[0]['ChEBI Name'], '(-)-epicatechin')

    def test_parse_xml(self):
        records = list(parse_xml(to_test_data("hmdb_record.xml")))
        self.assertEqual(records[0]['accession'], 'HMDB0000001')

    def test_parse_kegg_compound(self):
        records = list(parse_kegg_compound(to_test_data("kegg_record.txt")))
        self.assertEqual(records[0]['ENTRY'], 'C00001')

    def test_parse_delimited(self):
        records = list(parse_delimited(to_test_data("tab_delimited_record.txt"), "\t"))
        self.assertEqual(records[0]["Compound_id"], "2-METHYL-6-SOLANYL-14-BENZOQUINONE")


if __name__ == '__main__':
    unittest.main()