#!/bin/bash

beamspy group-features \
--peaklist ../tests/test_data/peaklist_lcms_pos_theoretical.txt \
--intensity-matrix ../tests/test_data/dataMatrix_lcms_theoretical.txt \
--gml graph.gml \
--db results.sqlite \
--max-rt-diff 5.0 \
--method pearson \
--coeff-threshold 0.7 \
--pvalue-threshold 0.01

beamspy annotate-peak-patterns \
--peaklist ../tests/test_data/peaklist_lcms_pos_theoretical.txt \
--intensity-matrix ../tests/test_data/dataMatrix_lcms_theoretical.txt \
--gml graph.gml \
--db results.sqlite \
--adducts \
--adducts-library ../beamspy/data/adducts.txt \
--isotopes \
--isotopes-library ../beamspy/data/isotopes.txt \
--ion-mode pos \
--ppm 5.0

beamspy annotate-compounds \
--peaklist ../tests/test_data/peaklist_lcms_pos_theoretical.txt \
--intensity-matrix ../tests/test_data/dataMatrix_lcms_theoretical.txt \
--db results.sqlite \
--db-name hmdb_full_v4_0_20200909_v1 \
--adducts-library ../beamspy/data/adducts.txt \
--ion-mode pos \
--ppm 3.0

beamspy summary-results \
--peaklist ../tests/test_data/peaklist_lcms_pos_theoretical.txt \
--intensity-matrix ../tests/test_data/dataMatrix_lcms_theoretical.txt \
--db results.sqlite \
--output summary.txt \
--sep tab
