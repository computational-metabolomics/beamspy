#!/bin/bash

beamspy group-features \
--peaklist ../tests/test_data/variableMetadata.txt \
--intensity-matrix ../tests/test_data/dataMatrix.txt \
--gml ../tests/test_results/graph.gml \
--db ../tests/test_results/results.sqlite \
--max-rt-diff 5.0 \
--method pearson \
--coeff-threshold 0.7 \
--pvalue-threshold 0.01

beamspy annotate-peak-patterns \
--peaklist ../tests/test_data/variableMetadata.txt \
--intensity-matrix ../tests/test_data/dataMatrix.txt \
--gml ../tests/test_results/graph.gml \
--db ../tests/test_results/results.sqlite \
--adducts \
--adducts-library ../beamspy/data/adducts.txt \
--isotopes \
--isotopes-library ../beamspy/data/isotopes.txt \
--ion-mode pos \
--ppm 5.0

beamspy annotate-mf \
--peaklist ../tests/test_data/variableMetadata.txt \
--intensity-matrix ../tests/test_data/dataMatrix.txt \
--db ../tests/test_results/results.sqlite \
--adducts-library ../beamspy/data/adducts.txt \
--ion-mode pos \
--ppm 3.0

beamspy annotate-compounds \
--peaklist ../tests/test_data/variableMetadata.txt \
--intensity-matrix ../tests/test_data/dataMatrix.txt \
--db ../tests/test_results/results.sqlite \
--db-name hmdb_full_v4_0_v1 \
--adducts-library ../beamspy/data/adducts.txt \
--ion-mode pos \
--ppm 3.0

beamspy summary-results \
--peaklist ../tests/test_data/variableMetadata.txt \
--intensity-matrix ../tests/test_data/dataMatrix.txt \
--db ../tests/test_results/results.sqlite \
--output ../tests/test_results/summary.txt \
--sep tab
