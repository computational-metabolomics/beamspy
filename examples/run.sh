#!/usr/bin/env bash

python -m beams group-features \
--peaklist tests/test_data/variableMetadata.txt \
--intensity-matrix tests/test_data/dataMatrix.txt \
--gml tests/test_results/graph.gml \
--db tests/test_results/results.sqlite \
--max-rt-diff 5.0 \
--method pearson \
--coeff-threshold 0.7 \
--pvalue-threshold 0.01

python -m beams annotate-peak-patterns \
--peaklist tests/test_data/variableMetadata.txt \
--intensity-matrix tests/test_data/dataMatrix.txt \
--gml tests/test_results/graph.gml \
--db tests/test_results/results.sqlite \
--adducts \
--adducts-library beams/data/adducts.txt \
--isotopes \
--isotopes-library beams/data/isotopes.txt \
--ion-mode pos \
--ppm 5.0

python -m beams annotate-mf \
--peaklist tests/test_data/variableMetadata.txt \
--intensity-matrix tests/test_data/dataMatrix.txt \
--db tests/test_results/results.sqlite \
--mf-db beams/data/db_mf.txt \
--adducts-library beams/data/adducts.txt \
--ion-mode pos \
--ppm 3.0

python -m beams annotate-compounds \
--peaklist tests/test_data/variableMetadata.txt \
--intensity-matrix tests/test_data/dataMatrix.txt \
--db tests/test_results/results.sqlite \
--db-compounds beams/data/db_compounds.txt \
--adducts-library beams/data/adducts.txt \
--ion-mode pos \
--ppm 3.0

python -m beams annotate-compounds \
--peaklist tests/test_data/variableMetadata.txt \
--intensity-matrix tests/test_data/dataMatrix.txt \
--db tests/test_results/results.sqlite \
--db-compounds beams/data/BEAMS_DB.sqlite \
--adducts-library beams/data/adducts.txt \
--db-name HMDB \
--ion-mode pos \
--ppm 3.0

python -m beams summary-results \
--peaklist tests/test_data/variableMetadata.txt \
--intensity-matrix tests/test_data/dataMatrix.txt \
--db tests/test_results/results.sqlite \
--output tests/test_results/summary.txt \
--sep tab
