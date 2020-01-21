beamspy group-features^
 --peaklist "C:\beams\tests\test_data\variableMetadata.txt"^
 --intensity-matrix "C:\beams\tests\test_data\dataMatrix.txt"^
 --gml "C:\beams\tests\test_results\graph.gml"^
 --db "C:\beams\tests\test_results\results.sqlite"^
 --max-rt-diff 5.0^
 --method pearson^
 --coeff-threshold 0.7^
 --pvalue-threshold 0.01

beamspy annotate-peak-patterns^
 --peaklist "C:\beams\tests\test_data\variableMetadata.txt"^
 --intensity-matrix "C:\beams\tests\test_data\dataMatrix.txt"^
 --gml "C:\beams\tests\test_results\graph.gml"^
 --db "C:\beams\tests\test_results\results.sqlite"^
 --adducts^
 --adducts-library "C:\beams\beams\data\adducts.txt"^
 --isotopes^
 --isotopes-library "C:\beams\beams\data\isotopes.txt"^
 --ion-mode pos^
 --ppm 5.0

beamspy annotate-mf^
 --peaklist "C:\beams\tests\test_data\variableMetadata.txt"^
 --intensity-matrix "C:\beams\tests\test_data\dataMatrix.txt"^
 --db "C:\beams\tests\test_results\results.sqlite"^
 --adducts-library "C:\beams\beams\data\adducts.txt"^
 --ion-mode pos^
 --ppm 3.0^
 --max-mz 700.0

beamspy annotate-compounds^
 --peaklist "C:\beams\tests\test_data\variableMetadata.txt"^
 --intensity-matrix "C:\beams\tests\test_data\dataMatrix.txt"^
 --db "C:\beams\tests\test_results\results.sqlite"^
 --db-name hmdb_full_v4_0_v1^
 --adducts-library "C:\beams\beams\data\adducts.txt"^
 --ion-mode pos^
 --ppm 3.0

beamspy summary-results^
 --peaklist "C:\beams\tests\test_data\variableMetadata.txt"^
 --intensity-matrix "C:\beams\tests\test_data\dataMatrix.txt"^
 --db "C:\beams\tests\test_results\results.sqlite"^
 --output "C:\beams\tests\test_results\summary.txt"^
 --sep tab
