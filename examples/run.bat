beams group-features^
 --peaklist "C:\beams_v0.1.0\tests\test_data\variableMetadata.txt"^
 --intensity-matrix "C:\beams_v0.1.0\tests\test_data\dataMatrix.txt"^
 --gml "C:\beams_v0.1.0\tests\test_results\graph.gml"^
 --db "C:\beams_v0.1.0\tests\test_results\results.sqlite"^
 --max-rt-diff 5.0^
 --method pearson^
 --coeff-threshold 0.7^
 --pvalue-threshold 0.01

beams annotate-peak-patterns^
 --peaklist "C:\beams_v0.1.0\tests\test_data\variableMetadata.txt"^
 --intensity-matrix "C:\beams_v0.1.0\tests\test_data\dataMatrix.txt"^
 --gml "C:\beams_v0.1.0\tests\test_results\graph.gml"^
 --db "C:\beams_v0.1.0\tests\test_results\results.sqlite"^
 --adducts^
 --adducts-library "C:\beams_v0.1.0\beams\data\adducts.txt"^
 --isotopes^
 --isotopes-library "C:\beams_v0.1.0\beams\data\isotopes.txt"^
 --ion-mode pos^
 --ppm 5.0

python -m beams annotate-mf^
 --peaklist "C:\beams_v0.1.0\tests\test_data\variableMetadata.txt"^
 --intensity-matrix "C:\beams_v0.1.0\tests\test_data\dataMatrix.txt"^
 --db "C:\beams_v0.1.0\tests\test_results\results.sqlite"^
 --db-mf "C:\beams_v0.1.0\beams\data\db_mf.txt"^
 --adducts-library "C:\beams_v0.1.0\beams\data\adducts.txt"^
 --ion-mode pos^
 --ppm 3.0^
 --max-mz 700.0

python -m beams annotate-compounds^
 --peaklist "C:\beams_v0.1.0\tests\test_data\variableMetadata.txt"^
 --intensity-matrix "C:\beams_v0.1.0\tests\test_data\dataMatrix.txt"^
 --db "C:\beams_v0.1.0\tests\test_results\results.sqlite"^
 --db-compounds "C:\beams_v0.1.0\beams\data\BEAMS_DB.sqlite"^
 --adducts-library "C:\beams_v0.1.0\beams\data\adducts.txt"^
 --db-name HMDB^
 --ion-mode pos^
 --ppm 3.0

python -m beams summary-results^
 --peaklist "C:\beams_v0.1.0\tests\test_data\variableMetadata.txt"^
 --intensity-matrix "C:\beams_v0.1.0\tests\test_data\dataMatrix.txt"^
 --db "C:\beams_v0.1.0\tests\test_results\results.sqlite"^
 --output "C:\beams_v0.1.0\tests\test_results\summary.txt"^
 --sep tab
