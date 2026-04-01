# CellBLASTer
A universal plant scRNA-seq annotation tool inspired by cellular BLAST strategies.
CellBlaster is a cross-species cell type identification and annotation tool designed specifically for plant single-cell transcriptome (scRNA-seq). Through cross-species orthogroup (OG) mapping, symbolic percentage encoding, and multi-round correction algorithms, it accurately maps the query dataset to the reference database, achieving high-confidence automatic cell type annotation.

## Installation
We recommend using conda to manage the installation of all dependencies. To do this, simply run:
Then, download this repo and install it.
```
git clone [repo_path]
cd CellBlaster-main
pip install .
```
The total installation time is around 1-2 mintunes. If error occuors, please upgrade pip and try again.




## Other Information
**Dicot Datasets Information**
| Species | Abbreviation | Count | SRA | Cell number |
| :--- | :--- | :---: | :--- | :---: |
| *Arabidopsis thaliana* | *A. thaliana* | 14 | SRP267870 | 188776 |
| | | | SRP235541 | 23519 |
| | | | SRP171040 | 33386 |
| | | | SRP182008 | 10353 |
| | | | SRP166333 | 16055 |
| | | | SRP285817 | 16572 |
| | | | SRP273996 | 9072 |
| | | | SRP330542 | 20193 |
| | | | SRP173393 | 10201 |
| | | | SRP169576 | 25338 |
| | | | SRP148288 | 14141 |
| | | | SRP332285 | 9837 |
| | | | SRP285040 | 973 |
| | | | SRP394711 | 310839 |
| *Glycine max* | *G. max* | 2 | CRA007122 | 14808 |
| | | | CRA008947 | 8021 |
| *Manihot esculenta* | *M. esculenta* | 1 | SRP406470 | 23765 |
| *Medicago truncatula* | *M. truncatula* | 1 | SRP390780 | 24978 |
| *Lotus japonicus* | *L. japonicus* | 1 | SRP376527 | 19921 |

**Monocot Datasets Information**
| Species | Abbreviation | Count | SRA | Cell number |
| :--- | :--- | :---: | :--- | :---: |
| *Oryza sativa* | *O. sativa* | 2 | SRP309176 | 33720 |
| | | | SRP250946 | 28918 |
| *Triticum aestivum* | *T. aestivum* | 3 | CRA008788 | 3756 |
| | | | GSE270342 | 3846 |
| | | | SRP543892 | 35907 |
| *Zea mays* | *Z. mays* | 3 | SRP145013 | 19669 |
| | | | GSE225118 | 7303 |
| | | | SRP335180 | 4251 |
| *Sorghum bicolor* | *S. bicolor* | 1 | SRP422815_1 | 12839 |
| *Setaria viridis* | *S. viridis* | 1 | SRP422815_2 | 25609 |
| *Phyllostachys edulis* | *P. edulis* | 1 | GSE229126 | 3573 |

**Isoform Datasets Information**
| Type | Parameter name | File name |
| :--- | :--- | :--- |
| Dicot | *A. thaliana* | *A. thaliana_isoform.fa* |
| Dicot | *G. max* | *G. max_isoform.fa* |
| Dicot | *L. japonicus* | *L. japonicus_isoform.fa* |
| Dicot | *M. esculenta* | *M. esculenta_isoform.fa* |
| Dicot | *M. truncatula* | *M. truncatula_isoform.fa* |
| Monocot | *O. sativa* | *O. sativa_isoform.fa* |
| Monocot | *P. edulis* | *P. edulis_isoform.fa* |
| Monocot | *S. bicolor* | *S. bicolor_isoform.fa* |
| Monocot | *S. viridis* | *S. viridis_isoform.fa* |
| Monocot | *T. aestivum* | *T. aestivum_isoform.A.fa* |
| Monocot | | *T. aestivum_isoform.B.fa* |
| Monocot | | *T. aestivum_isoform.D.fa* |
| Monocot | *Z. mays* | *Z. mays_isoform.fa* |
