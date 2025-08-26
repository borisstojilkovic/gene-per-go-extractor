# gene-per-go extractor

extract genes per GO term from RNA-seq result files, merge with annotations, and export:
- per-GO tables (rna-seq + annotations)
- a grouped table where each column is a GO term and rows are gene IDs

## features
- supports multiple input formats: `.tab` (tsv) or excel (`.xlsx/.xls`)
- species presets: tomato (`S`), tomato with RKN dataset (`SW`), arabidopsis (`A`)
- automatic output subfolders per input file
- de-duplicates by stable gene ID
- grouped matrix output for easy copy/paste or downstream analysis

## requirements
- python ≥ 3.9
- pandas
- openpyxl (for `.xlsx/.xls` files)

```bash
pip install pandas openpyxl
```

## folder layout
```
repo/
├─ extract_genes_per_GO_term.py
├─ input/                               # put your RNA-seq result files here
├─ output/                              # auto-created per input file
├─ annotations/
│  ├─ annotation_tom.xlsx
│  ├─ annotation_tom_with_RKN_Mi-Tomato.xlsx
│  └─ annotation_arab.xlsx
├─ GO term accession_mart_exportSL_3.0_with descr.tab
├─ GO term accession_mart_exportArabidopsis.tab
└─ Go_termnIDs_and_file_names.xlsx
```

## inputs
- **rna-seq files** (place in `input/`)
  - DESeq2 outputs as tab-separated (`.tab`) or excel (`.xlsx/.xls`)
  - must contain a column named **`GeneID`**
    - for tomato: `GeneID` like `SolycXXgXXXXX.X` (script derives `locus` by splitting at the first dot)
    - for arabidopsis: `GeneID` used directly as `locus`
- **go term list**: `Go_termnIDs_and_file_names.xlsx`
  - must include a column named **`names`** with accessions, e.g. `GO:0009535`
- **go legends**
  - tomato datasets: `GO term accession_mart_exportSL_3.0_with descr.tab`
  - arabidopsis: `GO term accession_mart_exportArabidopsis.tab`
  - must contain columns: `accession`, `Gene stable ID`, `GO term name`, `GO term definition`
- **annotation files** (in `annotations/`)
  - `annotation_tom.xlsx` (S)
  - `annotation_tom_with_RKN_Mi-Tomato.xlsx` (SW)
  - `annotation_arab.xlsx` (A)
  - each must contain a **`locus`** column and any additional annotation columns you want merged

## outputs
for each input file `<basename>`:
- `output/<basename>/<GO_XXXX term_name>.tab`  
  subset of the GO legend for that term
- `output/<basename>/<basename><GO_XXXX term_name>.tab`  
  merged table: rna-seq + annotations for genes in that GO term (`,` used as decimal in output)
- `output/<basename>/<basename>1_grouped.tab`  
  grouped matrix: each column is a GO term; rows are `GeneID` (padded to 10,000 rows)

## usage
1) place your rna-seq files in `input/`  
2) ensure `Go_termnIDs_and_file_names.xlsx` has a `names` column of GO accessions  
3) check the available go legend and annotation files for your species  
4) run the script:
```bash
python extract_genes_per_GO_term.py
```
5) choose species at the prompt:
- `S` → tomato
- `SW` → tomato with RKN (Mi–Tomato)
- `A` → arabidopsis

outputs will be written under `output/<file_basename>/`.

## notes & assumptions
- the script expects `GeneID` in your rna-seq files.
- for tomato (`S`/`SW`), `locus` is derived as text before the first dot in `GeneID`.
- duplicates are removed based on `Gene stable ID`.
- grouped output pads columns to 10,000 rows with empty strings for consistent table shape.
- output files use `,` as the decimal separator to match european spreadsheet defaults.

## common issues & troubleshooting
- **file not found (go legend / annotations)**  
  verify filenames and paths match those in this readme (case-sensitive on some systems).
- **missing `GeneID` column**  
  rename your gene identifier column to `GeneID` before running.
- **no entries found for a GO term**  
  the accession may not exist in the provided go legend; confirm the accession and legend source match your species.
- **encoding issues**  
  ensure files are UTF-8 (default for the script).
- **very large outputs**  
  consider filtering your rna-seq inputs (e.g., significant genes only) before running.

## example
- `input/DESeq2_example.tab`
- `Go_termnIDs_and_file_names.xlsx` with:
  ```
  names
  GO:0009535
  GO:0006355
  ```
run the script, select `S`, and inspect `output/DESeq2_example/` for per-GO tables and the grouped matrix.

## citation / acknowledgments
if you use this script in a publication, please cite 


