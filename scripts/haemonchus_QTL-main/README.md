# haemonchus_QTL
A repository containing files and a script to compare avermectin QTL between *C. elegans* and *H. contortus*.

The `print_interval_table.py` will parse the *C. elegans* and *H. contortus* GFF3 file (downloaded from WormBase and WormBase ParaSite; isoform filtered with AGAT 0.4.0), the single-copy orthogroups file (generated from OrthoFinder 2.4.0), and the locations of the respective QTL. The script will output (to STDOUT) a table (TSV) of names and coordinates of each gene in the *C. elegans* QTL, along the location of the *H. contortus* ortholog and whether it falls within the *H. contortus* QTL. 

The script should be run as follows: 

```
python print_interval_table.py c_elegans.PRJNA13758.WS273.annotations.gff3_longest_isoforms CELE_intervals.txt haemonchus_contortus.PRJEB506.WBPS15.annotations.gff3_longest_isoforms HCON_intervals.txt Orthogroups.singlecopy.txt >CELE_abamectin_QTL_gene_ortholog_coordinates.20112020.tsv
```
