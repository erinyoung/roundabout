# roundabout

The fundamental goal of this repo was to create something quick and reproducible to find large regions of similarity between highly related plasmids - like those in an outbreak.

roundabout uses blast to find regions of similarity (as opposed to groups of genes that would be required for synteny), and then assigns those regions a color. Those regions and colors are then visualized via circos.

## USAGE

```
nextflow run erinyoung/roundabout --fastas <directory to fastas>
```

There are not many editable parameters with this workflow, but the minimum length can be adjusted with `params.length` and the directory with the results can be adjusted with `params.outdir`.

This workflow uses singularity to use containers that run files through
dependencies:
- [blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi) : find similarities between fastas
- [circos](http://circos.ca/) : visualizing the end product
- [bedtools](https://bedtools.readthedocs.io/en/latest/) : combining regions of interest
- awk/sed/bash/python3 : a lot of file manipulation
- [samtools](http://www.htslib.org/) : to find genome size
- [prokka](http://www.htslib.org/) : to quickly annotate genes
- [amrfinder](https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/AMRFinder/) : for amr gene locations
- [plasmidfinder](https://bitbucket.org/genomicepidemiology/plasmidfinder/) : for plasmid identification

This was created for personal use with specific projects in mind, as opposed to general use. As such, other users may notice that the use case is highly specific and informal. Put in an issue if this is something that interests you and you need an additional feature.

<details>
   <summary>Final File Tree</summary>
```
roundabout/
├── bedtools
│   ├── example1.AC.bed
│   ├── example1.GC.bed
│   ├── example1.nuc.bed
│   ├── example1.skew.bed
│   ├── example1.windows.bed
│   ├── example2.AC.bed
│   ├── example2.GC.bed
│   ├── example2.nuc.bed
│   ├── example2.skew.bed
│   ├── example2.windows.bed
│   ├── example3.AC.bed
│   ├── example3.GC.bed
│   ├── example3.nuc.bed
│   ├── example3.skew.bed
│   ├── example3.windows.bed
│   ├── test_test1.AC.bed
│   ├── test_test1.GC.bed
│   ├── test_test1.nuc.bed
│   ├── test_test1.skew.bed
│   ├── test_test1.windows.bed
│   ├── test_test2.AC.bed
│   ├── test_test2.GC.bed
│   ├── test_test2.nuc.bed
│   ├── test_test2.skew.bed
│   └── test_test2.windows.bed
├── blastn
│   ├── example1.blast_hits.txt
│   ├── example1.divided.bed
│   ├── example1.sorted_divided.bed
│   ├── example1.starts_ends.txt
│   ├── example2.blast_hits.txt
│   ├── example2.divided.bed
│   ├── example2.sorted_divided.bed
│   ├── example2.starts_ends.txt
│   ├── example3.blast_hits.txt
│   ├── example3.divided.bed
│   ├── example3.sorted_divided.bed
│   ├── example3.starts_ends.txt
│   ├── test_test1.blast_hits.txt
│   ├── test_test1.divided.bed
│   ├── test_test1.sorted_divided.bed
│   ├── test_test1.starts_ends.txt
│   ├── test_test2.blast_hits.txt
│   ├── test_test2.divided.bed
│   ├── test_test2.sorted_divided.bed
│   └── test_test2.starts_ends.txt
├── circos
│   ├── example1
│   │   ├── data
│   │   │   ├── AC.bed
│   │   │   ├── amrfinder_text.txt
│   │   │   ├── amrfinder.txt
│   │   │   ├── GC.bed
│   │   │   ├── highlights_text.txt
│   │   │   ├── highlights.txt
│   │   │   ├── karyotype.txt
│   │   │   └── skew.bed
│   │   ├── example1_roundabout.png
│   │   └── example1_roundabout.svg
│   ├── example2
│   │   ├── data
│   │   │   ├── AC.bed
│   │   │   ├── amrfinder_text.txt
│   │   │   ├── amrfinder.txt
│   │   │   ├── GC.bed
│   │   │   ├── highlights_text.txt
│   │   │   ├── highlights.txt
│   │   │   ├── karyotype.txt
│   │   │   └── skew.bed
│   │   ├── example2_roundabout.png
│   │   └── example2_roundabout.svg
│   ├── example3
│   │   ├── data
│   │   │   ├── AC.bed
│   │   │   ├── amrfinder_text.txt
│   │   │   ├── amrfinder.txt
│   │   │   ├── GC.bed
│   │   │   ├── highlights_text.txt
│   │   │   ├── highlights.txt
│   │   │   ├── karyotype.txt
│   │   │   └── skew.bed
│   │   ├── example3_roundabout.png
│   │   └── example3_roundabout.svg
│   ├── test_test1
│   │   ├── data
│   │   │   ├── AC.bed
│   │   │   ├── amrfinder_text.txt
│   │   │   ├── amrfinder.txt
│   │   │   ├── GC.bed
│   │   │   ├── highlights_text.txt
│   │   │   ├── highlights.txt
│   │   │   ├── karyotype.txt
│   │   │   └── skew.bed
│   │   ├── test_test1_roundabout.png
│   │   └── test_test1_roundabout.svg
│   └── test_test2
│       ├── data
│       │   ├── AC.bed
│       │   ├── amrfinder_text.txt
│       │   ├── amrfinder.txt
│       │   ├── GC.bed
│       │   ├── highlights_text.txt
│       │   ├── highlights.txt
│       │   ├── karyotype.txt
│       │   └── skew.bed
│       ├── test_test2_roundabout.png
│       └── test_test2_roundabout.svg
├── highlight
│   ├── highlights_text.txt
│   └── highlights.txt
├── karyotype
│   ├── example1_karyotype_genes.bed
│   ├── example1_karyotype.txt
│   ├── example2_karyotype_genes.bed
│   ├── example2_karyotype.txt
│   ├── example3_karyotype_genes.bed
│   ├── example3_karyotype.txt
│   ├── test_test1_karyotype_genes.bed
│   ├── test_test1_karyotype.txt
│   ├── test_test2_karyotype_genes.bed
│   └── test_test2_karyotype.txt
├── ncbi-AMRFinderplus
│   ├── example1_amrfinder.bed
│   ├── example1_amrfinder_text.bed
│   ├── example1_amrfinder.txt
│   ├── example2_amrfinder.bed
│   ├── example2_amrfinder_text.bed
│   ├── example2_amrfinder.txt
│   ├── example3_amrfinder.bed
│   ├── example3_amrfinder_text.bed
│   ├── example3_amrfinder.txt
│   ├── test_test1_amrfinder.bed
│   ├── test_test1_amrfinder_text.bed
│   ├── test_test1_amrfinder.txt
│   ├── test_test2_amrfinder.bed
│   ├── test_test2_amrfinder_text.bed
│   └── test_test2_amrfinder.txt
├── plasmidfinder
│   ├── example1
│   │   ├── data.json
│   │   └── tmp
│   │       ├── out_enterobacteriaceae.xml
│   │       ├── out_Inc18.xml
│   │       ├── out_NT_Rep.xml
│   │       ├── out_Rep1.xml
│   │       ├── out_Rep2.xml
│   │       ├── out_Rep3.xml
│   │       ├── out_RepA_N.xml
│   │       ├── out_RepL.xml
│   │       └── out_Rep_trans.xml
│   ├── example2
│   │   ├── data.json
│   │   └── tmp
│   │       ├── out_enterobacteriaceae.xml
│   │       ├── out_Inc18.xml
│   │       ├── out_NT_Rep.xml
│   │       ├── out_Rep1.xml
│   │       ├── out_Rep2.xml
│   │       ├── out_Rep3.xml
│   │       ├── out_RepA_N.xml
│   │       ├── out_RepL.xml
│   │       └── out_Rep_trans.xml
│   ├── example3
│   │   ├── data.json
│   │   └── tmp
│   │       ├── out_enterobacteriaceae.xml
│   │       ├── out_Inc18.xml
│   │       ├── out_NT_Rep.xml
│   │       ├── out_Rep1.xml
│   │       ├── out_Rep2.xml
│   │       ├── out_Rep3.xml
│   │       ├── out_RepA_N.xml
│   │       ├── out_RepL.xml
│   │       └── out_Rep_trans.xml
│   ├── test_test1
│   │   ├── data.json
│   │   └── tmp
│   │       ├── out_enterobacteriaceae.xml
│   │       ├── out_Inc18.xml
│   │       ├── out_NT_Rep.xml
│   │       ├── out_Rep1.xml
│   │       ├── out_Rep2.xml
│   │       ├── out_Rep3.xml
│   │       ├── out_RepA_N.xml
│   │       ├── out_RepL.xml
│   │       └── out_Rep_trans.xml
│   └── test_test2
│       ├── data.json
│       └── tmp
│           ├── out_enterobacteriaceae.xml
│           ├── out_Inc18.xml
│           ├── out_NT_Rep.xml
│           ├── out_Rep1.xml
│           ├── out_Rep2.xml
│           ├── out_Rep3.xml
│           ├── out_RepA_N.xml
│           ├── out_RepL.xml
│           └── out_Rep_trans.xml
├── prep
│   ├── example1.fasta
│   ├── example2.fasta
│   ├── example3.fasta
│   ├── test_test1.fasta
│   └── test_test2.fasta
└── prokka
    ├── example1
    │   ├── example1.err
    │   ├── example1.faa
    │   ├── example1.ffn
    │   ├── example1.fna
    │   ├── example1.fsa
    │   ├── example1.gbk
    │   ├── example1.gff
    │   ├── example1.log
    │   ├── example1.sqn
    │   ├── example1.tbl
    │   ├── example1.tsv
    │   └── example1.txt
    ├── example2
    │   ├── example2.err
    │   ├── example2.faa
    │   ├── example2.ffn
    │   ├── example2.fna
    │   ├── example2.fsa
    │   ├── example2.gbk
    │   ├── example2.gff
    │   ├── example2.log
    │   ├── example2.sqn
    │   ├── example2.tbl
    │   ├── example2.tsv
    │   └── example2.txt
    ├── example3
    │   ├── example3.err
    │   ├── example3.faa
    │   ├── example3.ffn
    │   ├── example3.fna
    │   ├── example3.fsa
    │   ├── example3.gbk
    │   ├── example3.gff
    │   ├── example3.log
    │   ├── example3.sqn
    │   ├── example3.tbl
    │   ├── example3.tsv
    │   └── example3.txt
    ├── test_test1
    │   ├── test_test1.err
    │   ├── test_test1.faa
    │   ├── test_test1.ffn
    │   ├── test_test1.fna
    │   ├── test_test1.fsa
    │   ├── test_test1.gbk
    │   ├── test_test1.gff
    │   ├── test_test1.log
    │   ├── test_test1.sqn
    │   ├── test_test1.tbl
    │   ├── test_test1.tsv
    │   └── test_test1.txt
    └── test_test2
        ├── test_test2.err
        ├── test_test2.faa
        ├── test_test2.ffn
        ├── test_test2.fna
        ├── test_test2.fsa
        ├── test_test2.gbk
        ├── test_test2.gff
        ├── test_test2.log
        ├── test_test2.sqn
        ├── test_test2.tbl
        ├── test_test2.tsv
        └── test_test2.txt
```
</details>
