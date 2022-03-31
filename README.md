# roundabout

The fundamental goal of this repo was to create something quick and reproducible to find large regions of similarity between highly related plasmids - like those in an outbreak.

roundabout uses blast to find regions of similarity (as opposed to groups of genes that would be required for synteny), and then assigns those regions a color. Those regions and colors can then be visualized via circos.

dependencies:
- blast : find similarities
- circos : visualizing the end product (optional)
- bedtools : combining regions of interest
- awk/sed/bash : a lot of file manipulation
- samtools : to find genome size

INSTALL:

```
git clone https://github.com/erinyoung/roundabout.git

export PATH=$PATH:$(pwd)/roundabout/bin

# Using conda to install dependencies
conda create -n roundabout -c bioconda -c defaults samtools bedtools circos blast

# then activate the environment with
conda activate roundabout
```

USAGE:

- Put the completed/closed plasmids into a single directory (i.e. plasmids)
- (optional) Put AMRFinderPlus output in a directory (i.e. amrfinder)
- (optional) Put gff file from prokka/bakta/etc in a directory (i.e. gff)

Note : All files must have the same prefix. `plasmids/$sample.{fasta,fa,fna}`, `amrfinder/$sample*`, `gff/$sample*gff` respectively.

```
roundabout -d <directory with plasmid sequences>

# example
roundabout -d plasmids

# with amrfinder results
roundabout -a amrfinder

# with gff files
roundabout -g gff

# with both amrfinder results and gff files
roundabout -a amrfinder -g gff
```

This was created for personal use with specific projects in mind, as opposed to general use. As such, other users may notice that the use case is highly specific and informal. Put in an issue if this is something that interests you and you need an additional feature.

Future directions:
- removing the samtools dependency
- adding parallelization through either gnu parallel or snakemake
- creating a docker container
