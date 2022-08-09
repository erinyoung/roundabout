#!/usr/bin/env nextflow

params.config_directory   = workflow.projectDir + "/configs/"
params.outdir             = 'roundabout'
params.plasmid_directory  = 'plasmids' 
params.length             = ''

//# run plasmidfinder on each fasta - do I want to include the fasta backbone in the comparison?
//# visualize with pygenomeviz
//# contaienrs for all processes : karyotype, divide, highlight
//# choose colors

Channel
  .fromPath("${params.plasmid_directory}/*{.fa,.fasta,.fna}")
  .view { "input fasta file : " + it}
  .set { fastas }

circos_confs      = Channel.fromPath( workflow.projectDir + "/../conf", type:'dir' )
divide_script     = Channel.fromPath( workflow.projectDir + "/../bin/divide.sh", type: 'file' )
karyotype_script  = Channel.fromPath( workflow.projectDir + "/../bin/karyotype.py", type: 'file' )
highlights_script = Channel.fromPath( workflow.projectDir + "/../bin/groups.py")

include { amrfinder }                           from './modules/amrfinder'  addParams(outdir: params.outdir)
include { prokka }                              from './modules/prokka'     addParams(outdir: params.outdir)
include { blastn }                              from './modules/blast'      addParams(outdir: params.outdir)
include { divide; karyotype; highlight; prep }  from './modules/roundabout' addParams(outdir: params.outdir, length: params.length)
include { bedtools_nuc as nuc }                 from './modules/bedtools'   addParams(outdir: params.outdir)
include { circos }                              from './modules/circos'     addParams(outdir: params.outdir)

workflow {
  prep(fastas)
  prepped_fastas=prep.out.flatten()
  fasta_collection = prepped_fastas.collect()

  input = prepped_fastas.map{ it -> tuple(it.baseName, it) }
  amrfinder(input)
  prokka(input)
  karyotype(prokka.out.gff.combine(karyotype_script))

  blastn(input.map{ it -> it[0]}, fasta_collection )

  divide(blastn.out.hits.combine(divide_script))
  nuc(input)
  
  highlight(blastn.out.txt.collect(),divide.out.bed.collect(),highlights_script)

  karyotype.out.karyotype
    .join(amrfinder.out.bed)
    .join(nuc.out.skew)
    .combine(highlight.out.txt)
    .combine(circos_confs)
    .set{ circos_files}

  circos(circos_files)
}