#!/usr/bin/env nextflow

VERSION = "0.1.20220915"

params.outdir             = 'roundabout'
params.fastas             = 'plasmids' 
params.length             = ''

//# include the plasmidfinder fasta backbone in the comparison?
//# allow blast for specific gene(s) of interest

Channel
  .fromPath("${params.fastas}/*{.fa,.fasta,.fna}")
  .view { "input fasta file : " + it}
  .set { fastas }

circos_confs      = Channel.fromPath( workflow.projectDir + "/conf",              type: 'dir' )
divide_script     = Channel.fromPath( workflow.projectDir + "/bin/divide.py",     type: 'file')
karyotype_script  = Channel.fromPath( workflow.projectDir + "/bin/karyotype.py",  type: 'file')
highlights_script = Channel.fromPath( workflow.projectDir + "/bin/groups.py",     type: 'file')

include { amrfinder }                           from './modules/amrfinder'      addParams(outdir: params.outdir)
include { plasmidfinder }                       from './modules/plasmidfinder'  addParams(outdir: params.outdir)
include { prokka }                              from './modules/prokka'         addParams(outdir: params.outdir)
include { blastn }                              from './modules/blast'          addParams(outdir: params.outdir)
include { divide; karyotype; highlight; prep }  from './modules/roundabout'     addParams(outdir: params.outdir, length: params.length)
include { bedtools_nuc as nuc }                 from './modules/bedtools'       addParams(outdir: params.outdir)
include { circos }                              from './modules/circos'         addParams(outdir: params.outdir)
include { pmauve; mummer; mmseqs; dnadiff }     from './modules/pygenomeviz'    addParams(outdir: params.outdir)

workflow {
  prep(fastas)
  prepped_fastas=prep.out.flatten()
  fasta_collection = prepped_fastas.collect()
  pmauve(fasta_collection)

  input = prepped_fastas.map{ it -> tuple(it.baseName, it) }
  plasmidfinder(input)
  amrfinder(input)
  prokka(input)
  karyotype(prokka.out.gff.combine(karyotype_script))
  mummer(prokka.out.gbk.collect())
  mmseqs(prokka.out.gbk.collect())

  blastn(input.map{ it -> it[0]}, fasta_collection )
  dnadiff(input.map{ it -> it[0]}, fasta_collection )

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