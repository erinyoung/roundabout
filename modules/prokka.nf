process prokka {
  tag "${sample}"
  container 'staphb/prokka:latest'
  publishDir params.outdir, mode: 'copy'

  input:
  tuple val(sample), file(contigs)

  output:
  path "prokka/${sample}/*"                               
  tuple val(sample), file("prokka/${sample}/${sample}.gff"),  emit: gff
  path "prokka/${sample}/${sample}.gbk",                      emit: gbk

  shell:
  '''
    prokka \
      --cpu !{task.cpus} \
      --outdir prokka/!{sample} \
      --prefix !{sample} \
      --addgenes \
     !{contigs} 
  '''
}
