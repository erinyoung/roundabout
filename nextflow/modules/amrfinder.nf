process amrfinder {
  tag "${sample}"
  container 'staphb/ncbi-amrfinderplus:latest'
  publishDir params.outdir, mode: 'copy'

  input:
  tuple val(sample), file(contigs)

  output:
  path "ncbi-AMRFinderplus/${sample}_amrfinder.txt"
  tuple val(sample), file("ncbi-AMRFinderplus/${sample}_amrfinder*.bed"), emit: bed

  shell:
  '''
    mkdir -p ncbi-AMRFinderplus

    amrfinder \
      --nucleotide !{contigs} \
      --threads !{task.cpus} \
      --name !{sample} \
      --output ncbi-AMRFinderplus/!{sample}_amrfinder.txt

    awk '{print $3 " " $4 " " $5 }' ncbi-AMRFinderplus/!{sample}_amrfinder.txt | grep -v "identifier" >> ncbi-AMRFinderplus/!{sample}_amrfinder.bed
    awk '{print $3 " " $4 " " $5 " " $7 }' ncbi-AMRFinderplus/!{sample}_amrfinder.txt | grep -v "identifier" > ncbi-AMRFinderplus/!{sample}_amrfinder_text.bed
  '''
}
