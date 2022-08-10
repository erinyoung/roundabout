process plasmidfinder {
  tag "${sample}"
  container 'staphb/plasmidfinder:latest'
  publishDir params.outdir, mode: 'copy'

  input:
  tuple val(sample), file(file)

  output:
  path "plasmidfinder/${sample}/*", emit: files

  shell:
  '''
  mkdir -p plasmidfinder/!{sample}

  plasmidfinder.py \
    -i !{file} \
    -o plasmidfinder/!{sample} 
  '''
}
