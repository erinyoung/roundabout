process bedtools_nuc {
  tag "${sample}"
  container 'staphb/bedtools:latest'
  publishDir params.outdir, mode: 'copy'

  input:
  tuple val(sample), file(fasta)

  output:
  tuple val(sample), file("bedtools/${sample}.{AC,GC,skew}.bed"), emit: skew

  path "bedtools/${sample}*"

  shell:
  '''
    mkdir bedtools

    grep -v ">" !{fasta} | wc -c
    length=$(grep -v ">" !{fasta} | wc -c | awk '{print $1}')

    echo -e "#chr\\tstr\\tend" > bedtools/!{sample}.windows.bed

    for ((i=1;i<=length;i+=500))
    do
        if [ "$i" -lt "$((length - 499 ))" ]
        then
            echo -e "!{sample}\\t$i\\t$((i + 499 ))" >> bedtools/!{sample}.windows.bed
        else
            echo -e "!{sample}\\t$i\\t$length" >> bedtools/!{sample}.windows.bed
        fi
    done
    bedtools nuc -fi !{fasta} -bed bedtools/!{sample}.windows.bed > bedtools/!{sample}.nuc.bed

    grep -hv "#" bedtools/!{sample}.nuc.bed | awk '{print $1 " " $2 " " $3 " " $4 }' > bedtools/!{sample}.AC.bed
    grep -hv "#" bedtools/!{sample}.nuc.bed | awk '{print $1 " " $2 " " $3 " " $5 }' > bedtools/!{sample}.GC.bed
    #GC Skew is calculated as (G - C) / (G + C)
    grep -hv "#" bedtools/!{sample}.nuc.bed | awk '{print $1 " " $2 " " $3 " " $8 - $7 " " $7 + $8 }' | awk '{print $1 " " $2 " " $3 " " $4 / $5 }' > bedtools/!{sample}.skew.bed
  '''
}
