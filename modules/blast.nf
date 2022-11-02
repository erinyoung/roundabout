process blastn {
  tag "${sample}"
  container 'ncbi/blast:latest'
  publishDir params.outdir, mode: 'copy'

  input:
  val(sample)
  file(fasta)

  output:
  tuple val(sample), file("blastn/${sample}.blast_hits.txt"), emit: hits
  path "blastn/${sample}.blast_hits.txt", emit: txt
  path "blastn/${sample}*"

  shell:
  '''
    mkdir -p blastn
    for fasta in $(ls !{fasta} | grep -Fv "!{sample}." )
    do
      cat $fasta >> blastn/non_!{sample}.fasta
      echo "" >> blastn/non_!{sample}.fasta
    done

    blastn -query !{sample}.* -subject blastn/non_!{sample}.fasta -outfmt 6 -ungapped -out blastn/!{sample}.blast_hits.txt -evalue 1e-25
  '''
}

//-gapopen 10 -gapextend 6
//  -gapopen <Integer>
//    Cost to open a gap
//  -gapextend <Integer>
//    Cost to extend a gap
//  -penalty <Integer, <=0>
//    Penalty for a nucleotide mismatch
//  -reward <Integer, >=0>
//    Reward for a nucleotide match

