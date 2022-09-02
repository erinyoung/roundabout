process prep {
    tag "${fasta}"
    publishDir params.outdir, mode: 'copy'
    container "staphb/samtools:latest"

    input:
    file(fasta)

    output:
    path"prep/*fasta", emit: fastas

    shell:
    '''
    mkdir -p prep

    sample=$(echo !{fasta} | rev | cut -f 2- -d "." | rev)

    grep ">" !{fasta} | awk '{print $1}'
    headers=($(grep ">" !{fasta} | awk '{print $1}' | sed 's/>//g'))
    for header in ${headers[@]}
    do
        final_header=$(echo $header | sed 's/>$sample//g' )

        if [ "$final_header" == "$sample" ]
        then 
            final_header=$sample
        else
            final_header="${sample}_${header}"
        fi
        echo ">$final_header" > prep/$final_header.fasta
        samtools faidx !{fasta} $header | grep -v ">" | fold -w 75 >> prep/$final_header.fasta
        echo "" >> prep/$final_header.fasta
    done
    '''
}

process karyotype {
    tag "${sample}"
    publishDir params.outdir, mode: 'copy'
    container "python:3.10.6"

    input:
    tuple val(sample), file(gff), file(script)

    output:
    tuple val(sample), file("karyotype/${sample}_karyotype_genes.bed"), emit: genes
    tuple val(sample), file("karyotype/${sample}_karyotype.txt"),       emit: karyotype

    shell:
    '''
    mkdir -p karyotype
    python !{script} !{gff}
    cp karyotype_genes.bed karyotype/!{sample}_karyotype_genes.bed
    cp karyotype.txt karyotype/!{sample}_karyotype.txt
    '''
}

process divide {
    tag "${sample}"
    publishDir params.outdir, mode: 'copy'
    container "python:3.10.6"

    input:
    tuple val(sample), file(hits), file(script)

    output:
    path "blastn/${sample}.divided.bed", emit: bed

    shell:
    '''
    mkdir blastn
    python !{script} !{hits} !{params.length}
    mv !{sample}.divided.bed blastn/.
    '''
}

process highlight {
    tag "highlight"
    publishDir params.outdir, mode: 'copy'
    container "python:3.10.6"

    input:
    file(hits)
    file(beds)
    file(script)

    output:
    path "highlight/highlights*.txt", emit: txt

    shell:
    '''
    mkdir -p highlight combine
    cat *bed > combine/combined.bed
    cat *txt > combine/blast_hits.txt
    python !{script} !{params.length}
    mv highlights.txt highlight/highlights.txt
    mv highlights_text.txt highlight/highlights_text.txt
    '''
}

