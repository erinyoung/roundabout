process prep {
    tag "${fasta}"
    publishDir params.outdir, mode: 'copy'
    container = "staphb/samtools:latest"

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

    input:
    tuple val(sample), file(gff), file(script)

    output:
    tuple val(sample), file("karyotype/${sample}_karyotype_genes.bed"), emit: genes
    tuple val(sample), file("karyotype/${sample}_karyotype.txt"), emit: karyotype

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

    input:
    tuple val(sample), file(hits), file(script)

    output:
    path "blastn/${sample}.sorted_divided.bed", emit: bed
    path "blastn/${sample}.starts_ends.txt"
    path "blastn/${sample}.divided.bed"

    shell:
    '''
    bash !{script} !{hits} !{sample}
    '''
}

process highlight {
    tag "highlight"
    publishDir params.outdir, mode: 'copy'

    input:
    file(hits)
    file(beds)
    file(script)

    output:
    path "highlight/highlights.txt", emit: txt

    shell:
    '''
    mkdir -p highlight combine
    cat *bed > combine/combined.bed
    cat *txt > combine/blast_hits.txt
    python !{script} !{params.length}
    mv highlights.txt highlight/highlights.txt
    '''
}

process prep_blast {
    tag "${sample}"
    publishDir params.outdir, mode: 'copy'
    
    input:
    tuple val(sample), file(hits)

    output:
    tuple val(sample), file("blastn/!{sample}.divided.bed")
    path "whatever"

    shell:
    '''
    mkdir blastn
    awk '{if ($11 == "0.0") print $1 "\\t" $7 "\\t" $8 "\\t" $2 "\\t" $9 "\\t" $10 }' !{hits} | sort -k 1,1 -k2,2n -k3,3n > blastn/!{sample}.blast_hits.bed

    #cut -f 1-3 blastn/!{sample}.blast_hits.bed | sort -k 1,1 -k2,2n -k3,3n | uniq > blastn/!{sample}.uniq.bed

    awk '{if ($11 == "0.0") print $1 "\\t" $7 "\\t" $8 }' blastn/!{sample}.blast_hits.bed | awk '{ if ( $3 - $2 > 1000) print $1 "\\tstart:\t" $2 "\\n" $1 "\tend:\\t" $3 }' | sort -k 1,1 -k3,3n | uniq > blastn/!{sample}.starts_ends.txt
    
    while read line
    do
        line_parts=($(echo $line))
        divisions=($(grep "${line_parts[0]}" blastn/!{sample}.starts_ends.txt | awk -v start=${line_parts[1]} -v end=${line_parts[2]} '{if ($3 >= start && $3 <= end ) print $2 $3 }'))
        prior_division=""
        for division in ${divisions[@]}
        do
            print $division
            start_end=$(echo $division | cut -f 1 -d ":" )
            value=$(echo $division | cut -f 2 -d ":" )
            if [ -n "$prior_division" ]
            then
                if [ "$start_end" == "start" ]
                then
                    echo -e "${line_parts[0]}\\t$prior_division\\t$((value -1))" >>  blastn/!{sample}.divided.bed
                    prior_division=$value
                elif [ "$start_end" == "end" ]
                then
                    echo -e "${line_parts[0]}\\t$prior_division\\t$value" >>  blastn/!{sample}.divided.bed
                    prior_division=$((value +1))
                fi
            else
                prior_division=$value
            fi  
        done
    done < <(awk '{if ($11 == "0.0") print $1 "\\t" $7 "\\t" $8 }' !{hits} | awk '{if ($3 - $2 > 1000) print $0}' | sort -k 1,1 -k2,2n -k3,3n | uniq )

'''
}
