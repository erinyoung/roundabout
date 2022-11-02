process pmauve {
    tag "pgv-pmauve"
    publishDir params.outdir, mode: 'copy'
    container 'staphb/pygenomeviz:latest'
    errorStrategy 'ignore'

    input:
    file(files)

    output:
    path "pgv/pmauve/*"

    shell:
    '''
    mkdir -p pgv/pmauve

    pgv-pmauve --seq_files !{files} \
        --outdir pgv/pmauve \
        --cmap viridis \
        --tick_style axis \
        --curve \
        --normal_link_color '#323e4f' \
        --inverted_link_color '#3d6e70' \
        --dpi 900
    '''
}

process mummer {
    tag "pgv-mummer"
    publishDir params.outdir, mode: 'copy'
    container 'staphb/pygenomeviz:latest'
    errorStrategy 'ignore'

    input:
    file(files)

    output:
    path "pgv/mummer/*"

    shell:
    '''
    mkdir -p pgv/mummer

    pgv-mummer --gbk_resources !{files} \
        --outdir pgv/mummer \
        --tick_style axis \
        --curve \
        --feature_plotstyle arrow \
        --feature_color black \
        --inverted_link_color '#3d6e70' \
        --normal_link_color '#323e4f' \
        --dpi 900
    '''
}

process mmseqs {
    tag "pgv-mmseqs"
    publishDir params.outdir, mode: 'copy'
    container 'staphb/pygenomeviz:latest'
    cpus 8
    errorStrategy 'ignore'

    input:
    file(files)

    output:
    path "pgv/mmseqs/*"

    shell:
    '''
    mkdir -p pgv/mmseqs

    pgv-mmseqs --gbk_resources !{files} \
        --outdir pgv/mmseqs \
        --tick_style axis \
        --curve \
        --feature_plotstyle arrow \
        --feature_color black \
        --inverted_link_color "#3d6e70" \
        --normal_link_color "#323e4f" \
        --thread_num !{task.cpus} \
        --dpi 900
    '''
}

// yeah... I know I should probably use a new container, but mummer is in this one already
process dnadiff {
    tag "${sample}"
    publishDir params.outdir, mode: 'copy'
    container 'staphb/pygenomeviz:latest'

    input:
    val(sample)
    file(fasta)

    output:
    path "dnadiff/${sample}/*"

    shell:
    '''
    mkdir -p dnadiff/!{sample}

    for fasta in $(ls !{fasta} | grep -Fv "!{sample}.fasta" )
    do
      sample2=$(echo $fasta | rev | cut -f 2- -d "." | rev )
      dnadiff !{sample}.fasta $fasta --prefix !{sample}-$sample2

      mv !{sample}-$sample2* dnadiff/!{sample}/.
    done
    '''
}