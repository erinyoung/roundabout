process pmauve {
    tag "${sample}"
    publishDir params.outdir, mode: 'copy'

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
    tag "${sample}"
    publishDir params.outdir, mode: 'copy'

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
    tag "${sample}"
    publishDir params.outdir, mode: 'copy'

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

