process circos {
    tag "${sample}"
    container = 'quay.io/biocontainers/circos:0.69.8--hdfd78af_1'
    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(sample), file(karyotype), file(amr), file(skew), file(highlight), file(highlight_text), path(conf)

    output:
    path "circos/${sample}/data/*"
    path "circos/${sample}/*roundabout*"
    path "circos/*png"


    shell:
    '''
    mkdir -p circos/!{sample} data

    cat !{karyotype} | awk '{if ($5 != $6) print $0}' > data/karyotype.txt
    mv !{sample}_amrfinder_text.bed data/amrfinder_text.txt
    mv !{sample}_amrfinder.bed data/amrfinder.txt
    mv *AC.bed data/AC.bed
    mv *GC.bed data/GC.bed
    mv *skew.bed data/skew.bed
    mv highlights.txt data/highlights.txt
    mv highlights_text.txt data/highlights_text.txt

    circos -conf conf/template.conf
    mv *svg circos/!{sample}/!{sample}_roundabout.svg
    mv *png circos/!{sample}/!{sample}_roundabout.png
    cp circos/!{sample}/!{sample}_roundabout.png circos/.
    mv data circos/!{sample}/.
    '''
}