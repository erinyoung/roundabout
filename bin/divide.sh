#!/usr/bin/env bash

hits=$1
sample=$2

mkdir -p blastn

#awk '{if ($11 == "0.0") print $1 "\t" $7 "\t" $8 "\t" $2 "\t" $9 "\t" $10 }' $hits | sort -k 1,1 -k2,2n -k3,3n > blastn/$sample.blast_hits.bed

awk '{if ($11 == "0.0") print $1 "\t" $7 "\t" $8 }' $hits | awk '{ if ( $3 - $2 > 1000) print $1 "\tstart:\t" $2 "\n" $1 "\tend:\t" $3 }' | sort -k 1,1 -k3,3n | uniq > blastn/$sample.starts_ends.txt

while read line
do
    line_parts=($(echo $line))
    divisions=($(grep "${line_parts[0]}" blastn/$sample.starts_ends.txt | awk -v start=${line_parts[1]} -v end=${line_parts[2]} '{if ($3 >= start && $3 <= end ) print $2 $3 }'))
    prior_division=""
    for division in ${divisions[@]}
    do
        start_end=$(echo $division | cut -f 1 -d ":" )
        value=$(echo $division | cut -f 2 -d ":" )
        if [ -n "$prior_division" ]
        then
            if [ "$start_end" == "start" ]
            then
                echo -e "${line_parts[0]}\t$prior_division\t$((value -1))" >>  blastn/$sample.divided.bed
                prior_division=$value
            elif [ "$start_end" == "end" ]
            then
                    echo -e "${line_parts[0]}\t$prior_division\t$value" >>  blastn/$sample.divided.bed
                    prior_division=$((value +1))
            fi
        else
            prior_division=$value
        fi  
    done
done < <(awk '{if ($11 == "0.0") print $1 "\t" $7 "\t" $8 }' $hits | awk '{if ($3 - $2 > 1000) print $0}' | sort -k 1,1 -k2,2n -k3,3n | uniq )

cat blastn/$sample.divided.bed | sort -k 1,1 -k 2n,2 -k 3n,3 | uniq > blastn/$sample.sorted_divided.bed
