#!/usr/bin/python3

# This is just for creating the karyotype band colors. It's not a true karyotype and just for aesthetics and some visual matching
import sys

i=0
ranges={str(i): {'start':0, 'end':0, 'color':'black'}}

genes = open("karyotype_genes.bed", "w")

fasta_line=False
sequence=""
with open(sys.argv[1], 'r') as file :
    for line in file.readlines() :
        if "CDS" in line:
            chr=line.split()[0]
            start=int(line.split()[3])
            end=int(line.split()[4])
            strand=line.split()[6]
            info=line.split("\t")[8]
            info_dict = dict(inf.split("=") for inf in info.split(";"))

            # assinging colors
            color='black'
            if 'hypothetical protein' in info_dict['product']:
                color='white'
            elif 'transposase' in info_dict['product']:
                color='gpos100'
            elif 'recombinase' in info_dict['product']:
                color='gpos75'
            elif 'ase' in info_dict['product'] :
                color='gpos50'
            else:
                color='gpos25'

            if 'Name' in info_dict.keys():
                name=info_dict['Name']
            else:
                name="hypothetical protein"
            
            genes.write(chr + "\t" + str(start) + "\t" + str(end) + "\t" + strand + "\t" + name + "\n")
            
            # merging overlapping ranges (in likely the worst way possible)
            if start >= ranges[str(i)]['start'] and start < ranges[str(i)]['end']:
                if ranges[str(i)]['end'] < end:
                    ranges[str(i)]['end'] = end
                    if color != ranges[str(i)]['color']:
                        if color == 'gpos25' and (ranges[str(i)]['color'] == 'white' ):
                            ranges[str(i)]['color'] = color
                        elif color == 'gpos50' and (ranges[str(i)]['color'] == 'white' or ranges[str(i)]['color'] == 'gpos25'):
                            ranges[str(i)]['color'] = color
                        elif color == 'gpos75' and (ranges[str(i)]['color'] == 'white' or ranges[str(i)]['color'] == 'gpos25' or ranges[str(i)]['color'] == 'gpos50'):
                            ranges[str(i)]['color'] = color
                        elif color == 'gpos100':
                            ranges[str(i)]['color'] = color
            else:
                if (ranges[str(i)]['end']) < (start):
                    ranges[str(i) + "_between"] = {"start": ranges[str(i)]['end'], "end": start , "color": 'gneg' }
                i += 1
                ranges[str(i)] = {"start": start, "end": end, "color": color }
        if fasta_line:
            sequence += line
        if ">" in line:
            fasta_line=True
genes.close()

# adding the final area
if ranges[str(i)]['end'] < len(sequence):
    ranges['fin'] = {"start": ranges[str(i)]['end'],"end": len(sequence), "color": "gneg"}

del ranges["0"]

file = open("karyotype.txt", "w")
file.write("chr - " + chr + " " + chr + " 0 " + str(len(sequence)) + " black\n")

for key in ranges.keys():
    file.write("band " + chr + " " + key + " " + key + " " + str(ranges[key]['start']) + " " + str(ranges[key]['end']) + " " + ranges[key]['color'] + "\n")
file.close()