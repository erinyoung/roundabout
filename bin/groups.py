#!/usr/bin/python3

import sys

try:
    length = int(sys.argv[1])
except:
    length = 500

divided={}
with open("combine/combined.bed", 'r') as file :
    for line in file.readlines() :
        chr, start, end=line.split()
        if (int(end) - int(start) >= length):
            if chr not in divided.keys():
                divided[chr]={}
            divided[chr][start] = end

i=0
group = {str(i) : {'list': [], 'length' : 0}}
lengths = []

with open("combine/blast_hits.txt", 'r') as file :
    for line in file.readlines() :
        chr1 = line.split()[0]
        chr2 = line.split()[1]
        start1 = line.split()[6]
        end1 = line.split()[7]
        start2 = line.split()[8]
        end2 = line.split()[9]
        for key in divided[chr1].keys():
            if int(start1) <= int(key) and int(end1) >= int(divided[chr1][key]):
                start_diff = int(key) - int(start1)
                end_diff = int(end1) - int(divided[chr1][key])
                saved_line= chr1 + ":" + key + ':' + divided[chr1][key]
                match_length = int(divided[chr1][key]) - int(key)             
                if int(start2) < int(end2):
                    match_start = int(start2) + start_diff
                    match_end = int(end2) - end_diff
                elif (int(start2) > int(end2)):
                    match_start = int(start2) - int(start_diff)
                    match_end = int(end2) + int(end_diff)
                
                if match_start > match_end:
                    temp_start=match_start
                    temp_end=match_end
                    match_start=temp_end
                    match_end=temp_start

                saved_line2= chr2 + ":" + str(match_start) + ":" + str(match_end)
                new_group = True
                keep_looking = True
                for numbers in group.keys():
                    if saved_line in group[numbers]['list'] and keep_looking:
                        if saved_line2 not in group[numbers]['list']:
                            group[numbers]['list'].append(saved_line2)
                        new_group = False
                        keep_looking = False
                    elif saved_line2 in group[numbers]['list'] and keep_looking:
                        group[numbers]['list'].append(saved_line)
                        new_group = False
                        keep_looking = False

                if new_group:
                    i += 1
                    group[str(i)] = { 'list': [saved_line, saved_line2], "length" : match_length }
                    lengths.append(match_length)
del group["0"]

if i > 25:
    minimum=sorted(lengths, reverse=True)[:26][-1]
    for key in list(group.keys()):
        if group[key]['length'] <= minimum:
            group.pop(key)

file = open("highlights.txt", "w")
for key in group.keys():
    for region in group[key]['list']:
        chr, start, end = region.split(":") 
        file.write(chr + "\t" + start +  "\t" + end + "\tfill_color=color" + key + "\n")
file.close()
