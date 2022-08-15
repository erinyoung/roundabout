#!/usr/bin/python3

import sys

try:
    minlength = int(sys.argv[1])
except:
    minlength = 500

divided = {}
with open("combine/combined.bed", 'r') as file :
    for line in file.readlines() :
        chr, start, end=line.split()
        if (int(end) - int(start) >= minlength):
            if chr not in divided.keys():
                divided[chr]={}
            divided[chr][start] = end

i=0
group = {str(i) : {'list': [], 'length' : 0}}
lengths = []
matches = []
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
                    match_start = int(end2) + int(end_diff)
                    match_end = int(start2) - int(start_diff)

                #divided[chr][start]
                if str(match_start) not in divided[chr2].keys():
                    for potential_start in divided[chr2].keys():
                        if (int(potential_start) - 20) <= match_start <= (int(potential_start) + 20):
                            match_start=int(potential_start)

                # sometimes there are small indels
                if str(match_start) in divided[chr2].keys():
                    if int(divided[chr2][str(match_start)]) != int(match_end):
                        match_end=int(divided[chr2][str(match_start)])

                saved_line2= chr2 + ":" + str(match_start) + ":" + str(match_end)
                new_group = True

                for numbers in group.keys():
                    if saved_line in group[numbers]['list']:
                        new_group = False
                        if str(match_start) in divided[chr].keys():
                            if (saved_line2 not in group[numbers]['list']):
                                group[numbers]['list'].append(saved_line2)      
                    elif saved_line2 in group[numbers]['list']:
                        new_group = False
                        if saved_line not in group[numbers]['list']:
                            group[numbers]['list'].append(saved_line)

                if new_group:
                    i += 1
                    group[str(i)] = { 'list': [saved_line, saved_line2], "length" : match_length }
                    lengths.append(match_length)
del group["0"]

# adjusting duplicate groups
delete = []
k = i
while k > 0 :
    keep_looking = True
    for region in group[str(k)]['list']:
        if keep_looking:
            for grp in list(group.keys()):
                if (region in group[grp]['list']) and (int(grp) < int(k)) and (keep_looking):
                    group[grp]['list'] = group[grp]['list'] + group[str(k)]['list']
                    delete.append(str(k))
                    keep_looking = False
                    group[grp]['list'] = list(set(group[grp]['list']))
    k = k - 1

for element in delete:
    del group[element]

# filtering for the top 25 groups by length
if i > 25:
    minimum=sorted(lengths, reverse=True)[:26][-1]
    for key in list(group.keys()):
        if group[key]['length'] <= minimum:
            group.pop(key)

# giving each group a number
j = 1
file = open("highlights.txt", "w")
text_file = open("highlights_text.txt", "w")
for key in group.keys():
    for region in group[key]['list']:
        chr, start, end = region.split(":") 
        file.write(chr + "\t" + start +  "\t" + end + "\tfill_color=color" + str(j) + "\n")
        text_file.write(chr + " " + start +  " " + end + " g" + str(j) + "\n")
    j +=1
file.close()
text_file.close()