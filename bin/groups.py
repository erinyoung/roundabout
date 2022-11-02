#!/usr/bin/python3

import sys
#from xmlrpc.server import MultiPathXMLRPCServer

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
            divided[chr][start] = int(end)

if not divided:
    print("No blast matches found!")
    exit(0)

# blast file columns
# 0: qseqid, 1: sseqid, 2: pident, 3: length, 4: mismatch, 5: gapopen, 6: qstart, 7: qend, 8: sstart, 9: send, 10: evalue, 11: bitscore

i=0
group = {str(i) : {'list': [], 'length' : 0}}
lengths = []
matches = []
with open("combine/blast_hits.txt", 'r') as file :
    for line in file.readlines() :
        qseqid      = line.split()[0]
        sseqid      = line.split()[1]
        qstart      = int(line.split()[6])
        qend        = int(line.split()[7])
        sstart      = int(line.split()[8])
        send        = int(line.split()[9])

        for start in divided[qseqid].keys():
            saved_line          = ''
            start_diff          = 0
            end_diff            = 0

            # looking to see if potential groups in divided apply to this line
            if int(start) >= qstart and divided[qseqid][start] <= qend:
                start_diff      = int(start) - qstart
                end_diff        = qend - divided[qseqid][start]                
                saved_line      = qseqid + ":" + start + ':' + str(divided[qseqid][start])
                match_length    = divided[qseqid][start] - int(start)

                if sstart < send:
                    match_start = sstart + start_diff
                    match_end   = send   - end_diff
                elif sstart > send:
                    match_start = send   + end_diff
                    match_end   = sstart - start_diff

                match_found = False
                if str(match_start) not in divided[sseqid].keys():
                    for start_fix in divided[sseqid].keys():
                        if int(start_fix) - 50 <= int(match_start) <= int(start_fix) + 50 :
                            match_start     = start_fix
                            if divided[sseqid][start_fix] - 50 <= match_end <= divided[sseqid][start_fix] + 50 :
                                match_end   = divided[sseqid][start_fix]
                                match_found = True
                            else :
                                print("no match was found for 2 : " + line)
                else:
                    if divided[sseqid][str(match_start)] - 50 <= match_end <= divided[sseqid][str(match_start)] + 50 :
                        match_end   = divided[sseqid][str(match_start)]
                        match_found = True

                saved_line2 = sseqid + ":" + str(match_start) + ":" + str(match_end)

                if match_found :
                    new_group = True
                    for saved in group.keys():
                        if saved_line in group[saved]['list']:
                            new_group = False
                            if (saved_line2 not in group[saved]['list']):
                                group[saved]['list'].append(saved_line2)      
                        elif saved_line2 in group[saved]['list']:
                            new_group = False
                            if saved_line not in group[saved]['list']:
                                group[saved]['list'].append(saved_line)

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