#!/usr/bin/python3

import sys

hits=sys.argv[1]
try:
    minlength = int(sys.argv[2])
except:
    minlength = 500

# blast file columns
# 0: qseqid, 1: sseqid, 2: pident, 3: length, 4: mismatch, 5: gapopen, 6: qstart, 7: qend, 8: sstart, 9: send, 10: evalue, 11: bitscore

starts = []
ends   = []
chr    = ''
with open(hits, 'r') as file :
    for line in file.readlines() :
        if int(line.split()[3]) >= minlength and float(line.split()[2]) >= 90.0:
            chr   = line.split()[0]
            start = int(line.split()[6])
            end   = int(line.split()[7])
            if start not in starts:
                starts.append(start)
            if end not in ends:
                ends.append(end)

starts.sort()
ends.sort()

divisions = []
with open(hits, 'r') as file :
    for line in file.readlines() :
        if int(line.split()[3]) >= minlength and float(line.split()[2]) >= 90.0:
            start = int(line.split()[6])
            end   = int(line.split()[7])
            mids  = []
            for mid in starts:
                if start <= mid < end:
                    mids.append(mid)
            for mid in ends:
                if start < mid < end:
                    mids.append(mid+1)
            mids.sort()
            if len(mids) == 1:
                if start - end >= minlength:
                    div=str(start) + "-" + str(end)
                    if div not in divisions:
                        divisions.append(div)
            else:
                for i in range(1, len(mids)):
                    if mids[i]-1 - mids [i-1] >= minlength:
                        div=str(mids[i-1]) + "-" + str(mids[i]-1)
                        if div not in divisions:
                            divisions.append(div)
                if end - mids[-1] >= minlength:
                    div=str(mids[-1]) + "-" + str(end)
                    if div not in divisions:
                        divisions.append(div)

file = open(chr + ".divided.bed", "w")
for div in divisions:
    start,end = div.split("-")
    file.write(chr + "\t" + start + "\t" + end + "\n")
