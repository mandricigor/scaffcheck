#!/usr/bin/env python3

import argparse
import os
import copy
import shutil
import sys
from fastaq import sequences, utils, intervals, tasks


from Bio import SeqIO
import sys


# load contigs into memory
sw = {}
handle = open("assembly_contigs.fa", "rU")
for record in SeqIO.parse(handle, "fasta"):
    sw[record.id] = record.seq.tostring()
handle.close()



# check required nucmer programs are in path
progs = ['nucmer', 'delta-filter', 'show-coords']
not_in_path = [p for p in progs if shutil.which(p) is None]
    
if len(not_in_path):
    print('Error! Need these programs to be in your path:', file=sys.stderr)
    print('\n'.join(not_in_path), file=sys.stderr)
    print('Cannot continue', file=sys.stderr)
    sys.exit(1)


def nucmer_file_reader(fname):
    f = utils.open_file_read(fname)
    in_header = True

    for line in f:
        if in_header:
            if line.startswith('['):
                in_header = False
            continue
        yield NucmerHit(line)

    utils.close(f)


class NucmerHit:
    def __init__(self, line):
        # [S1]  [E1]    [S2]    [E2]    [LEN 1] [LEN 2] [% IDY] [LEN R] [LEN Q] [FRM]   [TAGS]
        #1162    25768   24536   4   24607   24533   99.32   640851  24536   1   -1  MAL1    NODE_25757_length_24482_cov_18.920391   [CONTAINS]
        try:
            l = line.rstrip().split('\t')
            self.ref_start = int(l[0])
            self.ref_end = int(l[1])
            self.qry_start = int(l[2])
            self.qry_end = int(l[3])
            self.hit_length_ref = int(l[4])
            self.hit_length_qry = int(l[5])
            self.percent_identity = float(l[6])
            self.ref_length = int(l[7])
            self.qry_length = int(l[8])
            self.frame = int(l[9])
            self.strand = int(l[10])
            self.ref_name = l[11]
            self.qry_name = l[12]

            if len(l) == 14:
                self.tag = l[13][1:-1]
            else:
                self.tag = None
        except:
            print('Error reading this nucmer line:\n' + line, file=sys.stderr)

    def __str__(self):
         return str(self.qry_name) + " : " + str(self.ref_start) + " - " + str(self.ref_end)



parser = argparse.ArgumentParser(
    description = 'Takes contigs and a reference sequence. Makes a new fasta file of the contigs, but they are now perfect sequences by using the reference instead',
    usage = '%(prog)s [options] <contigs.fa> <reference.fa> <outprefix>')
parser.add_argument('--min_seq_length', type=int, help='Minimum length of contig to output [%(default)s]', default=200)
parser.add_argument('--nucmer_options', help='Options when running nucmer [%(default)s]', default='')
parser.add_argument('contigs_fa', help='Name of contigs fasta file', metavar='contigs.fa')
parser.add_argument('ref_fa', help='Name of reference fasta file', metavar='reference.fa')
parser.add_argument('outprefix', help='Prefix of output files')
options = parser.parse_args()

ref_seqs = {}
tasks.file_to_dict(options.ref_fa, ref_seqs)

nucmer_out_prefix = options.outprefix + '.nucmer'
nucmer_out_delta = nucmer_out_prefix + '.delta'
nucmer_out_filter = nucmer_out_prefix + '.delta-filter'
nucmer_out_coords = nucmer_out_filter + '.coords'

# run nucmer of contigs vs ref
utils.syscall(' '.join(['nucmer', options.nucmer_options, '-p', nucmer_out_prefix, options.ref_fa, options.contigs_fa, '--maxmatch']))
utils.syscall(' '.join(['delta-filter', '-i 97 -l 180 -r', nucmer_out_delta, '>', nucmer_out_filter]))
utils.syscall(' '.join(['show-coords', '-dTlro', nucmer_out_filter, '>', nucmer_out_coords]))

# load hits into hash. key=ref_name, value=another hash with key=qry_name, value=list of hit positions in that ref seq
nucmer_hits = {}
contigs_to_print = {}




nucmer_reader = nucmer_file_reader(nucmer_out_coords)

# this dictionary contains hits for each contig per reference
hit_ref_dict = {}

for hit in nucmer_reader:
    hit_dict = hit_ref_dict.setdefault(hit.qry_name, {})
    ref_list = hit_dict.setdefault(hit.ref_name, [])
    ref_list.append(hit)



qry_segments = {}
for qry_name, qry_dict in hit_ref_dict.items():
    # determine if we have a full match, if so -> just add it to the list of contigs
    contains = False
    good_match = None
    for ref_name, ref_list in qry_dict.items():
        contains_contigs = [hit for hit in ref_list if hit.tag == "CONTAINS"]
        if contains_contigs:
            contains = True
            good_match = contains_contigs[0]
            break
    if contains:
        # if we have a hit for a contig which is completely contained 
        # then we just consider it -> we want to take biggest hit
        # and ignore the rest of them
        a, b = (good_match.qry_start, good_match.qry_end)
        if a > b:
            qry_segments[qry_name] = [(b, a)]
        else:
            qry_segments[qry_name] = [(a, b)]
    else:
        segments = []
        for ref_name, ref_list in qry_dict.items():
            hits = sorted(ref_list, key=lambda z: z.ref_start)
            # merge all overlapping hits for each reference
            merging = [[0]]
            for i in range(1, len(hits)):
                if hits[i].ref_start <= hits[merging[-1][-1]].ref_end + 1:
                    if hits[i].ref_end > hits[merging[-1][-1]].ref_end:
                        merging[-1].append(i)
                else:
                    merging.append([i])
            for merge_group in merging:
                a, b = (hits[merge_group[0]].qry_start, hits[merge_group[-1]].qry_end)
                if a > b:
                    segments.append((b, a))
                else:
                    segments.append((a, b))
        qry_segments[qry_name] = segments



###################################################

# print the final perfect contigs
with open(options.outprefix + '.fa', "w") as f:
    for x, y in qry_segments.items():
        if len(y) > 1:
            counter = 1
            for u, v in y:
                f.write(">%s\n" % (x + "_" + str(counter)))
                f.write("%s\n" % sw[x][(u - 1):(v - 1)])
                counter += 1
        else:
            f.write(">%s\n" % x)
            f.write("%s\n" % sw[x][(y[0][0] - 1):(y[0][1] - 1)])


