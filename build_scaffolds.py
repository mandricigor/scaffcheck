
from Bio import SeqIO
import sys


# load into memory the contigs
sw = {}
handle = open("IGOR.fa", "rU")
for record in SeqIO.parse(handle, "fasta"):
    sw[record.id] = record.seq.tostring()
handle.close()





with open("nucmerhit.coords") as f:
    a = f.readlines()

a = a[4:] # drop the header lines from nucmer hit file
a = map(lambda x: x.strip().split(), a)


# groups is the list of groups by reference
groups = [[a[0]]]
for i in range(1, len(a)):
    if a[i][11] == groups[-1][-1][11]:
        groups[-1].append(a[i])
    else:
        groups.append([a[i]])


# because we may have that some contigs 
allgroups = []
for group in groups:
    group2 = [[group[0]]]
    for j in range(1, len(group)):
        if group[j][12] == group2[-1][-1][12]:
            group2[-1].append(group[j])
        else:
            group2.append([group[j]])
    allgroups.append(group2)



contig_map = {}
global_counter = 0

scaffolds = []
for group in allgroups:
    scaf = []
    for group2 in group:
        if len(group2) == 1 and len(group2[0]) == 12:
            # if the contig does not match completely -> skip it
            # we do not need it anymore
            continue
        elif len(group2) == 1 and len(group2[0]) == 13:
            # it has probably a full match
            if group2[0][12] == "[CONTAINS]":
                if group2[0][12] in contig_map:
                    name = contig_map[group2[0][12]]
                else:
                    global_counter += 1
                    name = "contig_%s" % global_counter
                    contig_map[group2[0][12]] = name
                scaf.append(name + ":::" + group2[0][10])
        else:
            # here we may have contigs corresponding to circular genomes
            # if it has both [BEGINS] and [ENDS] -> we just split it
            # or we still have mis-assembly - just merge all hits
            if int(group2[-1][1]) - int(group2[0][0]) > 0.97 * int(group2[0][8]):
                if group2[0][12] in contig_map:
                    name = contig_map[group2[0][12]]
                else:
                    global_counter += 1
                    name = "contig_%s" % global_counter
                    contig_map[group2[0][12]] = name
                scaf.append(name + ":::" + group2[0][10])
    scaffolds.append(scaf)


rev_contig_map = {}
for x, y in contig_map.items():
    rev_contig_map[y] = x


with open("scaffolds.igor.scaf", "w") as f:
    counter = 1
    for scaffold in scaffolds:
        f.write(">scaffold_%s\n" % counter)
        for scaf in scaffold:
            f.write("%s\n" % scaf)
        counter += 1


with open("scaffolds.igor.fa", "w") as f:
    all_contigs = set()
    for scaffold in scaffolds:
        for contig in scaffold:
            all_contigs.add(contig.split(":::")[0])
    scontigs = []
    for sc in all_contigs:
        scontigs.append((sc, int(sc.split("_")[1])))
    for contig, order in sorted(scontigs, key=lambda z: z[1]):
        f.write(">%s\n" % contig)
        f.write("%s\n" % sw[rev_contig_map[contig]])









