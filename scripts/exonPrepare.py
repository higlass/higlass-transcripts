from __future__ import print_function

__author__ = "Alaleh Azhir,Peter Kerpedjiev,Alexander Veit"

import collections as col
import sys
import argparse
from gtfparse import read_gtf


class GeneInfo:
    def __init__(self):
        pass


def merge_gene_info(gene_infos, gene_info):
    """
    Add a new gene_info. If it's txStart and txEnd overlap with a previous entry for this
    gene, combine them.
    """
    merged = False

    for existing_gene_info in gene_infos[gene_info.geneId]:
        if (
            existing_gene_info.chrName == gene_info.chrName
            and existing_gene_info.txEnd > gene_info.txStart
            and gene_info.txEnd > existing_gene_info.txStart
        ):

            # overlapping genes, merge the exons of the second into the first
            existing_gene_info.txStart = min(
                existing_gene_info.txStart, gene_info.txStart
            )
            existing_gene_info.txEnd = max(existing_gene_info.txEnd, gene_info.txEnd)

            for (exon_start, exon_end) in gene_info.exonUnions:
                existing_gene_info.exonUnions.add((exon_start, exon_end))

            merged = True

    if not merged:
        gene_infos[gene_info.geneId].append(gene_info)

    return gene_infos


def main():
    parser = argparse.ArgumentParser(
        description="""

    python ExonUnion.py Calculate the union of the exons of a list
    of transcript.

    chr10   27035524        27150016        ABI1    76      -       NM_001178120    10006   protein-coding  abl-interactor 1        27037498        27149792        10      27035524,27040526,27047990,27054146,27057780,27059173,27060003,27065993,27112066,27149675,      27037674,27040712,27048164,27054247,27057921,27059274,27060018,27066170,27112234,27150016,
"""
    )

    parser.add_argument("--gtf_file")
    parser.add_argument("--transcript_bed")

    args = parser.parse_args()

    # returns GTF with essential columns such as "feature", "seqname", "start", "end"
    # alongside the names of any optional keys which appeared in the attribute column
    df = read_gtf(args.gtf_file)
    
    

    inputFile = open(args.transcript_bed, "r")

    gene_infos = col.defaultdict(list)

    for line in inputFile:
        words = line.strip().split("\t")

        gene_info = GeneInfo()

        try:
            gene_info.chrName = words[0]
            gene_info.txStart = words[1]
            gene_info.txEnd = words[2]
            gene_info.geneName = words[3]
            gene_info.score = words[4]
            gene_info.strand = words[5]
            gene_info.refseqId = words[6]
            gene_info.geneId = words[7]
            gene_info.geneType = words[8]
            gene_info.geneDesc = words[9]
            gene_info.cdsStart = words[10]
            gene_info.cdsEnd = words[11]
            gene_info.exonStarts = words[12]
            gene_info.exonEnds = words[13]
        except:
            print("ERROR: line:", line, file=sys.stderr)
            continue

        # for some reason, exon starts and ends have trailing commas
        gene_info.exonStartParts = gene_info.exonStarts.strip(",").split(",")
        gene_info.exonEndParts = gene_info.exonEnds.strip(",").split(",")
        gene_info.exonUnions = set(
            [
                (int(s), int(e))
                for (s, e) in zip(gene_info.exonStartParts, gene_info.exonEndParts)
            ]
        )

        # add this gene info by checking whether it overlaps with any existing ones
        gene_infos = merge_gene_info(gene_infos, gene_info)

    for gene_id in gene_infos:
        for contig in gene_infos[gene_id]:

            # rows = df[df["transcript_id"].str.contains(gene_id)] contig.geneName
            rows = df[df["transcript_name"] == contig.geneName]
            start_codon = rows[rows["feature"] == "start_codon"]
            stop_codon = rows[rows["feature"] == "stop_codon"]
            cds = rows[rows["feature"] == "CDS"]

            #cds_cur = cds[cds["transcript_id"].str.contains(gene_id)]
            cds_cur = cds

            if len(cds_cur.index) == 0:
                start_codon_start = "."
                stop_codon_start = "."
                protein_codng = False
            else: 
                protein_codng = True
                #start_codon = start_codons[start_codons["transcript_id"].str.contains(gene_id)]

                if len(start_codon.index) == 0:
                    cds_min = cds_cur[cds_cur.exon_number == cds_cur.exon_number.min()]
                    start_codon_start = cds_min.iloc[0]['start']
                    #print("No start codon in gencode file for ID ", gene_id)
                else:
                    start_codon_start = start_codon.iloc[0]['start']
                

                #stop_codon = stop_codons[stop_codons["transcript_id"].str.contains(gene_id)]
                

                if len(stop_codon.index) == 0:
                    cds_max = cds_cur[cds_cur.exon_number == cds_cur.exon_number.max()]
                    stop_codon_start = cds_max.iloc[0]['start']
                    #print("No stop codon in gencode file for ID ", gene_id)
                else:
                    stop_codon_start = stop_codon.iloc[0]['start']



            output = "\t".join(
                map(
                    str,
                    [
                        contig.chrName,
                        contig.txStart,
                        contig.txEnd,
                        contig.geneName,
                        contig.score,
                        contig.strand,
                        "union_" + gene_id,
                        gene_id,
                        "protein_coding" if protein_codng else "non_coding",
                        contig.geneDesc,
                        contig.cdsStart,
                        contig.cdsEnd,
                        ",".join([str(e[0]) for e in sorted(contig.exonUnions)]),
                        ",".join([str(e[1]) for e in sorted(contig.exonUnions)]),
                        start_codon_start,
                        stop_codon_start,
                    ],
                )
            )
            print(output)


if __name__ == "__main__":
    main()
