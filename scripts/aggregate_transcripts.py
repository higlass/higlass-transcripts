import gzip
import csv
import random
import negspy.coordinates as nc
import sqlite3
import numpy as np
import os
import os.path as op
import slugid
import math
import collections as col
import json

filepath = "extracted_transcripts_20200814.txt"
#filepath = "gene_table_v2_transcripts_names_new.txt"
output_file = "transcripts_20200814.beddb"
#output_file = "transcripts_20200723_3.beddb"
importance_column = 5
has_header = False
chromosome = None
max_transcripts_per_tile = 5
tile_size = 1024
delimiter = '\t'
chromsizes_filename = 'hg38_full.txt'
offset = 0


def load_chromsizes(chromsizes_filename, assembly=None):
    """
    Load a set of chromosomes from a file or using an assembly
    identifier. If using just an assembly identifier the chromsizes
    will be loaded from the negspy repository.

    Parameters:
    -----------
    chromsizes_filename: string
        The file containing the tab-delimited chromosome sizes
    assembly: string
        Assembly name (e.g. 'hg19'). Not necessary if a chromsizes_filename is passed in
    """
    if chromsizes_filename is not None:
        chrom_info = nc.get_chrominfo_from_file(chromsizes_filename)
        chrom_names = chrom_info.chrom_order
        chrom_sizes = [chrom_info.chrom_lengths[c] for c in chrom_info.chrom_order]
    else:
        if assembly is None:
            raise ValueError("No assembly or chromsizes specified")

        chrom_info = nc.get_chrominfo(assembly)
        chrom_names = nc.get_chromorder(assembly)
        chrom_sizes = nc.get_chromsizes(assembly)

    return (chrom_info, chrom_names, chrom_sizes)

def store_meta_data(
    cursor,
    zoom_step,
    max_length,
    assembly,
    chrom_names,
    chrom_sizes,
    tile_size,
    max_zoom,
    max_width,
    version,
    header=[],
):
    cursor.execute(
        """
        CREATE TABLE tileset_info
        (
            zoom_step INT,
            max_length INT,
            assembly text,
            chrom_names text,
            chrom_sizes text,
            tile_size REAL,
            max_zoom INT,
            max_width REAL,
            header text,
            version text
        )
        """
    )

    cursor.execute(
        "INSERT INTO tileset_info VALUES (?,?,?,?,?,?,?,?,?,?)",
        (
            zoom_step,
            max_length,
            assembly,
            "\t".join(chrom_names),
            "\t".join(map(str, chrom_sizes)),
            tile_size,
            max_zoom,
            max_width,
            "\t".join(header),
            version,
        ),
    )

    cursor.commit()

    pass


def aggregate_bedfile(
    filepath,
    output_file,
    importance_column,
    has_header,
    chromosome,
    max_transcripts_per_tile,
    tile_size,
    delimiter,
    chromsizes_filename,
    offset,
):
    BEDDB_VERSION = 3

    assembly = None

    if output_file is None:
        output_file = filepath + ".beddb"
    else:
        output_file = output_file

    if op.exists(output_file):
        os.remove(output_file)

    if filepath.endswith(".gz"):
        import gzip

        bed_file = gzip.open(filepath, "rt")
    else:
        bed_file = open(filepath, "r")

    try:
        (chrom_info, chrom_names, chrom_sizes) = load_chromsizes(
            chromsizes_filename, None
        )
    except FileNotFoundError:
        print(
                "Chromsizes filename not found:", chromsizes_filename, file=sys.stderr
            )
        return None

    rand = random.Random(3)

    
    def line_to_np_array(line):
        """
        Convert a bed file line to a numpy array which can later
        be used as an entry in an h5py file.
        """
        try:
            start = int(line[1])
            stop = int(line[2])
        except ValueError:
            raise ValueError("Error parsing the position, line: {}".format(line))

        chrom = line[0]
        gene_id = line[6]

        if importance_column is None:
            # assume a random importance when no aggregation strategy is given
            importance = rand.random()
        elif importance_column == "size":
            importance = stop - start
        elif importance_column == "random":
            importance = rand.random()
        else:
            importance = float(line[int(importance_column) - 1])

        if stop < start:
            print("WARNING: stop < start:", line, file=sys.stderr)

            start, stop = stop, start

        if len(line) > 3:
            bedline_name = line[3]
        else:
            bedline_name = ""
        # convert chromosome coordinates to genome coordinates
        genome_start = chrom_info.cum_chrom_lengths[chrom] + start + offset
        genome_end = chrom_info.cum_chrom_lengths[chrom] + stop + offset

        pos_offset = genome_start - start
        parts = {
            "startPos": genome_start,
            "endPos": genome_end,
            "uid": slugid.nice(),
            "name": bedline_name,
            "chrOffset": pos_offset,
            "geneId": gene_id,
            "fields": "\t".join(line),
            "importance": importance,
            "chromosome": str(chrom),
        }

        return parts

    dset = []

    print("delimiter:", delimiter)
    if has_header:
        line = bed_file.readline()
        header = line.strip().split(delimiter)
    else:
        line = bed_file.readline().strip()
        line_parts = line.strip().split(delimiter)
        try:
            dset += [line_to_np_array(line_parts)]
        except KeyError:
            print(
                f"Unable to find {line_parts[0]} in the list of chromosome sizes. "
                "Please make sure the correct assembly or chromsizes filename "
                "is passed in as a parameter",
                file=sys.stderr,
            )
            return None
        except IndexError:
            print("Invalid line:", line)
        header = map(str, list(range(1, len(line.strip().split(delimiter)) + 1)))

    for line in bed_file:
        line_parts = line.strip().split(delimiter)
        try:
            dset += [line_to_np_array(line_parts)]
        except IndexError:
            print("Invalid line:", line)

    if chromosome is not None:
        dset = [d for d in dset if d["chromosome"] == chromosome]

    # We neeed chromosome information as well as the assembly size to properly
    # tile this data
    tile_size = tile_size

    assembly_size = chrom_info.total_length + 1
    """
    else:
        try:
            assembly_size = chrom_info.chrom_lengths[chromosome]
        except KeyError:
            print(
                "ERROR: Chromosome {} not found in assembly {}.".format(
                    chromosome, assembly
                ),
                file=sys.stderr
            )
            return 1
    """

    max_zoom = int(math.ceil(math.log(assembly_size / tile_size) / math.log(2)))
    """
    if max_zoom is not None and max_zoom < max_zoom:
        max_zoom = max_zoom
    """

    # this script stores data in a sqlite database
    import sqlite3

    sqlite3.register_adapter(np.int64, lambda val: int(val))
    print("output_file:", output_file, "header:", header)
    conn = sqlite3.connect(output_file)

    # store some meta data
    store_meta_data(
        conn,
        1,
        max_length=assembly_size,
        assembly=assembly,
        chrom_names=chrom_names,
        chrom_sizes=chrom_sizes,
        tile_size=tile_size,
        max_zoom=max_zoom,
        max_width=tile_size * 2 ** max_zoom,
        header=header,
        version=BEDDB_VERSION,
    )

    # max_width = tile_size * 2 ** max_zoom
    uid_to_entry = {}


    gene_intervals = dict()

    # store each bed file entry as an interval
    #print("dset", json.dumps(dset[:10], indent = 4))
    for d in dset:
        uid = d["uid"]
        uid_to_entry[uid] = d

        if d["geneId"] in gene_intervals:
            gene_intervals[d["geneId"]] += [(d["startPos"], d["endPos"], d["importance"], uid)]
        else:
            gene_intervals[d["geneId"]] = [(d["startPos"], d["endPos"], d["importance"], uid)]


    #print('gene_intervals:',json.dumps(gene_intervals, indent = 4))

    gene_intervals_min_max = dict()

    # We assume that each transcript is equally important (gene importance)
    for key, gene_interval in gene_intervals.items():
        
        for transcript_interval in gene_interval:
            importance = transcript_interval[2]
            if key in gene_intervals_min_max:
                cur_min = gene_intervals_min_max[key][0]
                cur_max = gene_intervals_min_max[key][1]
                importance = transcript_interval[2]
                gene_intervals_min_max[key] = [ min(transcript_interval[0], cur_min), max(transcript_interval[1], cur_max), importance, key]
            else:
                #print(gene_interval)
                gene_intervals_min_max[key] = [transcript_interval[0], transcript_interval[1], importance, key]


    #print('gene_intervals_min_max:',json.dumps(gene_intervals_min_max, indent = 4))

    gene_intervals_min_max_arr = [] 

    for key, gene_interval in gene_intervals_min_max.items():
        gene_intervals_min_max_arr += [(gene_interval[0], gene_interval[1], gene_interval[2], gene_interval[3])]


    #print('gene_intervals_min_max_arr:',json.dumps(gene_intervals_min_max_arr[:2], indent = 4))
    # for gene_id in gene_ids:
    #     print(gene_intervals["s"] == None)

    tile_width = tile_size

    c = conn.cursor()
    c.execute(
        """
        CREATE TABLE intervals
        (
            id int PRIMARY KEY,
            zoomLevel int,
            importance real,
            startPos int,
            endPos int,
            chrOffset int,
            uid text,
            name text,
            fields text
        )
        """
    )

    c.execute(
        """
        CREATE VIRTUAL TABLE position_index USING rtree(
            id,
            rStartZoomLevel, rEndZoomLevel, rStartPos, rEndPos
        )
        """
    )

    curr_zoom = 0
    counter = 0

    max_viewable_zoom = max_zoom

    if max_zoom is not None and max_zoom < max_zoom:
        max_viewable_zoom = max_zoom

    # sorted_intervals = sorted(
    #     intervals, key=lambda x: -uid_to_entry[x[-1]]["importance"]
    # )

    sorted_gene_intervals = sorted(
        gene_intervals_min_max_arr, key=lambda x: -x[-2]
    )

    

    #print('si:',json.dumps(sorted_intervals[:10], indent = 4))
    #print('si:',json.dumps(sorted_gene_intervals, indent = 4))
    print("max_transcripts_per_tile:", max_transcripts_per_tile)

    tile_counts = col.defaultdict(int)

    for gene_interval in sorted_gene_intervals:
        # go through each interval from most important to least
        while curr_zoom <= max_viewable_zoom:
            # try to place it in the highest zoom level and go down from there
            tile_width = tile_size * 2 ** (max_zoom - curr_zoom)

            curr_pos = gene_interval[0]
            space_available = True

            # check if there's space at this zoom level
            while curr_pos < gene_interval[1]:
                curr_tile = math.floor(curr_pos / tile_width)
                tile_id = "{}.{}".format(curr_zoom, curr_tile)


                # print(tile_id, "tile_counts[tile_id]", tile_counts[tile_id])
                if tile_counts[tile_id] >= max_transcripts_per_tile:
                    space_available = False
                    break

                curr_pos += tile_width

            # if there is, then fill it up
            if space_available:
                curr_pos = gene_interval[0]
                while curr_pos < gene_interval[1]:
                    curr_tile = math.floor(curr_pos / tile_width)
                    tile_id = "{}.{}".format(curr_zoom, curr_tile)

                    # This counts the number of transcript blocks that we are adding
                    tile_counts[tile_id] += 1

                    """
                    # increment tile counts for lower level tiles
                    higher_zoom = curr_zoom + 1
                    higher_tile = math.floor(higher_zoom / 2)

                    while higher_zoom <= max_viewable_zoom:
                        new_tile_id = '{}.{}'.format(higher_zoom, higher_tile)
                        higher_zoom += 1
                        higher_tile = math.floor(higher_tile / 2)
                        tile_counts[new_tile_id] += 1
                    """

                    curr_pos += tile_width

            if space_available:
                # get all transcripts for that gene
                transcripts = gene_intervals[gene_interval[3]]

                for transcript in transcripts:
                    value = uid_to_entry[transcript[-1]]

                    # one extra question mark for the primary key
                    exec_statement = "INSERT INTO intervals VALUES (?,?,?,?,?,?,?,?,?)"

                    c.execute(
                        exec_statement,
                        # primary key, zoomLevel, startPos, endPos, chrOffset, line
                        (
                            counter,
                            curr_zoom,
                            value["importance"],
                            value["startPos"],
                            value["endPos"],
                            value["chrOffset"],
                            value["uid"],
                            value["name"],
                            value["fields"],
                        ),
                    )

                    if counter % 1000 == 0:
                        print("counter:", counter, value["endPos"] - value["startPos"])

                    exec_statement = "INSERT INTO position_index VALUES (?,?,?,?,?)"
                    c.execute(
                        exec_statement,
                        # add counter as a primary key
                        (counter, curr_zoom, curr_zoom, value["startPos"], value["endPos"]),
                    )

                    counter += 1
                break

            curr_zoom += 1

        curr_zoom = 0
    conn.commit()
    return True



aggregate_bedfile(
    filepath,
    output_file,
    importance_column,
    has_header,
    chromosome,
    max_transcripts_per_tile,
    tile_size,
    delimiter,
    chromsizes_filename,
    offset,
)