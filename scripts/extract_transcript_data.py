from gtfparse import read_gtf
import gzip
import csv
import random

# Input/Output file names (need to be in same folder)
# Gencode file
gencode_file = 'gencode.v29.annotation.gtf'
# Chromosome file
chr_file = 'hg38_full.txt'
# Output file
output_file = 'extracted_transcripts_20200814.txt'


df_orig = read_gtf(gencode_file)
#df = df_orig.head(20000)
df = df_orig
total_entries = len(df)
print("Length of dataframe: ", total_entries)

data = {}

# We use random ints for the importance column. Could be replaces with publication counts
# pub_count = {}
# with gzip.open('gencode.v29.metadata.Pubmed_id.gz', 'rb') as f:
#     for item in f:
#         line = item.decode().strip().split('\t')
#         if line[0] not in pub_count:
#             pub_count[line[0]] = 0
#         pub_count[line[0]] = pub_count[line[0]] + 1

with open(chr_file, 'r') as opf:
    chrsizes = opf.readlines()

chrms = [i.split('\t')[0] for i in chrsizes]

for i, v in df.iterrows():
    if i % 5000 == 0:
        print("Progress: ", str(i),"/", str(total_entries))
    if v['feature'] == 'transcript' or v['feature'] == 'exon' or v['feature'] == 'CDS' or v['feature'] == 'start_codon' or v['feature'] == 'stop_codon':
        if v['gene_id'] not in data:
            data[v['gene_id']] = {}

        if v['transcript_id'] not in data[v['gene_id']]:
            data[v['gene_id']][v['transcript_id']] = {
                'chr': '', 
                'start': '', 
                'end': '', 
                'transcript_name': '', 
                'citationCount': 1,  
                'strand': '', 
                'transcript_id': v['transcript_id'], 
                'gene_id': v['gene_id'], 
                'gene_type': '', 
                'CDSStarts': '', 
                'CDSEnds': '', 
                'ExonStarts': '', 
                'ExonEnds': '',
                'StartCodonStart': '.',
                'StopCodonStart': '.'}

        if v['feature'] == 'transcript':
            data[v['gene_id']][v['transcript_id']]['chr'] = v['seqname']
            data[v['gene_id']][v['transcript_id']]['start'] = v['start']
            data[v['gene_id']][v['transcript_id']]['end'] = v['end']
            data[v['gene_id']][v['transcript_id']]['transcript_name'] = v['transcript_name']
            data[v['gene_id']][v['transcript_id']]['strand'] = v['strand']
            data[v['gene_id']][v['transcript_id']]['gene_type'] = v['gene_type']

        if v['feature'] == 'exon':
            data[v['gene_id']][v['transcript_id']]['ExonStarts'] += str(v['start']) + ','
            data[v['gene_id']][v['transcript_id']]['ExonEnds'] += str(v['end']) + ','

        if v['feature'] == 'CDS':
            data[v['gene_id']][v['transcript_id']]['CDSStarts'] += str(v['start']) + ','
            data[v['gene_id']][v['transcript_id']]['CDSEnds'] += str(v['end']) + ','
        
        if v['feature'] == 'start_codon':
            data[v['gene_id']][v['transcript_id']]['StartCodonStart'] = str(v['start']) 

        if v['feature'] == 'stop_codon':
            data[v['gene_id']][v['transcript_id']]['StopCodonStart'] = str(v['start']) 
            
data = dict(sorted(data.items()))

for gene_id, transcripts in data.items():
    # Each transcript of the same gene gets the same importance value (could be changed)
    importance = random.randint(1,100)

    for transcript_id, info in transcripts.items():
        if info['gene_type'] == 'protein_coding' or info['gene_type'] == 'miRNA':
            if info['chr'] in chrms:
                exons_start = info['CDSStarts'].split(',')[:-1]
                exons_start_formated = [int(i) for i in exons_start]
                cds_start = min(exons_start_formated) if len(exons_start_formated)>0 else "." 
                exons_end = info['CDSEnds'].split(',')[:-1]
                exons_end_formated = [int(i) for i in exons_end]
                cds_end = max(exons_end_formated) if len(exons_end_formated)>0 else "."

                if info['gene_type'] == 'miRNA':
                    info['StartCodonStart'] = '.'
                    info['StopCodonStart'] = '.'
                
                if info['StartCodonStart'] != '.' and info['StopCodonStart'] == '.':
                    info['StopCodonStart'] = cds_end

                # if it is protein coding but we don't know the start codon, display it as non coding
                if info['StartCodonStart'] == '.' and info['StopCodonStart'] != '.':
                    info['StartCodonStart'] = '.'
                    info['StopCodonStart'] = '.'

                exons_start = info['ExonStarts'].split(',')[:-1]
                exons_start_formated = sorted([int(i) for i in exons_start])
                info['ExonStarts'] = ",".join([str(e) for e in exons_start_formated])

                exons_end = info['ExonEnds'].split(',')[:-1]
                exons_end_formated = sorted([int(i) for i in exons_end])
                info['ExonEnds'] = ",".join([str(e) for e in exons_end_formated])

                info['citationCount'] = importance
                # if transcript_id in pub_count:
                #     info['citationCount'] = pub_count[transcript_id]

# remove unnecessary fields
data_clean = {}
output = []
for gene_id, transcripts in data.items():
    for transcript_id, info in transcripts.items():
        
        if info['gene_type'] == 'protein_coding' or info['gene_type'] == 'miRNA':
            #print("TRANSCRIPT_ITEM",transcript_id, info)
            if info['chr'] in chrms:
                if info['gene_id'] not in data_clean:
                    data_clean[info['gene_id']] = {}

                if info['transcript_id'] not in data_clean[info['gene_id']]:
                    data_clean[info['gene_id']][info['transcript_id']] = {
                        'chr': info['chr'], 
                        'start': info['start'], 
                        'end': info['end'], 
                        'transcript_name': info['transcript_name'], 
                        'citationCount': info['citationCount'],  
                        'strand': info['strand'], 
                        'gene_id': info['gene_id'], 
                        'transcript_id': info['transcript_id'], 
                        'gene_type': info['gene_type'], 
                        'ExonStarts': info['ExonStarts'], 
                        'ExonEnds': info['ExonEnds'],
                        'StartCodonStart': info['StartCodonStart'],
                        'StopCodonStart': info['StopCodonStart']}

                    output.append(data_clean[info['gene_id']][info['transcript_id']])

headers = output[0].keys()
with open(output_file, 'w') as opf:
    myWriter = csv.DictWriter(opf, delimiter='\t', fieldnames=headers)
    for row in output:
        myWriter.writerow(row)