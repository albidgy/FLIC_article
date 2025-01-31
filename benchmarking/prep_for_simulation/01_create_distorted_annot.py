import argparse
import re
import random


def parser_args():
    parser = argparse.ArgumentParser(description='This script creates a distorted annotation')
    parser.add_argument('--inf_path', required=True, 
                        help='Path to the input annotation file')
    parser.add_argument('--transcript_modes_ouf_path', required=True, 
                        help='Path to the source transcripts mode file')
    parser.add_argument('--ouf_path', required=True, 
                        help='Path to the output distorted annotation file')

    return parser.parse_args()


def create_bins_for_transcripts(inf_path, ouf_path):
    '''
    0 - control
    1 - 5' shorter
    2 - 3' shorter
    3 - 5' longer
    4 - 3' longer
    5 - 5' and 3' shorter
    6 - 5' and 3' longer
    '''
    d_transcripts_mode = {}

    with open(inf_path) as inf:
        for line in inf:
            line_l = line.strip('\n').split('\t')
            feat_type = line_l[2]
            if feat_type == 'transcript':
                transcript_id = re.findall(r'transcript_id "(.+?)"', line_l[8])[0]
                d_transcripts_mode[transcript_id] = random.randint(0, 6)
    
    with open(ouf_path, 'w') as ouf:
        for transcript_id in d_transcripts_mode:
            mode = str(d_transcripts_mode[transcript_id])
            ouf.write(f'{transcript_id}\t{mode}\n')

    return d_transcripts_mode


def modify_transcripts(start, stop, orientation, type_mod):
    modifications_by_groups = {0: (0, 0), 1: (100, 0), 2: (0, -100), 3: (-100, 0), 
                               4: (0, 100), 5: (100, -100), 6: (-100, 100)}
    
    start_add, stop_add = modifications_by_groups[type_mod]
    if orientation == '+':
        start_mod = start + start_add
        stop_mod = stop + stop_add
    else:
        start_mod = start - stop_add
        stop_mod = stop - start_add
    
    if stop_mod - start_mod > 0:
        return start_mod, stop_mod
    else:
        return start, stop
    
    
def modify_exons(transcript_borders, exon_start, exon_end):
    transcript_start, transcript_end = transcript_borders
    if abs(transcript_start - exon_start) == 100:
        exon_start = transcript_start
    if abs(transcript_end - exon_end) == 100:
        exon_end = transcript_end
    
    return exon_start, exon_end


def change_borders_r1(inf_path, d_transcripts_mode):
    d_new_genes_coords = {}
    d_new_transcript_coords = {}
    d_exons = {}

    with open(inf_path) as inf:
        for line in inf:
            line_l = line.strip('\n').split('\t')
            feat_type = line_l[2]
            orientation = line_l[6]
            start = int(line_l[3])
            stop = int(line_l[4])

            if feat_type == 'transcript':
                gene_id = re.findall(r'gene_id "(.+?)"', line_l[8])[0]
                if gene_id not in d_new_genes_coords.keys():
                    d_new_genes_coords[gene_id] = [float('inf'), float('-inf')]

                transcript_id = re.findall(r'transcript_id "(.+?)"', line_l[8])[0]
                start, stop = modify_transcripts(start, stop, orientation, 
                                                          d_transcripts_mode[transcript_id])
                d_new_transcript_coords[transcript_id] = [start, stop]

                d_new_genes_coords[gene_id][0] = min(d_new_genes_coords[gene_id][0], start)
                d_new_genes_coords[gene_id][1] = max(d_new_genes_coords[gene_id][1], stop)

                d_exons[transcript_id] = []

            elif feat_type == 'exon':
                transcript_id = re.findall(r'transcript_id "(.+?)"', line_l[8])[0]
                d_exons[transcript_id].append(line_l)
    
    d_exons = {k: sorted(v, key=lambda x: x[3]) for k, v in d_exons.items()}
    return d_new_genes_coords, d_new_transcript_coords, d_exons


def correct_exons_borders(exons_l, start_transcript, stop_transcript):
    is_find_cor_start = False
    is_find_cor_stop = True
    new_exons_l = []

    for cur_exon in exons_l:
        if int(cur_exon[4]) > start_transcript and int(cur_exon[3]) < stop_transcript:
            new_exons_l.append(cur_exon)
    
    if new_exons_l == []:
        new_exons_l.append(cur_exon)
    
    new_exons_l.sort(key=lambda x: int(x[3]))
    new_exons_l[0][3] = start_transcript
    new_exons_l[-1][4] = stop_transcript
    return new_exons_l


def get_correct_d_exons(d_exons, d_new_transcript_coords):
    d_corr_exons = {}
    for transcript, exons_l in d_exons.items():
        start_transcript, stop_transcript = d_new_transcript_coords[transcript]
        d_corr_exons[transcript] = correct_exons_borders(exons_l, start_transcript, stop_transcript)
    
    return d_corr_exons


def write_final_res(inf_path, ouf_path, d_new_genes_coords, d_new_transcript_coords, d_corr_exons):
    with open(inf_path) as inf, open(ouf_path, 'w') as ouf:
        for line in inf:
            line_l = line.strip('\n').split('\t')
            feat_type = line_l[2]
            if feat_type == 'gene':
                gene_id = re.findall(r'gene_id "(.+?)"', line_l[8])[0]
                line_l[3] = str(d_new_genes_coords[gene_id][0])
                line_l[4] = str(d_new_genes_coords[gene_id][1])

                with open(ouf_path, 'a') as ouf:
                    ouf.write('\t'.join(line_l) + '\n')
            elif feat_type == 'transcript':
                transcript_id = re.findall(r'transcript_id "(.+?)"', line_l[8])[0]
                line_l[3] = str(d_new_transcript_coords[transcript_id][0])
                line_l[4] = str(d_new_transcript_coords[transcript_id][1])

                exons_l = d_corr_exons[transcript_id]
                with open(ouf_path, 'a') as ouf:
                    ouf.write('\t'.join(line_l) + '\n')
                    for cur_exon in exons_l:
                        ouf.write('\t'.join(list(map(str, cur_exon))) + '\n')


def main(inf_path, transcript_modes, ouf_path):
    d_transcripts_mode = create_bins_for_transcripts(inf_path, transcript_modes)
    d_new_genes_coords, d_new_transcript_coords, d_exons = change_borders_r1(inf_path, d_transcripts_mode)
    d_corr_exons = get_correct_d_exons(d_exons, d_new_transcript_coords)
    write_final_res(inf_path, ouf_path, d_new_genes_coords, d_new_transcript_coords, d_corr_exons)


if __name__ == '__main__':
    args = parser_args()
    main(args.inf_path, args.transcript_modes_ouf_path, args.ouf_path)
