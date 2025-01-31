import argparse
import numpy as np


def parser_args():
    parser = argparse.ArgumentParser(description='This script calculates statistics for FLIC isoforms')
    parser.add_argument('--inf_path', required=True, help='Path to the input file with isoforms structure')
    parser.add_argument('--ouf_path', required=True, help='Path to the output file')

    return parser.parse_args()


def prepare_introns_l(introns_str):
    correct_introns_l = []
    if introns_str == '':
        return correct_introns_l
    introns_l = introns_str.split(';')
    for elem in introns_l:
        elem = tuple(map(int, elem.split('-')))
        correct_introns_l.append(elem)
    return correct_introns_l


def get_exons_l(introns_l, start_transcript, end_transcript):
    exons_l = []
    next_exon_start = start_transcript
    for intron in introns_l:
        exons_l.append((next_exon_start, intron[0] - 1))
        next_exon_start = intron[1] + 1
    exons_l.append((next_exon_start, end_transcript))
    return exons_l


def calc_mean_introns_exons_len(introns_l, return_sum=False):
    l_introns_len = []
    if introns_l == []:
        return 'NA'
    for intron in introns_l:
        l_introns_len.append(intron[1] - intron[0] + 1)
    
    if return_sum:
        return round(np.mean(l_introns_len), 1), sum(l_introns_len)
    return round(np.mean(l_introns_len), 1)


def calc_stat(inf_path, ouf_path):
    with open(ouf_path, 'w') as ouf:
        ouf.write('#Chromosome\tstrand\tisoform_id\tintrons number\tstart len\tend len\tisoform len\tmean introns len\tmean exons len\tsum exons len\n')

    with open(inf_path) as iso:
        for line in iso:
            if line[0] == '#':
                continue
            line_l = line.strip('\n').split('\t')
            chrom = line_l[0]
            orientation = line_l[1]
            introns_l = prepare_introns_l(line_l[3])
            transcript_id = line_l[-1]

            if orientation == '+':
                start = tuple(map(int, line_l[2].split('-')))
                end = tuple(map(int, line_l[4].split('-')))
                exons_l = get_exons_l(introns_l, start[0], end[1])
                iso_len = end[1] - start[0] + 1
            elif orientation == '-':
                start = tuple(map(int, line_l[4].split('-')))
                end = tuple(map(int, line_l[2].split('-')))
                exons_l = get_exons_l(introns_l, end[0], start[1])
                iso_len = start[1] - end[0] + 1
            
            # statistics
            n_introns = len(introns_l)
            start_len = start[1] - start[0] + 1
            end_len = end[1] - end[0] + 1
            mean_introns_len = calc_mean_introns_exons_len(introns_l)
            mean_exons_len, sum_exons_len = list(map(str, calc_mean_introns_exons_len(exons_l, return_sum=True)))

            with open(ouf_path, 'a') as ouf:
                ouf.write(f'{chrom}\t{orientation}\t{transcript_id}\t{n_introns}\t{start_len}\t{end_len}\t{iso_len}\t{mean_introns_len}\t{mean_exons_len}\t{sum_exons_len}\n')


if __name__ == '__main__':
    args = parser_args()
    calc_stat(args.inf_path, args.ouf_path)
       