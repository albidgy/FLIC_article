import argparse

from collections import defaultdict


def parser_args():
    parser = argparse.ArgumentParser(description='This script calculates gene statistics from FLIC isoforms')
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


def read_final_iso_file(fpath):
    d_of_iso_gene_comb = {}
    d_sorted = defaultdict(dict)

    with open(fpath) as iso:
        for line in iso:
            if line[0] == '#':
                continue
            line_l = line.strip('\n').split('\t')
            *line_l, iso_id = line.strip('\n').split('\t')
            gene_id = iso_id.split('.')[0] + line_l[1]
            start = tuple(map(int, line_l[2].split('-')))
            end = tuple(map(int, line_l[4].split('-')))
            introns_l = prepare_introns_l(line_l[3])

            d_sorted[gene_id][iso_id] = (len(introns_l), end[1] - start[0])

            if gene_id not in d_of_iso_gene_comb.keys():
                d_of_iso_gene_comb[gene_id] = {}
            d_of_iso_gene_comb[gene_id][iso_id] = [start, introns_l, end]

    return d_of_iso_gene_comb, d_sorted


def is_intersected_introns(start_major, end_major, start_comp, end_comp):
    max_start = max(start_major, start_comp)
    min_end = min(end_major, end_comp)
    if min_end - max_start > 0:
        return True
    return False


def compare_introns(major_struct, compared_struct):
    d_introns_compare_for1_iso = {'aib5': 0, 'aib3': 0, 'introns_retention': 0,
                                  'exon_skip': 0, 'exon_extra': 0}

    major_introns = set(major_struct[1])
    compared_start = int(compared_struct[0][0])
    compared_introns = set(compared_struct[1])
    compared_end = int(compared_struct[2][1])
    diff_major_introns = major_introns - compared_introns
    diff_compared_introns = compared_introns - major_introns

    if compared_introns != set():
        min_start_comp_intron = min(map(lambda x: x[0], compared_introns))
        max_start_comp_intron = max(map(lambda x: x[1], compared_introns))
    else:
        min_start_comp_intron = compared_end
        max_start_comp_intron = compared_start
    s_intersected_compared_introns = set()
    d_n_intersected_major_introns = defaultdict(int)

    for major_intron in diff_major_introns:
        s_is_intersected = set()
        start_major, end_major = major_intron
        if end_major < min_start_comp_intron:
            if start_major > compared_start:
                d_introns_compare_for1_iso['introns_retention'] += 1
            continue
        if start_major > max_start_comp_intron:
            if end_major < compared_end:
                d_introns_compare_for1_iso['introns_retention'] += 1
            continue

        for compared_intron in diff_compared_introns:
            start_comp, end_comp = compared_intron
            is_intersected = is_intersected_introns(start_major, end_major, start_comp, end_comp)
            s_is_intersected.add(is_intersected)
            if is_intersected:
                d_n_intersected_major_introns[major_intron] += 1
                if compared_intron in s_intersected_compared_introns:
                    d_introns_compare_for1_iso['exon_skip'] += 1

                if start_major != start_comp:
                    d_introns_compare_for1_iso['aib5'] += 1
                if end_major != end_comp:
                    d_introns_compare_for1_iso['aib3'] += 1

                s_intersected_compared_introns.add(compared_intron)

        if True not in s_is_intersected:
            d_introns_compare_for1_iso['introns_retention'] += 1

    d_introns_compare_for1_iso['exon_extra'] = len([x for x in d_n_intersected_major_introns.values() if x > 1])
    d_introns_compare_for1_iso['aib5'] -= d_introns_compare_for1_iso['exon_skip'] + d_introns_compare_for1_iso['exon_extra']
    d_introns_compare_for1_iso['aib3'] -= d_introns_compare_for1_iso['exon_skip'] + d_introns_compare_for1_iso['exon_extra']

    return d_introns_compare_for1_iso


def main(inf_path, ouf_path):
    with open(ouf_path, 'w') as ouf:
        ouf.write('gene_id\tn_iso\tn_starts\tn_ends\texon_skip\texon_extra\talt_introns_5\talt_introns_3'
                  '\tintrons_retention\tn_exons\n')

    d_of_iso_gene_comb, d_sorted = read_final_iso_file(inf_path)
    for gene_id, isoforms_d in d_sorted.items():
        d_of_genes_stat = {'n_iso': len(isoforms_d), 'n_starts': 0, 'n_ends': 0, 'exon_skip': 0,
                           'exon_extra': 0, 'aib5': 0, 'aib3': 0, 'introns_retention': 0,
                           'n_exons_longest_iso': 0}

        major_isoid, *sorted_isoids = list({k: v for k, v in sorted(isoforms_d.items(),
                                                                  key=lambda item: item[1], reverse=True)}.keys())
        d_of_genes_stat['n_exons_longest_iso'] = len(d_of_iso_gene_comb[gene_id][major_isoid][1]) + 1
        strand = gene_id[-1]
        if strand == '+':
            d_of_genes_stat['n_starts'] = len(set(isoform[0] for isoform in d_of_iso_gene_comb[gene_id].values()))
            d_of_genes_stat['n_ends'] = len(set(isoform[-1] for isoform in d_of_iso_gene_comb[gene_id].values()))
        else:
            d_of_genes_stat['n_starts'] = len(set(isoform[-1] for isoform in d_of_iso_gene_comb[gene_id].values()))
            d_of_genes_stat['n_ends'] = len(set(isoform[0] for isoform in d_of_iso_gene_comb[gene_id].values()))

        for compared_isoid in sorted_isoids:
            d_introns_compare_for1_iso = compare_introns(d_of_iso_gene_comb[gene_id][major_isoid],
                                                         d_of_iso_gene_comb[gene_id][compared_isoid])

            d_of_genes_stat['exon_skip'] += d_introns_compare_for1_iso['exon_skip']
            d_of_genes_stat['introns_retention'] += d_introns_compare_for1_iso['introns_retention']
            d_of_genes_stat['exon_extra'] += d_introns_compare_for1_iso['exon_extra']
            if strand == '+':
                d_of_genes_stat['aib5'] += d_introns_compare_for1_iso['aib5']
                d_of_genes_stat['aib3'] += d_introns_compare_for1_iso['aib3']
            else:
                d_of_genes_stat['aib5'] += d_introns_compare_for1_iso['aib3']
                d_of_genes_stat['aib3'] += d_introns_compare_for1_iso['aib5']

        with open(ouf_path, 'a') as ouf:
            res_l = [gene_id[:-1], d_of_genes_stat['n_iso'],
                     d_of_genes_stat['n_starts'], d_of_genes_stat['n_ends'],
                     d_of_genes_stat['exon_skip'], d_of_genes_stat['exon_extra'],
                     d_of_genes_stat['aib5'], d_of_genes_stat['aib3'],
                     d_of_genes_stat['introns_retention'], d_of_genes_stat['n_exons_longest_iso']]
            res_str = '\t'.join(list(map(str, res_l)))
            ouf.write(res_str + '\n')


if __name__ == '__main__':
    args = parser_args()
    main(args.inf_path, args.ouf_path)
    # inf_path = '../flic_final_res/RES_filtered/filt_isoforms.tsv'
    # ouf_path = './statistics/FINAL_genes_stat.tsv'
    # main(inf_path, ouf_path)
