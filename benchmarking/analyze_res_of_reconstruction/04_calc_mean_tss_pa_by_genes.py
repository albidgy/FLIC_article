import argparse
import re

from collections import defaultdict


def parser_args():
    parser = argparse.ArgumentParser(description='This script calculates peak width by genes based on provided annotation and CAGE data')
    parser.add_argument('--ref_annot_fpath', required=True, 
                        help='Path to the reference GTF annotation file')
    parser.add_argument('--cagefightr_tss', required=True, 
                        help='Path to the CAGEfigthR TSSs')
    parser.add_argument('--cagefightr_pa', required=True, 
                        help='Path to the CAGEfigthR PAs')
    parser.add_argument('--ouf_path', required=True, 
                        help='Path to the output file for peak width by genes')

    return parser.parse_args()


def read_annot_file(ref_annot_fpath):
    d_gene_borders = defaultdict(dict)
    with open(ref_annot_fpath) as inf:
        for line in inf:
            if line[0] == '#':
                continue
            
            line_l = line.strip('\n').split('\t')
            if line_l[2] == 'gene':
                chrom_and_orientation = f'{line_l[0]}*{line_l[6]}'
                start = int(line_l[3])
                stop = int(line_l[4])
                gene_id = re.findall(r'gene_id "(.+?)"', line_l[8])[0]
                d_gene_borders[chrom_and_orientation][(start, stop)] = gene_id
    
    return d_gene_borders


def intersection_with_genes(gene_coords, d_of_ref):
    start_gene, stop_gene = gene_coords
    
    for key, gene_id in d_of_ref.items():
        start_ref, stop_ref = key

        max_start = max(start_ref, start_gene)
        min_end = min(stop_ref, stop_gene)
        prop_intersection = min_end - max_start
        
        if prop_intersection > 0:
            return gene_id   
    return 'unassigned_gene'


def find_mean_peak_width_by_gene(cage_fpath, d_gene_borders):
    d_peaks_len_by_genes = {}
    n_unassigned_peaks = 0

    with open(cage_fpath) as inf:
        for line in inf:
            line_l = line.strip('\n').split('\t')
            chrom_and_orientation = f'{line_l[0]}*{line_l[5]}'
            peak_start = int(line_l[1])
            peak_stop = int(line_l[2])
            gene_id = intersection_with_genes((peak_start, peak_stop), d_gene_borders[chrom_and_orientation])

            if gene_id == 'unassigned_gene':
                n_unassigned_peaks += 1
            else:
                if gene_id not in d_peaks_len_by_genes.keys():
                    d_peaks_len_by_genes[gene_id] = []
                d_peaks_len_by_genes[gene_id].append(peak_stop - peak_start)

    return {key: max(values) for key, values in d_peaks_len_by_genes.items()}


def main(ref_annot_fpath, cage_start_fpath, cage_pa_fpath, ouf_path):
    d_gene_borders = read_annot_file(ref_annot_fpath)
    d_max_tss_by_genes = find_mean_peak_width_by_gene(cage_start_fpath, d_gene_borders)
    d_max_pa_by_genes = find_mean_peak_width_by_gene(cage_pa_fpath, d_gene_borders)
    common_genes = set(d_max_tss_by_genes.keys()) & set(d_max_pa_by_genes.keys())

    with open(ouf_path, 'w') as ouf:
        ouf.write('#gene_id\tTSS_width\tPA_width\n')
        for gene_id in common_genes:
            ouf.write(f'{gene_id}\t{d_max_tss_by_genes[gene_id]}\t{d_max_pa_by_genes[gene_id]}\n')


if __name__ == '__main__':
    args = parser_args()
    main(args.ref_annot_fpath, args.cagefightr_tss, args.cagefightr_pa, args.ouf_path)
