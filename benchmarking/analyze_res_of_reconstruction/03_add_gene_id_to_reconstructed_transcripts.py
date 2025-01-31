import argparse
import os
import re

from collections import defaultdict


def parser_args():
    parser = argparse.ArgumentParser(description='This script adds gene IDs to transcripts based on the provided GTF annotation')

    parser.add_argument('--inp_dir', required=True, 
                        help='Path to the input directory containing reconstructed isoforms')
    parser.add_argument('--annot_fpath', required=True, 
                        help='Path to the reference GTF annotation file')
    parser.add_argument('--out_dir', required=True, 
                        help='Path to the output directory')

    return parser.parse_args()


def read_annot_file(annot_file):
    d_of_annot_genes = defaultdict(dict)
    with open(annot_file) as inf:
        for line in inf:
            if line[0] == '#':
                continue
            line_l = line.strip('\n').split('\t')
            feat_id = line_l[2]
            if feat_id == 'gene':
                chrom_and_orientation = f'{line_l[0]}*{line_l[6]}'
                start = int(line_l[3])
                stop = int(line_l[4])
                gene_id = re.findall(r'gene_id "(.+?)"', line_l[8])[0]
                d_of_annot_genes[chrom_and_orientation][(start, stop)] = gene_id
                
    return d_of_annot_genes


def intersection_with_genes(gene_coords, d_of_ref):
    start_gene, stop_gene = gene_coords
    best_gene = 'unassigned_gene'
    best_prop_intersection = 0
    
    for key, gene_id in d_of_ref.items():
        start_ref, stop_ref = key

        max_start = max(start_ref, start_gene)
        min_end = min(stop_ref, stop_gene)
        min_length = min(stop_ref - start_ref,
                         stop_gene - start_gene)
        prop_intersection = max(min_end - max_start, 0) / min_length
        
        if prop_intersection >= 0.1:
            if best_prop_intersection < prop_intersection:
                best_gene = gene_id
                best_prop_intersection = prop_intersection
            
            if prop_intersection >= 0.8:
                return gene_id
    return best_gene


def add_gene_ids(inf_path, d_of_annot_genes, ouf_path):
    n_uniassigned_genes = 0
    
    with open(inf_path) as inf, open(ouf_path, 'w') as ouf:
        for line in inf:
            line_l = line.strip('\n').split('\t')
            line_l = line_l[:5]
            chrom_and_orientation = f'{line_l[0]}*{line_l[1]}'
            if '-' in line_l[2]:
                start = int(line_l[2].split('-')[0])
                stop = int(line_l[4].split('-')[1])
            else:
                start = int(line_l[2])
                stop = int(line_l[4])
                
            gene_id = intersection_with_genes((start, stop), d_of_annot_genes[chrom_and_orientation])
            line_l.append(gene_id)

            ouf.write('\t'.join(line_l) + '\n')

            
def main(inp_dir, annot_fpath, out_dir):
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    d_of_annot_genes = read_annot_file(annot_fpath)
    
    for file in os.listdir(inp_dir):
        add_gene_ids(os.path.join(inp_dir, file), d_of_annot_genes, os.path.join(out_dir, file))


if __name__ == '__main__':
    args = parser_args()
    main(args.inp_dir, args.annot_fpath, args.out_dir)
