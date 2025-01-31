import argparse
import os


def parser_args():
    parser = argparse.ArgumentParser(description='This script adds peak width information to isoform')
    parser.add_argument('--peak_width', required=True, 
                        help='Path to the file with peak width by genes')
    parser.add_argument('--inp_dir', required=True, 
                        help='Path to the input directory containing isoform structures')
    parser.add_argument('--out_dir', required=True, 
                        help='Path to the output directory')

    return parser.parse_args()


def read_peak_width(peak_width_fpath):
    d_peak_width_by_genes = {}
    with open(peak_width_fpath) as inf:
        for line in inf:
            if line[0] == '#':
                continue
            gene_id, tss, pa = line.strip('\n').split('\t')
            d_peak_width_by_genes[gene_id] = (int(round(int(tss) / 2, 0)), int(round(int(pa) / 2, 0)))
    
    return d_peak_width_by_genes


def create_peaks(inf_path, out_dir, d_peak_width_by_genes):
    ouf_path = os.path.join(out_dir, os.path.basename(inf_path))

    with open(inf_path) as inf, open(ouf_path, 'w') as ouf:
        for line in inf:
            line_l = line.strip('\n').split('\t')
            gene_id = line_l[-1]
            if gene_id in d_peak_width_by_genes.keys():
                if not '-' in line_l[2]:
                    strand = line_l[1]
                    if strand == '+':
                        start_start = int(line_l[2]) - d_peak_width_by_genes[gene_id][0]
                        start_end = int(line_l[2]) + d_peak_width_by_genes[gene_id][0]
                        end_start = int(line_l[4]) - d_peak_width_by_genes[gene_id][1]
                        end_end = int(line_l[4]) + d_peak_width_by_genes[gene_id][1]
                    else:
                        start_start = int(line_l[2]) - d_peak_width_by_genes[gene_id][1]
                        start_end = int(line_l[2]) + d_peak_width_by_genes[gene_id][1]
                        end_start = int(line_l[4]) - d_peak_width_by_genes[gene_id][0]
                        end_end = int(line_l[4]) + d_peak_width_by_genes[gene_id][0]

                    line_l[2] = f'{start_start}-{start_end}'
                    line_l[4] = f'{end_start}-{end_end}'
                ouf.write('\t'.join(line_l) + '\n')


def main(peak_width_fpath, inp_dir, out_dir):
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    
    d_peak_width_by_genes = read_peak_width(peak_width_fpath)
    for file in os.listdir(inp_dir):
        inf_path = os.path.join(inp_dir, file)
        if os.path.isfile(inf_path):
            create_peaks(inf_path, out_dir, d_peak_width_by_genes)


if __name__ == '__main__':
    args = parser_args()
    main(args.peak_width, args.inp_dir, args.out_dir)
