import argparse


def parser_args():
    parser = argparse.ArgumentParser(description='This script filters the reference isoforms, leaving only those isoforms for which peaks of TSSs and PAs were found')
    parser.add_argument('--inf_path', required=True, 
                        help='Path to the input file with reference isoforms structure')
    parser.add_argument('--peak_width', required=True, 
                        help='Path to the file with peak width by genes')
    parser.add_argument('--ouf_path', required=True, 
                        help='Path to the output file')

    return parser.parse_args()


def get_s_good_genes(peak_width_fpath):
    s_good_genes = set()
    with open(peak_width_fpath) as inf:
        for line in inf:
            line_l = line.strip('\n').split('\t')
            s_good_genes.add(line_l[0])
    return s_good_genes


def filt_ref_iso(inf_path, peak_width_fpath, ouf_path):
    s_good_genes = get_s_good_genes(peak_width_fpath)

    with open(inf_path) as inf, open(ouf_path, 'w') as ouf:
        for line in inf:
            line_l = line.strip('\n').split('\t')
            transcript_id = line_l[-1]
            gene_id = '.'.join(transcript_id.split('.')[:-1])

            if gene_id in s_good_genes:
                ouf.write(line)


if __name__ == '__main__':
    args = parser_args()
    filt_ref_iso(args.inf_path, args.peak_width, args.ouf_path)
