import argparse
import os
import shutil


def parser_args():
    parser = argparse.ArgumentParser(description='This script represents isoforms from the GTF annotation as Start -> Splice sites -> End and keeps only those isoforms that are present in at least 2 repetitions for IsoQuant and StringTie results')

    # Add arguments
    parser.add_argument('--inp_dir', required=True, 
                        help='Path to the input directory containing annotations in GTF format')
    parser.add_argument('--ouf_path', required=True, 
                        help='Path to the output reconstructed isoforms file')

    # Parse the arguments
    return parser.parse_args()


def get_iso_struct(inf_path, ouf_path):
    introns_s = set()
    prev_start_end = tuple()
    prev_exon_end = False

    with open(ouf_path, 'w') as ouf:
        ouf.write('')

    with open(inf_path) as inf:
        for line in inf:
            if line[0] == '#':
                continue
            line_l = line.strip('\n').split('\t')
            feat_type = line_l[2]
            if feat_type == 'transcript':
                if prev_start_end != tuple():
                    with open(ouf_path, 'a') as ouf:
                        introns_s = list(introns_s)
                        chrom, orientation, _, _ = prev_start_end
                        ready_iso = []
                        for elem in introns_s:
                            _, _, start, stop = elem
                            ready_iso.append((start, stop))
                        ready_iso = sorted(list(ready_iso))
                        ready_iso = ';'.join(['%s-%s' % pos for pos in ready_iso])
                        ouf.write(f'{chrom}\t{orientation}\t{prev_start_end[2]}\t{ready_iso}\t{prev_start_end[3]}\n')
                chrom_transcript = line_l[0]
                orientation_transcript = line_l[6]
                prev_start_end = (chrom_transcript, orientation_transcript, line_l[3], line_l[4])

                prev_exon_end = False
                introns_s = set()
            if feat_type == 'exon':
                chrom = line_l[0]
                orientation = line_l[6]
                start = int(line_l[3])
                end = int(line_l[4])
                if orientation == '+':
                    if not prev_exon_end:
                        prev_exon_end = end
                    else:
                        introns_s.add((chrom, orientation, prev_exon_end + 1, start - 1))
                        prev_exon_end = end
                elif orientation == '-':
                    if not prev_exon_end:
                        prev_exon_end = start
                    else:
                        introns_s.add((chrom, orientation, end + 1, prev_exon_end - 1))
                        prev_exon_end = start

    with open(ouf_path, 'a') as ouf:
        introns_s = list(introns_s)
        chrom, orientation, _, _ = prev_start_end
        ready_iso = []
        for elem in introns_s:
            _, _, start, stop = elem
            ready_iso.append((start, stop))
        ready_iso = sorted(list(ready_iso))
        ready_iso = ';'.join(['%s-%s' % pos for pos in ready_iso])
        ouf.write(f'{chrom}\t{orientation}\t{prev_start_end[2]}\t{ready_iso}\t{prev_start_end[3]}\n')
        
        
def comb_by_splice_sites(fpath, d_of_iso, idx):
    with open(fpath) as inf:
        for line in inf:
            line = line.replace(' ', '\t')
            iso_splice_sites_only = line
            
            if iso_splice_sites_only not in d_of_iso.keys():
                d_of_iso[iso_splice_sites_only] = [0, 0, 0]
            d_of_iso[iso_splice_sites_only][idx] += 1
    
    return d_of_iso


def keep_iso_in_2reps_min(d_of_iso, ouf_path):
    with open(ouf_path, 'w') as ouf:
        for key, val in d_of_iso.items():
            if sum([x >= 1 for x in val]) > 1:
                ouf.write(key)

                
def main(inp_dir, ouf_path):
    tmp_dir = './tmp/'
    os.mkdir(tmp_dir)
    
    d_of_iso = {}
    
    for idx, file in enumerate(sorted(os.listdir(inp_dir))):
        inf_path = os.path.join(inp_dir, file)
        tmp_fpath = os.path.join(tmp_dir, file.replace('.gtf', '.tsv'))
        get_iso_struct(inf_path, tmp_fpath)
        
        d_of_iso = comb_by_splice_sites(tmp_fpath, d_of_iso, idx)
        
    
    keep_iso_in_2reps_min(d_of_iso, ouf_path)
    shutil.rmtree(tmp_dir)


if __name__ == '__main__':
    args = parser_args()
    main(args.inp_dir, args.ouf_path)
