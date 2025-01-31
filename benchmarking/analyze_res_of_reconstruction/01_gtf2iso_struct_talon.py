import argparse


def parser_args():
    parser = argparse.ArgumentParser(description='This script processes the GTF annotation file and generates isoform structures for the reconstruction')
    parser.add_argument('--inf_path', required=True, 
                        help='Path to the input GTF file')
    parser.add_argument('--ouf_path', required=True, 
                        help='Path to the output reconstructed isoforms file')

    return parser.parse_args()


def get_iso_struct(inf_path, ouf_path):
    prev_exon_end = False
    introns_s = set()
    prev_start_end = tuple()

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


if __name__ == '__main__':
    args = parser_args()
    get_iso_struct(args.inf_path, args.ouf_path)
