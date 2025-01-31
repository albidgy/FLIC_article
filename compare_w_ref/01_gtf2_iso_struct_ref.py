import argparse
import re


def create_transcript_struct(inf_path, ouf_path):
    prev_exon_end = False
    prev_transcript_start_stop_chrom_orientation = ()
    introns_s = set()
    prev_transcript_id = None

    with open(ouf_path, 'w') as ouf:
        ouf.write('')

    with open(inf_path) as inf:
        for line in inf:
            if line[0] == '#':
                continue
            line_l = line.strip('\n').split('\t')
            chrom = line_l[0]

            feat_type = line_l[2]
            if feat_type == 'transcript':
                if prev_transcript_id != None:
                    with open(ouf_path, 'a') as ouf:
                        introns_s = list(introns_s)
                        ready_iso = []
                        if introns_s == []:
                                ready_iso = []
                        else:
                            for elem in introns_s:
                                start, stop = elem
                                ready_iso.append((start, stop))
                            ready_iso = sorted(list(ready_iso))
                        ready_iso = ';'.join(['%s-%s' % pos for pos in ready_iso])
                        ouf.write(f'{prev_transcript_start_stop[0]}\t{prev_transcript_start_stop[1]}\t{prev_transcript_start_stop[2]}\t{ready_iso}\t{prev_transcript_start_stop[3]}\t{prev_transcript_id}\n')


                prev_exon_end = False
                introns_s = set()
                prev_transcript_start_stop = (line_l[0], line_l[6], line_l[3], line_l[4])
                prev_transcript_id = re.findall(r'gene_id "(.+?)"', line)[0] + f'.{str(counter_iso)}'
                counter_iso += 1
            if feat_type == 'exon':
                orientation = line_l[6]
                start = int(line_l[3])
                end = int(line_l[4])
                if orientation == '+':
                    if not prev_exon_end:
                        prev_exon_end = end
                    else:
                        introns_s.add((prev_exon_end + 1, start - 1))
                        prev_exon_end = end
                elif orientation == '-':
                    if not prev_exon_end:
                        prev_exon_end = start
                    else:
                        introns_s.add((end + 1, prev_exon_end - 1))
                        prev_exon_end = start
            if feat_type == 'gene':
                counter_iso = 1

    with open(ouf_path, 'a') as ouf:
        introns_s = list(introns_s)
        ready_iso = []
        if introns_s == []:
                ready_iso = []
        else:
            for elem in introns_s:
                start, stop = elem
                ready_iso.append((start, stop))
            ready_iso = sorted(list(ready_iso))
        ready_iso = ';'.join(['%s-%s' % pos for pos in ready_iso])
        ouf.write(f'{prev_transcript_start_stop[0]}\t{prev_transcript_start_stop[1]}\t{prev_transcript_start_stop[2]}\t{ready_iso}\t{prev_transcript_start_stop[3]}\t{prev_transcript_id}\n')


def parser_args():
    parser = argparse.ArgumentParser(description='This script processes the GTF annotation file and generates isoform structures')
    parser.add_argument('--annot_fpath', required=True, help='Path to the input reference GTF file')
    parser.add_argument('--ouf_path', required=True, help='Path to the output isoform structure file')

    return parser.parse_args()


if __name__ == '__main__':
    args = parser_args()
    create_transcript_struct(args.annot_fpath, args.ouf_path)
    # annot_fpath = '../data/assembly/GCF_000001735.4_TAIR10.1_genomic.gtf'
    # ouf_path = '../tmp/ref_isoforms_w_starts_ends.tsv'
    # create_transcript_struct(inf_path, ouf_path)
