import argparse
import os


def parser_args():
    parser = argparse.ArgumentParser(description="This script converts TSS and PA peaks of reconstructed FLIC isoforms into points, where point represents the summit of the peak")
    parser.add_argument('--cagefightr_tss', required=True, 
                        help='Path to the CAGEfigthR TSSs')
    parser.add_argument('--cagefightr_pa', required=True, 
                        help='Path to the CAGEfigthR PAs')
    parser.add_argument('--isoform_fpath', required=True, 
                        help='Path to the file containing reconstructed FLIC isoforms')
    parser.add_argument('--ouf_path', required=True, 
                        help='Path to the output file')

    return parser.parse_args()


def get_max_tss_pa_peaks(inf_fpath, d_tss_w_max, d_polya_w_max):
    with open(inf_fpath) as inf:
        for line in inf:
            line_l = line.strip('\n').split('\t')
            chrom = line_l[0]
            orientation = line_l[5]
            start_peak = int(line_l[1]) + 1
            end_peak = int(line_l[2])
            max_start_peak = int(line_l[6]) + 1
            max_end_peak = int(line_l[7])
            max_peak = int(round((max_end_peak - max_start_peak) / 2 + max_start_peak, 0))

            if orientation == '+':
                d_tss_w_max[f'{chrom}*{orientation}*{start_peak}-{end_peak}'] = max_peak
            else:
                d_polya_w_max[f'{chrom}*{orientation}*{start_peak}-{end_peak}'] = max_peak
    return d_tss_w_max, d_polya_w_max


def write_max_peaks_by_iso(d_tss_w_max, d_polya_w_max, isoform_fpath, ouf_path):
    with open(ouf_path, 'w') as ouf:
        ouf.write('')
        
    d_isoforms_w_max_peaks = {}
    with open(isoform_fpath) as inf:
        for line in inf:
            if line[0] == '#':
                continue
            *line_l, _ = line.strip('\n').split('\t')
            iso_start_info = f'{line_l[0]}*{line_l[1]}*{line_l[2]}'
            iso_end_info = f'{line_l[0]}*{line_l[1]}*{line_l[4]}'
            line_l[2] = str(d_tss_w_max[iso_start_info])
            line_l[4] = str(d_polya_w_max[iso_end_info])
            
            with open(ouf_path, 'a') as ouf:
                ouf.write('\t'.join(line_l) + '\n')
   
            
def main(inf_fpath_start, inf_fpath_end, isoform_fpath, ouf_path):
    d_tss_w_max = {}
    d_polya_w_max = {}
    d_tss_w_max, d_polya_w_max = get_max_tss_pa_peaks(inf_fpath_start, d_tss_w_max, d_polya_w_max)
    d_polya_w_max, d_tss_w_max = get_max_tss_pa_peaks(inf_fpath_end, d_polya_w_max, d_tss_w_max)
    write_max_peaks_by_iso(d_tss_w_max, d_polya_w_max, isoform_fpath, ouf_path)


if __name__ == '__main__':
    args = parser_args()
    main(args.cagefightr_tss, args.cagefightr_pa, args.isoform_fpath, args.ouf_path)
