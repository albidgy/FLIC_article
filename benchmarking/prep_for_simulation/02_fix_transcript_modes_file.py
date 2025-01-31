import argparse
import re


def parser_args():
    parser = argparse.ArgumentParser(description='This script corrects the transcript group information file based on the distorted annotation')
    parser.add_argument('--real_annot', required=True, 
                        help='Path to the real annotation file')
    parser.add_argument('--bad_annot', required=True, 
                        help='Path to the distorted annotation file')
    parser.add_argument('--source_modes', required=True, 
                        help='Path to the source transcripts mode file')
    parser.add_argument('--ouf_name', required=True, 
                        help='Path to the fixed transcripts mode file')

    return parser.parse_args()


def read_annot_file(inf_path):
    transcript_coords = {}
    with open(inf_path) as inf:
        for line in inf:
            line_l = line.strip('\n').split('\t')
            feat_type = line_l[2]
            if feat_type == 'transcript':
                transcript_id = re.findall(r'transcript_id "(.+?)"', line_l[8])[0]
                transcript_coords[transcript_id] = (int(line_l[3]), int(line_l[4]))
    return transcript_coords


def read_mod_info(mod_info_file):
    d_modes = {}
    with open(mod_info_file) as inf:
        for line in inf:
            transcript_id, mode = line.strip('\n').split('\t')
            d_modes[transcript_id] = mode
    return d_modes


def modify_transcripts(start, stop, orientation, type_mod):
    modifications_by_groups = {0: (0, 0), 1: (100, 0), 2: (0, -100), 3: (-100, 0), 
                               4: (0, 100), 5: (100, -100), 6: (-100, 100)}
    
    start_add, stop_add = modifications_by_groups[type_mod]
    if orientation == '+':
        start_mod = start + start_add
        stop_mod = stop + stop_add
    else:
        start_mod = start - stop_add
        stop_mod = stop - start_add
    
    if stop_mod - start_mod > 0:
        return start_mod, stop_mod
    else:
        return start, stop

    
def main(real_annot, bad_annot, source_modes, ouf_name):
    real_transcript_coords = read_annot_file(real_annot)
    bad_transcript_coords = read_annot_file(bad_annot)
    d_modes = read_mod_info(source_modes)

    for transcript_id, ref_borders in real_transcript_coords.items():
        bad_borders = bad_transcript_coords[transcript_id]
        if ref_borders == bad_borders:
            d_modes[transcript_id] = '0'

    with open(ouf_name, 'w') as ouf:
        for key, val in d_modes.items():
            ouf.write(f'{key}\t{val}\n')
    
    return d_modes


if __name__ == '__main__':
    args = parser_args()
    main(args.real_annot, args.bad_annot, args.source_modes, args.ouf_name)
