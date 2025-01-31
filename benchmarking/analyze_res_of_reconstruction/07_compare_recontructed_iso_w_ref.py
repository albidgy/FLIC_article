import argparse
import os


def parser_args():
    parser = argparse.ArgumentParser(description='This script calculates statistics on isoform reconstruction results')

    parser.add_argument('--ref_iso', required=True, 
                        help='Path to the reference isoforms file')
    parser.add_argument('--sim_transcripts_dir', required=True, 
                        help='Path to the simulated transcripts directory')
    parser.add_argument('--modes_info', required=True, help='Path to transcripts mode information file')
    parser.add_argument('--reconstructed_iso_dir', required=True, 
                        help='Path to the directory containing reconstructed isoforms')
    parser.add_argument('--out_dir', required=True, 
                        help='Path to the output directory for statistics')

    return parser.parse_args()


def extract_real_transcripts_struct(ref_transcripts_fpath):
    d_ref_transcripts_all = {}
    with open(ref_transcripts_fpath) as inf:
        for line in inf:
            line_l = line.strip('\n').split('\t')
            transcript_id = line_l[-1]
            d_ref_transcripts_all[transcript_id] = '\t'.join(line_l[:-1])

    return d_ref_transcripts_all


def read_simulated_transcripts(inf_path, d_ref_transcripts_all, d_transcripts_cov, idx):
    with open(inf_path) as inf:
        for line in inf:
            line_l = line.strip('\n').split('\t')
            transcript_id = line_l[1]

            if transcript_id in d_ref_transcripts_all.keys():
                transcript_id_w_struct = f'{transcript_id} {d_ref_transcripts_all[transcript_id]}'
                if transcript_id_w_struct not in d_transcripts_cov.keys():
                    d_transcripts_cov[transcript_id_w_struct] = [0, 0, 0]
                d_transcripts_cov[transcript_id_w_struct][idx] += int(line_l[-1])
    return d_transcripts_cov


def split_data_by_expr(d_transcripts_cov):
    d_good_expr_transcripts = {}
    d_bad_expr_transcripts = {}
    for key, val_l in d_transcripts_cov.items():
        transcript_id, iso_struct = key.split(' ')
        gene_id = '.'.join(transcript_id.split('.')[:-1])

        if sum([x >= 1 for x in val_l]) > 1 and sum([x >= 5 for x in val_l]) > 0:
            if gene_id not in d_good_expr_transcripts.keys():
                d_good_expr_transcripts[gene_id] = []
            d_good_expr_transcripts[gene_id].append((transcript_id, iso_struct))
        else:
            if gene_id not in d_bad_expr_transcripts.keys():
                d_bad_expr_transcripts[gene_id] = []
            d_bad_expr_transcripts[gene_id].append((transcript_id, iso_struct))

    return d_good_expr_transcripts, d_bad_expr_transcripts


def prep_ref_data(ref_iso_fpath, sim_transcripts_path):
    d_ref_transcripts_all = extract_real_transcripts_struct(ref_iso_fpath)

    d_transcripts_cov = {}
    for idx, file in enumerate(os.listdir(sim_transcripts_path)):
        d_transcripts_cov = read_simulated_transcripts(os.path.join(sim_transcripts_path, file),
                                                       d_ref_transcripts_all,
                                                       d_transcripts_cov,
                                                       idx)
    ref_pos, ref_neg = split_data_by_expr(d_transcripts_cov)
    return ref_pos, ref_neg, len(d_ref_transcripts_all)


def is_find_ref_iso(s_ref_iso, peak_start, peak_end, reconstructed_splice_sites):
    for ref_iso_w_id in s_ref_iso:
        transcript_id, ref_iso = ref_iso_w_id
        ref_iso_l = ref_iso.split('\t')
        ref_start = int(ref_iso_l[2])
        ref_end = int(ref_iso_l[4])
        ref_splice_sites = ref_iso_l[3]
        if (reconstructed_splice_sites == ref_splice_sites and
                peak_start[0] <= ref_start <= peak_start[1] and
                peak_end[0] <= ref_end <= peak_end[1]):
            return transcript_id
    return 'None'


def calc_n_matched_iso(reconstructed_iso_fpath, ref_isoforms):
    s_matches = set()
    n_mismatches = 0
    with open(reconstructed_iso_fpath) as inf:
        for line in inf:
            *iso_struct, gene_id = line.strip('\n').split('\t')
            peak_start = tuple(map(int, iso_struct[2].split('-')))
            peak_end = tuple(map(int, iso_struct[4].split('-')))
            reconstructed_splice_sites = iso_struct[3]

            if gene_id in ref_isoforms.keys():
                transcript_id = is_find_ref_iso(ref_isoforms[gene_id], peak_start,
                                                peak_end, reconstructed_splice_sites)
                if transcript_id != 'None':
                    s_matches.add(transcript_id)
                else:
                    n_mismatches += 1

    return s_matches, n_mismatches


def write_stat(ouf_path, tp, fp, fn, tn):
    precision = len(tp) / (len(tp) + fp)
    recall = len(tp) / (len(tp) + fn)
    f1_score = 2 * precision * recall / (precision + recall)

    with open(ouf_path, 'w') as ouf:
        ouf.write(f'Precision\t{precision:.4f}\n')
        ouf.write(f'Recall\t{recall:.4f}\n')
        ouf.write(f'f1-score\t{f1_score:.4f}\n')

        ouf.write(f'TP\t{len(tp)}\n')
        ouf.write(f'FP\t{fp}\n')
        ouf.write(f'FN\t{fn}\n')
        ouf.write(f'TN\t{tn}\n')


def read_mod_info(mod_info_file):
    d_mods = {}
    with open(mod_info_file) as inf:
        for line in inf:
            transcript_id, mode = line.strip('\n').split('\t')
            d_mods[transcript_id] = mode
    return d_mods


def write_tp_by_modes(mod_info_file, tp, ouf_path):
    d_preds_split_by_modes = dict.fromkeys(list(map(str, range(0, 7))), 0)
    d_mods = read_mod_info(mod_info_file)
    for elem in tp:
        if elem in d_mods.keys():
            d_preds_split_by_modes[d_mods[elem]] += 1

    with open(ouf_path, 'a') as ouf:
        ouf.write('Stat by modes:\n')
        for key, val in d_preds_split_by_modes.items():
            ouf.write(f'{key}\t{val}\n')


def main(ref_iso_fpath, sim_transcripts_path, reconstructed_iso_fpath, out_dir, mod_info_file):
    ouf_path = os.path.join(out_dir, os.path.basename(reconstructed_iso_fpath))
    ref_pos, ref_neg, n_transcripts = prep_ref_data(ref_iso_fpath, sim_transcripts_path)
    tp, fp = calc_n_matched_iso(reconstructed_iso_fpath, ref_pos)
    _, tn = calc_n_matched_iso(reconstructed_iso_fpath, ref_neg)
    fn = n_transcripts - len(tp)
    write_stat(ouf_path, tp, fp, fn, tn)
    write_tp_by_modes(mod_info_file, tp, ouf_path)


if __name__ == '__main__':
    args = parser_args()

    for file in os.listdir(args.reconstructed_iso_dir):
        if os.path.isfile(os.path.join(args.reconstructed_iso_dir, file)):
            reconstructed_iso_fpath = os.path.join(args.reconstructed_iso_dir, file)
            main(args.ref_iso, args.sim_transcripts_dir, reconstructed_iso_fpath, args.out_dir, args.modes_info)
