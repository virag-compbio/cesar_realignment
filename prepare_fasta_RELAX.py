import sys, os
import argparse
import utils

def flip_stop_codon(nt_seq):
    ''' Given a nucleotide sequence find if it contains a stop codon '''
    stop_codons = ['TAA', 'TAG', 'TGA']

    codon_positions = list(range(0, len(nt_seq), 3))
    codon_seq = [nt_seq[i:i + 3] for i in codon_positions]
    return ''.join([codon if codon not in stop_codons else 'NNN' for codon in codon_seq])

if __name__ == '__main__':

    ## argument parsing business
    parser = argparse.ArgumentParser(description = "Convert a CESAR generated multi-fasta to a file that can be used as an input to RELAX. Flips stop codons to NNN and removes insertions that disrupt the reading frame")
    parser.add_argument("input", help = 'The input fasta', type = str)
    parser.add_argument("output", help = 'The output fasta (RELAX-ready)', type = str)
    
    args = parser.parse_args()
    input_file   = args.input
    output_file  = args.output

    seqs_dict = utils.fasta_to_dict(input_file)
    ref_seq = seqs_dict["hg38"]
    ref_seq_gaps = utils.get_indels(ref_seq)

    ## get only the frame-disrupting indels
    ref_seq_gaps_fd = {gap_pos: gap_length for gap_pos, gap_length in ref_seq_gaps.items() if gap_length%3 != 0}

    for seq_id, aligned_seq in seqs_dict.items():
        for gap_pos, gap_length in ref_seq_gaps_fd.items():
            nt_seq_updated = aligned_seq[:gap_pos] + aligned_seq[gap_pos + gap_length:]
            seqs_dict[seq_id] = nt_seq_updated

    ## after removing all inserts that wreck the reading frame, flip any potential stop-codons to NNN and write to the output file
    with open(output_file, 'w') as FO:
        for seq_id, aligned_seq in seqs_dict.items():
            aligned_seq = flip_stop_codon(aligned_seq)
            FO.write(f'>{seq_id}\n{aligned_seq}\n')


    