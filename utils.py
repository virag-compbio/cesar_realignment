import subprocess, os

def run_subprocess(command, return_output = False):
    ''' Given a command or a call to a process, execute it '''
    command_return = subprocess.run(command.split(), capture_output = True, text = True)
    assert command_return.returncode == 0, f'Aborting, the call {command} did not run successfully'
    if return_output:
        return command_return.stdout.rstrip('\n')

def fasta_to_dict(fasta_file):
    ''' Given a fasta file, convert this to a dictionary where the key is the sequence header while the value is the sequence '''
    assert os.path.isfile(fasta_file), f'Aborting, the fasta file {fasta_file} could not be found'

    fasta_file_dict = dict()
    seq, seq_id = '', ''

    with open(fasta_file, 'r') as FI:
        for line in FI:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if seq != '':
                    fasta_file_dict[seq_id] = seq
                    seq, seq_id = '', ''
                seq_id = line.replace('>','')
            else:
                seq += line
    fasta_file_dict[seq_id] = seq
    return fasta_file_dict

def reverse_comp(sequence):
    ''' Return the reverse complement of a given sequence '''
    return sequence.translate({65: 84, 97: 116, 84: 65, 116: 97, 99: 103, 103: 99, 71: 67, 67: 71})[::-1]

def get_sline_data(line):
    ''' parse the slines from the maf and return the elements associated with the sline'''
    pos_strands = ["+", "-"]
    assert line.startswith('s'), f'Error in get_sline_data. The sline {line} does not start with s'

    line_list = line.split()
    assert len(line_list) == 7, f'Error in get_sline_data. Cannot parse 7 elements from {line}. The resulting list is {line_list}'
    _, src, start, size, strand, src_size, seq = line.split()
    ## sanity checks. Convert start, size and src_size into integers. This would fail if these values cannot be converted
    start = int(start)
    size = int(size)
    src_size = int(src_size)
    assert strand in pos_strands, f'Error in get_sline_data. Unidentified value of strand {strand} in the line {line}'

    species, *_ = src.split('.')
    chromosome = src[len(species) + 1:]
    return species, chromosome, start, size, strand, src_size, seq

def run_maf_extract(maf_index, output_file, chromosome, cds):
    ''' Given a maf_index i.e. .bb file and a coordinate i.e. chromosome and genomic coordinates,
    run mafextract and return the output'''
    mafex_call = f'mafExtract {maf_index} -region={chromosome}:{cds} {output_file}'
    run_subprocess(mafex_call)

def filter_maf_file(input_maf, output_maf, species_list):
    ''' A poor man's implementation of mafSpeciesSubset. Only filters slines for the species of interest'''

    with open(input_maf, 'r') as FI, open(output_maf, 'w') as FO:
        for line in FI:
            if line.startswith("s"):
                species, *_ = get_sline_data(line)
                if species in species_list:
                    FO.write(f'{line}')

def create_temp_devshm():
    return(run_subprocess("mktemp /dev/shm/XXXXX", True))

def get_sequence_from_2bit(cds_start, cds_stop, query_chr, query_2bit_file):
    two_bit_to_fa_call = f'twoBitToFa {query_2bit_file}:{query_chr}:{cds_start}-{cds_stop} stdout'
    query_seq = run_subprocess(two_bit_to_fa_call, True)
    return ''.join([seq_line for seq_line in query_seq.split('\n') if '>' not in seq_line])

def find_stop_codon(nt_seq):
    ''' Given a nucleotide sequence find if it contains a stop codon '''
    stop_codons = ['TAA', 'TAG', 'TGA']
    if len(nt_seq)%3 == 1:
        nt_seq = nt_seq + 'NN'
    elif len(nt_seq)%3 == 2:
        nt_seq = nt_seq + 'N'

    codon_positions = list(range(0, len(nt_seq), 3))
    codon_seq = [nt_seq[i:i + 3] for i in codon_positions]

    stop_present = [codon for codon in codon_seq if codon in stop_codons]
    if stop_present:
        return True
    return False

def get_indels(aligned_seq):
    ''' Given a nucleotide sequence, report the gaps i.e --- in this sequence '''
    indices = [ind for ind, ele in enumerate(aligned_seq) if ele == '-']
    ## from https://stackoverflow.com/questions/38740424/merge-adjacent-number-in-a-list-in-python: guillaume.deslandes
    if len(indices) == 0:
        return dict()
    gaps_merged = [[indices[0]]]
    for e in indices[1:]:
        if gaps_merged[-1][-1] == e - 1:
            gaps_merged[-1].append(e)
        else:
            gaps_merged.append([e])

    gap_pos_dict = dict()
    for element in gaps_merged:
        gap_pos = element[0]
        gap_pos_dict[gap_pos] = len(element)

    return gap_pos_dict


