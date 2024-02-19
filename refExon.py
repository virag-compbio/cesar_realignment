import sys, os
import utils

class ReferenceExon:

    def __init__(self, exons_list, chromosome, strand, two_bit_file, verbosity):
        self.exons_list = exons_list
        self.chromosome = chromosome
        self.strand = strand
        self.two_bit_file = two_bit_file
        self.verbosity = verbosity

    def extract_seqs_splice_sites(self):
        ''' wrapper function which has calls to other functions that do the actual work '''

        seqs_dict = self.__get_reference_exon_seqs()
        splice_sites_dict = self.__get_reference_splice_sites(seqs_dict)
        exon_phases_list = self.__get_exon_phases()
        exon_seqs_formatted = self.__get_reference_exon_seqs_formatted(seqs_dict, exon_phases_list)

        if self.verbosity > 1:
            print('Here are the input sequences for CESAR\n')
            for exon_ct, exon_seq in exon_seqs_formatted.items():
                print(f'{exon_ct}\t{exon_seq}')

        return exon_seqs_formatted, splice_sites_dict, exon_phases_list

    def __get_reference_exon_seqs(self):
        ''' Given a two-bit file, a list of exonic coordinates, chromosome and strand,
        extract the exonic sequences together with the splice sites '''

        temp_bed  = utils.create_temp_devshm()
        temp_fasta = utils.create_temp_devshm()

        with open(temp_bed, 'w') as FO:
            for exon_ct, exon in enumerate(self.exons_list):
                start_ss, stop_ss = exon.split('-')
                start_ss = int(start_ss) - 2
                stop_ss = int(stop_ss) + 2
                FO.write(f"{self.chromosome}\t{start_ss}\t{stop_ss}\texon_{exon_ct}\n")

        two_bit_to_fa_call = f'twoBitToFa {self.two_bit_file} {temp_fasta} -bed={temp_bed}'
        utils.run_subprocess(two_bit_to_fa_call)
        seqs_dict = utils.fasta_to_dict(temp_fasta)

        if self.strand == '-':
            seqs_dict = {exon_ct:utils.reverse_comp(exonic_seq) for exon_ct, exonic_seq in seqs_dict.items()}

        os.remove(temp_bed)
        os.remove(temp_fasta)    
        return seqs_dict

    def __get_reference_splice_sites(self, seqs_dict):
        ''' The input is a dictionary of sequences where the key is the exonic index and the value is the exonic sequence together with the splice sites.
        The output is a 2 dimensional dictionary - one dimension is the type of splice site i.e. acceptor or donor. The second dimension is the exonic
        index. The value is the splice site dinucleotide itself '''

        splice_sites_dict = dict()
        splice_sites_dict["acc"] = dict()
        splice_sites_dict["donor"] = dict()

        for exon_ct, exon_seq in seqs_dict.items():
            splice_sites_dict["acc"][exon_ct] = exon_seq[0:2]
            splice_sites_dict["donor"][exon_ct] = exon_seq[-2:]
        return splice_sites_dict

    def __get_exon_phases(self):

        exons_list = self.exons_list
        
        length_cum = 0
        exons_length  = list()
        for exon in exons_list:
            start, stop = exon.split('-')
            start = int(start)
            stop = int(stop)
            exon_len = abs(start - stop)
            length_cum += exon_len
            exons_length.append(length_cum)

        exons_phase_list = list()
        for exon_ct, exon in enumerate(exons_length):
            phase_5prime, phase_3prime = '', ''

            if exon_ct != 0:
                phase_5prime = 3 - exons_length[exon_ct -1]%3
                phase_3prime = exons_length[exon_ct]%3
            else:
                phase_5prime = 0
                phase_3prime = exons_length[exon_ct]%3

            phase_5prime = 0 if phase_5prime == 3 else phase_5prime
            phase_3prime = 0 if phase_3prime == 3 else phase_3prime
            exons_phase_list.append(f'{phase_5prime}-{phase_3prime}')
        
        return exons_phase_list

    def __get_reference_exon_seqs_formatted(self, seqs_dict, exons_phase_list):
        ''' The input is a dictionary of sequences where the key is the exonic index and the value is the exonic sequence together with splice sites.
        The function first strips the sequence of the splice sites and then formats the exonic sequence according to the exonic phases (another input provided
        in a list), so that CESAR can accept it '''

        seqs_dict_formatted = dict()
        i = 0
        for exon_ct, sequence in seqs_dict.items():
            phase_5prime, phase_3prime = exons_phase_list[i].split('-')
            phase_5prime = int(phase_5prime)
            phase_3prime = int(phase_3prime)

            sequence = sequence[2:-2].upper() ## stripping the splice sites
            sequence = "|" + sequence if phase_5prime == 0 else sequence[0:phase_5prime].lower() + "|" + sequence[phase_5prime:]
            sequence = sequence + "|" if phase_3prime == 0 else sequence[:-phase_3prime] + "|" + sequence[-phase_3prime:].lower()

            seqs_dict_formatted[exon_ct] = sequence
            i += 1
        return seqs_dict_formatted
