import utils
import indel_functions

class AlignedSeqPair:

    def __init__(self, ref_sequence, query_sequence, query_header, phase_string, ref_acc, ref_donor, exon_index, exon_ct, ref_strand, verbosity):
        self.ref_sequence = ref_sequence
        self.query_sequence = query_sequence
        self.query_header = query_header
        self.phase_string = phase_string
        self.ref_acc = ref_acc
        self.ref_donor = ref_donor
        self.exon_index = exon_index
        self.exon_ct = exon_ct
        self.ref_strand = ref_strand
        self.verbosity = verbosity
        
    def check_exon_intactness(self):
        ''' The wrapper method that calls other methods which do the actual work/heavy lifting '''

        reference = self.ref_sequence
        query     = self.query_sequence
        phase_5prime, phase_3prime = self.phase_string.split('-')
        phase_5prime, phase_3prime = int(phase_5prime), int(phase_3prime)

        alignment_start, alignment_stop = self.__get_alignment_boundaries(reference, query)
        ref_exonic_seq, query_exonic_seq = self.__extract_aligned_exonic_seq(alignment_start, alignment_stop, phase_5prime, phase_3prime)
        codon_pos_dict_ref = self.__get_codon_pos_dict(ref_exonic_seq)
        codon_pos_dict_query = self.__get_codon_pos_dict(query_exonic_seq)

        inframe_stops_count = self.__find_inframe_stops(ref_exonic_seq, query_exonic_seq, codon_pos_dict_ref)
        ss_mutations_count  = self.__get_splice_site_mutations(query, alignment_start, alignment_stop, self.ref_acc, self.ref_donor, self.exon_index, self.exon_ct)
        inserts_with_stop, frameshifts_count   = self.__get_frameshifts(ref_exonic_seq, query_exonic_seq, codon_pos_dict_ref, codon_pos_dict_query)

        if self.verbosity > 1:
            print(f'The number of inframe stops is {inframe_stops_count}')
            print(f'The number of splice site mutations is {ss_mutations_count}')
            print(f'The number of frameshifts is {frameshifts_count}')
            print(f'The number of inserts with stop codons is {inserts_with_stop}')

        if inframe_stops_count + ss_mutations_count + frameshifts_count + inserts_with_stop > 1:
            return "NI"
        else:
            return (self.__get_intact_exon_coordinates(self.query_header, alignment_start, alignment_stop, self.exon_index, len(ref_exonic_seq), self.query_sequence, self.ref_strand))

    def __get_alignment_boundaries(self, reference, query):

        ref_seq_list = reference.split()
        alignment_start = reference.index(ref_seq_list[0])
        alignment_end   = alignment_start + len(ref_seq_list[0])
        return alignment_start, alignment_end

    def __extract_aligned_exonic_seq(self, alignment_start, alignment_end, phase_5prime, phase_3prime):

        ref_exonic_seq   = self.ref_sequence[alignment_start:alignment_end]
        query_exonic_seq = self.query_sequence[alignment_start:alignment_end]

        nt_add_dict = {0 : '', 1 : 'NN', 2 : 'N'}
        nt_add_5p = nt_add_dict[phase_5prime]
        nt_add_3p = nt_add_dict[phase_3prime]

        ref_exonic_seq   = (nt_add_5p + ref_exonic_seq + nt_add_3p).upper()
        query_exonic_seq = (nt_add_5p + query_exonic_seq + nt_add_3p).upper()
        return ref_exonic_seq, query_exonic_seq

    def __get_codon_pos_dict(self, nt_seq):
        codon_pos_dict = dict()

        ct, base_ungapped = 0, 0
        for nt in nt_seq:
            if nt != '-':
                base_ungapped += 1
                codon_pos = base_ungapped%3
                if codon_pos == 0:
                    codon_pos = 3

                codon_pos_dict[ct] = codon_pos
            else:
                codon_pos_dict[ct] = "-"

            ct += 1

        return codon_pos_dict

    def __find_inframe_stops(self, ref_sequence, query_sequence, codon_pos_dict_ref):

        stop_codons = ['TAA', 'TAG', 'TGA']
        stop_codon_count = 0
        i = 0
        while i < len(query_sequence):
            codon_query = query_sequence[i:i+3]
            codon_ref   = ref_sequence[i:i+3]

            if codon_query in stop_codons and codon_ref not in stop_codons:
                codon_pos_string = str(codon_pos_dict_ref[i]) + str(codon_pos_dict_ref[i+1]) + str(codon_pos_dict_ref[i+2])
                if codon_pos_string == '123':
                    stop_codon_count += 1

            i += 3

        return stop_codon_count

    def __get_splice_site_mutations(self, query_aligned, alignment_start, alignment_stop, ref_acc, ref_donor, exon_index, exon_ct):

        acc_query = query_aligned[alignment_start - 2:alignment_start].upper()
        donor_query = query_aligned[alignment_stop:alignment_stop + 2].upper() 

        splice_site_muts_count = 0
        if acc_query != 'AG' and acc_query != ref_acc and exon_index != 0:
            splice_site_muts_count += 1
        if donor_query != 'GT' and donor_query != 'GC' and donor_query != ref_donor and exon_index != exon_ct - 1:
            splice_site_muts_count += 1
        
        return splice_site_muts_count

    def __get_frameshifts(self, ref_exonic_sequence, query_exonic_sequence, codon_dict_ref, codon_dict_query):
        insertions = utils.get_indels(ref_exonic_sequence)
        deletions  = utils.get_indels(query_exonic_sequence)

        ## first deal with framepreserving insertions as they may contain an in-frame stop codon, this needs to be reported
        fp_insertions = {pos:length for pos, length in insertions.items() if length%3 == 0}
        inserts_with_stop = 0
        for insert, insert_length in fp_insertions.items():
            seq_insert = query_exonic_sequence[insert:insert + insert_length]
            if utils.find_stop_codon(seq_insert):
                inserts_with_stop += 1
        ## that's it. Now deal with frameshifts

        ## filter insertions or deletions that are frame-disrupting
        fs_insertions = {pos:length for pos, length in insertions.items() if length%3 != 0}
        fs_deletions  = {pos:0-length for pos, length in deletions.items() if length%3 != 0}
        all_frameshifts_dict = {**fs_insertions, **fs_deletions}
        fs_count = len(all_frameshifts_dict)

        ## check for the possibility of compensating frameshifts
        if fs_count >= 2:
            all_frameshifts_dict = dict(sorted(all_frameshifts_dict.items()))
            excluded_insertions, excluded_deletions = dict(), dict()

            compensated_events = indel_functions.get_compensated_events(all_frameshifts_dict, codon_dict_ref, codon_dict_query, query_exonic_sequence)
            ignore_compensated_events = list()

            for comp_event_key, comp_events_list in compensated_events.items():
                ignore_compensated_events.extend(comp_events_list)

            all_frameshifts_dict_non_compensated = [fs for fs in all_frameshifts_dict.keys() if fs not in ignore_compensated_events]
            fs_count = len(all_frameshifts_dict_non_compensated)

        return fs_count, inserts_with_stop

    def __get_intact_exon_coordinates(self, header, alignment_start, alignment_stop, exon_index, length_aln, query_sequence_full, ref_strand):

        species, start, stop, query_strand, query_chr  = header.split('#')
        
        if ref_strand != query_strand:
            alignment_start = len(query_sequence_full) - alignment_stop

        exon_start = int(start) + alignment_start
        exon_end   = exon_start + length_aln

        query_strand_final = '-' if ref_strand != query_strand else '+'
        bed_line = f'{query_chr}\t{exon_start}\t{exon_end}\t{query_strand_final}\t{exon_index}'
        return bed_line
