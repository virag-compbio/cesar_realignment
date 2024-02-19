import utils
import os, sys

class Alignment:

    def __init__(self, cesar_output, output_fasta, ref_species, mode = 'both'):
        self.cesar_output = cesar_output
        self.output_fasta = output_fasta
        self.ref_species  = ref_species
        self.mode = mode

    def process_split_alignments(self):
    ## the master function that acts as wrapper to other functions

        alignments_dict = self.__read_cesar_output()    
        alignments_split_exon_wise_full, alignments_split_exon_wise, exons_ct = self.split_alignments(alignments_dict)
        
        output_mafs_dict, species_list_final = self.__combine_pairwise_alignments(alignments_split_exon_wise, exons_ct)
        
        ## create a mullti-fasta
        self.__generate_multi_fasta(output_mafs_dict, self.output_fasta, species_list_final)
        
        ## return alignments that can be checked for exon-intactness
        return alignments_split_exon_wise_full

    def __read_cesar_output(self):
        ''' Reads the CESAR output and returns a dictionary where the key is the species while the value is a \n separated
        string containing the aligned reference, the query header and the aligned query sequence '''
        
        cesar_alignment = self.cesar_output.split('\n')
        aligned_seqs_list = [line.rstrip('\n') for line in cesar_alignment]
        
        alignments_dict = dict()
        i = 0
        while i < len(aligned_seqs_list):
            if aligned_seqs_list[i].startswith('>referenceExon'):

                query_seq = aligned_seqs_list[i+2]
                query_species, *_ = query_seq.split('#')
                query_species = query_species.replace('>','')

                alignments_dict[query_species] = f'{aligned_seqs_list[i+1]}\n{aligned_seqs_list[i+2]}\n{aligned_seqs_list[i+3]}'
                i += 3
            else:
                i += 1

        return alignments_dict

    def split_alignments(self, alignments_dict):
        ''' Takes the pairwise alignments and split them exon wise. Returns a dictionary of dictionaries. The inner-level dict has the exon as the key and the alignment 
        for that exon as the value. The outer level dict has the species as the key and the inner level dict as the value'''

        alignments_exon_wise_dict_full = dict()
        alignments_exon_wise_dict = dict()

        for species, alignment in alignments_dict.items():
            ref_seq, query_header, query_seq = alignment.split('\n')    

            alignments_dict_full, alignments_dict = dict(), dict()
            alignment_seqs_list = ref_seq.split()
            start_index = 0
            
            for exon_ct, aligned_exon in enumerate(alignment_seqs_list):
                exon_start = alignment.find(aligned_exon, start_index)
                exon_end   = exon_start + len(aligned_exon)
                
                aligned_ref = exon_start * ' ' + ref_seq[exon_start:exon_end]
                alignments_dict_full[exon_ct] = f'{aligned_ref}\n{query_header}\n{query_seq}'
                alignments_dict[exon_ct] = f'{ref_seq[exon_start:exon_end]}\n{query_seq[exon_start:exon_end]}'
                

            alignments_exon_wise_dict[species] = alignments_dict
            alignments_exon_wise_dict_full[species] = alignments_dict_full

        return alignments_exon_wise_dict_full, alignments_exon_wise_dict, exon_ct + 1


    def __combine_pairwise_alignments(self, alignments_exon_wise_dict, exons_ct):
        species_list_final = list(alignments_exon_wise_dict.keys())

        tmp_files_dict   = {species : utils.create_temp_devshm() for species in species_list_final}
        output_mafs_dict = {exon_index : '' for exon_index in list(range(exons_ct))}

        i = 0
        while i < exons_ct:

            len_seq = lambda input_seq: len(input_seq.replace('-','')) ## lambda functions that returns the length of the sequence minus the gaps     
            mafs_list = ''

            for species in species_list_final:
                ref_seq, query_exon = alignments_exon_wise_dict[species][i].rstrip('\n').split('\n')
                tmp_file = tmp_files_dict[species]

                with open(tmp_file, 'w') as FO:
                    FO.write(f"s {self.ref_species}.chr16 100000 {len_seq(ref_seq)} + 90338345 {ref_seq}\n")
                    FO.write(f"s {species}.chr16 100000 {len_seq(query_exon)} + 90338345 {query_exon}\n")
                mafs_list += f' {tmp_file}'
                    
            maf_join_call = f"maf-join {mafs_list}"
            maf_output = utils.run_subprocess(maf_join_call, True)
            output_mafs_dict[i] = maf_output

            i += 1

        for f in (list(tmp_files_dict.values())):
            os.remove(f)

        return output_mafs_dict, species_list_final
        
    def __generate_multi_fasta(self, output_mafs_dict, output_fasta, species_list_final):
    
        multiple_aligned_seqs_dict = {species : '' for species in species_list_final}    
        multiple_aligned_seqs_dict[self.ref_species] = ''

        for exon_index, combined_maf in output_mafs_dict.items():
            aligned_seqs_dict = dict()
            
            for line in combined_maf.split('\n'):
                if line.startswith('s'):  ## "s hg38.chr16       47429 72 + 90338345 CTAATCCCTCCAGCGGTGTCCACACTGAGCATTGCAGCACTTGTAGAAGGTGGTCATCGGCTCATCTGCAGA"
                    species, *_, aligned_seq = utils.get_sline_data(line)
                    multiple_aligned_seqs_dict[species] += aligned_seq

        with open(output_fasta, 'w') as FO:
            for species, aligned_seq in multiple_aligned_seqs_dict.items():
                FO.write(f'>{species}\n{aligned_seq}\n')        