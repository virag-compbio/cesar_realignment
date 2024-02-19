import sys, os
import utils

class QueryExonSeqs:

    def __init__(self, ref_species, exons_list, chromosome, ref_strand, two_bit_dir, query_species_list, maf_index, verbosity):
        self.ref_species = ref_species
        self.exons_list = exons_list
        self.chromosome = chromosome
        self.ref_strand = ref_strand
        self.two_bit_dir = two_bit_dir
        self.query_species_list = query_species_list
        self.maf_index = maf_index
        self.verbosity = verbosity

    def get_query_seqs(self): ## this is the wrapper function that contains calls to different functions which do the actual work

        temp_maf = utils.create_temp_devshm()
        temp_maf_filtered = utils.create_temp_devshm()

        start_dict, stop_dict, chr_dict, strand_dict = dict(), dict(), dict(), dict()
        exon_indices_list = list(range(0, len(self.exons_list)))

        for species in self.query_species_list:
            start_dict[species] = {exon_index:None for exon_index in exon_indices_list}
            stop_dict[species]  = {exon_index:None for exon_index in exon_indices_list}
            chr_dict[species]   = {exon_index:None for exon_index in exon_indices_list}
            strand_dict[species] = {exon_index:None for exon_index in exon_indices_list}
        
        if self.verbosity > 1:
            print(f'Species under consideration are {self.query_species_list}')
    
        for exon_index, exon in enumerate(self.exons_list):
            if self.verbosity > 1:
                print(f'The following exon is under consideration: {self.chromosome}-{exon}')

            utils.run_maf_extract(self.maf_index, temp_maf, self.chromosome, exon)
            utils.filter_maf_file(temp_maf, temp_maf_filtered, self.query_species_list)

            cds_start_query, cds_stop_query, strand_dict_query, chr_dict_query = self.__get_coordinates_from_maf(temp_maf_filtered)
            os.remove(temp_maf)
            os.remove(temp_maf_filtered)

            for species in self.query_species_list:
                if species in cds_start_query:
                    start_dict[species][exon_index]  = cds_start_query[species]
                    stop_dict[species][exon_index]   = cds_stop_query[species]
                    chr_dict[species][exon_index]    = chr_dict_query[species]
                    strand_dict[species][exon_index] = strand_dict_query[species]            

        ## now check if all the exons for a transcript are colinear. If not then exclude these
        ## colinear means - one strand, one scaffold/chromosome. all the start_cds should be linearly arranged. all the stop_cds should be linearly arranged
        query_species_cds = dict()
        valid_species = 0
        species_ignore = list()

        for species in self.query_species_list:
            cds_starts = list(start_dict[species].values())
            cds_stops  = list(stop_dict[species].values())
            strand_set = set(strand_dict[species].values())
            chromosome_set = set(chr_dict[species].values())

            ## allow for those cases where one of the exons in the middle has a None value for cds_starts/cds_stop but do not allow for the first or the last exon to have a None value
            if cds_starts[0] == None or cds_starts[-1] == None:
                species_ignore.append(species)

                if self.verbosity > 1:
                    print(f'Omitting species {species} because the cds_starts have a None value at the beginning/end. Here is the list {cds_starts}')
                continue
            
            ########################################################################
            ## Remove None values from these lists
            if None in strand_set: 
               chromosome_set.remove(None)
               strand_set.remove(None)
               cds_starts = [start for start in cds_starts if start != None]
               cds_stops  = [stop for stop in cds_stops if stop != None]
            ########################################################################  

            if len(strand_set)  == 1 and len(chromosome_set) == 1:
                strand = list(strand_set)[0]
                chromosome = list(chromosome_set)[0]

                if strand == self.ref_strand:
                    if cds_starts == sorted(cds_starts) and cds_stops == sorted(cds_stops):
                        query_species_cds[species] = f'{species}#{cds_starts[0]}#{cds_stops[-1]}#{strand}#{chromosome}'
                        valid_species += 1
                else:
                    if cds_starts == sorted(cds_starts)[::-1] and cds_stops == sorted(cds_stops)[::-1]:
                        query_species_cds[species] = f'{species}#{cds_stops[-1]}#{cds_starts[0]}#{strand}#{chromosome}'
                        valid_species += 1
        
        return query_species_cds, valid_species

    def __get_coordinates_from_maf(self, input_maf):

        ref_species = self.ref_species
        species_list = self.query_species_list

        with open(input_maf, 'r') as FI:
            maf_extract_result = FI.readlines()
        maf_extract_result = [line.rstrip('\n') for line in maf_extract_result]

        cds_list_dict = dict()
        chr_tracker_dict, strand_tracker_dict = dict(), dict()
        for species in species_list:
            cds_list_dict[species] = list()
            chr_tracker_dict[species] = set()
            strand_tracker_dict[species] = set()

        for line in maf_extract_result:
            if line.startswith('s'):
                species, chromosome, start, stop, strand, src_size, seq = utils.get_sline_data(line)

                if species == ref_species:
                    continue

                cds_start_query = cds_stop_query = '', ''
                if strand == '-':
                    cds_start_query = src_size - start - stop
                    cds_stop_query = src_size - start
                else:
                    cds_start_query = start
                    cds_stop_query = start + stop

                cds_string = f'{cds_start_query}-{cds_stop_query}'

                cds_list_dict[species].append(cds_string)
                chr_tracker_dict[species].add(chromosome)
                strand_tracker_dict[species].add(strand)

        ## Now identify those species for which the aligning sequence comes from more than one scaffold/chromosome or strand
        species_exclude = [species for species in species_list if len(strand_tracker_dict[species]) != 1 or len(chr_tracker_dict[species]) != 1]
        species_list = [species for species in species_list if species not in species_exclude] ## and update the species_list

        strand_query_dict = {species:list(element)[0] for species, element in strand_tracker_dict.items() if len(element) == 1}
        chr_query_dict    = {species:list(element)[0] for species, element in chr_tracker_dict.items() if len(element) == 1}

        cds_start_query_dict, cds_stop_query_dict = dict(), dict()
        for species in species_list:
            cds_list = cds_list_dict[species]

            if len(cds_list) > 0:
                if strand_query_dict[species] == '-':
                    cds_list = cds_list[::-1]

                first_element = cds_list[0]
                last_element = cds_list[-1]
                cds_start_query_dict[species], _ = first_element.split('-')
                _, cds_stop_query_dict[species]  = last_element.split('-')
            else:
                cds_start_query_dict[species] = None
                cds_stop_query_dict[species] = None

        if self.verbosity > 1:
            print("The __get_coordinates_from_maf function returns the following values:")
            print(cds_start_query_dict)
            print(cds_stop_query_dict)
            print(strand_query_dict)
            print(chr_query_dict)
            print("\n")

        return cds_start_query_dict, cds_stop_query_dict, strand_query_dict, chr_query_dict
