import sys, os
import argparse
import utils

from refExon import ReferenceExon
from queryExon import QueryExonSeqs
from alignmentClass import Alignment
from alignedSeqs import AlignedSeqPair

if __name__ == "__main__":

    ## argument parsing business
    parser = argparse.ArgumentParser(description = "CESAR multiexon annotation engine")
    parser.add_argument("input", help = 'Transcript/Gene identifier', type = str)
    parser.add_argument("maf_index", help = 'The path to the mafIndex file', type = str)
    parser.add_argument("gene_pred", help = 'The master genepred file that contains information about all the genes/transcripts', type = str)
    parser.add_argument("reference", help = 'Reference species', type = str)
    parser.add_argument("species_list", help = 'A file containing a list of query species', type = str)
    parser.add_argument("two_bit_dir", help = 'Path to two bit diectory that contains two bit files for every query species', type = str)
    parser.add_argument("cesar_dir", help = 'Path to the directory that contains CESAR binary', type = str)
    parser.add_argument("output_fasta", help = 'Path to the directory that contains CESAR generated alignments for every input transcript/gene', type = str)
    parser.add_argument("output_annotation", help = 'Path to the directory that contains the annotation output', type = str)
    parser.add_argument("-l", "--list", help = 'If the input is a list', action = 'store_true') ## if the input is a list of genes
    parser.add_argument('--verbosity', type = int, default = 0,  choices = [0, 1, 2], help=('Level of details to print out along with the progress.'))

    args = parser.parse_args()
    input_gene       = args.input
    maf_index        = args.maf_index
    gene_pred        = args.gene_pred
    reference        = args.reference
    species_list     = args.species_list
    two_bit_dir      = args.two_bit_dir
    cesar_dir        = args.cesar_dir
    output_fasta_dir = args.output_fasta
    output_anno_dir  = args.output_annotation
    input_list       = args.list
    verbosity        = args.verbosity

    clade = 'human'
######################################################################
##############           Do some sanity checks          ##############
######################################################################

    assert os.path.isfile(maf_index), f'Aborting!! The maf_index file {maf_index} does not exist'
    assert os.path.isfile(gene_pred), f'Aborting!! The gene_pred file {gene_pred} does not exist'
    assert os.path.isdir(two_bit_dir), f'Aborting!! The two_bit_dir {two_bit_dir} does not exist'
    assert os.path.isdir(cesar_dir), f'Aborting!! The CESAR binary directory {cesar_dir} does not exist'

    os.makedirs(output_fasta_dir, exist_ok = True)
    os.makedirs(output_anno_dir, exist_ok = True)

##############################################################################################################################
    ## check if two bit files are present for every query species:
    with open(species_list, 'r') as FI:
        query_species_list = FI.readlines()
    query_species_list = [species.rstrip('\n') for species in query_species_list]

    for species in query_species_list:
        two_bit_file = os.path.join(two_bit_dir, species, f'{species}.2bit')
        assert os.path.isfile(two_bit_file), f'Aborting, the 2bit file for {species} i.e {two_bit_file} is missing'

    two_bit_files_dict = {species:os.path.join(two_bit_dir, species, f'{species}.2bit') for species in query_species_list}    
    two_bit_file_ref = os.path.join(two_bit_dir, reference, f'{reference}.2bit')
##############################################################################################################################

    ## if the input is a file that contains a list of transcripts, read this file into a list
    genes_list = list()
    if input_list:
        with open(input_gene, 'r') as FI:
            genes_list = [gene.rstrip('\n') for gene in FI.readlines()]
    else:
        genes_list.append(input_gene)

    transcripts_info_dict = dict()
    with open(gene_pred, 'r') as FI:
        for line in FI:
            line_list = line.rstrip('\n').split('\t')
            transcript_id = line_list.pop(0)
            if transcript_id in genes_list:
                transcripts_info_dict[transcript_id] = '\t'.join(line_list)

    if len([gene for gene in genes_list if gene in transcripts_info_dict]) != len(genes_list):
        print(f'Aborting. Not all the transcripts from the input list are found in the input gene pred file {gene_pred}')
        for gene in genes_list:
            if gene not in transcripts_info_dict:
                print(f'The following transcript {gene} is not present in the master annotation file {gene_pred}')
        sys.exit()
        
    if verbosity > 1:
        print("Sanity checks done")
        print(f"The following gene(s) will be subjected to this realignment procedure {genes_list}.")
        print(f"The gene info dictionary looks like the following: {transcripts_info_dict}")
        print(f"The two bit directory is {two_bit_dir}, the reference species is {reference}")

######################################################################
############## sanity checks passed. Time for business  ##############
######################################################################

    ## create all temp files upfront and store them in temp_files_list
    maf_file     = utils.create_temp_devshm()
    cesar_input  = utils.create_temp_devshm()
    temp_files_list = [maf_file, cesar_input]
######################################################################

    for transcript_id, transcript_info in transcripts_info_dict.items():

        output_fasta     = os.path.join(output_fasta_dir, f'{transcript_id}.fa')
        output_anno_file = os.path.join(output_anno_dir, f'{transcript_id}.txt')

        chromosome, strand, *_, exons_list_str = transcript_info.split('\t')    
        exons_list = exons_list_str.split(',')
        number_of_exons = len(exons_list)
        if strand == '-':
            exons_list = exons_list[::-1]

        ref_exon = ReferenceExon(exons_list, chromosome, strand, os.path.join(two_bit_dir, reference, f'{reference}.2bit'), verbosity)
        exon_seqs, ss_dict, exon_phases_list = ref_exon.extract_seqs_splice_sites()

        if verbosity > 1:
            print("The following object is created for the reference_exon: BEGIN #####\n")
            print(f'{vars(ref_exon)}\n{exon_seqs}\n{ss_dict}')
            print("###### END\n\n")

        query_obj = QueryExonSeqs(reference, exons_list, chromosome, strand, two_bit_dir, query_species_list, maf_index, verbosity)
        query_cds_dict, valid_seqs_ct = query_obj.get_query_seqs()

        if verbosity > 1:
            print("The following is the information for the query_sequence_object: BEGIN #####\n")
            print(f'{vars(query_obj)}\n{query_cds_dict}')
            print("###### END\n\n")

        if valid_seqs_ct == 0:
            continue

        ######################################################################
        ##############     Prepare the input file for CESAR     ##############
        ######################################################################
        if verbosity > 0:
            print("Preparing the input for CESAR : running mafExtract and then twoBitToFa on the query species\n")

        with open(cesar_input, 'w') as FO:
            for exon_id, exon_seq in exon_seqs.items():
                _, exon_ct = exon_id.split('_')
                exon_ct = int(exon_ct)

                acc_profile_def = os.path.join("extra",  "tables", f'{clade}', "acc_profile.txt")
                do_profile_def  = os.path.join("extra",  "tables", f'{clade}', "do_profile.txt")
                acc_profile = os.path.join("extra",  "tables", f'{clade}', "firstCodon_profile.txt") if exon_ct == 0 else acc_profile_def
                do_profile = os.path.join("extra",  "tables", f'{clade}', "lastCodon_profile.txt") if exon_ct == number_of_exons - 1 else do_profile_def

                if ss_dict["acc"][exon_id] == "AC" and exon_ct != 0:
                    acc_profile = os.path.join("extra",  "tables", f'{clade}', "u12_acc_profile.txt")
                if ss_dict["donor"][exon_id] == "AT" and exon_ct != len(exons_list) - 1:
                    do_profile = os.path.join("extra",  "tables", f'{clade}', "u12_donor_profile.txt")

                FO.write(f">{transcript_id}_{exon_ct}\t{acc_profile}\t{do_profile}\n{exon_seq}\n")
            ######################################################
            ####  now populate this file with query sequences ####
            ######################################################

            FO.write("####\n")
            species_for_alignment = 0
            for query_species, query_species_cds in query_cds_dict.items():
                _,  cds_start, cds_stop, query_strand, query_chr = query_species_cds.split('#') #5235109#5202543#-#JH863841

                cds_start = int(cds_start) - 500
                cds_stop  = int(cds_stop) + 500

                if cds_start > cds_stop:
                    if verbosity > 1:
                        print(f'Omitting the species {query_species} because the cds_start is greater than the cds_stop {cds_start} > {cds_stop}')
                        continue
                    
                query_sequence = utils.get_sequence_from_2bit(cds_start, cds_stop, query_chr, two_bit_files_dict[query_species])
                if query_strand != strand:
                    query_sequence = utils.reverse_comp(query_sequence)
                FO.write(f">{query_species_cds}\n{query_sequence}\n")
                species_for_alignment += 1

        ######################################################################
        ##############               Run CESAR                  ##############
        ######################################################################
        if species_for_alignment == 0:
            print("Nothing to run CESAR on!! All query species are problematic - either not on the same scaffold/chromosome or the query locus has some issues with the genome assembly.")
            continue

        if verbosity > 0:
            print('Running CESAR now!!!')
        cesar_call = f'cesar {cesar_input} -x 32'
        cesar_output = utils.run_subprocess(cesar_call, True)

        if verbosity > 0:
            print(f"CESAR run finished. This was the call {cesar_call}\n")

        ### process CESAR's output
        alignments_obj = Alignment(cesar_output, output_fasta, reference)
        alignments_split_exon_wise = alignments_obj.process_split_alignments()

        ######################################################################
        ###########           CESAR output processing              ###########
        ######################################################################
        ## multiple-alignment file for RELAX/Hyphy has already been generated and stored in output_fasta  ##
        ## now annotate exons or check for exon intactness ##

        print("Now looking into the annotation business")
        with open(output_anno_file, 'w') as FO:
            for species, alignments_dict in alignments_split_exon_wise.items():
                for exon_index, alignment_block in alignments_dict.items():
                    ref_acceptor = ss_dict["acc"][f'exon_{exon_index}']
                    ref_donor    = ss_dict["donor"][f'exon_{exon_index}']

                    aligned_ref, query_header, aligned_query = alignment_block.split('\n')
                    phase_string = exon_phases_list[exon_index]
                    aligned_ref_query_pair = AlignedSeqPair(aligned_ref, aligned_query, query_header, phase_string, ref_acceptor, ref_donor, exon_index, number_of_exons, strand, verbosity)
                    exon_cds = aligned_ref_query_pair.check_exon_intactness()
                    FO.write(f'{transcript_id}\t{species}\t{exon_index}\t{exon_cds}\n')
                                
    for f in temp_files_list:
        os.remove(f)
