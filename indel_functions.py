import utils
import sys

def get_ancestral_reading_frame(start_pos, codon_pos_dict_ref):

    species_rf = ''
    while start_pos >= 0: ## This is the position in the sequence from where I ancestralize the
        codon_pos = codon_pos_dict_ref[start_pos] ## reading frame of the species i.e. a positon just upstream of the indel

        if codon_pos != '-':       ## Get the line where there is a character in the Reference species. Remember we need
            species_rf = codon_pos ## the position where there is a character in the Reference, not in the query species
            break                  ## After all, the function returns the ancestral reading frame
        start_pos -=1

    if species_rf == '':
        return 0
    return species_rf

def check_for_stop_codon(gene_seq, reading_frame):

    gene_seq = gene_seq.replace('-','')
    
    if reading_frame == 2:
        gene_seq = "N" + gene_seq
    elif reading_frame == 3:
        gene_seq = "NN" + gene_seq

    stop_codon = utils.find_stop_codon(gene_seq)
    return stop_codon

def get_compensated_events(all_fs_dict, codon_pos_dict_ref, codon_pos_dict_query, query_sequence):

    all_fs_list = list(all_fs_dict.keys())
    compensation_dict, excluded_events = dict(), list()
    i = 0
    while i < len(all_fs_list):
        query_event = all_fs_list[i]
        if query_event in excluded_events:
            i += 1
            continue

        length_indel = all_fs_dict[query_event]
        length_cum_pot_comp_events = length_indel
        list_pot_comp_events    = list()

        j = i + 1
        while j < len(all_fs_list):
            pot_comp_event = all_fs_list[j]
            if pot_comp_event in excluded_events:
                j += 1
                continue

            length_pot_comp_event = all_fs_dict[pot_comp_event]
            length_cum_pot_comp_events += length_pot_comp_event
            list_pot_comp_events.append(pot_comp_event)

            if length_cum_pot_comp_events%3 == 0:
                start = query_event
                stop  = pot_comp_event + abs(length_pot_comp_event) - 1

                start_upstream = query_event - 1
                ancestral_rf = get_ancestral_reading_frame(start_upstream, codon_pos_dict_ref)

                if ancestral_rf == 0:
                    species_rf = 1
                else:
                    species_rf = ancestral_rf + 1
                    if species_rf == 4:
                        species_rf = 1


                stop_codon_present = check_for_stop_codon(query_sequence[start:stop], species_rf)
                if stop_codon_present == False:
                    excluded_events.extend(list_pot_comp_events)
                    compensation_dict[query_event] = list_pot_comp_events
                break
            j += 1
        i += 1

    ### Perform the "Return To Ancestral Reading Frame Test"
    for event, comp_events in dict(sorted(compensation_dict.items())).items():
        stop = comp_events[-1]

        start_upstream = event - 1
        species_rf = get_ancestral_reading_frame(start_upstream, codon_pos_dict_ref)
        if species_rf == 0: ## Only deletions at the beginning of the gene sequence are excluded from the "Return To Ancestral Reading Frame Test"
            continue

        start = event
        stop_pos = stop + abs(all_fs_dict[stop]) - 1
        reading_frame = ''

        while start <= stop_pos:
            if codon_pos_dict_ref[start] != '-':
                reading_frame = codon_pos_dict_ref[start]

            if codon_pos_dict_query[start] != '-':
                species_rf += 1
            start += 1

        species_rf_final = species_rf%3
        if species_rf_final == 0:
            species_rf_final = 3
        assert reading_frame == species_rf_final, f'Aborting!! error in the get_compensated_events function'

    return compensation_dict
