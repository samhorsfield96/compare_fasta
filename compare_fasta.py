from Bio import SeqIO

# used for analysing old ggCaller output (colours held in square brackets)
def compare_3prime_old(genome_fasta, ref_fasta, query_fasta, caller_type, min_size):
    genome_list = []
    genome_rec = {}
    ref_rec = {}
    query_rec = {}

    total_correct_query_records = 0
    total_ref_records = 0
    total_query_records = 0

    unmatched_query_list = []
    unmatched_ref_list = []
    incorrect_query_list = []

    # parse genome_fasta
    for rec in SeqIO.parse(genome_fasta, "fasta"):
        description = (rec.description).split("_")
        id = description[0]
        genome_list.append(id)
        genome_rec[id] = str(rec.seq)
        ref_rec[id] = {}
        query_rec[id] = {}

    # parse ref_fasta
    for rec in SeqIO.parse(ref_fasta, "fasta"):
        description = (rec.description).split("_")
        id = description[1]

        if len(str(rec.seq)) < min_size:
           continue


        # look for the 3prime index of the string
        start_index = genome_rec[id].find(str(rec.seq))
        if start_index == -1:
            prime3 = genome_rec[id].find(str(rec.seq.reverse_complement()))
        else:
            prime3 = start_index + (len(str(rec.seq)) - 1)
        #check that sequence is present in genome is says, and that the gene sequence is valid and no Ns present. Also check that no duplicate entries, as ggc can't call duplicate genes
        if prime3 != -1 and remove_invalid(str(rec.seq)) and "N" not in str(rec.seq) and prime3 not in ref_rec[id]:
            ref_rec[id][prime3] = str(rec.seq)
            total_ref_records += 1
            unmatched_ref_list.append((id, str(rec.seq)))

    if caller_type == "ggc":
        # parse query_fasta
        for rec in SeqIO.parse(query_fasta, "fasta"):
            colours = (((((rec.description.strip()).split("["))[4]).replace("]", "")).replace("'", "")).replace(", ", "")
            colours = list(colours)

            for index, col in enumerate(colours):
                if col == "1":
                    id = genome_list[index]
                    # look for the 3prime index of the string
                    start_index = genome_rec[id].find(str(rec.seq))
                    if start_index == -1:
                        prime3 = genome_rec[id].find(str(rec.seq.reverse_complement()))
                    else:
                        prime3 = start_index + (len(str(rec.seq)) - 1)

                    # determine if gene sequence is real, if not add to incorrect_query_list
                    if prime3 != -1:
                        query_rec[id][prime3] = str(rec.seq)
                    else:
                        incorrect_query_list.append((id, str(rec.seq)))
                    total_query_records += 1

    elif caller_type == "prod":
        # parse query_fasta
        for rec in SeqIO.parse(query_fasta, "fasta"):
            description = (rec.description).split("_")
            id = description[0]

            # look for the 3prime index of the string
            start_index = genome_rec[id].find(str(rec.seq))
            if start_index == -1:
                prime3 = genome_rec[id].find(str(rec.seq.reverse_complement()))
            else:
                prime3 = start_index + (len(str(rec.seq)) - 1)
            if remove_invalid(str(rec.seq)) and "N" not in str(rec.seq) and prime3 not in query_rec[id]:
                query_rec[id][prime3] = str(rec.seq)
            else:
                incorrect_query_list.append((id, str(rec.seq)))
            total_query_records += 1

    # iterate over query_rec, count number of times each 3prime match found
    for colour, prime3_dict in query_rec.items():
        for prime3_key in prime3_dict.keys():
            if prime3_key in ref_rec[colour]:
                total_correct_query_records += 1
                unmatched_ref_list.remove((colour, ref_rec[colour][prime3_key]))
            else:
                unmatched = (colour, prime3_dict[prime3_key])
                unmatched_query_list.append(unmatched)


    recall = total_correct_query_records / total_ref_records
    precision = total_correct_query_records / total_query_records

    print("Total ORFs: {}".format(total_query_records))
    print("Recall: {}".format(recall))
    print("Precision: {}".format(precision))
    return(unmatched_ref_list, unmatched_query_list, incorrect_query_list)

# used for analysing new ggCaller output (colours binary array)
def compare_3prime_new(genome_fasta, ref_fasta, query_fasta, caller_type, min_size):
    genome_list = []
    genome_rec = {}
    ref_rec = {}
    query_rec = {}

    total_correct_query_records = 0
    total_ref_records = 0
    total_query_records = 0

    unmatched_query_list = []
    unmatched_ref_list = []
    incorrect_query_list = []

    # parse genome_fasta
    for rec in SeqIO.parse(genome_fasta, "fasta"):
        description = (rec.description).split("_")
        id = description[0]
        genome_list.append(id)
        genome_rec[id] = str(rec.seq)
        ref_rec[id] = {}
        query_rec[id] = {}

    # parse ref_fasta
    for rec in SeqIO.parse(ref_fasta, "fasta"):
        description = (rec.description).split("_")
        id = description[1]

        if len(str(rec.seq)) < min_size:
           continue


        # look for the 3prime index of the string
        start_index = genome_rec[id].find(str(rec.seq))
        if start_index == -1:
            prime3 = genome_rec[id].find(str(rec.seq.reverse_complement()))
        else:
            prime3 = start_index + (len(str(rec.seq)) - 1)
        #check that sequence is present in genome is says, and that the gene sequence is valid and no Ns present. Also check that no duplicate entries, as ggc can't call duplicate genes
        if prime3 != -1 and remove_invalid(str(rec.seq)) and "N" not in str(rec.seq) and prime3 not in ref_rec[id]:
            ref_rec[id][prime3] = str(rec.seq)
            total_ref_records += 1
            unmatched_ref_list.append((id, str(rec.seq)))

    if caller_type == "ggc":
        # parse query_fasta
        for rec in SeqIO.parse(query_fasta, "fasta"):
            description = (rec.description).split("_")
            colours = description[1]

            for index, col in enumerate(colours):
                if col == "1":
                    id = genome_list[index]
                    # look for the 3prime index of the string
                    start_index = genome_rec[id].find(str(rec.seq))
                    if start_index == -1:
                        prime3 = genome_rec[id].find(str(rec.seq.reverse_complement()))
                    else:
                        prime3 = start_index + (len(str(rec.seq)) - 1)

                    # determine if gene sequence is real, if not add to incorrect_query_list
                    if prime3 != -1:
                        query_rec[id][prime3] = str(rec.seq)
                    else:
                        incorrect_query_list.append((id, str(rec.seq)))
                    total_query_records += 1

    elif caller_type == "prod":
        # parse query_fasta
        for rec in SeqIO.parse(query_fasta, "fasta"):
            description = (rec.description).split("_")
            id = description[0]

            # look for the 3prime index of the string
            start_index = genome_rec[id].find(str(rec.seq))
            if start_index == -1:
                prime3 = genome_rec[id].find(str(rec.seq.reverse_complement()))
            else:
                prime3 = start_index + (len(str(rec.seq)) - 1)
            if remove_invalid(str(rec.seq)) and "N" not in str(rec.seq) and prime3 not in query_rec[id]:
                query_rec[id][prime3] = str(rec.seq)
            else:
                incorrect_query_list.append((id, str(rec.seq)))
            total_query_records += 1

    # iterate over query_rec, count number of times each 3prime match found
    for colour, prime3_dict in query_rec.items():
        for prime3_key in prime3_dict.keys():
            if prime3_key in ref_rec[colour]:
                total_correct_query_records += 1
                unmatched_ref_list.remove((colour, ref_rec[colour][prime3_key]))
            else:
                unmatched = (colour, prime3_dict[prime3_key])
                unmatched_query_list.append(unmatched)


    recall = total_correct_query_records / total_ref_records
    precision = total_correct_query_records / total_query_records

    print("Total ORFs: {}".format(total_query_records))
    print("Recall: {}".format(recall))
    print("Precision: {}".format(precision))
    return(unmatched_ref_list, unmatched_query_list, incorrect_query_list)

def remove_invalid(query_seq):
    import re
    start_codons = ["ATG", "GTG", "TTG"]
    stop_codons = ["TAA", "TGA", "TAG"]
    # check if codons in right place, and in correct frame
    if query_seq[0:3] not in start_codons or query_seq[-3:] not in stop_codons or len(query_seq) % 3 != 0:
        return False

    # check if sequence contains a premature stop codon
    for stop in stop_codons:
        stop_indices = [m.start() for m in re.finditer(stop, query_seq)]
        # remove last stop as this is valid
        if query_seq[-3:] == stop:
            stop_indices = stop_indices[:-1]
        if any([i % 3 == 0 for i in stop_indices]):
            return False

    # If all tests come back fine, return true
    return True

def select_seq_length(infasta, outfasta, length):
    with open(outfasta, "w") as f:
        for rec in SeqIO.parse(infasta, "fasta"):
            if len(rec.seq) >= length:
                f.write(">" + str(rec.description) + "\n" + str(rec.seq) + "\n")


if __name__ == '__main__':
    from Bio import SeqIO
    #unmatched_ref, unmatched_query = compare_3prime("clique_556_seqs_all.fasta", "clique_556_CDS_all.fasta", "clique556_calls_4threads_post_unitigDict_alteration_3.fasta", "ggc", 90)
    unmatched_ref_post, unmatched_query_post = compare_3prime_old("clique_556_seqs_all.fasta", "clique_556_CDS_all.fasta",
                                                              "clique_556_list_chunking_4threads_filtered.fasta", "ggc",
                                                              90)
