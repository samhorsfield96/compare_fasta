from Bio import SeqIO

def compare_fasta(ref_fasta, query_fasta, type):
    total_ref_records = 0
    total_query_records = 0
    total_correct_query_records = 0

    ref_seq_list = []

    #create list of all reference genes
    for ref_rec in SeqIO.parse(ref_fasta, "fasta"):
        total_ref_records += 1
        ref_seq_list.append(str(ref_rec.seq))

    unmatched_seq_list = ref_seq_list[:]

    #iterate through query, checking if entry is in ref_seq_list, depending on whether running prodigal or ggCaller due to gene merging in ggCaller
    if type == 'ggc':
        for ggc_rec in SeqIO.parse(query_fasta, "fasta"):
            total_query_records += 1
            total_correct_query_records += ref_seq_list.count(str(ggc_rec.seq))
            unmatched_seq_list = list(filter(lambda a: a != str(ggc_rec.seq), unmatched_seq_list))
    elif type == 'prod':
        for prod_rec in SeqIO.parse(query_fasta, "fasta"):
            total_query_records += 1
            if str(prod_rec.seq) in unmatched_seq_list:
                total_correct_query_records += 1
                unmatched_seq_list.remove(str(prod_rec.seq))



    #calculate recall and precision
    recall = total_correct_query_records / total_ref_records
    precision = total_correct_query_records / total_query_records

    print(query_fasta)
    print("Recall: {}".format(recall))
    print("Precision: {}".format(precision))
    return(unmatched_seq_list)

def compare_fasta_exact(ref_fasta, query_fasta, type, group):
    total_ref_records = 0
    total_query_records = 0
    total_correct_query_records = 0

    ref_seq_list = []

    if group == "1":
        group_list = ["CR931645.1", "CR931648.1", "CR931646.1", "CR931647.1"]
    elif group == "2":
        group_list = ["CR931719.1", "CR931658.1", "CR931659.1", "CR931660.1", "CR931717.1"]
    elif group == "3":
        group_list = ["CR931662.1", "CR931663.1", "CR931664.1", "CR931665.1", "CR931666.1"]

    #create list of all reference genes
    for ref_rec in SeqIO.parse(ref_fasta, "fasta"):
        source_id = ((ref_rec.description.strip()).split("_"))[1]
        total_ref_records += 1
        ref_seq_list.append((source_id, str(ref_rec.seq)))

    unmatched_seq_list = ref_seq_list[:]

    #iterate through query, checking if entry is in ref_seq_list, depending on whether running prodigal or ggCaller due to gene merging in ggCaller
    if type == 'ggc':
        for ggc_rec in SeqIO.parse(query_fasta, "fasta"):
            ORF_colours = (((((ggc_rec.description.strip()).split("["))[4]).replace("]", "")).replace("'", "")).replace(", ", "")
            ORF_colours = list(ORF_colours)
            for i in range(0, len(ORF_colours)):
                if ORF_colours[i] == "1":
                    ORF_pair = (group_list[i], str(ggc_rec.seq))
                    total_query_records += 1
                    if ORF_pair in unmatched_seq_list:
                        total_correct_query_records += 1
                        unmatched_seq_list.remove(ORF_pair)
    elif type == 'prod':
        for prod_rec in SeqIO.parse(query_fasta, "fasta"):
            ORF_id = ((prod_rec.description.strip()).split("_"))[0] + ".1"
            ORF_pair = (ORF_id, str(prod_rec.seq))
            total_query_records += 1
            if ORF_pair in unmatched_seq_list:
                total_correct_query_records += 1
                unmatched_seq_list.remove(ORF_pair)



    #calculate recall and precision
    recall = total_correct_query_records / total_ref_records
    precision = total_correct_query_records / total_query_records

    print(query_fasta)
    print("Total ORFs: {}".format(total_query_records))
    print("Recall: {}".format(recall))
    print("Precision: {}".format(precision))
    return(unmatched_seq_list)

def compare_fasta_py_cpp(ref_fasta, query_fasta, type, group):
    total_ref_records = 0
    total_query_records = 0
    total_correct_query_records = 0

    ref_seq_list = []

    if group == "1":
        group_list = ["CR931645.1", "CR931648.1", "CR931646.1", "CR931647.1"]
    elif group == "2":
        group_list = ["CR931719.1", "CR931658.1", "CR931659.1", "CR931660.1", "CR931717.1"]
    elif group == "3":
        group_list = ["CR931662.1", "CR931663.1", "CR931664.1", "CR931665.1", "CR931666.1"]

    #create list of all reference genes
    for py_rec in SeqIO.parse(ref_fasta, "fasta"):
        ORF_colours = (((((py_rec.description.strip()).split("["))[4]).replace("]", "")).replace("'", "")).replace(
            ", ", "")
        ORF_colours = list(ORF_colours)
        for i in range(0, len(ORF_colours)):
            if ORF_colours[i] == "1":
                ORF_pair = (group_list[i], str(py_rec.seq))
                total_ref_records += 1
                ref_seq_list.append(ORF_pair)

    unmatched_seq_list = ref_seq_list[:]

    #iterate through query, checking if entry is in ref_seq_list, depending on whether running prodigal or ggCaller due to gene merging in ggCaller
    if type == 'ggc':
        for ggc_rec in SeqIO.parse(query_fasta, "fasta"):
            ORF_description = (ggc_rec.description).split("_")
            ORF_colours = list(ORF_description[2])
            for i in range(0, len(ORF_colours)):
                if ORF_colours[i] == "1":
                    ORF_pair = (group_list[i], str(ggc_rec.seq))
                    total_query_records += 1
                    if ORF_pair in unmatched_seq_list:
                        total_correct_query_records += 1
                        unmatched_seq_list.remove(ORF_pair)
    elif type == 'prod':
        for prod_rec in SeqIO.parse(query_fasta, "fasta"):
            ORF_id = ((prod_rec.description.strip()).split("_"))[0] + ".1"
            ORF_pair = (ORF_id, str(prod_rec.seq))
            total_query_records += 1
            if ORF_pair in unmatched_seq_list:
                total_correct_query_records += 1
                unmatched_seq_list.remove(ORF_pair)



    #calculate recall and precision
    recall = total_correct_query_records / total_ref_records
    precision = total_correct_query_records / total_query_records

    print(query_fasta)
    print("Total ORFs: {}".format(total_query_records))
    print("Recall: {}".format(recall))
    print("Precision: {}".format(precision))
    return(unmatched_seq_list)

def compare_fasta_cpp_cpp(ref_fasta, query_fasta, group):
    total_ref_records = 0
    total_query_records = 0
    total_correct_query_records = 0

    ref_seq_list = []

    if group == "1":
        group_list = ["CR931645.1", "CR931648.1", "CR931646.1", "CR931647.1"]
    elif group == "2":
        group_list = ["CR931719.1", "CR931658.1", "CR931659.1", "CR931660.1", "CR931717.1"]
    elif group == "3":
        group_list = ["CR931662.1", "CR931663.1", "CR931664.1", "CR931665.1", "CR931666.1"]

    #create list of all reference genes
    for ref_rec in SeqIO.parse(ref_fasta, "fasta"):
        ORF_description = (ref_rec.description).split("_")
        ORF_colours = list(ORF_description[2])
        for i in range(0, len(ORF_colours)):
            if ORF_colours[i] == "1":
                ORF_pair = (group_list[i], str(ref_rec.seq))
                total_ref_records += 1
                ref_seq_list.append(ORF_pair)

    unmatched_seq_list = ref_seq_list[:]

    #iterate through query, checking if entry is in ref_seq_list, depending on whether running prodigal or ggCaller due to gene merging in ggCaller
    for query_rec in SeqIO.parse(query_fasta, "fasta"):
        ORF_description = (query_rec.description).split("_")
        ORF_colours = list(ORF_description[2])
        for i in range(0, len(ORF_colours)):
            if ORF_colours[i] == "1":
                ORF_pair = (group_list[i], str(query_rec.seq))
                total_query_records += 1
                if ORF_pair in unmatched_seq_list:
                    total_correct_query_records += 1
                    unmatched_seq_list.remove(ORF_pair)

    #calculate recall and precision
    recall = total_correct_query_records / total_ref_records
    precision = total_correct_query_records / total_query_records

    print(query_fasta)
    print("Total ORFs: {}".format(total_query_records))
    print("Recall: {}".format(recall))
    print("Precision: {}".format(precision))
    return(unmatched_seq_list)

def compare_fasta_exact_cpp(ref_fasta, query_fasta, type, group):
    total_ref_records = 0
    total_query_records = 0
    total_correct_query_records = 0

    ref_seq_list = []

    if group == "1":
        group_list = ["CR931645.1", "CR931648.1", "CR931646.1", "CR931647.1"]
    elif group == "2":
        group_list = ["CR931719.1", "CR931658.1", "CR931659.1", "CR931660.1", "CR931717.1"]
    elif group == "3":
        group_list = ["CR931662.1", "CR931663.1", "CR931664.1", "CR931665.1", "CR931666.1"]

    #create list of all reference genes
    for ref_rec in SeqIO.parse(ref_fasta, "fasta"):
        source_id = ((ref_rec.description.strip()).split("_"))[1]
        total_ref_records += 1
        ref_seq_list.append((source_id, str(ref_rec.seq)))

    unmatched_seq_list = ref_seq_list[:]

    #iterate through query, checking if entry is in ref_seq_list, depending on whether running prodigal or ggCaller due to gene merging in ggCaller
    if type == 'ggc':
        for ggc_rec in SeqIO.parse(query_fasta, "fasta"):
            ORF_description = (ggc_rec.description).split("_")
            ORF_colours = list(ORF_description[2])
            for i in range(0, len(ORF_colours)):
                if ORF_colours[i] == "1":
                    ORF_pair = (group_list[i], str(ggc_rec.seq))
                    total_query_records += 1
                    if ORF_pair in unmatched_seq_list:
                        total_correct_query_records += 1
                        unmatched_seq_list.remove(ORF_pair)
    elif type == 'prod':
        for prod_rec in SeqIO.parse(query_fasta, "fasta"):
            ORF_id = ((prod_rec.description.strip()).split("_"))[0] + ".1"
            ORF_pair = (ORF_id, str(prod_rec.seq))
            total_query_records += 1
            if ORF_pair in unmatched_seq_list:
                total_correct_query_records += 1
                unmatched_seq_list.remove(ORF_pair)



    #calculate recall and precision
    recall = total_correct_query_records / total_ref_records
    precision = total_correct_query_records / total_query_records

    print(query_fasta)
    print("Total ORFs: {}".format(total_query_records))
    print("Recall: {}".format(recall))
    print("Precision: {}".format(precision))
    return(unmatched_seq_list)

def compare_3prime(genome_fasta, ref_fasta, query_fasta, caller_type, min_size):
    genome_list = []
    genome_rec = {}
    ref_rec = {}
    query_rec = {}

    total_correct_query_records = 0
    total_ref_records = 0
    total_query_records = 0

    unmatched_query_list = []
    unmatched_ref_list = []


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
        id = description[0]

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
                    query_rec[id][prime3] = str(rec.seq)
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
            query_rec[id][prime3] = str(rec.seq)
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
    return(unmatched_ref_list, unmatched_query_list)

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
    unmatched_ref, unmatched_query = compare_3prime("clique_556_seqs_all.fasta", "clique_556_CDS_all.fasta", "clique556_calls_4threads_post_unitigDict_alteration_3.fasta", "ggc", 90)
