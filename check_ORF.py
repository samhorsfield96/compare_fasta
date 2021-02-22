def check_ORF(ref_fasta_for, ref_fasta_rev, query_fasta, outfasta):
    from Bio import SeqIO
    ref_seq_list = []
    for ref_rec in SeqIO.parse(ref_fasta_for, "fasta"):
        ref_seq_list.append(str(ref_rec.seq))
    for ref_rec in SeqIO.parse(ref_fasta_rev, "fasta"):
        ref_seq_list.append(str(ref_rec.seq))
    with open(outfasta, "w") as f:
        for query_rec in SeqIO.parse(query_fasta, "fasta"):
            if not any(str(query_rec.seq) in ref_seq for ref_seq in ref_seq_list):
                f.write(">" + str(query_rec.description) + "\n" + str(query_rec.seq) + "\n")

def check_ORF_colours(ref_fasta_for, ref_fasta_rev, query_fasta, outfasta):
    from Bio import SeqIO
    ref_seq_list_pos = []
    ref_seq_list_neg = []
    for ref_rec in SeqIO.parse(ref_fasta_for, "fasta"):
        ref_seq_list_pos.append(str(ref_rec.seq))
    for ref_rec in SeqIO.parse(ref_fasta_rev, "fasta"):
        ref_seq_list_neg.append(str(ref_rec.seq))
    with open(outfasta, "w") as f:
        for query_rec in SeqIO.parse(query_fasta, "fasta"):
            pos_list = []
            neg_list = []
            for pos_ref_seq in ref_seq_list_pos:
                if str(query_rec.seq) in pos_ref_seq:
                    pos_list.append('1')
                else:
                    pos_list.append('0')
            for neg_ref_seq in ref_seq_list_neg:
                if str(query_rec.seq) in neg_ref_seq:
                    neg_list.append('1')
                else:
                    neg_list.append('0')
            f.write(str(query_rec.description) + " Positive matrix: " + str(pos_list) + " Negative matrix: " + str(neg_list) + "\n")


def check_kmer_pairs(ref_fasta, query_fasta, ksize, strand, outfile):
    from Bio import SeqIO
    mismatch_dict = {}
    for ref_rec in SeqIO.parse(ref_fasta, "fasta"):
        mismatch_dict[ref_rec.id] = {}
        kmers = []
        mismatch_dict[ref_rec.id]['kmer_pairs'] = []
        for i in range(len(str(ref_rec.seq)) - ksize + 1):
            kmer = (str(ref_rec.seq))[i:i + ksize]
            kmers.append(kmer)
        for i in range(len(kmers) - 1):
            pair = (kmers[i], kmers[i + 1])
            mismatch_dict[ref_rec.id]['kmer_pairs'].append(pair)

    for query_rec in SeqIO.parse(query_fasta, "fasta"):
        if strand in query_rec.description:
            kmers = []
            kmer_pairs = []
            for i in range(len(str(query_rec.seq)) - ksize + 1):
                kmer = (str(query_rec.seq))[i:i+ksize]
                kmers.append(kmer)

            for i in range(len(kmers) - 1):
                pair = (kmers[i], kmers[i + 1])
                kmer_pairs.append(pair)

            for ref in mismatch_dict:
                mismatch_dict[ref][query_rec.description] = []
                for pair in kmer_pairs:
                    if pair not in mismatch_dict[ref]['kmer_pairs']:
                        mismatch_dict[ref][query_rec.description].append(pair)

    for ref in mismatch_dict:
        del mismatch_dict[ref]['kmer_pairs']

    with open(outfile, "w") as f:
        for key, item in mismatch_dict.items():
            f.write(str(key) + "\n" + str(item) + "\n")

def check_ref_in_query(ref_fasta, query_fasta, outfile):
    from Bio import SeqIO
    query_list = []
    for query_rec in SeqIO.parse(query_fasta, "fasta"):
        query_list.append(str(query_rec.seq))

    with open(outfile, "w") as f:
        for ref_rec in SeqIO.parse(ref_fasta, "fasta"):
            if not any(str(ref_rec.seq) in query_rec for query_rec in query_list):
                f.write(">" + str(ref_rec.description) + "\n" + str(ref_rec.seq) + "\n")