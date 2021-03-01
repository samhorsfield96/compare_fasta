def search_valid(genome_fasta, query_file, outfile):
	from Bio import SeqIO
	remaining_valid = []
	genome_seq = []
	for rec in SeqIO.parse(genome_fasta, "fasta"):
		genome_seq.append(str(rec.seq))
	
	for line in query_file:
		found = false
		seq = line.strip("\n")
		for genome in genome_seq:
			index = genome.find(seq)
			if index == -1:
				rev_seq = str(Seq(seq).reverse_complement())
				index = genome.find(rev_seq)
			if index != -1:
				found = true
				break
		if found:		
			remaining_valid.append(seq)
	
	with open(outfile, "w") as f:
		for i in remaining_valid:
			f.write(i + "\n")	
