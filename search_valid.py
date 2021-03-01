def search_valid(genome_fasta, query_file, outfile):
	from Bio import SeqIO
	from Bio.Seq import Seq
	remaining_valid = []
	genome_seq = []
	for rec in SeqIO.parse(genome_fasta, "fasta"):
		genome_seq.append(str(rec.seq))
	
	with open(query_file, "r") as q:
		for line in q:
			found = False
			seq = line.strip("\n")
			for genome in genome_seq:
				index = genome.find(seq)
				if index == -1:
					rev_seq = str(Seq(seq).reverse_complement())
					index = genome.find(rev_seq)
				if index != -1:
					found = True
					break
			if found:
				remaining_valid.append(seq)
	
	with open(outfile, "w") as f:
		for i in remaining_valid:
			f.write(i + "\n")	


if __name__ == '__main__':
	search_valid("Pneumo_capsular_data/group3_forward.txt", "Pneumo_capsular_data/test.fasta", "test.txt")