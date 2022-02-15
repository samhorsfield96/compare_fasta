from Bio import SeqIO
import argparse

def get_options():
    description = "Compares pangenome fasta from ggCaller to known genes"
    parser = argparse.ArgumentParser(description=description,
                                     prog='python compare_fasta.py')

    IO = parser.add_argument_group('Input/options.out')
    IO.add_argument('--genes',
                    help='Reference gene panel (FASTA format)')
    IO.add_argument('--query',
                    help='Genes generated by gene caller to query (FASTA format)')
    return parser.parse_args()

def find_extra(ref_fasta, query_fasta, aa=False):
    ref_calls = set()
    # read in reference set of calls
    with open(ref_fasta, "r") as f:
        for record in SeqIO.parse(f, "fasta"):
            ref_calls.add(str(record.seq))

    # read in queries and translate
    query_dict = {}
    if aa == True:
        with open(query_fasta, "r") as f:
            for record in SeqIO.parse(f, "fasta"):
                query_dict[str(record.id)] = str((record.seq).translate())
    else:
        with open(query_fasta, "r") as f:
            for record in SeqIO.parse(f, "fasta"):
                query_dict[str(record.id)] = str(record.seq)

    # find calls not present in ref
    found = set()
    for id, seq in query_dict.items():
        if seq in ref_calls:
            found.add(id)

    # remove all found from query_dict
    for id in found:
        del query_dict[id]

    return query_dict.keys()

def main():
    options = get_options()
    output = find_extra(options.genes, options.query)
    return 0

if __name__ == '__main__':
    main()



