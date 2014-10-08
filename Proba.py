from Bio import SeqIO

for seq_record in SeqIO.parse("/home/mare/Preuzimanja/ERR031114.fa", "fasta"):
    print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record))