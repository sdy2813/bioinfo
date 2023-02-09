




def parse_fasta(file_path):
    with open(file_path, "r") as file:
        data = file.read().split(">")[1:]
    sequences = {}
    for sequence in data:
        lines = sequence.strip().split("\n")
        seqid = lines[0].split()[0]
        seq = "".join(lines[1:])
        sequences[seqid] = seq
    return sequences

  def extract_sequences(fasta_sequences, seqid_file_path):
    with open(seqid_file_path, "r") as file:
        seqids = file.read().strip().split("\n")
    extracted_sequences = {}
    for seqid in seqids:
        if seqid in fasta_sequences:
            extracted_sequences[seqid] = fasta_sequences[seqid]
    return extracted_sequences

 def write_sequences(sequences, output_file_path):
    with open(output_file_path, "w") as file:
        for seqid, seq in sequences.items():
            file.write(f">{seqid}\n{seq}\n")

fasta_sequences = parse_fasta("input.fasta")
extracted_sequences = extract_sequences(fasta_sequences, "seqids.txt")
write_sequences(extracted_sequences, "output.fasta")
