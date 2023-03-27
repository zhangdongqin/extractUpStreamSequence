import argparse
import os
from Bio import SeqIO
from collections import defaultdict

def parse_gff_columns(columns, gene_aliases):
    attr_str = columns[8]
    attrs = attr_str.split(";")
    gene_id = None
    for attr in attrs[:3]:
        if "=" in attr:
            key, value = attr.strip().split("=")
        elif " " in attr:
            key, value = attr.strip().split(" ", 1)
        else:
            continue

        if key in ["gene_id", "ID", "transcript_id"]:
            candidate_gene_id = value.strip('"')
            if candidate_gene_id in gene_aliases:
                gene_id = candidate_gene_id
                break
    return gene_id



def extract_upstream_sequences(genome_file, gene_file, gff_file, output_folder, window, feature_type):
    os.makedirs(output_folder, exist_ok=True)
    gene_seq_output = os.path.join(output_folder, "gene_upStream.fasta")
    gene_bed_output = os.path.join(output_folder, "gene_upStream.bed")

    
    genome = SeqIO.to_dict(SeqIO.parse(genome_file, "fasta"))

    
    gene_aliases = {}
    with open(gene_file, 'r') as f:
        for line in f:
            gene_id, alias = line.strip().split("\t")
            gene_aliases[gene_id] = alias
    
    gene_info = defaultdict(dict)
    found_gene = False

    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            columns = line.strip().split("\t")
            current_feature_type = columns[2]
            if current_feature_type == feature_type or (not found_gene and current_feature_type == "mRNA"):
                gene_id = parse_gff_columns(columns, gene_aliases)
                if gene_id:
                    found_gene = True
                    gene_info[gene_id] = {
                        "seq_id": columns[0],
                        "start": int(columns[3]),
                        "end": int(columns[4]),
                        "strand": columns[6]
                    }

    
    upstream_seqs = {}
    for gene_id, info in gene_info.items():
        seq_id, start, end, strand = info["seq_id"], info["start"], info["end"], info["strand"]
        if strand == "+":
            upstream_start = max(start - window, 1)
            upstream_end = start - 1
        else:
            upstream_start = end + 1
            upstream_end = min(end + window, len(genome[seq_id]))

        upstream_seq = genome[seq_id].seq[upstream_start-1:upstream_end]
        if strand == "-":
            upstream_seq = upstream_seq.reverse_complement()

        upstream_seqs[gene_id] = (seq_id, upstream_start, upstream_end, upstream_seq.lower())

    
    with open(gene_seq_output, 'w') as fasta_out, open(gene_bed_output, 'w') as bed_out:
        for gene_id, (seq_id, start, end, seq) in upstream_seqs.items():
            fasta_out.write(f">{gene_aliases[gene_id]}_{window}\n{seq}\n")
            bed_out.write(f"{seq_id}\t{start-1}\t{end}\t{gene_aliases[gene_id]}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract upstream sequences of given genes from a genome.")
    parser.add_argument("--genome", required=True, help="Genome sequence file in FASTA format.")
    parser.add_argument("--gene", required=True, help="Gene ID and alias list file, tab-separated.")
    parser.add_argument("--gff", required=True, help="Genome annotation file in GFF3 format.")
    parser.add_argument("--output", default="results", help="Output folder name. Default is 'results'.")
    parser.add_argument("--window", type=int, default=2000, help="Length of the upstream region to extract. Default is 2000bp.")
    parser.add_argument("--feature", default="mRNA", help="Feature type to extract from the GFF file. Default is 'mRNA'.")
    args = parser.parse_args()

    extract_upstream_sequences(args.genome, args.gene, args.gff, args.output, args.window, args.feature)
