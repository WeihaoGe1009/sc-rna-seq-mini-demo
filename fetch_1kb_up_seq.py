import requests
from Bio.Seq import Seq

#ensembl_ids = [
#    "ENSG00000105369", "ENSG00000156738", "ENSG00000247982", "ENSG00000007312", "ENSG00000100721",
#    "ENSG00000196735", "ENSG00000179344", "ENSG00000237541", "ENSG00000204257", "ENSG00000242574"
#]
ensembl_ids = []
with open('./data/sig_genes.csv', 'r') as fin:
    for lines in fin:
        ensembl_ids.append(lines.strip())

def get_gene_location(ensembl_id):
    server = "https://rest.ensembl.org"
    ext = f"/lookup/id/{ensembl_id}?expand=1"
    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
    if not r.ok:
        return None
    return r.json()

def get_upstream_sequence(chrom, start, end, strand, length=1000):
    server = "https://rest.ensembl.org"
    if strand == 1:
        region = f"{chrom}:{max(start - length, 1)}..{start - 1}:1"
    else:
        region = f"{chrom}:{end + 1}..{end + length}:1"
    ext = f"/sequence/region/human/{region}"
    r = requests.get(server+ext, headers={ "Content-Type" : "text/plain"})
    if not r.ok:
        return None
    seq = r.text.strip()
    return str(Seq(seq).reverse_complement()) if strand == -1 else seq

# Fetch sequences
upstream_sequences = {}
for gene_id in ensembl_ids:
    info = get_gene_location(gene_id)
    if info:
        seq = get_upstream_sequence(info['seq_region_name'], info['start'], info['end'], info['strand'])
        if seq:
            upstream_sequences[gene_id] = seq

# Print preview
for gene, seq in upstream_sequences.items():
    print(f">{gene}\n{seq[:100]}... (truncated)\n")

with open("./data/upstream_sequences.fasta", "w") as f:
    for gene, seq in upstream_sequences.items():
        f.write(f">{gene}\n{seq}\n")

