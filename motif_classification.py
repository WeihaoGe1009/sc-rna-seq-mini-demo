# an over-simplified 0-shot classifier to practice language model
from transformers import AutoTokenizer, AutoModel
from Bio import SeqIO
import torch
import numpy as np
import umap.umap_ as umap
import matplotlib.pyplot as plt

# Step 0: Load reference set
def load_fasta_with_label(filepath, label):
    return [(str(record.seq), label) for record in SeqIO.parse(filepath, "fasta")]

enhancers = load_fasta_with_label("./data/enh_sequences.fasta", "enhancer")
promoters = load_fasta_with_label("./data/tss_sequences.fasta", "promoter")

reference_set = enhancers + promoters


# Step 1: Load sequences
fasta_file = "./data/upstream_sequences.fasta"
sequences = {record.id: str(record.seq) for record in SeqIO.parse(fasta_file, "fasta")}

# Step 2: Load DNABERT-2 tokenizer/model (6-mer)
model_name = "zhihan1996/DNABERT-2-117M"
tokenizer = AutoTokenizer.from_pretrained(model_name, trust_remote_code=True)
model = AutoModel.from_pretrained(model_name, trust_remote_code=True)

# Helper: convert DNA to 6-mer tokens
def tokenize_dna(seq, k=6):
    return " ".join([seq[i:i+k] for i in range(len(seq)-k+1)])

# Step 3.1: Get embeddings
embedding_dict = {}
model.eval()
with torch.no_grad():
    for gene_id, seq in sequences.items():
        tokenized = tokenizer(tokenize_dna(seq), return_tensors="pt", truncation=True, max_length=512)
        outputs = model(**tokenized)
        # Use [CLS] token embedding (first token)
        #cls_embedding = outputs[0][:, 0, :].squeeze().numpy()
        cls_embedding = outputs[0][:, 0, :].detach().numpy().reshape(1, -1)
        embedding_dict[gene_id] = cls_embedding

# Step 3.2: Embed references 
ref_embeddings = []
ref_labels = []
for seq, label in reference_set:
    tokenized = tokenizer(seq, return_tensors="pt")
    output = model(**tokenized)
    #cls_embedding = output[0][:, 0, :].detach().numpy()
    cls_embedding = outputs[0][:, 0, :].detach().numpy().reshape(-1)
    ref_embeddings.append(cls_embedding)
    ref_labels.append(label)
ref_embeddings = np.vstack(ref_embeddings) 

# Step 4: 0-shot prediction based on cosine_similarity
# need to change the names and types of the variables to make them consistent
from sklearn.metrics.pairwise import cosine_similarity

for kk in embedding_dict.keys():
    similarities = cosine_similarity(embedding_dict[kk].reshape(1, -1), ref_embeddings)
    predicted_labels = [ref_labels[idx] for idx in similarities.argmax(axis=1)]
    print(predicted_labels)

