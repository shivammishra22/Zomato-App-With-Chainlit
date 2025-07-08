 #Step 2: Load the Data

import pandas as pd

# Load your dataset
df = pd.read_csv("your_file.csv")  # Or Excel or other formats
df = df[["PMID", "Title", "Abstract"]]  # Keep only necessary columns
df.dropna(subset=["Abstract"], inplace=True)  # Remove rows with empty abstracts
df.reset_index(drop=True, inplace=True)

# Step 3: Generate Embeddings Locally using Sentence Transformers

from sentence_transformers import SentenceTransformer

# Load local embedding model (you can replace with a better one like 'all-MiniLM-L6-v2' if offline)
model = SentenceTransformer('all-MiniLM-L6-v2')  # Downloads model first time only

# Create embeddings for all abstracts
abstracts = df["Abstract"].tolist()
embeddings = model.encode(abstracts, show_progress_bar=True)

# Step 4: Store Embeddings in FAISS

import faiss
import numpy as np

# Convert embeddings to float32 (FAISS requirement)
embedding_dim = embeddings[0].shape[0]
embedding_matrix = np.array(embeddings).astype("float32")

# Build FAISS index
index = faiss.IndexFlatL2(embedding_dim)
index.add(embedding_matrix)

# Step 5: Define Search Function

def search_abstracts(query, top_k=5):
    query_embedding = model.encode([query])[0].astype("float32")
    D, I = index.search(np.array([query_embedding]), top_k)
    
    # Return top_k results
    results = df.iloc[I[0]]
    return results

#  Step 6: Run a Query and See Top Matches

query = "pancreatic cancer treatment"
top_results = search_abstracts(query, top_k=5)
print(top_results[["PMID", "Title", "Abstract"]])


# Optional) Step 7: Use Ollama Locally to Summarize or Interact
# If you have Ollama running a model like llama3:

import subprocess

def ask_ollama(prompt, model="llama3"):
    result = subprocess.run(
        ["ollama", "run", model],
        input=prompt,
        capture_output=True,
        text=True
    )
    return result.stdout

# Example: Summarize the abstract
abstract_text = top_results.iloc[0]["Abstract"]
summary = ask_ollama(f"Summarize the following abstract:\n\n{abstract_text}")
print(summary)
