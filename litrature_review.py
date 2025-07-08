# part1_build_and_save_faiss.py

import pandas as pd
import numpy as np
import faiss
import pickle
from sentence_transformers import SentenceTransformer

# === Step 1: Load Data ===
df = pd.read_csv("your_file.csv")
df = df[["PMID", "Title", "Abstract"]].dropna(subset=["Abstract"]).reset_index(drop=True)

# === Step 2: Load Embedding Model ===
model = SentenceTransformer("all-MiniLM-L6-v2")

# === Step 3: Generate Embeddings ===
abstracts = df["Abstract"].tolist()
embeddings = model.encode(abstracts, show_progress_bar=True)
embedding_matrix = np.array(embeddings).astype("float32")

# === Step 4: Create and Save FAISS Index ===
dimension = embedding_matrix.shape[1]
index = faiss.IndexFlatL2(dimension)
index.add(embedding_matrix)

faiss.write_index(index, "faiss_index.bin")  # Save index to file

# === Step 5: Save Metadata (PMID, Title, Abstract) ===
df.to_pickle("metadata.pkl")  # Save dataframe as pickle




# part2_load_and_query_faiss.py

import numpy as np
import pandas as pd
import faiss
import pickle
from sentence_transformers import SentenceTransformer

# === Step 1: Load Metadata and FAISS Index ===
df = pd.read_pickle("metadata.pkl")
index = faiss.read_index("faiss_index.bin")

# === Step 2: Load Embedding Model ===
model = SentenceTransformer("all-MiniLM-L6-v2")

# === Step 3: Search Function ===
def search_local_faiss(query, top_k=5):
    query_embedding = model.encode([query])[0].astype("float32")
    D, I = index.search(np.array([query_embedding]), top_k)
    return df.iloc[I[0]].reset_index(drop=True)

# === Step 4: Run a Sample Query ===
query = "lung cancer treatment"
results = search_local_faiss(query, top_k=5)

# === Step 5: Print Results ===
print(results[["PMID", "Title", "Abstract"]])
