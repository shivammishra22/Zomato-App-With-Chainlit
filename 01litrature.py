import pandas as pd
import numpy as np
import faiss
import re
from sentence_transformers import SentenceTransformer

# === Step 1: Load CSV ===
df = pd.read_csv("your_file.csv")  # replace with your CSV file path

# Drop missing abstracts
df = df.dropna(subset=["Abstract"]).reset_index(drop=True)

# === Step 2: Text Preprocessing to Reduce Size ===
def preprocess(text):
    text = re.sub(r"\s+", " ", text)  # Remove extra whitespaces/newlines
    text = re.sub(r"[^\w\s.,;:?!-]", "", text)  # Remove strange characters
    return text.strip()

df["Processed_Abstract"] = df["Abstract"].apply(preprocess)

# === Step 3: Generate Embeddings Locally ===
model = SentenceTransformer("all-MiniLM-L6-v2")  # Fast and small

embeddings = model.encode(df["Processed_Abstract"].tolist(), show_progress_bar=True)
embeddings = np.array(embeddings).astype("float32")

# === Step 4: Create FAISS Index ===
dimension = embeddings.shape[1]
index = faiss.IndexFlatL2(dimension)
index.add(embeddings)

# === Step 5: Combined Semantic + Keyword Search Function ===
def search(query, keyword=None, top_k=5):
    query_vector = model.encode([query]).astype("float32")
    distances, indices = index.search(query_vector, k=top_k)

    match_df = df.iloc[indices[0]].copy()
    match_df["Similarity"] = distances[0]

    if keyword:
        match_df = match_df[match_df["Processed_Abstract"].str.contains(keyword, case=False)]

    return match_df[["PMID", "Title", "Abstract", "Similarity"]]

# === Example Usage ===
result = search(query="adverse reaction in pregnancy", keyword="pregnancy", top_k=10)
print(result)
