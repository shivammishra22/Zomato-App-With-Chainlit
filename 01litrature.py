import pandas as pd
import numpy as np
import faiss
import re
from concurrent.futures import ThreadPoolExecutor
from langchain_ollama import OllamaEmbeddings
from langchain.schema import Document

# === Load Data ===
df = pd.read_csv("your_file.csv")  # Replace with actual file path
df = df[["PMID", "Title", "Abstract"]].dropna(subset=["Abstract"])

# === Convert Abstracts to Embeddings ===
embedding_model = OllamaEmbeddings(model="mxbai-embed-large:latest")
abstracts = df["Abstract"].tolist()
abstract_embeddings = embedding_model.embed_documents(abstracts)
abstract_embeddings = np.array(abstract_embeddings)

# === Build FAISS Index ===
index = faiss.IndexFlatL2(abstract_embeddings.shape[1])
index.add(abstract_embeddings)

# === Semantic Search Query ===
query = "case of adverse reaction in paediatric or elderly population"
keywords = ["adverse reaction", "overdose", "medication error", "pregnancy", "efficacy", "clinical trial"]
threshold = 0.6  # Lower means more loose, higher means more strict

query_vector = embedding_model.embed_query(query)
query_vector = np.array([query_vector])
distances, indices = index.search(query_vector, k=20)  # Fetch top 20 similar

# === Filter Results by Distance Threshold and Keyword Match ===
similar_rows = []
for dist, idx in zip(distances[0], indices[0]):
    if idx >= len(df):
        continue
    abstract = df.iloc[idx]["Abstract"]
    if dist <= threshold and any(re.search(kw, abstract, re.IGNORECASE) for kw in keywords):
        similar_rows.append((idx, abstract))

# === Prompt-based Classification and Summary ===
INCLUSION_PROMPT = """
Check if the abstract discusses one or more of the following:
- Suspected adverse reactions in humans, including those from published abstracts, solicited reports, or manuscripts.
- Specific situations like pregnancy, paediatrics, elderly, off-label use, overdose, medication error, or misuse.
- Adverse reactions due to product quality, falsified medicine, or transmission of infection.
- Lack of therapeutic efficacy.
- Review of non-company-sponsored clinical trial outcomes.
- Aggregated adverse reaction data that could become a valid ICSR.
If any of these are present, classify it as INCLUSION.
"""

EXCLUSION_PROMPT = """
Classify the abstract as EXCLUSION if:
- No adverse event (AE) with company suspect product is discussed.
- It refers only to animal/preclinical/in-vitro/ex-vivo studies.
- There’s no or negative causality with company suspect product.
- Suspect product is from a non-company (different MAH).
 No identifiable ICSR or medical relevance. 
"""

import ollama
def classify_and_summarize(row):
    pmid, title, abstract = row["PMID"], row["Title"], row["Abstract"]

    # Classification prompt
    response_inc = ollama.generate(model="llama3", prompt=INCLUSION_PROMPT + f"\n\nAbstract:\n{abstract}")
    inclusion_result = response_inc["response"].strip()

    if "INCLUSION" in inclusion_result.upper():
        classification = "INCLUSION"
    else:
        response_exc = ollama.generate(model="llama3", prompt=EXCLUSION_PROMPT + f"\n\nAbstract:\n{abstract}")
        classification = "EXCLUSION" if "EXCLUSION" in response_exc["response"].upper() else "UNCERTAIN"

    # Summary prompt
    prompt_summary = f'''Generate a short and precise summary of this abstract: "{abstract}". Start with "The abstract reports..."'''
    response_sum = ollama.generate(model="llama3", prompt=prompt_summary)
    summary = response_sum["response"].strip()

    return {
        "PMID": pmid,
        "Title": title,
        "Classification": classification,
        "Summary": summary
    }

# === Run Classification and Summary in Parallel ===
from tqdm import tqdm
tqdm.pandas()

filtered_df = df.loc[[idx for idx, _ in similar_rows]]
results = []
with ThreadPoolExecutor(max_workers=4) as executor:
    results = list(executor.map(classify_and_summarize, [row for _, row in filtered_df.iterrows()]))

# === Print Results ===
for res in results:
    print(f"\nPMID: {res['PMID']}")
    print(f"Title: {res['Title']}")
    print(f"Classification: {res['Classification']}")
    print(f"Summary: {res['Summary']}")
    print("-" * 100)
