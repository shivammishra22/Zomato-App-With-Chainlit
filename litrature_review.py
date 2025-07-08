# part1_build_and_save_faiss.py

import pandas as pd
import numpy as np
import faiss
from sentence_transformers import SentenceTransformer

# === Step 1: Load Data ===
df = pd.read_csv("your_file.csv")  # Change to your real path
df = df[["PMID", "Title", "Abstract"]].dropna(subset=["Abstract"]).reset_index(drop=True)

# === Step 2: Load Local Embedding Model ===
model = SentenceTransformer("all-MiniLM-L6-v2")

# === Step 3: Generate Embeddings ===
embeddings = model.encode(df["Abstract"].tolist(), show_progress_bar=True)
embedding_matrix = np.array(embeddings).astype("float32")

# === Step 4: Create and Save FAISS Index ===
index = faiss.IndexFlatL2(embedding_matrix.shape[1])
index.add(embedding_matrix)
faiss.write_index(index, "faiss_index.bin")  # Save index

# === Step 5: Save Metadata (PMID, Title, Abstract) ===
df.to_pickle("metadata.pkl")
#################################################


# part2_load_and_classify_with_ollama.py

import pandas as pd
import numpy as np
import faiss
from sentence_transformers import SentenceTransformer
import subprocess

# === PROMPTS ===
INCLUSION_PROMPT = """
Check if the abstract discusses one or more of the following:
1. Suspected adverse reactions in humans, including those from published abstracts, solicited reports, or manuscripts.
2. Specific situations like pregnancy, paediatrics, elderly, off-label use, overdose, medication error, or misuse.
3. Adverse reactions due to product quality, falsified medicine, or transmission of infection.
4. Lack of therapeutic efficacy.
5. Review of non-company-sponsored clinical trial outcomes.
6. Aggregated adverse reaction data that could become a valid ICSR.
If any of these are present, classify it as INCLUSION.
"""

EXCLUSION_PROMPT = """
Classify the abstract as EXCLUSION if:
1. No adverse event (AE) with company suspect product is discussed.
2. It refers only to animal/preclinical/in-vitro/ex-vivo studies.
3. There’s no or negative causality with company suspect product.
4. Suspect product is from a non-company (different MAH).
5. No identifiable ICSR or medical relevance.
"""

FULL_PROMPT_TEMPLATE = f"""{INCLUSION_PROMPT}\n{EXCLUSION_PROMPT}
Now read the abstract below and respond only with the word INCLUSION or EXCLUSION.\nAbstract:\n{{abstract}}
"""

# === Step 1: Load Metadata & FAISS Index ===
df = pd.read_pickle("metadata.pkl")
index = faiss.read_index("faiss_index.bin")

# === Step 2: Load Embedding Model ===
model = SentenceTransformer("all-MiniLM-L6-v2")

# === Step 3: Ollama Function ===
def query_ollama(text, model_name="llama3"):
    prompt = FULL_PROMPT_TEMPLATE.format(abstract=text)
    result = subprocess.run(
        ["ollama", "run", model_name],
        input=prompt,
        capture_output=True,
        text=True
    )
    return result.stdout.strip().upper()

# === Step 4: (Optional) Filter FAISS Top-K ===
def semantic_filter(query, top_k=500):
    query_vector = model.encode([query])[0].astype("float32")
    _, I = index.search(np.array([query_vector]), top_k)
    return df.iloc[I[0]]

# === Step 5: Run Classification ===
filtered_df = semantic_filter("adverse drug reaction", top_k=500)

inclusion_rows = []
exclusion_rows = []

for _, row in filtered_df.iterrows():
    abstract = row["Abstract"]
    decision = query_ollama(abstract)

    if "INCLUSION" in decision:
        inclusion_rows.append(row)
    elif "EXCLUSION" in decision:
        exclusion_rows.append(row)
    else:
        print("⚠️ Ambiguous result:", decision)

# === Step 6: Save Results ===
inclusion_df = pd.DataFrame(inclusion_rows)
exclusion_df = pd.DataFrame(exclusion_rows)

inclusion_df.to_csv("relevant_df.csv", index=False)
exclusion_df.to_csv("irrelevant_df.csv", index=False)

print("✅ Classification completed and saved to CSV.")


###############################################################




