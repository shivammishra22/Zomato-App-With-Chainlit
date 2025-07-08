import streamlit as st
import pandas as pd
import numpy as np
import faiss
import subprocess
from sentence_transformers import SentenceTransformer

# === Load FAISS and Metadata ===
@st.cache_resource
def load_resources():
    df = pd.read_pickle("metadata.pkl")
    index = faiss.read_index("faiss_index.bin")
    model = SentenceTransformer("all-MiniLM-L6-v2")
    return df, index, model

df, index, model = load_resources()

# === Prompts ===
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

PROMPT_TEMPLATE = f"{INCLUSION_PROMPT}\n{EXCLUSION_PROMPT}\n\nAbstract:\n{{abstract}}\n\nRespond with only INCLUSION or EXCLUSION."

# === Query Ollama LLM Locally ===
def query_ollama(abstract_text, model_name="llama3"):
    prompt = PROMPT_TEMPLATE.format(abstract=abstract_text)
    try:
        result = subprocess.run(
            ["ollama", "run", model_name],
            input=prompt,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="ignore"
        )
        return result.stdout.strip().upper()
    except Exception as e:
        return f"ERROR: {e}"

# === Semantic Filter using FAISS ===
def semantic_filter(query, top_k):
    query_vec = model.encode([query])[0].astype("float32")
    _, I = index.search(np.array([query_vec]), top_k)
    return df.iloc[I[0]]

# === Streamlit App Interface ===
st.title("📑 ICSR Abstract Classifier (Local Ollama + FAISS)")

top_k = st.number_input("How many abstracts (Top K) to classify?", min_value=1, max_value=1000, value=100)
run_button = st.button("🔍 Run Classification")

if run_button:
    st.info("Running semantic search and classification...")
    filtered_df = semantic_filter("adverse drug reaction", top_k=top_k)

    inclusion_results = []

    for _, row in filtered_df.iterrows():
        decision = query_ollama(row["Abstract"])
        if "INCLUSION" in decision:
            inclusion_results.append(f"🔹 **PMID**: {row['PMID']}\n\n{row['Abstract']}\n")
        elif "EXCLUSION" not in decision:
            st.warning(f"⚠️ Ambiguous or error response: {decision[:100]}...")

    if inclusion_results:
        st.subheader("✅ Included Abstracts:")
        for item in inclusion_results:
            st.markdown(item)
    else:
        st.error("❌ No relevant abstracts classified as INCLUSION.")
