import streamlit as st
import pandas as pd
import numpy as np
import faiss
import subprocess
from sentence_transformers import SentenceTransformer

# === Load FAISS index, metadata, and embedding model ===
@st.cache_resource
def load_resources():
    df = pd.read_pickle("metadata.pkl")
    index = faiss.read_index("faiss_index.bin")
    model = SentenceTransformer("all-MiniLM-L6-v2")
    return df, index, model

df, index, model = load_resources()

# === Prompts for LLM ===
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

# === Local Ollama query function ===
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

# === Semantic FAISS filter ===
def semantic_filter(query, top_k):
    query_vec = model.encode([query])[0].astype("float32")
    _, I = index.search(np.array([query_vec]), top_k)
    return df.iloc[I[0]]

# === Streamlit UI ===
st.title("📑 ICSR Abstract Classifier (FAISS + Ollama Local LLM)")

user_query = st.text_input("🔍 Enter your medical search query:", value="adverse reaction")
top_k = st.number_input("📊 How many top similar abstracts to check?", min_value=1, max_value=1000, value=50)
run_button = st.button("🚀 Classify Abstracts")

if run_button and user_query:
    st.info("Running semantic search and classification...")

    filtered_df = semantic_filter(user_query, top_k=top_k)

    inclusions, exclusions, errors = [], [], []

    for _, row in filtered_df.iterrows():
        abstract = row.get("Abstract", "")
        pmid = row.get("PMID", "N/A")

        decision = query_ollama(abstract)
        if "INCLUSION" in decision:
            inclusions.append(f"🔹 **PMID**: {pmid}\n\n{abstract}")
        elif "EXCLUSION" in decision:
            exclusions.append(f"🔸 **PMID**: {pmid}\n\n{abstract}")
        else:
            errors.append(f"⚠️ PMID: {pmid} - Unexpected response: {decision}")

    # === Display Results ===
    st.success("✅ Classification completed.")

    col1, col2 = st.columns(2)

    with col1:
        st.subheader("✅ INCLUSION Abstracts")
        if inclusions:
            for item in inclusions:
                st.markdown(item)
                st.markdown("---")
        else:
            st.info("No relevant abstracts classified as INCLUSION.")

    with col2:
        st.subheader("❌ EXCLUSION Abstracts")
        if exclusions:
            for item in exclusions:
                st.markdown(item)
                st.markdown("---")
        else:
            st.info("No abstracts classified as EXCLUSION.")

    if errors:
        st.warning("⚠️ Some responses were unclear or failed to classify.")
        for err in errors:
            st.text(err)
