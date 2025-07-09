import pandas as pd
import json
import subprocess
import numpy as np
from sentence_transformers import SentenceTransformer, util
import torch

# === Step 1: Load Data ===
csv_path = r"C:\Users\shivam.mishra2\Downloads\literature_data.csv"
df = pd.read_csv(csv_path)

# === Step 2: Print Missing Abstracts ===
missing_abstracts = df[df['Abstract'].isnull() | df['Abstract'].str.strip().eq("")]
print("\n=== Missing Abstracts ===")
print(missing_abstracts[['PMID', 'Title']])

# Filter valid abstracts
df = df.dropna(subset=["Abstract"])
df = df[df["Abstract"].str.strip() != ""]

# === Step 3: Embedding Setup ===
model = SentenceTransformer("all-MiniLM-L6-v2")

# Encode all abstracts
print("🔄 Encoding abstracts...")
abstract_texts = df["Abstract"].tolist()
abstract_embeddings = model.encode(abstract_texts, convert_to_tensor=True)

# === Step 4: Semantic Search Setup ===
inclusive_query = """
Suspected adverse reactions in humans; pregnancy, pediatrics, overdose, misuse,
lack of efficacy, quality defects, falsified medicines, infectious transmission, ICSR, etc.
"""
inclusive_keywords = [
    "adverse reaction", "ICSR", "overdose", "lack of efficacy", "pregnancy", "elderly",
    "misuse", "medication error", "paediatrics", "falsified", "infectious"
]

query_embedding = model.encode(inclusive_query, convert_to_tensor=True)
keyword_embeddings = model.encode(inclusive_keywords, convert_to_tensor=True)

similarity_threshold = 0.5

# === Step 5: Ollama Utilities ===
def classify_with_ollama(abstract):
    prompt = f"""
Classify the abstract as INCLUSION or EXCLUSION based on these criteria:

INCLUSION:
• Suspected adverse reactions in humans.
• Specific situations: pregnancy, overdose, misuse, elderly, pediatric, etc.
• Lack of efficacy, falsified drugs, quality issues, etc.
• Clinical trial results (non-company).
• Aggregated data with potential for ICSR.

EXCLUSION:
• No adverse event with company product.
• Preclinical/animal/in-vitro studies.
• No/negative causality.
• Non-company suspect product.

Abstract:
{abstract}

Answer in one word only: INCLUSION or EXCLUSION
"""
    try:
        result = subprocess.run(
            ["ollama", "run", "llama3", prompt],
            capture_output=True,
            text=True,
            check=True
        )
        return result.stdout.strip().upper()
    except subprocess.CalledProcessError as e:
        print("❌ Classification error:", e.stderr)
        return "UNKNOWN"

def summarize_with_ollama(abstract, abstract_type):
    prompt = f"""
You are a medical reviewer. Provide a concise summary of this {abstract_type} abstract.
Highlight adverse reactions, population groups, study type, and significance.

Abstract:
{abstract}
"""
    try:
        result = subprocess.run(
            ["ollama", "run", "llama3", prompt],
            capture_output=True,
            text=True,
            check=True
        )
        return result.stdout.strip()
    except subprocess.CalledProcessError as e:
        print("❌ Summary error:", e.stderr)
        return ""

# === Step 6: Semantic Search + Ollama Loop ===
inclusive_abstracts = []
exclusive_abstracts = []
inclusive_summaries = []
exclusive_summaries = []

print("🔍 Starting semantic filtering + classification...")
for idx, (pmid, title, abstract) in df[["PMID", "Title", "Abstract"]].iterrows():
    emb = abstract_embeddings[idx]

    sim_query = util.cos_sim(emb, query_embedding).item()
    sim_keywords = util.cos_sim(emb, keyword_embeddings).squeeze()
    max_sim_keyword = torch.max(sim_keywords).item()

    if sim_query >= similarity_threshold or max_sim_keyword >= similarity_threshold:
        decision = classify_with_ollama(abstract)
        record = {"PMID": pmid, "Title": title, "Abstract": abstract}

        if decision == "INCLUSION":
            inclusive_abstracts.append(record)
            summary = summarize_with_ollama(abstract, "INCLUSIVE")
            inclusive_summaries.append((pmid, title, summary))
        elif decision == "EXCLUSION":
            exclusive_abstracts.append(record)
            summary = summarize_with_ollama(abstract, "EXCLUSIVE")
            exclusive_summaries.append((pmid, title, summary))
        else:
            print(f"⚠️ UNKNOWN classification for PMID: {pmid}")
    else:
        print(f"⏭️ Skipped PMID {pmid} (similarity too low: {max(sim_query, max_sim_keyword):.2f})")

# === Step 7: Save Outputs ===
with open("inclusive_abstracts.json", "w", encoding="utf-8") as f:
    json.dump(inclusive_abstracts, f, indent=2, ensure_ascii=False)

with open("exclusive_abstracts.json", "w", encoding="utf-8") as f:
    json.dump(exclusive_abstracts, f, indent=2, ensure_ascii=False)

with open("inclusive_summary_with_pmid.txt", "w", encoding="utf-8") as f:
    for pmid, title, summary in inclusive_summaries:
        f.write(f"PMID: {pmid}\nTitle: {title}\nSummary: {summary}\n\n{'='*100}\n\n")

with open("exclusive_summary_with_pmid.txt", "w", encoding="utf-8") as f:
    for pmid, title, summary in exclusive_summaries:
        f.write(f"PMID: {pmid}\nTitle: {title}\nSummary: {summary}\n\n{'='*100}\n\n")

print("\n✅ DONE: Embedding → Filtering → Classification → Summarization completed.")
