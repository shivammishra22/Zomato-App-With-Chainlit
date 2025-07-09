import pandas as pd
import json
import subprocess

# === Load Data ===
csv_path = r"C:\Users\shivam.mishra2\Downloads\literature_data.csv"
df = pd.read_csv(csv_path)

# === Step 1: Print Missing Abstracts Only ===
missing_abstracts = df[df['Abstract'].isnull() | df['Abstract'].str.strip().eq("")]
print("\n=== Missing Abstracts (PMID, Title) ===")
print(missing_abstracts[['PMID', 'Title']])

# Remove rows with missing Abstracts
df = df.dropna(subset=["Abstract"])
df = df[df["Abstract"].str.strip() != ""]

# === Step 2: Prompts for Classification ===
INCLUSION_PROMPT = """
Inclusion Criteria for Literature article:
• Suspected adverse reactions in humans (abstracts, manuscripts, solicited reports).
• Product brand not specified assumed to be Jubilant.
• Specific Situations: pregnancy, paediatrics, elderly, overdose, misuse, off-label use.
• Lack of efficacy, quality defects, falsified products, transmission of infection.
• Clinical trial reviews (non-company sponsored).
• Aggregate ADR data with potential to become valid ICSR.
"""

EXCLUSION_PROMPT = """
Exclusion Criteria for Literature article:
• No adverse event with company suspect product.
• Preclinical/Animal/In-vitro/In-vivo/Ex-vivo study.
• No or negative causality with company suspect product.
• Suspect product not from company.
• Other: no ICSR(s) (explain briefly).
"""

# === Step 3: Classification Function Using Ollama ===
def classify_with_ollama(abstract):
    prompt = f"""
Classify the abstract as INCLUSION or EXCLUSION based on these criteria:

INCLUSION:
{INCLUSION_PROMPT}

EXCLUSION:
{EXCLUSION_PROMPT}

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

# === Step 4: Summary Function Using Ollama ===
def summarize_with_ollama(abstract, abstract_type):
    prompt = f"""
You are a medical reviewer. Provide a concise summary of this {abstract_type} abstract.
Highlight relevant findings: adverse reactions, population groups, study type, and significance.

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
        print("❌ Summarization error:", e.stderr)
        return ""

# === Step 5: Classify and Summarize ===
inclusive_abstracts = []
exclusive_abstracts = []
inclusive_summaries = []
exclusive_summaries = []

print("\n=== Processing Abstracts ===")
for _, row in df.iterrows():
    pmid = row["PMID"]
    title = row["Title"]
    abstract = row["Abstract"]

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

# === Step 6: Save JSON Files ===
with open("inclusive_abstracts.json", "w", encoding="utf-8") as f:
    json.dump(inclusive_abstracts, f, indent=2, ensure_ascii=False)

with open("exclusive_abstracts.json", "w", encoding="utf-8") as f:
    json.dump(exclusive_abstracts, f, indent=2, ensure_ascii=False)

# === Step 7: Save Summary Text Files ===
with open("inclusive_summary_with_pmid.txt", "w", encoding="utf-8") as f:
    for pmid, title, summary in inclusive_summaries:
        f.write(f"PMID: {pmid}\nTitle: {title}\nSummary: {summary}\n\n{'='*100}\n\n")

with open("exclusive_summary_with_pmid.txt", "w", encoding="utf-8") as f:
    for pmid, title, summary in exclusive_summaries:
        f.write(f"PMID: {pmid}\nTitle: {title}\nSummary: {summary}\n\n{'='*100}\n\n")

print("\n✅ All processing complete. JSON and summary TXT files generated.")
