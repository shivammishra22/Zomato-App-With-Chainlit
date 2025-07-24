#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import re
from tqdm.auto import tqdm
from langchain_core.prompts import ChatPromptTemplate
from langchain_community.chat_models import ChatOllama

# --- 1. Load your data ---
# Replace "input.csv" with your actual file
df = pd.read_csv("input.csv", dtype=str)

# --- 2. Get the indication from the user ---
indication = input("Enter the therapeutic indication to evaluate against the abstracts: ").strip()

# --- 3. Build the prompt template with reason request ---
template = ChatPromptTemplate.from_messages([
    ("system", """You are a pharmacovigilance expert generating the PSUR sub‑section titled “Characterisation of Benefits”.

For each abstract, evaluate:
1) Strength of evidence (comparators, effect size, stats, consistency, methodology)  
2) Clinical relevance of effect size  
3) Generalisability  
4) Dose‑response  
5) Duration of effect  
6) Comparative efficacy

Also:
- Does the abstract support the indication: "{Indication}"?  Provide “Yes” or “No” **and** a brief reason for your choice.
- List any *other* indications mentioned.

Respond exactly in this format (one value per line):
Relevance: Relevant or Not Relevant  
Indication Match: Yes or No  
Reason: [brief explanation]  
Other Indications: [comma‑separated list or 'None']"""),
    ("user", "{Abstract}")
])

# --- 4. Initialize the model & chain ---
llm = ChatOllama(model="gemma3:4b", temperature=0.1, num_ctx=1000)
chain = template | llm

# --- 5. Prepare regexes for parsing ---
pat_rel    = re.compile(r"^Relevance:\s*(Relevant|Not\s+Relevant)", re.IGNORECASE)
pat_match  = re.compile(r"^Indication\s+Match:\s*(Yes|No)", re.IGNORECASE)
pat_reason = re.compile(r"^Reason:\s*(.+)", re.IGNORECASE)
pat_other  = re.compile(r"^Other\s+Indications:\s*(.*)", re.IGNORECASE)

# --- 6. Run through abstracts with tqdm ---
results = []
for abstract in tqdm(df["Abstract"], desc="Evaluating abstracts"):
    resp = chain.invoke({"Abstract": abstract, "Indication": indication}).content.strip()

    # default fallbacks
    rel    = "Not Available"
    match  = "Not Available"
    reason = ""
    other  = ""

    for line in resp.splitlines():
        if m := pat_rel.match(line):
            rel = m.group(1)
        elif m := pat_match.match(line):
            match = m.group(1)
        elif m := pat_reason.match(line):
            reason = m.group(1).strip()
        elif m := pat_other.match(line):
            other = m.group(1).strip()
            if other.lower() in ("none", "n/a", "not applicable", ""):
                other = ""

    results.append((rel, match, reason, other))

# --- 7. Assemble the output DataFrame ---
out_df = df[["PMID", "Title", "Abstract", "Publication_Year"]].copy()
out_df[["Relevance", "Indication Match", "Reason", "Other Indications"]] = pd.DataFrame(
    results, index=out_df.index
)

# --- 8. Save a clean CSV ---
out_file = "results_with_reasons.csv"
out_df.to_csv(out_file, index=False)

print(f"\n✅ Done! Results saved to '{out_file}'\n")
print(out_df.head(10).to_markdown(index=False))
