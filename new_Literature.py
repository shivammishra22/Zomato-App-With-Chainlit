#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import re
from tqdm.auto import tqdm
from langchain_core.prompts import ChatPromptTemplate
from langchain_community.chat_models import ChatOllama

# --- 1. Load input CSV ---
df = pd.read_csv("input.csv", dtype=str)  # Ensure this file exists with 'PMID', 'Title', 'Abstract', 'Publication_Year'

# --- 2. Get indication from the user ---
indication = input("Enter the therapeutic indication to evaluate against the abstracts: ").strip().lower()

# --- 3. Build structured evaluation prompt ---
template = ChatPromptTemplate.from_messages([
    ("system", f"""You are a pharmacovigilance expert preparing the “Characterisation of Benefits” section for a PSUR.

For each abstract provided, perform a structured evaluation based on the following criteria:

1. Strength of evidence (comparators, effect size, statistics, methodology)
2. Clinical relevance of the effect size
3. Generalisability of the results to wider populations
4. Presence of a dose–response relationship
5. Duration and consistency of observed effects
6. Comparative efficacy with existing therapies

Then, answer the following based on the content:

- Relevance: Is the abstract relevant to the benefit assessment? (Respond with "Relevant" or "Not Relevant")
- Indication Match: Does the study support the therapeutic indication: **"{indication}"**? (Respond with "Yes" or "No")
- Reason: Provide a short explanation for your answer to Indication Match.
- Other Indications: Mention any additional therapeutic uses or benefits **different from "{indication}"** that are supported in the abstract. If none, respond with "None".

⚠️ Return your output in exactly the following format:
Relevance: Relevant or Not Relevant  
Indication Match: Yes or No  
Reason: [One-sentence justification]  
Other Indications: [Comma-separated list or 'None']"""),
    ("user", "{Abstract}")
])

# --- 4. Initialize local LLM ---
llm = ChatOllama(model="gemma3", temperature=0.1, num_ctx=1000)
chain = template | llm

# --- 5. Regex for response extraction ---
pat_rel    = re.compile(r"^Relevance:\s*(Relevant|Not\s+Relevant)", re.IGNORECASE)
pat_match  = re.compile(r"^Indication\s+Match:\s*(Yes|No)", re.IGNORECASE)
pat_reason = re.compile(r"^Reason:\s*(.+)", re.IGNORECASE)
pat_other  = re.compile(r"^Other\s+Indications:\s*(.*)", re.IGNORECASE)

# --- 6. Process each abstract ---
results = []
for abstract in tqdm(df["Abstract"], desc="Evaluating abstracts"):
    try:
        resp = chain.invoke({"Abstract": abstract, "Indication": indication}).content.strip()
    except Exception as e:
        print(f"⚠️ LLM error: {e}")
        resp = ""

    # Fallback values
    rel = "Not Available"
    match = "Not Available"
    reason = "No reason extracted"
    other = "None"

    for line in resp.splitlines():
        if m := pat_rel.match(line):
            rel = m.group(1).strip()
        elif m := pat_match.match(line):
            match = m.group(1).strip()
        elif m := pat_reason.match(line):
            reason = m.group(1).strip()
        elif m := pat_other.match(line):
            other_raw = m.group(1).strip()
            if other_raw.lower() in ("none", "n/a", "not applicable", ""):
                other = "None"
            else:
                other_list = [ind.strip() for ind in other_raw.split(",")]
                filtered = [ind for ind in other_list if ind.lower() != indication]
                other = ", ".join(filtered) if filtered else "None"

    results.append((rel, match, reason, other))

# --- 7. Create output DataFrame ---
out_df = df[["PMID", "Title", "Abstract", "Publication_Year"]].copy()
out_df[["Relevance", "Indication Match", "Reason", "Other Indications"]] = pd.DataFrame(results, index=out_df.index)

# --- 8. Save to CSV ---
out_file = "hello.csv"
out_df.to_csv(out_file, index=False)

print(f"\n✅ Done! Results saved to '{out_file}'\n")
print(out_df.head(10).to_markdown(index=False))
