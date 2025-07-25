#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import re
from tqdm.auto import tqdm
from langchain_core.prompts import ChatPromptTemplate
from langchain_community.chat_models import ChatOllama

# --- 1. Load CSV ---
df = pd.read_csv("input.csv", dtype=str)

# --- 2. Check required columns ---
required_cols = {"PMID", "Title", "Abstract", "Authors", "Publication_Date"}
missing_cols = required_cols - set(df.columns)
if missing_cols:
    raise ValueError(f"❌ Missing columns in input file: {', '.join(missing_cols)}")

# --- 3. Get indication from user ---
indication = input("Enter the therapeutic indication to evaluate against the abstracts: ").strip().lower()

# --- 4. Unified Prompt ---
unified_prompt = ChatPromptTemplate.from_messages([
    ("system", f"""You are a pharmacovigilance expert preparing the “Characterisation of Benefits” section of a PSUR.

Your task is to:
1. Evaluate the abstract for benefit characteristics such as:
   - Strength of evidence (comparators, effect size, statistics, methodology)
   - Clinical relevance of the effect size
   - Generalisability
   - Dose–response
   - Duration and consistency
   - Comparative efficacy

2. Determine if it provides evidence for the therapeutic indication: "{indication}"

3. Also, identify if any other therapeutic indications are mentioned in the abstract apart from the primary indication "{indication}". If other indications are found, list them.

Respond **only in the following format**:

Relevance: Relevant or Not Relevant  
Reason: <brief reason for benefit evaluation>  
Indication Match: Yes or No  
Indication Reason: <brief justification regarding indication match or mismatch>  
Other Indications: <list of other indications or 'None'>"""),
    ("user", "{Abstract}")
]).partial(indication=indication)

# --- 5. Initialize Ollama LLM ---
llm = ChatOllama(model="gemma3", temperature=0.1, num_ctx=1000)
chain = unified_prompt | llm

# --- 6. Evaluation Loop ---
rel_benefit_list = []
res_benefit_list = []
ind_match_list = []
ind_reason_list = []
other_indications_list = []

for _, row in tqdm(df.iterrows(), total=len(df), desc="Evaluating abstracts"):
    abstract = row.get("Abstract", "")
    if not isinstance(abstract, str) or not abstract.strip():
        rel_benefit_list.append("Not Available")
        res_benefit_list.append("Abstract is empty or missing")
        ind_match_list.append("Not Available")
        ind_reason_list.append("Abstract is empty or missing")
        other_indications_list.append("Not Available")
        continue

    # --- Unified LLM Call ---
    response = chain.invoke({"Abstract": abstract})
    result_text = response.content if hasattr(response, "content") else str(response)

    # --- Parse LLM Response ---
    match_rel = re.search(r"Relevance:\s*(Relevant|Not Relevant)", result_text, re.IGNORECASE)
    match_reason = re.search(r"Reason:\s*(.+)", result_text, re.IGNORECASE)
    match_ind = re.search(r"Indication Match:\s*(Yes|No)", result_text, re.IGNORECASE)
    match_ind_reason = re.search(r"Indication Reason:\s*(.+)", result_text, re.IGNORECASE)
    match_other_ind = re.search(r"Other Indications:\s*(.+)", result_text, re.IGNORECASE)

    # Store results
    rel_benefit_list.append(match_rel.group(1).strip() if match_rel else "Unclear")
    res_benefit_list.append(match_reason.group(1).strip() if match_reason else "Could not extract reason")
    ind_match_list.append(match_ind.group(1).strip() if match_ind else "Unclear")
    ind_reason_list.append(match_ind_reason.group(1).strip() if match_ind_reason else "Could not extract indication reason")
    other_indications_list.append(match_other_ind.group(1).strip() if match_other_ind else "None")

# --- 7. Assign Output Columns ---
df["Relavance_benifit"] = rel_benefit_list
df["Reson_Benifit"] = res_benefit_list
df["Indication Match"] = ind_match_list
df["Indication Reason"] = ind_reason_list
df["Other Indications"] = other_indications_list

# --- 8. Select Final Columns ---
final_columns = [
    "PMID",
    "Title",
    "Abstract",
    "Authors",
    "Publication_Date",
    "Relavance_benifit",
    "Reson_Benifit",
    "Indication Match",
    "Indication Reason",
    "Other Indications"
]
df_final = df[final_columns]

# --- 9. Export to CSV ---
output_path = "psur_evaluation_output.csv"
df_final.to_csv(output_path, index=False)
print(f"✅ Evaluation complete. Output saved to '{output_path}'")
