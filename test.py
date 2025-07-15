import pandas as pd
from tqdm import tqdm
from langchain_community.document_loaders.csv_loader import CSVLoader
from langchain_community.llms import Ollama
from langchain_core.documents import Document

# === Step 1: Load CSV ===
loader = CSVLoader(
    file_path=r"C:\Users\assis\Downloads\Pubmed_cleaned_csv.csv",
    encoding="utf-8",
    source_column="Abstract",
    metadata_columns=["PMID", "Title"],
    csv_args={
        "delimiter": ",",
        "quotechar": '"',
        "skipinitialspace": True
    }
)
docs = loader.load()

# === Step 2: Initialize Local Ollama LLM ===
llm = Ollama(model="mistral")  # Replace with your local model

# === Step 3: Define Prompt Template ===
def build_classification_prompt(abstract: str) -> str:
    return f"""
You are a biomedical expert.

Please read the following abstract and do the following:
1. Determine if the abstract discusses only the drug "Valsartan" (not sacubitril/valsartan or any other drug).
2. If it mentions only "Valsartan", classify it as INCLUDE.
3. If it mentions any other drug (even along with Valsartan), classify it as EXCLUDE.
4. If it is INCLUDE, also summarize the main finding or focus in 2-3 lines.

Return your answer in this format exactly:
Class: <INCLUDE/EXCLUDE>
Summary: <brief summary or 'NA' if excluded>

Abstract:
\"\"\"
{abstract}
\"\"\"
"""

# === Step 4: Process Abstracts with Ollama ===
results = []
print("🔍 Processing abstracts...")

for doc in tqdm(docs, desc="Evaluating"):
    abstract = doc.page_content
    pmid = doc.metadata.get("PMID", "")
    title = doc.metadata.get("Title", "")
    
    try:
        prompt = build_classification_prompt(abstract)
        response = llm.invoke(prompt)

        # Parse the response
        lines = response.strip().splitlines()
        class_line = next((line for line in lines if line.startswith("Class:")), "Class: EXCLUDE")
        summary_line = next((line for line in lines if line.startswith("Summary:")), "Summary: NA")

        include_status = class_line.replace("Class:", "").strip()
        summary = summary_line.replace("Summary:", "").strip()

        results.append({
            "PMID": pmid,
            "Title": title,
            "Abstract": abstract,
            "Classification": include_status,
            "Summary": summary
        })
    except Exception as e:
        print(f"⚠️ Error processing PMID {pmid}: {e}")
        results.append({
            "PMID": pmid,
            "Title": title,
            "Abstract": abstract,
            "Classification": "ERROR",
            "Summary": str(e)
        })

# === Step 5: Save Results to DataFrame ===
df_result = pd.DataFrame(results)

# Optional: Save to CSV
df_result.to_csv("filtered_valsartan_results.csv", index=False)

print("✅ Processing complete. Output saved to 'filtered_valsartan_results.csv'")
