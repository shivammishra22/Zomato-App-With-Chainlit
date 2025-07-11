import pandas as pd
import numpy as np
from langchain_community.document_loaders.csv_loader import CSVLoader
from langchain_community.embeddings import OllamaEmbeddings
from langchain_community.llms import Ollama
from langchain_community.vectorstores import FAISS
from langchain.prompts import PromptTemplate
from langchain.chains import LLMChain
from langchain_core.documents import Document
from tqdm import tqdm

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

# === Step 2: Initialize Embedding Model ===
embedding_model = OllamaEmbeddings(model="mxbai-embed-large")

# === Step 3: Build FAISS Vector Store ===
db = FAISS.from_documents(docs, embedding_model)

# === Step 4: Define Evaluation Prompt ===
evaluation_prompt = PromptTemplate.from_template("""
You are a pharmacovigilance expert. Evaluate the abstract below and answer:
"Does this abstract provide information relevant to the 'Characterisation of Benefit Data' section of a PSUR based on the criteria listed?"

Criteria:
1. Dose-response characterisation
2. Duration of effect
3. Comparative efficacy
4. Strengths and limitations of evidence
5. Adequacy of dose-response characterisation
6. Effect size, comparators, statistical rigor, consistency across studies
7. Clinical relevance
8. Generalisability across populations

Return only "Yes" or "No".

Abstract:
{abstract}
""")

llm = Ollama(model="llama3")  # Local LLM

# === Step 5: Create LLM Evaluation Chain ===
llm_chain = LLMChain(llm=llm, prompt=evaluation_prompt)

# === Step 6: Evaluate Each Abstract ===
results = []
for doc in tqdm(docs, desc="Evaluating abstracts with LLM..."):
    abstract = doc.page_content.strip()
    metadata = doc.metadata
    if not abstract:
        continue

    try:
        llm_response = llm_chain.run({"abstract": abstract}).strip()
        if "yes" in llm_response.lower():
            decision = "Yes"
        elif "no" in llm_response.lower():
            decision = "No"
        else:
            decision = "Unclear"
    except Exception as e:
        decision = f"Error: {str(e)}"

    results.append({
        "PMID": metadata.get("PMID", ""),
        "Title": metadata.get("Title", ""),
        "Abstract": abstract,
        "LLM_Evaluation": decision
    })

# === Step 7: Convert to DataFrame ===
df_result = pd.DataFrame(results)

# === Step 8: Display Table ===
from ace_tools import display_dataframe_to_user
display_dataframe_to_user("PSUR Abstract Evaluation", df_result)
