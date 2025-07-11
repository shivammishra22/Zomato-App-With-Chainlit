import pandas as pd
import numpy as np
from langchain_community.document_loaders.csv_loader import CSVLoader
from langchain_community.vectorstores import FAISS
from langchain_community.embeddings import OllamaEmbeddings
from langchain_community.llms import Ollama
from langchain.chains import RetrievalQA
from langchain.prompts import PromptTemplate

# === Step 1: Load CSV Data ===
loader = CSVLoader(
    file_path=r"C:\Users\assis\Downloads\Pubmed_cleaned_csv.csv",
    encoding="utf-8",
    source_column="Abstract",
    metadata_columns=["PMID"],
    csv_args={
        "delimiter": ",",
        "quotechar": '"',
        "skipinitialspace": True
    }
)

documents = loader.load()

# === Step 2: Convert to Embeddings ===
embedding_model = OllamaEmbeddings(model="mxbai-embed-large")  # or try "nomic-embed-text"

# === Step 3: Create FAISS Vector Store ===
vectorstore = FAISS.from_documents(documents, embedding_model)

# === Step 4: Define Prompt ===
template = """
You are the pharmacovigilance expert generating the PSUR Sub-Section titled "Characterisation of Benefit Data".
Integrate the baseline benefit information with any new benefit data provided in the input.

Find the output from the given abstract based on the following instructions:
1. Dose-response characterisation.
2. Duration of effect.
3. Comparative efficacy.
4. Provide a concise but critical evaluation of the strengths and limitations of the evidence.
5. Adequacy of dose-response characterisation.
6. Strength of evidence of benefit, including comparator(s), effect size, statistical rigor, methodological strengths and deficiencies, and consistency across studies.
7. Clinical relevance of the effect size.
8. Generalisability of treatment response across the indicated patient population, including sub-populations.

Be concise but include key findings. Use bullet points where clarity is needed. Use professional pharmacovigilance language suitable for regulatory reporting.

Abstract:
{context}
"""

prompt = PromptTemplate.from_template(template)

# === Step 5: Setup LLM Locally via Ollama ===
llm = Ollama(model="llama3")  # you can use any supported local model

# === Step 6: Create Retrieval QA Chain ===
retriever = vectorstore.as_retriever(search_kwargs={"k": 5})

qa_chain = RetrievalQA.from_chain_type(
    llm=llm,
    retriever=retriever,
    chain_type="stuff",
    chain_type_kwargs={"prompt": prompt}
)

# === Step 7: Ask a Question ===
query = "Summarize the benefit data for regulatory submission using the given criteria."
result = qa_chain.run(query)

# === Step 8: Print Final Output ===
print("========== Final Summary ==========")
print(result)
