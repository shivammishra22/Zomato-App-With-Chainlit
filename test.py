import ssl
ssl._create_default_https_context = ssl._create_unverified_context

import streamlit as st
import pandas as pd
import io
import time
import re
from datetime import datetime
from typing import List, Dict, Optional, Tuple
from tqdm import tqdm

from langchain_community.utilities import SerpAPIWrapper
from langchain_community.document_loaders import DataFrameLoader
from langchain_ollama import ChatOllama
from langchain_huggingface import HuggingFaceEmbeddings
from langchain_core.prompts import ChatPromptTemplate, HumanMessagePromptTemplate
from langchain.chains import LLMChain

try:
    from Bio import Entrez
except ImportError:
    st.error("Biopython not installed. Please run `pip install biopython`")
    st.stop()

# -- PubMed Extractor Class --
class PubMedExtractor:
    def __init__(self, email: str, api_key: Optional[str] = None):
        Entrez.email = email
        if api_key:
            Entrez.api_key = api_key

    def search_pubmed(self, search_term: str, retmax: int = 20) -> List[str]:
        handle = Entrez.esearch(db="pubmed", term=search_term, retmax=retmax)
        search_results = Entrez.read(handle)
        handle.close()
        return search_results.get("IdList", [])

    def fetch_abstracts(self, pmids: List[str]) -> List[Dict]:
        if not pmids:
            return []
        handle = Entrez.efetch(db="pubmed", id=pmids, rettype="medline", retmode="xml")
        records = Entrez.read(handle)
        handle.close()
        articles = []
        for record in records.get("PubmedArticle", []):
            medline = record.get("MedlineCitation", {})
            article = medline.get("Article", {})
            pmid = medline.get("PMID", "")
            title = article.get("ArticleTitle", "N/A")
            abstract = article.get("Abstract", {}).get("AbstractText", ["N/A"])
            if isinstance(abstract, list):
                abstract = " ".join(abstract)
            authors = article.get("AuthorList", [])
            author_names = "; ".join([a.get("LastName", "") + ", " + a.get("ForeName", "") for a in authors if "LastName" in a and "ForeName" in a])
            journal = article.get("Journal", {}).get("Title", "N/A")
            year = article.get("Journal", {}).get("JournalIssue", {}).get("PubDate", {}).get("Year", "N/A")
            articles.append({
                "PMID": pmid,
                "Title": title,
                "Abstract": abstract,
                "Authors": author_names,
                "Journal": journal,
                "Publication_Year": year
            })
        return articles

# -- Streamlit App --
def main():
    st.set_page_config(page_title="Pharmacovigilance Search Tool", layout="wide")
    st.title("🔎 Pharmacovigilance Abstract Analyzer")

    with st.sidebar:
        st.header("🔧 Configuration")
        email = "your_email@example.com"  # Default email to avoid manual entry
        api_key = "6f5df4899c545b65d2b584c22e70ec181608"  # Automatically used
        serpapi_key = "adb5d6da4a13ced8ad8f6f0d7b41804ae6df887f43d142ecfedaaa3c223eeebe"  # Automatically used
        model_name = st.selectbox("LLM Model", ["gemma3:4b", "qwen3:4b"], index=0)

    # SerpAPI Search
    st.subheader("🧠 Step 1: Drug Discovery using SerpAPI")
    base_query = st.text_input("Query", value="what are the other phosphodiesterase inhibitor medicine similar to Tadalafil ?")
    if st.button("🔍 Find Related Drugs"):
        serpapi = SerpAPIWrapper(serpapi_api_key=serpapi_key)
        serp_results = serpapi.run(base_query)
        llm = ChatOllama(model=model_name, temperature=0.1, num_ctx=2000)
        llm_prompt = f"""
        From the following search result text, extract and list all phosphodiesterase inhibitor drugs similar to Tadalafil. Just give names:

        {serp_results}
        """
        llm_response = llm([HumanMessagePromptTemplate.from_template(llm_prompt).format()])
        extracted_drugs = llm_response.content
        st.success("✅ Extracted Drugs:")
        st.code(extracted_drugs)
        st.session_state["drug_names"] = extracted_drugs

    # PubMed Search
    st.subheader("📚 Step 2: PubMed Abstract Extraction")
    drug_input = st.text_area("Enter Drug Names to Search", value=st.session_state.get("drug_names", "Tadalafil"))
    if st.button("🔬 Fetch PubMed Abstracts"):
        search_term = f"({drug_input}) AND (efficacy OR effectiveness OR treatment OR outcome) AND humans"
        extractor = PubMedExtractor(email, api_key)
        pmids = extractor.search_pubmed(search_term)
        st.info(f"Found {len(pmids)} articles, fetching abstracts...")
        articles = extractor.fetch_abstracts(pmids)
        df = pd.DataFrame(articles)
        st.session_state["abstract_df"] = df
        st.dataframe(df)

    # Summarization
    if "abstract_df" in st.session_state:
        st.subheader("🧾 Step 3: Summarization of Abstracts")
        df = st.session_state["abstract_df"]
        summarizer = ChatOllama(model=model_name, temperature=0.1, num_ctx=2000)
        prompt = ChatPromptTemplate.from_messages([
            ("system", "You are a pharmacovigilance expert. Summarize the following medical abstract clearly."),
            ("human", "{Abstract}")
        ])
        chain = LLMChain(llm=summarizer, prompt=prompt)
        summaries = []
        for abstract in tqdm(df['Abstract'], desc="Summarizing"):
            result = chain.invoke({"Abstract": abstract})
            summaries.append(result['text'])
        df['Summary'] = summaries
        st.dataframe(df[['PMID', 'Title', 'Summary']])

    # Step 4: Select Relevant Summaries
    st.subheader("🎯 Step 4: Select Relevant Summaries")
    if "abstract_df" in st.session_state:
        df = st.session_state["abstract_df"]

    # Combine summaries with PMIDs
    combined_summaries = "\n\n".join(
        f"PMID: {pmid}\nAbstract: {summary}" for pmid, summary in zip(df['PMID'], df['Summary'])
    )

    # Define the LLM and prompt
    filter_llm = ChatOllama(model="qwen3:4b", temperature=0.1, num_ctx=15000)
    human_prompt = HumanMessagePromptTemplate.from_template(
        "Select the top 3 abstracts from the following list that best match the indication: Treatment of schizophrenia or psychosis.\n\n{summaries}"
    )
    chat_prompt = ChatPromptTemplate.from_messages([human_prompt])
    select_chain = LLMChain(llm=filter_llm, prompt=chat_prompt)

    # Run the chain
    response = select_chain.invoke({"summaries": combined_summaries})

    # Extract PMIDs using regex
    selected_pmids = re.findall(r'PMID:\s*(\d+)', response['text'])

    # Filter the DataFrame
    filtered_df = df[df['PMID'].astype(str).isin(selected_pmids)]

    # Display results
    st.success("Top 3 Selected Abstracts:")
    st.dataframe(filtered_df)

    # Download option
    csv_buffer = io.StringIO()
    filtered_df.to_csv(csv_buffer, index=False)
    st.download_button("📥 Download Selected", csv_buffer.getvalue(), file_name="top3_selected.csv", mime="text/csv")


if __name__ == "__main__":
    main()


  C:\Users\shivam.mishra2\Downloads\Literature\Section18\final.py:166 in <module>

    163 
    164 
    165 if __name__ == "__main__":
  ❱ 166 │   main()
    167 

  C:\Users\shivam.mishra2\Downloads\Literature\Section18\final.py:135 in main

    132 │   
    133 │   # Combine summaries with PMIDs
    134 │   combined_summaries = "\n\n".join(
  ❱ 135 │   │   f"PMID: {pmid}\nAbstract: {summary}" for pmid, summary in zip(df['PMID
    136 │   )
    137 │   
    138 │   # Define the LLM and prompt
────────────────────────────────────────────────────────────────────────────────────────
UnboundLocalError: cannot access local variable 'df' where it is not associated with a
value
C:\Users\shivam.mishra2\Downloads\Literature\Section18\final.py:92: LangChainDeprecationWarning: The method `BaseChatModel.__call__` was deprecated in langchain-core 0.1.7 and will be removed in 1.0. Use :meth:`~invoke` instead.
  llm_response = llm([HumanMessagePromptTemplate.from_template(llm_prompt).format()])
────────────────────────── Traceback (most recent call last) ───────────────────────────
  C:\Users\shivam.mishra2\Downloads\Literature\myenv\Lib\site-packages\streamlit\runti
  me\scriptrunner\exec_code.py:128 in exec_func_with_error_handling

  C:\Users\shivam.mishra2\Downloads\Literature\myenv\Lib\site-packages\streamlit\runti
  me\scriptrunner\script_runner.py:669 in code_to_exec

  C:\Users\shivam.mishra2\Downloads\Literature\Section18\final.py:166 in <module>

    163 
    164 
    165 if __name__ == "__main__":
  ❱ 166 │   main()
    167 

  C:\Users\shivam.mishra2\Downloads\Literature\Section18\final.py:135 in main

    132 │   
    133 │   # Combine summaries with PMIDs
    134 │   combined_summaries = "\n\n".join(
  ❱ 135 │   │   f"PMID: {pmid}\nAbstract: {summary}" for pmid, summary in zip(df['PMID
    136 │   )
    137 │   
    138 │   # Define the LLM and prompt
────────────────────────────────────────────────────────────────────────────────────────
────────────────────────────────────────────────────────────────────────────────────────
UnboundLocalError: cannot access local variable 'df' where it is not associated with a
value
────────────────────────── Traceback (most recent call last) ───────────────────────────
  C:\Users\shivam.mishra2\Downloads\Literature\myenv\Lib\site-packages\streamlit\runti  
  me\scriptrunner\exec_code.py:128 in exec_func_with_error_handling

  C:\Users\shivam.mishra2\Downloads\Literature\myenv\Lib\site-packages\streamlit\runti  
  me\scriptrunner\script_runner.py:669 in code_to_exec

  C:\Users\shivam.mishra2\Downloads\Literature\Section18\final.py:166 in <module>       

    163 
    164 
    165 if __name__ == "__main__":
  ❱ 166 │   main()
    167 

  C:\Users\shivam.mishra2\Downloads\Literature\Section18\final.py:135 in main

    132 │   
    133 │   # Combine summaries with PMIDs
    134 │   combined_summaries = "\n\n".join(
  ❱ 135 │   │   f"PMID: {pmid}\nAbstract: {summary}" for pmid, summary in zip(df['PMID
    136 │   )
    137 │   
    138 │   # Define the LLM and prompt
