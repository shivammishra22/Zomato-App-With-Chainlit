#!/usr/bin/env python3
import ssl
ssl._create_default_https_context = ssl._create_unverified_context

import streamlit as st
import pandas as pd
import time
from datetime import datetime, date
from typing import List, Optional

from Bio import Entrez

# LLM functionality: Ollama chat wrapper (local inference)
try:
    from langchain_ollama import ChatOllama
except ImportError:
    class ChatOllama:
        def __init__(self, model: str, temperature: float = 0, num_ctx: int = 4096):
            self.model = model
            self.temperature = temperature
            self.num_ctx = num_ctx

        def invoke(self, prompt: str) -> str:
            return "Relevance: Yes\nSummary: Mock summary text for testing."

class PubMedExtractor:
    def __init__(self):
        self.email = "adityachannadelhi@gmail.com"
        self.api_key = "45287b47125a7f74898e42afd44238d85a08"
        Entrez.email = self.email
        Entrez.api_key = self.api_key

    def search_pubmed(self, search_term: str, start_date: str, end_date: str) -> List[str]:
        query = search_term + f' AND ("{start_date}"[Date - Publication] : "{end_date}"[Date - Publication])'
        handle = Entrez.esearch(db="pubmed", term=query, retmax=0)
        results = Entrez.read(handle)
        handle.close()
        count = int(results["Count"])
        if count == 0:
            return []
        max_ret = min(count, 10000)
        handle = Entrez.esearch(db="pubmed", term=query, retmax=max_ret, sort="relevance")
        results = Entrez.read(handle)
        handle.close()
        return results.get("IdList", [])

    def fetch_abstracts(self, pmids: List[str]) -> List[dict]:
        if not pmids:
            return []
        articles = []
        BATCH = 50
        for i in range(0, len(pmids), BATCH):
            batch = pmids[i: i + BATCH]
            handle = Entrez.efetch(db="pubmed", id=batch, rettype="medline", retmode="xml")
            records = Entrez.read(handle)
            handle.close()
            for rec in records["PubmedArticle"]:
                art = self._parse_article(rec)
                if art:
                    articles.append(art)
            time.sleep(0.1)
        return [a for a in articles if a.get("Abstract") and a["Abstract"].strip()]

    def _parse_article(self, record) -> Optional[dict]:
        try:
            art = record["MedlineCitation"]["Article"]
            pmid = record["MedlineCitation"]["PMID"]
            title = art.get("ArticleTitle", "N/A")
            abs_segs = art.get("Abstract", {}).get("AbstractText", [])
            abstract = " ".join(str(s) for s in abs_segs) if isinstance(abs_segs, list) else str(abs_segs)
            return {"PMID": str(pmid), "Title": title, "Abstract": abstract}
        except Exception:
            return None

INCLUSIVE_PROMPT = """
Evaluate the abstract based on the following:

1. If there is **no new benefit information**, characterise the baseline benefit data only.
2. If **new benefit data is present**, critically integrate and evaluate it with the baseline data.
3. Consider the following elements when available:
   - Strength of evidence: comparator(s), effect size, statistical significance, consistency, methodology.
   - Validity of endpoints used (especially surrogate endpoints).
   - Clinical relevance of effect size.
   - Generalisability to full indicated population (e.g., any sub-population not showing benefit).
   - Dose-response characterisation.
   - Duration of effect.
   - Comparative efficacy.
   - Whether trial data is applicable to real-world practice.

Abstract:
{abstract}

Respond in this exact format:
Relevance: Yes/No
Summary: <Concise one-sentence evaluation>
""".strip()

EXCLUSION_PROMPT = """
Classify the abstract as EXCLUSION if:
1. No adverse event (AE) with company suspect product is discussed.
2. It refers only to animal/preclinical/in-vitro/ex-vivo studies.
3. There’s no or negative causality with company suspect product.
4. Suspect product is from a non-company (different MAH).
5. No identifiable ICSR or medical relevance.

Abstract:
{abstract}

Respond in this exact format:
Relevance: Yes/No
Summary: <Concise reason for exclusion>
""".strip()

def main():
    st.set_page_config(page_title="PubMed + LLM Inclusive/Exclusive", layout="wide")
    st.title("🔬 PubMed → LLM Relevance (Inclusive/Exclusive)")

    drug_name = st.text_input("Drug name", value="Valsartan")
    col1, col2 = st.columns(2)
    with col1:
        start_date = st.date_input("Start date", value=date(2020, 1, 1))
    with col2:
        end_date = st.date_input("End date", value=date.today())
    search_btn = st.button("🔍 Search & Analyze")

    if search_btn:
        if not drug_name.strip() or start_date > end_date:
            st.error("Check drug name and date inputs.")
            return

        query = (
            f"({drug_name}/drug therapy[mh] OR "
            f"{drug_name}/therapeutic use[mh] OR "
            f"{drug_name}/treatment outcome[mh]) AND "
            "(efficacy[tiab] OR effectiveness[tiab] OR clinical benefit[tiab]) "
            "AND humans[mesh]"
        )

        extractor = PubMedExtractor()
        pmids = extractor.search_pubmed(query, start_date.strftime("%Y/%m/%d"), end_date.strftime("%Y/%m/%d"))
        if not pmids:
            st.warning("No articles found.")
            return

        articles = extractor.fetch_abstracts(pmids)
        if not articles:
            st.warning("No abstracts to analyze.")
            return

        df = pd.DataFrame(articles)[["PMID", "Title", "Abstract"]]
        st.subheader(f"📝 {len(df)} Articles with Abstracts")
        st.dataframe(df, use_container_width=True)

        llm = ChatOllama(model="gemma3:4b", temperature=0, num_ctx=4096)
        relevance, summary = [], []

        for idx, row in df.iterrows():
            abs_text = row["Abstract"]
            inclusive_prompt = INCLUSIVE_PROMPT.format(abstract=abs_text)
            response = llm.invoke(inclusive_prompt)
            if hasattr(response, "content"):
                response = response.content
            elif not isinstance(response, str):
                response = str(response)

            lines = response.split("\n")
            rel = next((l.split(":", 1)[1].strip() for l in lines if l.lower().startswith("relevance:")), "No")
            summ = next((l.split(":", 1)[1].strip() for l in lines if l.lower().startswith("summary:")), "")

            if rel.lower() == "no":
                excl_prompt = EXCLUSION_PROMPT.format(abstract=abs_text)
                excl_resp = llm.invoke(excl_prompt)
                if hasattr(excl_resp, "content"):
                    excl_resp = excl_resp.content
                elif not isinstance(excl_resp, str):
                    excl_resp = str(excl_resp)
                lines = excl_resp.split("\n")
                rel = next((l.split(":", 1)[1].strip() for l in lines if l.lower().startswith("relevance:")), "No")
                summ = next((l.split(":", 1)[1].strip() for l in lines if l.lower().startswith("summary:")), "")

            relevance.append(rel)
            summary.append(summ)

        df["Relevance"] = relevance
        df["Summary"] = summary

        df_rel = df[df["Relevance"].str.lower() == "yes"].reset_index(drop=True)
        st.subheader(f"✅ {len(df_rel)} Relevant Articles")
        st.dataframe(df_rel, use_container_width=True)

        if not df_rel.empty:
            combined = "\n\n".join(df_rel["Abstract"].tolist())
            overall_prompt = f"Provide a concise overall summary of these abstracts:\n\n{combined}"
            ov = llm.invoke(overall_prompt)
            if hasattr(ov, "content"):
                ov = ov.content
            elif not isinstance(ov, str):
                ov = str(ov)
            st.subheader("🗒️ Overall Summary")
            st.write(ov)

            csv_rel = df_rel.to_csv(index=False)
            st.download_button("Download Relevant Articles as CSV", data=csv_rel, file_name=f"relevant_{drug_name}_{datetime.now():%Y%m%d}.csv", mime="text/csv")
            st.download_button("Download Overall Summary", data=ov, file_name=f"summary_{drug_name}_{datetime.now():%Y%m%d}.txt", mime="text/plain")

if __name__ == "__main__":
    main()
