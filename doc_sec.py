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
            # Fallback mock response
            return (
                "Relevance: Yes\n"
                "Summary: This article shows strong dose–response data with clear effect sizes "
                "and methodological rigor supporting efficacy."
            )

class PubMedExtractor:
    def __init__(self):
        # Replace with your own email & API key if desired
        self.email = "adityachannadelhi@gmail.com"
        self.api_key = "45287b47125a7f74898e42afd44238d85a08"
        Entrez.email = self.email
        Entrez.api_key = self.api_key

    def search_pubmed(
        self,
        search_term: str,
        start_date: str,
        end_date: str
    ) -> List[str]:
        query = search_term + f' AND ("{start_date}"[Date - Publication] : "{end_date}"[Date - Publication])'
        st.info(f"🔍 Searching PubMed: {query}")
        handle = Entrez.esearch(db="pubmed", term=query, retmax=0)
        results = Entrez.read(handle); handle.close()
        count = int(results["Count"])
        st.info(f"Found {count} articles")
        if count == 0:
            return []
        max_ret = min(count, 10000)
        handle = Entrez.esearch(
            db="pubmed",
            term=query,
            retmax=max_ret,
            sort="relevance"
        )
        results = Entrez.read(handle); handle.close()
        return results.get("IdList", [])

    def fetch_abstracts(self, pmids: List[str]) -> List[dict]:
        if not pmids:
            return []
        articles = []
        BATCH = 50
        pbar = st.progress(0)
        status = st.empty()
        for i in range(0, len(pmids), BATCH):
            batch = pmids[i : i + BATCH]
            status.text(f"Fetching {i+1}–{min(i+BATCH, len(pmids))}...")
            handle = Entrez.efetch(
                db="pubmed",
                id=batch,
                rettype="medline",
                retmode="xml"
            )
            records = Entrez.read(handle); handle.close()
            for rec in records["PubmedArticle"]:
                art = self._parse_article(rec)
                if art:
                    articles.append(art)
            time.sleep(0.1)
            pbar.progress(min((i+BATCH)/len(pmids), 1.0))
        # drop articles without abstracts
        filtered = [
            a for a in articles
            if a.get("Abstract") and a["Abstract"].strip()
        ]
        removed = len(articles) - len(filtered)
        if removed:
            st.info(f"Removed {removed} articles with empty abstracts")
        return filtered

    def _parse_article(self, record) -> Optional[dict]:
        try:
            art = record["MedlineCitation"]["Article"]
            pmid = record["MedlineCitation"]["PMID"]
            title = art.get("ArticleTitle", "N/A")
            abs_segs = art.get("Abstract", {}).get("AbstractText", [])
            if isinstance(abs_segs, list):
                abstract = " ".join(str(s) for s in abs_segs)
            else:
                abstract = str(abs_segs) if abs_segs else ""
            return {"PMID": str(pmid), "Title": title, "Abstract": abstract}
        except Exception:
            return None

# Inclusive prompt template
INCLUSIVE_PROMPT = """
Find the output from the given abstract based on the following instructions:
a) Dose–response characterisation.
b) Duration of effect.
c) Comparative efficacy.
d) Provide a concise but critical evaluation of the strengths and limitations of the evidence.
e) Adequacy of dose–response characterisation.
f) Strength of evidence of benefit, including comparator(s), effect size, statistical rigor, methodological strengths and deficiencies, and consistency across studies.
g) Clinical relevance of the effect size.
h) Generalisability of treatment response across the indicated patient population, including sub‑populations.

Abstract:
{abstract}

Respond EXACTLY in this format:
Relevance: Yes/No
Summary: <one‑sentence critical evaluation>
""".strip()

def main():
    st.set_page_config(page_title="PubMed + LLM Analysis", layout="wide")
    st.title("🔬 PubMed → LLM Relevance & Summary")

    # Input controls
    drug_name = st.text_input(
        "Drug name",
        value="Valsartan",
        help="Enter the drug name to build the PubMed query"
    )
    col1, col2 = st.columns(2)
    with col1:
        start_date = st.date_input("Start date", value=date(2020,1,1))
    with col2:
        end_date = st.date_input("End date", value=date.today())
    search_btn = st.button("🔍 Search & Analyze")

    if search_btn:
        if not drug_name.strip():
            st.error("Please enter a drug name.")
            return
        if start_date > end_date:
            st.error("Start date must be before end date.")
            return

        # Build dynamic query
        q = (
            f"({drug_name}/drug therapy[mh] OR "
            f"{drug_name}/therapeutic use[mh] OR "
            f"{drug_name}/treatment outcome[mh]) AND "
            "(efficacy[tiab] OR effectiveness[tiab] OR clinical benefit[tiab]) "
            "AND humans[mesh]"
        )

        extractor = PubMedExtractor()
        pmids = extractor.search_pubmed(
            q,
            start_date.strftime("%Y/%m/%d"),
            end_date.strftime("%Y/%m/%d")
        )
        if not pmids:
            st.warning("No articles found.")
            return

        articles = extractor.fetch_abstracts(pmids)
        if not articles:
            st.warning("No abstracts to analyze.")
            return

        # DataFrame & drop empties
        df = pd.DataFrame(articles)[["PMID","Title","Abstract"]]
        df = df.dropna(subset=["Abstract"])
        df = df[df["Abstract"].str.strip() != ""]
        st.subheader(f"📝 {len(df)} Articles with Abstracts")
        st.dataframe(df, use_container_width=True)

        # LLM relevance + summary
        llm = ChatOllama(model="gemma3:4b", temperature=0, num_ctx=4096)
        relevance, summary = [], []
        for idx, row in df.iterrows():
            prompt = INCLUSIVE_PROMPT.format(abstract=row["Abstract"])
            with st.spinner(f"Evaluating {idx+1}/{len(df)}..."):
                resp = llm.invoke(prompt)
                if hasattr(resp, "content"):
                    resp = resp.content
                elif not isinstance(resp, str):
                    resp = str(resp)
            lines = resp.split("\n")
            r = next((l.split(":",1)[1].strip() for l in lines if l.lower().startswith("relevance:")), "No")
            s = next((l.split(":",1)[1].strip() for l in lines if l.lower().startswith("summary:")), "")
            relevance.append(r)
            summary.append(s)

        df["Relevance"] = relevance
        df["Summary"]   = summary

        # Show relevant only
        df_rel = df[df["Relevance"].str.lower() == "yes"].reset_index(drop=True)
        st.subheader(f"✅ {len(df_rel)} Relevant Articles")
        st.dataframe(df_rel, use_container_width=True)

        # Overall summary of relevant abstracts
        if not df_rel.empty:
            combined = "\n\n".join(df_rel["Abstract"].tolist())
            overall_prompt = f"Provide a concise overall summary of these abstracts:\n\n{combined}"
            with st.spinner("Generating overall summary..."):
                ov = llm.invoke(overall_prompt)
                if hasattr(ov, "content"):
                    ov = ov.content
                elif not isinstance(ov, str):
                    ov = str(ov)
            st.subheader("🗒️ Overall Summary")
            st.write(ov)

            # Download buttons
            csv_rel = df_rel.to_csv(index=False)
            st.download_button(
                "Download Relevant Articles as CSV",
                data=csv_rel,
                file_name=f"relevant_{drug_name}_{datetime.now():%Y%m%d}.csv",
                mime="text/csv"
            )
            st.download_button(
                "Download Overall Summary",
                data=ov,
                file_name=f"summary_{drug_name}_{datetime.now():%Y%m%d}.txt",
                mime="text/plain"
            )

if __name__ == "__main__":
    main()
