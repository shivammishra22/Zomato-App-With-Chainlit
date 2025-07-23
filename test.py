#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import ssl
ssl._create_default_https_context = ssl._create_unverified_context  # <-- bypass SSL cert verify

import streamlit as st
import io
import time
import pandas as pd
from datetime import datetime, date
from typing import List, Dict, Optional, Tuple

try:
    from Bio import Entrez
except ImportError:
    st.error("Error: Biopython not installed. Run:  pip install biopython")
    st.stop()


class PubMedExtractor:
    def __init__(self, email: str, api_key: Optional[str] = None):
        """Initialize PubMed extractor"""
        Entrez.email = email
        if api_key:
            Entrez.api_key = api_key

    def search_pubmed(
        self,
        search_term: str,
        start_date: Optional[str] = None,
        end_date: Optional[str] = None,
        journal: Optional[str] = None,
        species: Optional[str] = None,
        publication_type: Optional[str] = None,
        language: Optional[str] = None,
        age_group: Optional[str] = None,
        free_full_text: bool = False,
        author: Optional[str] = None,
        mesh_terms: Optional[str] = None,
    ) -> Tuple[List[str], str, int]:
        """Search PubMed and return list of PMIDs, final query, and total count"""

        query_parts = [search_term]

        # Date range filter
        if start_date and end_date:
            query_parts.append(f'("{start_date}"[Date - Publication] : "{end_date}"[Date - Publication])')
        elif start_date:
            query_parts.append(f'"{start_date}"[Date - Publication] : 3000[Date - Publication]')
        elif end_date:
            query_parts.append(f'1900[Date - Publication] : "{end_date}"[Date - Publication]')

        # Journal filter
        if journal:
            query_parts.append(f'"{journal}"[Journal]')

        # Species filter
        if species and species.lower() != "all":
            if species.lower() == "humans":
                query_parts.append('"humans"[MeSH Terms]')
            elif species.lower() == "animals":
                query_parts.append('"animals"[MeSH Terms]')

        # Publication type filter
        if publication_type and publication_type != "All":
            query_parts.append(f'"{publication_type}"[Publication Type]')

        # Language filter
        if language and language != "All":
            query_parts.append(f'"{language}"[Language]')

        # Age group filter
        if age_group and age_group != "All":
            query_parts.append(f'"{age_group}"[MeSH Terms]')

        # Free full text filter
        if free_full_text:
            query_parts.append('"free full text"[sb]')

        # Author filter
        if author:
            query_parts.append(f'"{author}"[Author]')

        # MeSH terms filter
        if mesh_terms:
            terms = [term.strip() for term in mesh_terms.split(",")]
            for term in terms:
                query_parts.append(f'"{term}"[MeSH Terms]')

        final_query = " AND ".join(query_parts)

        try:
            # First call just for count
            handle = Entrez.esearch(db="pubmed", term=final_query, retmax=0)
            search_results = Entrez.read(handle)
            handle.close()

            total_count = int(search_results["Count"])
            if total_count == 0:
                return [], final_query, 0

            max_retrievable = min(total_count, 10000)

            # Second call to actually get ids
            handle = Entrez.esearch(
                db="pubmed",
                term=final_query,
                retmax=max_retrievable,
                sort="relevance"
            )
            search_results = Entrez.read(handle)
            handle.close()

            pmids = search_results.get("IdList", [])
            return pmids, final_query, total_count

        except Exception as e:
            st.error(f"Error searching PubMed: {e}")
            return [], "", 0

    def fetch_abstracts(self, pmids: List[str]) -> List[Dict]:
        """Fetch detailed info including abstracts for given PMIDs"""
        if not pmids:
            return []

        articles = []
        batch_size = 500

        progress_bar = st.progress(0)
        status_text = st.empty()

        for i in range(0, len(pmids), batch_size):
            batch_pmids = pmids[i:i + batch_size]
            progress = (i + batch_size) / max(len(pmids), 1)
            progress_bar.progress(min(progress, 1.0))
            status_text.text(f"Fetching articles {i + 1} to {min(i + batch_size, len(pmids))} of {len(pmids)}...")

            try:
                handle = Entrez.efetch(
                    db="pubmed",
                    id=batch_pmids,
                    rettype="medline",
                    retmode="xml"
                )
                records = Entrez.read(handle)
                handle.close()

                for record in records.get("PubmedArticle", []):
                    article_info = self._parse_article(record)
                    if article_info:
                        articles.append(article_info)

                time.sleep(0.1)

            except Exception as e:
                st.warning(f"Error fetching batch {i // batch_size + 1}: {e}")
                continue

        progress_bar.progress(1.0)
        status_text.text(f"Successfully extracted information for {len(articles)} articles")

        return articles

    def _parse_article(self, record) -> Optional[Dict]:
        """Parse individual PubMed article record"""
        try:
            medline = record.get("MedlineCitation", {})
            article = medline.get("Article", {})
            pmid = medline.get("PMID", "")

            # Title
            title = article.get("ArticleTitle", "N/A")

            # Abstract
            abstract_sections = article.get("Abstract", {}).get("AbstractText", [])
            if isinstance(abstract_sections, list):
                abstract = " ".join([str(section) for section in abstract_sections])
            else:
                abstract = str(abstract_sections) if abstract_sections else "N/A"

            # Authors
            authors = []
            for auth in article.get("AuthorList", []):
                if "LastName" in auth and "ForeName" in auth:
                    authors.append(f"{auth['LastName']}, {auth['ForeName']}")
                elif "CollectiveName" in auth:
                    authors.append(auth["CollectiveName"])
            authors_str = "; ".join(authors) if authors else "N/A"

            # Journal info
            journal_info = article.get("Journal", {})
            journal_title = journal_info.get("Title", "N/A")

            # Publication date
            pub_date = "N/A"
            pub_year = "N/A"
            if "JournalIssue" in journal_info and "PubDate" in journal_info["JournalIssue"]:
                date_info = journal_info["JournalIssue"]["PubDate"]
                year = str(date_info.get("Year", "")).strip()
                month = str(date_info.get("Month", "")).strip()
                day = str(date_info.get("Day", "")).strip()
                # Normalize month/day if missing
                pub_date = "-".join([x for x in [year, month, day] if x])
                pub_year = year if year else "N/A"

            # DOI
            doi = "N/A"
            article_ids = record.get("PubmedData", {}).get("ArticleIdList", [])
            for article_id in article_ids:
                if article_id.attributes.get("IdType") == "doi":
                    doi = str(article_id)
                    break

            # Keywords
            keywords = []
            for kw_group in medline.get("KeywordList", []):
                for kw in kw_group:
                    keywords.append(str(kw))
            keywords_str = "; ".join(keywords) if keywords else "N/A"

            return {
                "PMID": str(pmid),
                "Title": title,
                "Abstract": abstract,
                "Authors": authors_str,
                "Journal": journal_title,
                "Publication_Date": pub_date,
                "Publication_Year": pub_year,   # <-- added
                "DOI": doi,
                "Keywords": keywords_str,
                "PubMed_URL": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
            }

        except Exception:
            return None


def main():
    st.set_page_config(page_title="PubMed Search Tool", page_icon="🔬", layout="wide")

    st.title("🔬 PubMed Abstract Extractor")
    st.markdown("Search PubMed and extract abstracts with comprehensive filtering options.")

    # Sidebar
    with st.sidebar:
        st.header("⚙️ Configuration")
        email = st.text_input("Email Address*", help="Required by NCBI for API access")
        api_key = st.text_input("NCBI API Key", type="password", help="Optional: For higher rate limits")

        if not email:
            st.warning("Please enter your email address to proceed")
            st.stop()

    # Main input
    col1, col2 = st.columns([2, 1])

    with col1:
        st.header("🔍 Search Query")
        search_term = st.text_area(
            "Search Terms*",
            placeholder="e.g., (diabetes) AND (treatment) NOT (case report)",
            help="Use PubMed syntax: AND, OR, NOT",
        )

    with col2:
        st.header("📅 Date Range")
        c1, c2 = st.columns(2)
        with c1:
            start_date = st.date_input("Start Date", value=None)
        with c2:
            end_date = st.date_input("End Date", value=None)

    # Advanced filters
    with st.expander("🔧 Advanced Filters", expanded=False):
        c1, c2, c3 = st.columns(3)

        with c1:
            journal = st.text_input("Journal Name", placeholder="e.g., Nature")
            author = st.text_input("Author Name", placeholder="e.g., Smith J")
            mesh_terms = st.text_input("MeSH Terms", placeholder="e.g., Diabetes Mellitus, Hypertension")

        with c2:
            species = st.selectbox("Species", ["All", "Humans", "Animals"])
            publication_type = st.selectbox(
                "Publication Type",
                [
                    "All",
                    "Clinical Trial",
                    "Randomized Controlled Trial",
                    "Review",
                    "Meta-Analysis",
                    "Case Reports",
                    "Observational Study",
                    "Letter",
                    "Editorial",
                ],
            )
            language = st.selectbox(
                "Language",
                ["All", "English", "Spanish", "French", "German", "Italian", "Japanese"],
            )

        with c3:
            age_group = st.selectbox(
                "Age Group",
                ["All", "Infant", "Child", "Adolescent", "Adult", "Middle Aged", "Aged"],
            )
            free_full_text = st.checkbox("Free Full Text Only")

            st.subheader("Result Limits")
            return_all = st.checkbox(
                "Return All Available Articles",
                help="Get all articles found (up to PubMed's 10,000 limit)",
                value=False,
            )
            if return_all:
                max_results = 10000
                st.info("Will return all available articles (up to 10,000)")
            else:
                max_results = st.number_input("Max Results", 1, 10000, 1000)

        # Extra options
        remove_empty_abstracts = st.checkbox("Remove rows with empty 'Abstract'", value=True)

    # Preview search
    if search_term:
        st.subheader("🔍 Search Preview")
        preview_parts = [search_term]

        if start_date and end_date:
            preview_parts.append(
                f'("{start_date.strftime("%Y/%m/%d")}"[Date - Publication] : "{end_date.strftime("%Y/%m/%d")}"[Date - Publication])'
            )
        elif start_date:
            preview_parts.append(f'"{start_date.strftime("%Y/%m/%d")}"[Date - Publication] : 3000[Date - Publication]')
        elif end_date:
            preview_parts.append(f'1900[Date - Publication] : "{end_date.strftime("%Y/%m/%d")}"[Date - Publication]')

        if journal:
            preview_parts.append(f'"{journal}"[Journal]')
        if species and species != "All":
            preview_parts.append(f'"{species.lower()}"[MeSH Terms]')
        if publication_type and publication_type != "All":
            preview_parts.append(f'"{publication_type}"[Publication Type]')
        if language and language != "All":
            preview_parts.append(f'"{language}"[Language]')
        if age_group and age_group != "All":
            preview_parts.append(f'"{age_group}"[MeSH Terms]')
        if free_full_text:
            preview_parts.append('"free full text"[sb]')
        if author:
            preview_parts.append(f'"{author}"[Author]')
        if mesh_terms:
            for t in [m.strip() for m in mesh_terms.split(",")]:
                preview_parts.append(f'"{t}"[MeSH Terms]')

        st.code(" AND ".join(preview_parts), language="text")

        if return_all:
            st.info("📊 Will retrieve all available articles (up to 10,000)")
        else:
            st.info(f"📊 Will retrieve up to {max_results} articles")

    # Search button
    disabled = not (search_term and email)
    if st.button("🔍 Search PubMed", type="primary", disabled=disabled):
        extractor = PubMedExtractor(email, api_key)

        start_date_str = start_date.strftime("%Y/%m/%d") if start_date else None
        end_date_str = end_date.strftime("%Y/%m/%d") if end_date else None

        with st.spinner("Searching PubMed..."):
            pmids, query, total_count = extractor.search_pubmed(
                search_term=search_term,
                start_date=start_date_str,
                end_date=end_date_str,
                journal=journal or None,
                species=species if species != "All" else None,
                publication_type=publication_type if publication_type != "All" else None,
                language=language if language != "All" else None,
                age_group=age_group if age_group != "All" else None,
                free_full_text=free_full_text,
                author=author or None,
                mesh_terms=mesh_terms or None,
            )

        if not pmids:
            st.warning("No articles found matching your criteria.")
            st.info(f"**Search Query Used:** {query}")
            st.stop()

        # Summary
        st.success(f"Found {total_count} articles (retrieving up to {min(len(pmids), max_results)})")
        st.info(f"**Search Query Used:** {query}")

        original_count = len(pmids)
        if return_all:
            if total_count > 10000:
                st.warning(
                    f"⚠️ Found {total_count} articles, but PubMed limits results to 10,000. Retrieving available {len(pmids)}."
                )
            else:
                st.info(f"✅ Retrieving all {len(pmids)} articles found.")
        else:
            if len(pmids) > max_results:
                pmids = pmids[:max_results]
                st.info(f"📊 Limited to first {max_results} results (out of {original_count} available)")
            else:
                st.info(f"✅ Retrieving all {len(pmids)} articles (less than your {max_results} limit)")

        with st.spinner("Fetching article details..."):
            articles = extractor.fetch_abstracts(pmids)

        if not articles:
            st.error("No article details could be fetched.")
            st.stop()

        # DataFrame
        df = pd.DataFrame(articles)

        if remove_empty_abstracts:
            before = len(df)
            df = df[df["Abstract"].str.strip().ne("N/A") & df["Abstract"].str.strip().ne("")]
            st.info(f"🧹 Removed {before - len(df)} rows with empty abstracts.")

        # Metrics
        st.header("📊 Results")
        c1, c2, c3, c4, c5 = st.columns(5)
        with c1:
            st.metric("Total Articles", len(articles))
        with c2:
            st.metric("After Clean", len(df))
        with c3:
            st.metric("Unique Journals", df["Journal"].nunique())
        with c4:
            st.metric("With Abstracts", (df["Abstract"].str.strip().ne("N/A")).sum())
        with c5:
            st.metric("With DOI", (df["DOI"].str.strip().ne("N/A")).sum())

        # Display table (shortened)
        st.subheader("📋 Article Details (truncated)")
        display_df = df.copy()
        display_df["Title"] = display_df["Title"].str.slice(0, 100) + "..."
        display_df["Abstract"] = display_df["Abstract"].str.slice(0, 150) + "..."
        st.dataframe(
            display_df[
                ["PMID", "Title", "Authors", "Publication_Year", "Journal", "Publication_Date"]
            ],
            use_container_width=True,
            hide_index=True
        )

        # Download CSV
        st.subheader("💾 Download Results")
        csv_buffer = io.StringIO()
        df.to_csv(csv_buffer, index=False)
        csv_string = csv_buffer.getvalue()

        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        result_type = "all" if return_all else f"top{max_results}"
        filename = f"pubmed_results_{result_type}_{timestamp}.csv"

        st.download_button(
            label="📥 Download CSV",
            data=csv_string,
            file_name=filename,
            mime="text/csv",
            type="primary"
        )

        with st.expander("👀 Preview Full DataFrame"):
            st.dataframe(df, use_container_width=True)


if __name__ == "__main__":
    main()
            
