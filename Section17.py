
#!/usr/bin/env python3
"""
PubMed Abstract Extractor Interface
Streamlit web interface for extracting PubMed abstracts with caching
"""
import ssl
ssl._create_default_https_context = ssl._create_unverified_context
import streamlit as st
import pandas as pd
import csv
import time
import sys
from datetime import datetime, date
from typing import List, Dict, Optional

# LLM functionality - replace with actual implementation or use Ollama
try:
    from langchain_ollama import ChatOllama
except ImportError:
    # Fallback ChatOllama implementation when langchain_ollama is not available
    class ChatOllama:
        def __init__(self, model: str, temperature: float = 0, num_ctx: int = 4096):
            self.model = model
            self.temperature = temperature
            self.num_ctx = num_ctx
        
        def invoke(self, prompt: str) -> str:
            # Mock response for demonstration when Ollama is not available
            return "Relevance: Yes\nSummary: This article contains new safety findings that require pharmacovigilance assessment.\nThe authors express concerns about potential safety implications and recommend continued monitoring.\nThese findings are highly relevant for PSUR reporting and drug safety evaluation."

try:
    from Bio import Entrez
except ImportError:
    st.error("Error: Biopython not installed. Install with: pip install biopython")
    st.stop()


class PubMedExtractor:
    def __init__(self):
        """
        Initialize PubMed extractor with hardcoded credentials
        """
        # Hardcoded credentials
        self.email = "adityachannadelhi@gmail.com"
        self.api_key = "45287b47125a7f74898e42afd44238d85a08"
        
        Entrez.email = self.email
        Entrez.api_key = self.api_key
    
    def search_pubmed(self, 
                     search_term: str,
                     start_date: str,
                     end_date: str,
                     journal: Optional[str] = None,
                     species: Optional[str] = None) -> List[str]:
        """
        Search PubMed and return list of ALL PMIDs matching criteria
        """
        # Build search query
        query_parts = [search_term]
        
        # Add date range filter
        query_parts.append(f'("{start_date}"[Date - Publication] : "{end_date}"[Date - Publication])')
        
        # Add journal filter
        if journal and journal.strip():
            query_parts.append(f'"{journal}"[Journal]')
        
        # Add species filter
        if species and species != "All":
            if species.lower() == 'humans':
                query_parts.append('"humans"[MeSH Terms]')
            elif species.lower() == 'animals':
                query_parts.append('"animals"[MeSH Terms]')
        
        # Combine all query parts
        final_query = ' AND '.join(query_parts)
        
        st.info(f"Searching PubMed with query: {final_query}")
        
        try:
            # First, get the total count
            handle = Entrez.esearch(db="pubmed", 
                                  term=final_query, 
                                  retmax=0)
            search_results = Entrez.read(handle)
            handle.close()
            
            total_count = int(search_results["Count"])
            st.info(f"Found {total_count} total articles")
            
            if total_count == 0:
                return []
            
            # PubMed has a limit of 10,000 results per search
            max_retrievable = min(total_count, 10000)
            
            if total_count > 10000:
                st.warning(f"PubMed limits results to 10,000. Retrieving first {max_retrievable} articles.")
            
            # Get all available PMIDs
            handle = Entrez.esearch(db="pubmed", 
                                  term=final_query, 
                                  retmax=max_retrievable,
                                  sort="relevance")
            search_results = Entrez.read(handle)
            handle.close()
            
            pmids = search_results["IdList"]
            st.success(f"Retrieved {len(pmids)} article IDs")
            return pmids
            
        except Exception as e:
            st.error(f"Error searching PubMed: {e}")
            return []
    
    def fetch_abstracts(self, pmids: List[str]) -> List[Dict]:
        """
        Fetch detailed information including abstracts for given PMIDs
        """
        if not pmids:
            return []
        
        articles = []
        batch_size = 50 # Process in batches to avoid overwhelming the server
        
        progress_bar = st.progress(0)
        status_text = st.empty()
        
        for i in range(0, len(pmids), batch_size):
            batch_pmids = pmids[i:i+batch_size]
            progress = (i + batch_size) / len(pmids)
            progress_bar.progress(min(progress, 1.0))
            status_text.text(f"Fetching articles {i+1} to {min(i+batch_size, len(pmids))}...")
            
            try:
                # Fetch detailed records
                handle = Entrez.efetch(db="pubmed", 
                                     id=batch_pmids, 
                                     rettype="medline", 
                                     retmode="xml")
                records = Entrez.read(handle)
                handle.close()
                
                # Parse each record
                for record in records['PubmedArticle']:
                    article_info = self._parse_article(record)
                    if article_info:
                        articles.append(article_info)
                
                # Be nice to NCBI servers
                time.sleep(0.1)
                
            except Exception as e:
                st.error(f"Error fetching batch {i//batch_size + 1}: {e}")
                continue
        
        progress_bar.progress(1.0)
        status_text.text(f"Successfully extracted information for {len(articles)} articles")
        
        # Filter out articles with no abstracts
        articles_with_abstracts = [article for article in articles if article.get('Abstract', 'N/A') != 'N/A' and article.get('Abstract', '').strip()]
        if len(articles_with_abstracts) < len(articles):
            removed_count = len(articles) - len(articles_with_abstracts)
            st.info(f"Removed {removed_count} articles with no abstracts. {len(articles_with_abstracts)} articles remaining.")
        
        return articles_with_abstracts
    
    def _parse_article(self, record) -> Optional[Dict]:
        """
        Parse individual PubMed article record
        """
        try:
            article = record['MedlineCitation']['Article']
            pmid = record['MedlineCitation']['PMID']
            
            # Extract basic information
            title = article.get('ArticleTitle', 'N/A')
            
            # Extract abstract
            abstract_sections = article.get('Abstract', {}).get('AbstractText', [])
            if isinstance(abstract_sections, list):
                abstract = ' '.join([str(section) for section in abstract_sections])
            else:
                abstract = str(abstract_sections) if abstract_sections else 'N/A'
            
            # Extract authors
            authors = []
            author_list = article.get('AuthorList', [])
            for author in author_list:
                if 'LastName' in author and 'ForeName' in author:
                    authors.append(f"{author['LastName']}, {author['ForeName']}")
                elif 'CollectiveName' in author:
                    authors.append(author['CollectiveName'])
            authors_str = '; '.join(authors) if authors else 'N/A'
            
            # Extract journal information
            journal_info = article.get('Journal', {})
            journal_title = journal_info.get('Title', 'N/A')
            
            # Extract publication date
            pub_date = 'N/A'
            if 'JournalIssue' in journal_info and 'PubDate' in journal_info['JournalIssue']:
                date_info = journal_info['JournalIssue']['PubDate']
                year = date_info.get('Year', '')
                month = date_info.get('Month', '')
                day = date_info.get('Day', '')
                pub_date = f"{year}-{month}-{day}".strip('-')
            
            # Extract DOI if available
            doi = 'N/A'
            article_ids = record.get('PubmedData', {}).get('ArticleIdList', [])
            for article_id in article_ids:
                if article_id.attributes.get('IdType') == 'doi':
                    doi = str(article_id)
                    break
            
            # Extract keywords
            keywords = []
            keyword_list = record['MedlineCitation'].get('KeywordList', [])
            for kw_group in keyword_list:
                for keyword in kw_group:
                    keywords.append(str(keyword))
            keywords_str = '; '.join(keywords) if keywords else 'N/A'
            
            return {
                'PMID': str(pmid),
                'Title': title,
                'Abstract': abstract,
                'Authors': authors_str,
                'Journal': journal_title,
                'Publication_Date': pub_date,
                'DOI': doi,
                'Keywords': keywords_str,
                'PubMed_URL': f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
            }
            
        except Exception as e:
            st.error(f"Error parsing article: {e}")
            return None

def main():
    st.set_page_config(
        page_title="PubMed Abstract Extractor",
        page_icon="🔬",
        layout="wide"
    )
    
    st.title("🔬 PubMed Abstract Extractor")
    st.markdown("Extract abstracts from PubMed with advanced filtering and predefined drug queries")
    
    # Initialize extractor
    if 'extractor' not in st.session_state:
        st.session_state.extractor = PubMedExtractor()
    
    # Center the search parameters
    st.markdown("---")
    st.header("Search Parameters")
    
    # Create centered columns for better layout
    col1, col2, col3 = st.columns([1, 2, 1])
    
    with col2:
        # Drug name input for predefined queries
        drug_name = st.text_input(
            "Drug Name (for predefined queries)",
            value="TADALAFIL",
            help="Enter drug name to generate predefined search queries automatically"
        )
        
        # Toggle between manual and predefined queries
        query_mode = st.radio(
            "Query Mode",
            ["Predefined Drug Queries", "Manual Query"],
            help="Choose to use predefined queries based on drug name or enter manual query"
        )
        
        if query_mode == "Manual Query":
            # Manual search term input
            search_term = st.text_area(
                "Search Term/Query",
                value='("TADALAFIL") NOT (ICSR OR CASE REPORT)',
                help="Enter your PubMed search query. Use quotes for exact phrases and Boolean operators (AND, OR, NOT).",
                height=100
            )
        else:
            # Define questions and their corresponding queries
            questions = [
                "Literature articles on pregnancy outcomes (including termination) with/ without adverse outcomes",
                "Literature articles with use in paediatric populations", 
                "Literature articles with compassionate supply, named patient use",
                "Literature articles on lack of efficacy",
                "Literature articles on asymptomatic overdose, abuse or misuse",
                "Literature articles with medication error where no adverse events occurred",
                "Literature articles with important non-clinical safety results",
                "Other relevant literature articles"
            ]
            
            predefined_queries = [
                f'{drug_name} AND pregnan* AND humans[MeSH Terms]',
                f'{drug_name} AND paediatric AND humans[MeSH Terms]',
                f'{drug_name} AND (Compassionate Use Trials OR supply)',
                f'{drug_name} AND Treatment Failure AND humans[MeSH Terms]',
                f'{drug_name} AND (overdose OR abuse OR misuse) AND humans[MeSH Terms]',
                f'{drug_name} AND medication error AND humans[MeSH Terms]',
                f'{drug_name} AND In Vitro',
                f'{drug_name} AND (Pharmacogenetics OR drug interactions OR Pharmacovigilance OR ("Risk Factors") OR ("Adverse Events") OR ("Adverse effects")) AND humans[MeSH Terms]'
            ]
            
            # Allow users to select which questions to search
            st.subheader("Select Research Questions to Execute:")
            
            # Option to select all
            select_all = st.checkbox("Select All Questions", value=True)
            
            selected_questions = []
            if select_all:
                selected_questions = list(range(len(questions)))
                st.success("All 8 research questions will be executed.")
            else:
                for i, question in enumerate(questions):
                    if st.checkbox(f"Q{i+1}: {question}", key=f"q_{i}"):
                        selected_questions.append(i)
                
                # Only show queries if not all are selected
                if selected_questions:
                    with st.expander("Show Selected Query Details", expanded=False):
                        st.info("The following queries will be executed:")
                        for i in selected_questions:
                            st.text(f"Q{i+1}: {predefined_queries[i]}")
                            st.caption(f"Question: {questions[i]}")
                            st.markdown("---")
                else:
                    st.warning("Please select at least one question to search.")
            
            search_term = None  # Will be handled differently for predefined queries
        
        # Date range inputs
        st.subheader("Date Range")
        date_col1, date_col2 = st.columns(2)
        
        with date_col1:
            start_date = st.date_input(
                "Start Date",
                value=date(2017, 7, 27),
                min_value=date(1900, 1, 1),
                max_value=date.today()
            )
        
        with date_col2:
            end_date = st.date_input(
                "End Date",
                value=date.today(),
                min_value=date(1900, 1, 1),
                max_value=date.today()
            )
        
        # Additional filters
        st.subheader("Additional Filters")
        
        filter_col1, filter_col2 = st.columns(2)
        
        with filter_col1:
            species = st.selectbox(
                "Species Filter",
                options=["All", "Humans", "Animals"],
                index=1
            )
        
        with filter_col2:
            journal = st.text_input(
                "Specific Journal (Optional)",
                help="Enter a specific journal name to filter results"
            )
        
        # Center the search button
        st.markdown("</br>", unsafe_allow_html=True)
        search_col1, search_col2, search_col3 = st.columns([1, 1, 1])
        with search_col2:
            search_button = st.button("🔍 Search PubMed", type="primary", use_container_width=True)
    
    # Search logic
    if search_button:
        if query_mode == "Manual Query" and not search_term.strip():
            st.error("Please enter a search term")
            return
        elif query_mode == "Predefined Drug Queries" and not drug_name.strip():
            st.error("Please enter a drug name")
            return
        elif query_mode == "Predefined Drug Queries" and 'selected_questions' not in locals():
            st.error("Please select at least one question to search")
            return
        
        if start_date > end_date:
            st.error("Start date must be before end date")
            return
        
        # Convert dates to string format
        start_date_str = start_date.strftime("%Y/%m/%d")
        end_date_str = end_date.strftime("%Y/%m/%d")
        
        if query_mode == "Manual Query":
            # Single query search
            with st.spinner("Searching PubMed..."):
                pmids = st.session_state.extractor.search_pubmed(
                    search_term=search_term,
                    start_date=start_date_str,
                    end_date=end_date_str,
                    journal=journal if journal.strip() else None,
                    species=species if species != "All" else None
                )
                
                if not pmids:
                    st.warning("No articles found matching your criteria")
                    return
                
                # Fetch abstracts
                with st.spinner("Fetching article details..."):
                    articles = st.session_state.extractor.fetch_abstracts(pmids)
                    
                    if articles:
                        # Convert to DataFrame with origin tracking
                        df_results = pd.DataFrame(articles)
                        df_results = df_results[['PMID', 'Title', 'Abstract']].copy()
                        df_results['Query_Origin'] = 'Manual Query'
                        st.session_state.articles = df_results.to_dict('records')
                        st.session_state.combined_df = df_results
                        st.session_state.current_drug_name = drug_name  # Store drug name
                    else:
                        st.error("Failed to fetch article details")
                        return
        else:
            # Predefined queries search
            questions = [
                "Literature articles on pregnancy outcomes (including termination) with/ without adverse outcomes",
                "Literature articles with use in paediatric populations", 
                "Literature articles with compassionate supply, named patient use",
                "Literature articles on lack of efficacy",
                "Literature articles on asymptomatic overdose, abuse or misuse",
                "Literature articles with medication error where no adverse events occurred",
                "Literature articles with important non-clinical safety results",
                "Other relevant literature articles"
            ]
            
            predefined_queries = [
                f'{drug_name} AND pregnan* AND humans[MeSH Terms]',
                f'{drug_name} AND paediatric AND humans[MeSH Terms]',
                f'{drug_name} AND (Compassionate Use)',
                f'{drug_name} AND Treatment failure AND humans[MeSH Terms]',
                f'{drug_name} AND (overdose OR abuse OR misuse) AND humans[MeSH Terms]',
                f'{drug_name} AND medication error AND humans[MeSH Terms]',
                f'{drug_name} AND In Vitro',
                f'{drug_name} AND (Pharmacogenetics OR drug interactions OR Pharmacovigilance OR ("Risk Factors") OR ("Adverse Events") OR ("Adverse effects")) AND humans[MeSH Terms]'
            ]
            
            query_labels = [
                'Pregnancy Outcomes',
                'Paediatric Use',
                'Compassionate Use/Supply',
                'efficacy',
                'Overdose/Abuse/Misuse',
                'Medication Error',
                'In Vitro Studies',
                'others'
            ]
            
            all_results = []
            combined_df = pd.DataFrame()
            
            # Only process selected questions
            selected_queries = [(predefined_queries[i], query_labels[i], i) for i in selected_questions]
            progress_bar = st.progress(0)
            
            for idx, (query, label, original_idx) in enumerate(selected_queries):
                progress_bar.progress((idx + 1) / len(selected_queries))
                
                with st.spinner(f"Searching: {label}..."):
                    pmids = st.session_state.extractor.search_pubmed(
                        search_term=query,
                        start_date=start_date_str,
                        end_date=end_date_str,
                        journal=journal if journal.strip() else None,
                        species=species if species != "All" else None
                    )
                    
                    if pmids:
                        articles = st.session_state.extractor.fetch_abstracts(pmids)
                        if articles:
                            # Convert to DataFrame and add origin tracking
                            df_query = pd.DataFrame(articles)
                            df_query = df_query[['PMID', 'Title', 'Abstract']].copy()
                            df_query['Query_Origin'] = label
                            df_query['Query_Number'] = original_idx + 1
                            
                            all_results.extend(articles)
                            combined_df = pd.concat([combined_df, df_query], ignore_index=True)
            
            if not combined_df.empty:
                # Create origin tracking for duplicates
                pmid_origins = combined_df.groupby('PMID')['Query_Origin'].apply(
                    lambda x: ' | '.join(sorted(set(x)))
                ).to_dict()
                
                pmid_numbers = combined_df.groupby('PMID')['Query_Number'].apply(
                    lambda x: ', '.join(map(str, sorted(set(x))))
                ).to_dict()
                
                # Create final combined DataFrame without duplicates
                final_df = combined_df.drop_duplicates(subset=['PMID']).copy()
                final_df['Query_Origin'] = final_df['PMID'].map(pmid_origins)
                final_df['Query_Numbers'] = final_df['PMID'].map(pmid_numbers)
                final_df = final_df.drop('Query_Number', axis=1)
                
                st.session_state.articles = final_df.to_dict('records')
                st.session_state.combined_df = final_df
                st.session_state.current_drug_name = drug_name  # Store drug name
                st.session_state.individual_query_results = {
                    query_labels[i]: combined_df[combined_df['Query_Origin'] == query_labels[i]].drop(['Query_Number'], axis=1)
                    for i in selected_questions if query_labels[i] in combined_df['Query_Origin'].values
                }
            else:
                st.warning("No articles found for any of the predefined queries")
                return
    
    # Display results
    if 'articles' in st.session_state and st.session_state.articles:
        articles = st.session_state.articles
        
        st.header(f"📊 Results ({len(articles)} articles)")
        
        # Show individual query results if predefined queries were used
        if 'individual_query_results' in st.session_state:
            st.subheader("📋 Results by Query")
            
            # Create tabs for each query
            query_labels = list(st.session_state.individual_query_results.keys())
            tabs = st.tabs([f"Q{i+1}: {label[:20]}..." if len(label) > 20 else f"Q{i+1}: {label}" 
                           for i, label in enumerate(query_labels)])
            
            for i, (label, tab) in enumerate(zip(query_labels, tabs)):
                with tab:
                    query_df = st.session_state.individual_query_results[label]
                    st.write(f"**Query:** {label}")
                    st.write(f"**Results:** {len(query_df)} articles")
                    
                    if not query_df.empty:
                        st.dataframe(query_df, use_container_width=True)
                        
                        # Download button for individual query
                        csv_data = query_df.to_csv(index=False)
                        st.download_button(
                            label=f"📄 Download Query {i+1} as CSV",
                            data=csv_data,
                            file_name=f"query_{i+1}_{label.replace('/', '_')}.csv",
                            mime="text/csv",
                            key=f"download_query_{i}"
                        )
                    else:
                        st.info("No results for this query")
            
            st.markdown("---")
            st.subheader("📊 Combined Results (Deduplicated)")
        
        # Summary statistics
        col1, col2, col3, col4 = st.columns(4)
        with col1:
            st.metric("Total Articles", len(articles))
        with col2:
            if 'combined_df' in st.session_state and 'Query_Origin' in st.session_state.combined_df.columns:
                unique_origins = st.session_state.combined_df['Query_Origin'].nunique()
                st.metric("Query Sources", unique_origins)
            else:
                st.metric("Query Sources", 1)
        with col3:
            abstracts_with_content = [article for article in articles if article.get('Abstract', 'N/A') != 'N/A']
            st.metric("With Abstracts", len(abstracts_with_content))
        with col4:
            if 'combined_df' in st.session_state:
                duplicates = len(st.session_state.combined_df) - len(articles)
                st.metric("Duplicates Removed", duplicates)
            else:
                st.metric("Duplicates Removed", 0)
        
        # Download options
        st.subheader("📥 Download Options")
        
        # Prepare CSV data
        df = pd.DataFrame(articles)
        csv_data = df.to_csv(index=False)
        
        col1, col2 = st.columns(2)
        with col1:
            st.download_button(
                label="📄 Download as CSV",
                data=csv_data,
                file_name=f"pubmed_results_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv",
                mime="text/csv"
            )
        
        with col2:
            # JSON download
            json_data = df.to_json(orient='records', indent=2)
            st.download_button(
                label="📋 Download as JSON",
                data=json_data,
                file_name=f"pubmed_results_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json",
                mime="application/json"
            )
        
        # Display articles in expandable format
        st.subheader("📚 Article Details")
        
        # Search within results
        search_filter = st.text_input("🔍 Search within results (title/abstract):")
        
        filtered_articles = articles
        if search_filter:
            filtered_articles = [
                article for article in articles
                if search_filter.lower() in article['Title'].lower() or 
                   search_filter.lower() in article['Abstract'].lower()
            ]
            st.info(f"Showing {len(filtered_articles)} articles matching '{search_filter}'")
        
        # Pagination
        articles_per_page = 10
        total_pages = (len(filtered_articles) + articles_per_page - 1) // articles_per_page
        
        if total_pages > 1:
            page = st.selectbox("Page", range(1, total_pages + 1))
            start_idx = (page - 1) * articles_per_page
            end_idx = start_idx + articles_per_page
            page_articles = filtered_articles[start_idx:end_idx]
        else:
            page_articles = filtered_articles
        
        # Display articles
        for i, article in enumerate(page_articles):
            # Handle both old and new article formats
            title = article.get('Title', 'N/A')
            pub_date = article.get('Publication_Date', 'N/A')
            pmid = article.get('PMID', 'N/A')
            abstract = article.get('Abstract', 'N/A')
            query_origin = article.get('Query_Origin', 'N/A')
            
            with st.expander(f"**{title}** ({pub_date})"):
                col1, col2 = st.columns([2, 1])
                
                with col1:
                    if 'Authors' in article:
                        st.write(f"**Authors:** {article['Authors']}")
                    if 'Journal' in article:
                        st.write(f"**Journal:** {article['Journal']}")
                    st.write(f"**PMID:** {pmid}")
                    if 'DOI' in article and article['DOI'] != 'N/A':
                        st.write(f"**DOI:** {article['DOI']}")
                    if 'Keywords' in article and article['Keywords'] != 'N/A':
                        st.write(f"**Keywords:** {article['Keywords']}")
                    
                    # Show query origin if available
                    if query_origin != 'N/A':
                        st.write(f"**Query Origin:** {query_origin}")
                        if 'Query_Numbers' in article:
                            st.write(f"**Query Numbers:** {article['Query_Numbers']}")
                
                with col2:
                    if 'PubMed_URL' in article:
                        st.link_button("🔗 View on PubMed", article['PubMed_URL'])
                    else:
                        st.link_button("🔗 View on PubMed", f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/")
                
                if abstract != 'N/A':
                    st.write("**Abstract:**")
                    st.write(abstract)
                else:
                    st.info("No abstract available")
    
    # LLM-based row-wise relevance judgment section - Individual Question Analysis
    st.markdown("---")
    st.header("🤖 LLM-based Individual Question Analysis")
    
    # Get available query results from session state
    available_queries = []
    if 'individual_query_results' in st.session_state:
        available_queries = list(st.session_state.individual_query_results.keys())
    
    if available_queries:
        llm_analysis_option = st.checkbox("Perform LLM relevance analysis on individual question results")
        
        if llm_analysis_option:
            # Get drug name from session state
            drug_name = st.session_state.get('current_drug_name', 'Unknown Drug')
            
            # Create selectbox for query selection
            selected_query = st.selectbox(
                "Select Query for LLM Analysis:",
                options=available_queries,
                help="Choose which individual query results you want to evaluate with LLM"
            )
            
            # Get the dataframe for selected query
            query_df = st.session_state.individual_query_results[selected_query].copy()
            
            st.write(f"**LLM analysis for:** {selected_query}")
            st.write(f"**Drug:** {drug_name}")
            st.write(f"**Articles to evaluate:** {len(query_df)}")
            
            if not query_df.empty:
                # Safety check for large datasets
                if len(query_df) > 20:
                    st.warning(f"This will process {len(query_df)} rows. This may take a long time. Consider filtering your results first.")
                    proceed = st.button("Proceed with LLM analysis")
                    if not proceed:
                        st.stop()
                
                # Define enhanced PSUR prompts with safety focus, sentiment analysis, and explicit consideration of author suspicion, recommendations, and negative findings
                psur_prompts = {
                    'Pregnancy Outcomes': '''You are evaluating this article for PSUR Section 11.1 - Literature articles on pregnancy outcomes (including termination) with/without adverse outcomes.

INCLUSION CRITERIA:
I. Articles with specific safety implications/safety issue with company suspect drug use in pregnancy:
- Articles describing teratogenic effects, fetal toxicity, or other maternal complications related to the use of the company suspect drug.
- Articles discussing mechanisms of harm or regulatory warnings related to pregnancy and the drug.

II. Articles with no specific safety implications/safety issue:
- Articles describing termination/abortion of pregnancy which assess maternal complications during termination.
- Articles describing clinical factors associated with subsequent surgical intervention in women undergoing medical termination of viable or non-viable pregnancies.

Article about {drug_name}:
Title: {title}
Abstract: {abstract}

INSTRUCTIONS: Mark as RELEVANT if the article fits any of the above criteria. Consider as relevant even if the main result is negative or inconclusive, if there is author suspicion, concern, or recommendation for further study.

Respond EXACTLY in this format:
Relevance: Yes/No
Summary: [Line 1: State which inclusion criteria the article fits and summarize the main finding.]
[Line 2: Explain the evidence level and author sentiment (positive/negative findings, suspicions, recommendations).]
[Line 3: State the safety relevance and implications for PSUR reporting and pharmacovigilance.]''',

                    'Paediatric Use': '''You are evaluating this article for PSUR Section 11.2 - Literature articles with use in paediatric populations.

INCLUSION CRITERIA:
- Articles including paediatric populations for the age range for which the company suspect product is approved and indicated.
- Articles focusing on safety findings of the drug in paediatric populations.

Article about {drug_name}:
Title: {title}
Abstract: {abstract}

INSTRUCTIONS: Mark as RELEVANT if the article fits any of the above criteria. Consider as relevant even if the main result is negative or inconclusive, if there is author suspicion, concern, or recommendation for further study.

Respond EXACTLY in this format:
Relevance: Yes/No
Summary: [Line 1: State if the article includes the approved paediatric age range and summarize safety findings.]
[Line 2: Explain the evidence level and author sentiment (positive/negative findings, suspicions, recommendations).]
[Line 3: State the safety relevance and implications for PSUR reporting and pharmacovigilance.]''',

                    'Compassionate Use/Supply': '''You are evaluating this article for PSUR Section 11.3 - Literature articles with compassionate supply, named patient use.

INCLUSION CRITERIA:
- Articles investigating the relationship between compassionate use and safety monitoring, such as identifying and managing safety concerns during drug use.
- Articles describing named-patient basis treatments (doctor obtains medicine directly from manufacturer before authorisation/approval, on individual basis, under doctor's responsibility).
- Note: Compassionate use is a treatment option for unauthorised medicines under strict conditions for patients with no satisfactory authorised therapies and who cannot enter clinical trials. Do NOT confuse with named-patient basis.

Article about {drug_name}:
Title: {title}
Abstract: {abstract}

INSTRUCTIONS: Mark as RELEVANT if the article fits any of the above criteria. Consider as relevant even if the main result is negative or inconclusive, if there is author suspicion, concern, or recommendation for further study.

Respond EXACTLY in this format:
Relevance: Yes/No
Summary: [Line 1: State if the article is about compassionate use or named-patient basis and summarize safety findings.]
[Line 2: Explain the evidence level and author sentiment (positive/negative findings, suspicions, recommendations).]
[Line 3: State the safety relevance and implications for PSUR reporting and pharmacovigilance.]''',

                    'efficacy': '''You are evaluating this article for PSUR Section 11.4 - Literature articles on lack of efficacy.

INCLUSION CRITERIA:
- Articles clearly defining treatment failure, lack of effect, inefficacy, or partial/incomplete effect of the drug.
- Articles that rule out non-compliance or incorrect dosing as a cause of lack of efficacy.

Article about {drug_name}:
Title: {title}
Abstract: {abstract}

INSTRUCTIONS: Mark as RELEVANT if the article fits any of the above criteria. Consider as relevant even if the main result is negative or inconclusive, if there is author suspicion, concern, or recommendation for further study.

Respond EXACTLY in this format:
Relevance: Yes/No
Summary: [Line 1: State if the article defines treatment failure/lack of efficacy and if non-compliance/incorrect dosing is ruled out.]
[Line 2: Explain the evidence level and author sentiment (positive/negative findings, suspicions, recommendations).]
[Line 3: State the safety relevance and implications for PSUR reporting and pharmacovigilance.]''',

                    'Overdose/Abuse/Misuse': '''You are evaluating this article for PSUR Section 11.5 - Literature articles on asymptomatic overdose, abuse or misuse.

INCLUSION CRITERIA:
- Articles describing instances of overdose, abuse, or misuse of the company suspect drug, with details on dose, duration, and patient outcome.
- Articles specifying whether overdose was intentional or accidental.
- Articles should be considered even if no symptoms occurred.

Article about {drug_name}:
Title: {title}
Abstract: {abstract}

INSTRUCTIONS: Mark as RELEVANT if the article fits any of the above criteria. Consider as relevant even if the main result is negative or inconclusive, if there is author suspicion, concern, or recommendation for further study.

Respond EXACTLY in this format:
Relevance: Yes/No
Summary: [Line 1: State if the article describes overdose/abuse/misuse, and if intentional/accidental, and if symptoms occurred.]
[Line 2: Explain the evidence level and author sentiment (positive/negative findings, suspicions, recommendations).]
[Line 3: State the safety relevance and implications for PSUR reporting and pharmacovigilance.]''',

                    'Medication Error': '''You are evaluating this article for PSUR Section 11.6 - Literature articles with medication error where no adverse events occurred.

INCLUSION CRITERIA:
- Articles of medication error with company suspect product, including type of error (prescribing, dispensing, administration).
- Articles describing contributing factors responsible for the medication error.

Article about {drug_name}:
Title: {title}
Abstract: {abstract}

INSTRUCTIONS: Mark as RELEVANT if the article fits any of the above criteria. Consider as relevant even if the main result is negative or inconclusive, if there is author suspicion, concern, or recommendation for further study.

Respond EXACTLY in this format:
Relevance: Yes/No
Summary: [Line 1: State the type of medication error and any contributing factors described.]
[Line 2: Explain the evidence level and author sentiment (positive/negative findings, suspicions, recommendations).]
[Line 3: State the safety relevance and implications for PSUR reporting and pharmacovigilance.]''',

                    'In Vitro Studies': '''You are evaluating this article for PSUR Section 11.7 - Literature articles with important non-clinical safety results.

INCLUSION CRITERIA:
- Animal studies showing safety issues and toxicity related to company suspect product.
- In vitro studies revealing potential safety concerns (e.g., genotoxicity, carcinogenicity; includes cell cultures, tissues, or other in-vitro models).
- Articles reflecting mechanistic insights relevant to human safety (how the product works/interacts with biological systems to produce potential safety concerns/adverse events).
- Articles focusing on findings with potential clinical relevance.

Article about {drug_name}:
Title: {title}
Abstract: {abstract}

INSTRUCTIONS: Mark as RELEVANT if the article fits any of the above criteria. Consider as relevant even if the main result is negative or inconclusive, if there is author suspicion, concern, or recommendation for further study.

Respond EXACTLY in this format:
Relevance: Yes/No
Summary: [Line 1: State if the article is animal/in vitro/mechanistic/clinically relevant and summarize the safety finding.]
[Line 2: Explain the evidence level and author sentiment (positive/negative findings, suspicions, recommendations).]
[Line 3: State the safety relevance and implications for PSUR reporting and pharmacovigilance.]''',

                    'others': '''You are evaluating this article for PSUR Section 11.8 - Other relevant literature articles.

INCLUSION CRITERIA:
- Articles discussing safety concerns or risks associated with the company suspect product that do not fit into the above sections, including but not limited to:
  - Drug interactions
  - Real-world evidence studies
  - Pharmacogenomics (relationship between genetic variants and drug safety/risk)
  - Regulatory updates or safety communications

Article about {drug_name}:
Title: {title}
Abstract: {abstract}

INSTRUCTIONS: Mark as RELEVANT if the article fits any of the above criteria. Consider as relevant even if the main result is negative or inconclusive, if there is author suspicion, concern, or recommendation for further study.

Respond EXACTLY in this format:
Relevance: Yes/No
Summary: [Line 1: State which inclusion criteria the article fits and summarize the main finding.]
[Line 2: Explain the evidence level and author sentiment (positive/negative findings, suspicions, recommendations).]
[Line 3: State the safety relevance and implications for PSUR reporting and pharmacovigilance.]'''
                }
                
                # Get appropriate prompt for selected query
                selected_prompt = psur_prompts.get(selected_query, psur_prompts['others'])
                
                # Initialize LLM
                try:
                    rowwise_llm = ChatOllama(
                        model="gemma3:4b",
                        temperature=0,
                        num_ctx=4096
                    )
                    llm_available = True
                except Exception as e:
                    st.warning(f"Ollama not available: {e}. Using mock responses for demonstration.")
                    llm_available = False
                
                # Process each article
                relevance_results = []
                summary_results = []
                progress_bar = st.progress(0)
                total_rows = len(query_df)
                
                for i, (idx, row) in enumerate(query_df.iterrows()):
                    progress_bar.progress((i + 1) / total_rows)
                    
                    # Format prompt
                    prompt = selected_prompt.format(
                        title=row['Title'],
                        abstract=row['Abstract'],
                        drug_name=drug_name
                    )
                    
                    with st.spinner(f"LLM evaluating article {i+1} of {total_rows}..."):
                        if llm_available:
                            try:
                                response = rowwise_llm.invoke(prompt)
                                if hasattr(response, 'content'):
                                    response = response.content
                                elif not isinstance(response, str):
                                    response = str(response)
                            except Exception as e:
                                st.error(f"LLM evaluation failed for row {i+1}: {str(e)}")
                                response = "Relevance: No\nSummary: Error in evaluation\nUnable to process\nPlease review manually"
                        else:
                            # Mock response for demonstration with safety focus
                            response = f"Relevance: Yes\nSummary: This article reports new safety findings related to {drug_name} in the context of {selected_query}.\nThe authors express concerns about potential safety implications and recommend further monitoring and investigation.\nThese findings are highly relevant for PSUR reporting and require pharmacovigilance attention for risk assessment."
                    
                    # Parse response to extract relevance and summary
                    lines = response.strip().split('\n')
                    relevance = "No"
                    summary = "Unable to parse summary"
                    
                    for line in lines:
                        if line.startswith('Relevance:'):
                            relevance = line.replace('Relevance:', '').strip()
                        elif line.startswith('Summary:'):
                            summary_start = lines.index(line)
                            summary_lines = lines[summary_start:summary_start+4]  # Get Summary line + 3 lines
                            summary = '\n'.join([l.replace('Summary:', '').strip() if l.startswith('Summary:') else l.strip() for l in summary_lines if l.strip()])
                            break
                    
                    relevance_results.append(relevance)
                    summary_results.append(summary)
                
                # Add results to dataframe
                query_df['Relevance'] = relevance_results
                query_df['Summary'] = summary_results
                
                # Create final output dataframe with required columns
                output_df = query_df[['PMID', 'Title', 'Abstract', 'Relevance', 'Summary']].copy()
                
                # Display results
                st.subheader("LLM Analysis Results")
                st.dataframe(output_df, use_container_width=True)
                
                # Download option with proper filename
                csv_data = output_df.to_csv(index=False)
                st.download_button(
                    label=f"📄 Download LLM Analysis Results",
                    data=csv_data,
                    file_name=f"llm_analysis_{selected_query.replace('/', '_')}_{drug_name.replace(' ', '_')}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv",
                    mime="text/csv"
                )
                
                # Summary statistics
                st.subheader("📊 Analysis Summary")
                col1, col2 = st.columns(2)
                
                relevant_count = len([r for r in relevance_results if r.lower() == 'yes'])
                not_relevant_count = len([r for r in relevance_results if r.lower() == 'no'])
                
                with col1:
                    st.metric("Relevant Articles", relevant_count)
                
                with col2:
                    st.metric("Not Relevant Articles", not_relevant_count)
                
                # Show sample of relevant articles
                if relevant_count > 0:
                    st.subheader("📋 Sample Relevant Articles")
                    relevant_articles = output_df[output_df['Relevance'].str.lower() == 'yes'].head(3)
                    for idx, article in relevant_articles.iterrows():
                        with st.expander(f"PMID: {article['PMID']} - {article['Title'][:100]}..."):
                            st.write(f"**Relevance:** {article['Relevance']}")
                            st.write(f"**Summary:**")
                            st.write(article['Summary'])
            else:
                st.info("No articles found for the selected query.")
    else:
        st.info("No individual query results available. Please perform searches first.")

if __name__ == "__main__":
    main()
