import streamlit as st
import csv
import time
import io
import pandas as pd
from datetime import datetime, date
from typing import List, Dict, Optional
import sys

try:
    from Bio import Entrez
except ImportError:
    st.error("Error: Biopython not installed. Install with: pip install biopython")
    st.stop()

class PubMedExtractor:
    def __init__(self, email: str, api_key: Optional[str] = None):
        """Initialize PubMed extractor"""
        Entrez.email = email
        if api_key:
            Entrez.api_key = api_key
        
    def search_pubmed(self, 
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
                     mesh_terms: Optional[str] = None) -> List[str]:
        """Search PubMed and return list of PMIDs matching criteria"""
        
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
        if species and species != "All":
            if species.lower() == 'humans':
                query_parts.append('"humans"[MeSH Terms]')
            elif species.lower() == 'animals':
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
            terms = [term.strip() for term in mesh_terms.split(',')]
            for term in terms:
                query_parts.append(f'"{term}"[MeSH Terms]')
        
        final_query = ' AND '.join(query_parts)
        
        try:
            # Get total count
            handle = Entrez.esearch(db="pubmed", term=final_query, retmax=0)
            search_results = Entrez.read(handle)
            handle.close()
            
            total_count = int(search_results["Count"])
            
            if total_count == 0:
                return [], final_query, 0
            
            max_retrievable = min(total_count, 10000)
            
            # Get PMIDs
            handle = Entrez.esearch(db="pubmed", 
                                  term=final_query, 
                                  retmax=max_retrievable,
                                  sort="relevance")
            search_results = Entrez.read(handle)
            handle.close()
            
            pmids = search_results["IdList"]
            return pmids, final_query, total_count
            
        except Exception as e:
            st.error(f"Error searching PubMed: {e}")
            return [], "", 0
    
    def fetch_abstracts(self, pmids: List[str]) -> List[Dict]:
        """Fetch detailed information including abstracts for given PMIDs"""
        if not pmids:
            return []
        
        articles = []
        batch_size = 500
        
        progress_bar = st.progress(0)
        status_text = st.empty()
        
        for i in range(0, len(pmids), batch_size):
            batch_pmids = pmids[i:i+batch_size]
            progress = (i + batch_size) / len(pmids)
            progress_bar.progress(min(progress, 1.0))
            status_text.text(f"Fetching articles {i+1} to {min(i+batch_size, len(pmids))} of {len(pmids)}...")
            
            try:
                handle = Entrez.efetch(db="pubmed", 
                                     id=batch_pmids, 
                                     rettype="medline", 
                                     retmode="xml")
                records = Entrez.read(handle)
                handle.close()
                
                for record in records['PubmedArticle']:
                    article_info = self._parse_article(record)
                    if article_info:
                        articles.append(article_info)
                
                time.sleep(0.1)
                
            except Exception as e:
                st.warning(f"Error fetching batch {i//batch_size + 1}: {e}")
                continue
        
        progress_bar.progress(1.0)
        status_text.text(f"Successfully extracted information for {len(articles)} articles")
        
        return articles
    
    def _parse_article(self, record) -> Optional[Dict]:
        """Parse individual PubMed article record"""
        try:
            article = record['MedlineCitation']['Article']
            pmid = record['MedlineCitation']['PMID']
            
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
            
            # Extract DOI
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
            return None

def main():
    st.set_page_config(
        page_title="PubMed Search Tool",
        page_icon="🔬",
        layout="wide"
    )
    
    st.title("🔬 PubMed Abstract Extractor")
    st.markdown("Search PubMed and extract abstracts with comprehensive filtering options")
    
    # Sidebar for configuration
    with st.sidebar:
        st.header("⚙️ Configuration")
        email = st.text_input(
            "Email Address*", 
            help="Required by NCBI for API access",
            key="email_input"
        )
        api_key = st.text_input(
            "NCBI API Key", 
            type="password", 
            help="Optional: For higher rate limits",
            key="api_key_input"
        )
        
        if not email:
            st.warning("Please enter your email address to proceed")
            st.stop()
    
    # Main search interface
    col1, col2 = st.columns([2, 1])
    
    with col1:
        st.header("🔍 Search Query")
        search_term = st.text_area(
            "Search Terms*",
            placeholder="e.g., (diabetes) AND (treatment) NOT (case report)",
            help="Use PubMed search syntax. AND, OR, NOT operators supported",
            key="search_term_input"
        )
    
    with col2:
        st.header("📅 Date Range")
        col_start, col_end = st.columns(2)
        with col_start:
            start_date = st.date_input("Start Date", value=None, key="start_date_input")
        with col_end:
            end_date = st.date_input("End Date", value=None, key="end_date_input")
    
    # Advanced filters
    with st.expander("🔧 Advanced Filters", expanded=False):
        col1, col2, col3 = st.columns(3)
        
        with col1:
            journal = st.text_input(
                "Journal Name", 
                placeholder="e.g., Nature, Science",
                key="journal_input"
            )
            author = st.text_input(
                "Author Name", 
                placeholder="e.g., Smith J",
                key="author_input"
            )
            mesh_terms = st.text_input(
                "MeSH Terms", 
                placeholder="e.g., Diabetes Mellitus, Hypertension",
                key="mesh_input"
            )
        
        with col2:
            species = st.selectbox(
                "Species",
                ["All", "Humans", "Animals"],
                key="species_input"
            )
            
            publication_type = st.selectbox(
                "Publication Type",
                ["All", "Clinical Trial", "Randomized Controlled Trial", "Review", 
                 "Meta-Analysis", "Case Reports", "Observational Study", "Letter", "Editorial"],
                key="pub_type_input"
            )
            
            language = st.selectbox(
                "Language",
                ["All", "English", "Spanish", "French", "German", "Italian", "Japanese"],
                key="language_input"
            )
        
        with col3:
            age_group = st.selectbox(
                "Age Group",
                ["All", "Infant", "Child", "Adolescent", "Adult", "Middle Aged", "Aged"],
                key="age_group_input"
            )
            
            free_full_text = st.checkbox("Free Full Text Only", key="free_text_input")
            
            # Max results options
            st.subheader("Result Limits")
            return_all = st.checkbox("Return All Available Articles", 
                                   help="Get all articles found (up to PubMed's 10,000 limit)",
                                   key="return_all_input")
            
            if not return_all:
                max_results = st.number_input("Max Results", min_value=1, max_value=10000, value=1000, key="max_results_input")
            else:
                max_results = 10000  # PubMed's maximum limit
                st.info("Will return all available articles (up to 10,000)")
    
    # Show search preview when terms are entered
    if search_term:
        st.subheader("🔍 Search Preview")
        
        # Build preview query (similar to actual search)
        preview_parts = [search_term]
        
        if start_date and end_date:
            start_str = start_date.strftime("%Y/%m/%d")
            end_str = end_date.strftime("%Y/%m/%d") 
            preview_parts.append(f'("{start_str}"[Date - Publication] : "{end_str}"[Date - Publication])')
        elif start_date:
            start_str = start_date.strftime("%Y/%m/%d")
            preview_parts.append(f'"{start_str}"[Date - Publication] : 3000[Date - Publication]')
        elif end_date:
            end_str = end_date.strftime("%Y/%m/%d")
            preview_parts.append(f'1900[Date - Publication] : "{end_str}"[Date - Publication]')
        
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
            terms = [term.strip() for term in mesh_terms.split(',')]
            for term in terms:
                preview_parts.append(f'"{term}"[MeSH Terms]')
        
        preview_query = ' AND '.join(preview_parts)
        st.code(preview_query, language='text')
        
        # Show result limit info
        if return_all:
            st.info("📊 Will retrieve all available articles (up to 10,000)")
        else:
            st.info(f"📊 Will retrieve up to {max_results} articles")
    
    # Search button
    search_disabled = not search_term or not email
    if search_disabled and not email:
        st.warning("⚠️ Please enter your email address in the sidebar to enable searching")
    elif search_disabled and not search_term:
        st.info("💡 Enter search terms above to begin")
    
    if st.button("🔍 Search PubMed", type="primary", disabled=search_disabled):
        if not email:
            st.error("Please enter your email address in the sidebar")
            st.stop()
        
        with st.spinner("Searching PubMed..."):
            extractor = PubMedExtractor(email, api_key)
            
            # Convert dates to string format
            start_date_str = start_date.strftime("%Y/%m/%d") if start_date else None
            end_date_str = end_date.strftime("%Y/%m/%d") if end_date else None
            
            # Search PubMed
            pmids, query, total_count = extractor.search_pubmed(
                search_term=search_term,
                start_date=start_date_str,
                end_date=end_date_str,
                journal=journal if journal else None,
                species=species if species != "All" else None,
                publication_type=publication_type if publication_type != "All" else None,
                language=language if language != "All" else None,
                age_group=age_group if age_group != "All" else None,
                free_full_text=free_full_text,
                author=author if author else None,
                mesh_terms=mesh_terms if mesh_terms else None
            )
            
            if not pmids:
                st.warning("No articles found matching your criteria")
                st.info(f"**Search Query Used:** {query}")
                st.stop()
            
            # Display search results summary
            st.success(f"Found {total_count} articles (retrieving up to {min(len(pmids), max_results)})")
            st.info(f"**Search Query Used:** {query}")
            
            # Handle result limiting
            original_count = len(pmids)
            if return_all:
                if total_count > 10000:
                    st.warning(f"⚠️ Found {total_count} articles, but PubMed limits results to 10,000. Retrieving all available {len(pmids)} articles.")
                else:
                    st.info(f"✅ Retrieving all {len(pmids)} articles found.")
                # Use all available PMIDs
            else:
                if len(pmids) > max_results:
                    pmids = pmids[:max_results]
                    st.info(f"📊 Limited to first {max_results} results as requested (out of {original_count} available)")
                else:
                    st.info(f"✅ Retrieving all {len(pmids)} articles (less than your {max_results} limit)")
            
            # Fetch abstracts
            with st.spinner("Fetching article details..."):
                articles = extractor.fetch_abstracts(pmids)
            
            if articles:
                # Display results
                st.header("📊 Results")
                
                # Convert to DataFrame for display
                df = pd.DataFrame(articles)
                
                # Display summary stats
                col1, col2, col3, col4 = st.columns(4)
                with col1:
                    st.metric("Total Articles", len(articles))
                with col2:
                    unique_journals = df['Journal'].nunique()
                    st.metric("Unique Journals", unique_journals)
                with col3:
                    articles_with_abstracts = len(df[df['Abstract'] != 'N/A'])
                    st.metric("With Abstracts", articles_with_abstracts)
                with col4:
                    articles_with_doi = len(df[df['DOI'] != 'N/A'])
                    st.metric("With DOI", articles_with_doi)
                
                # Display articles table
                st.subheader("📋 Article Details")
                
                # Make the table more readable
                display_df = df.copy()
                display_df['Title'] = display_df['Title'].str[:100] + "..."
                display_df['Abstract'] = display_df['Abstract'].str[:150] + "..."
                
                st.dataframe(
                    display_df[['PMID', 'Title', 'Authors', 'Journal', 'Publication_Date']],
                    use_container_width=True
                )
                
                # Download button
                st.subheader("💾 Download Results")
                
                # Create CSV
                csv_buffer = io.StringIO()
                df.to_csv(csv_buffer, index=False)
                csv_string = csv_buffer.getvalue()
                
                # Generate filename
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
                
                # Show sample of full data
                with st.expander("👀 Preview Full Data", expanded=False):
                    st.dataframe(df, use_container_width=True)
            
            else:
                st.error("Failed to retrieve article details")
    
    # Footer
    st.markdown("---")
    st.markdown("""
    **Tips for better searches:**
    - 🔄 **All inputs update automatically** - no need to press Enter
    - Use Boolean operators: AND, OR, NOT
    - Use quotes for exact phrases: "machine learning"
    - Use parentheses to group terms: (diabetes OR "diabetes mellitus") AND treatment
    - Use MeSH terms for more precise results
    - Combine multiple filters for refined searches
    - Check "Return All Available Articles" to get comprehensive results (up to 10,000)
    - For large datasets, be patient as processing may take several minutes
    """)

if __name__ == "__main__":
    main()
