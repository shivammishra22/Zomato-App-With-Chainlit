import streamlit as st
import pandas as pd
from langchain_ollama import ChatOllama
from langchain_core.prompts import ChatPromptTemplate
from langchain_huggingface import HuggingFaceEmbeddings
from langchain.document_loaders import DataFrameLoader
from langchain.vectorstores import FAISS
import time

# --- Streamlit UI ---
st.title("🧠 Medical Abstract Summarizer & Semantic Search")
uploaded_file = st.file_uploader("📄 Upload CSV file with abstracts", type=["csv"])
query = st.text_input("🔍 Enter your query (e.g., treatment of schizophrenia using olanzapine):")
run_button = st.button("Run")

# --- Define Model and Prompt Once ---
@st.cache_resource
def load_model():
    return ChatOllama(model="qwen3:0.6b", reasoning=False, temperature=0.1, max_tokens=500, num_ctx=2048)

@st.cache_resource
def load_prompt():
    return ChatPromptTemplate.from_messages([
        ("system", "You are a helpful assistant that summarizes medical abstracts. Just provide a concise summary that captures the main points of the abstract. Don't include any additional information or context."),
        ("human", "Summarize the following abstract: {abstract}")
    ])

@st.cache_resource
def get_embedding_model():
    return HuggingFaceEmbeddings(model_name='NeuML/pubmedbert-base-embeddings')

# --- Main Logic ---
if run_button and uploaded_file and query:
    try:
        df = pd.read_csv(uploaded_file)
        if 'Abstract' not in df.columns:
            st.error("❌ 'Abstract' column not found in the uploaded CSV.")
            st.stop()

        df.dropna(subset=['Abstract'], inplace=True)
        df['Abstract_Summary_Qwen'] = ''

        model = load_model()
        prompt = load_prompt()

        st.info("⏳ Summarizing abstracts...")
        progress_bar = st.progress(0)
        for i, (idx, row) in enumerate(df.iterrows()):
            abstract = row['Abstract']
            if pd.notna(abstract) and abstract.strip() and abstract != 'N/A':
                try:
                    summary = model.invoke(prompt.invoke({"abstract": abstract}))
                    df.at[idx, 'Abstract_Summary_Qwen'] = summary.content
                except Exception as e:
                    df.at[idx, 'Abstract_Summary_Qwen'] = f"Error: {str(e)[:100]}"
            else:
                df.at[idx, 'Abstract_Summary_Qwen'] = "No abstract available"
            progress_bar.progress((i + 1) / len(df))

        st.success("✅ Summarization completed.")
        df_display = df[["PMID", "Title", "Abstract_Summary_Qwen"]]
        st.dataframe(df_display.head(10))

        st.info("🔎 Performing semantic search...")

        # Prepare documents and FAISS index
        df_for_search = df[["PMID", "Title", "Abstract_Summary_Qwen"]].copy()
        loader = DataFrameLoader(df_for_search, page_content_column="Abstract_Summary_Qwen")
        docs = loader.load()

        embedding_model = get_embedding_model()
        faiss_index = FAISS.from_documents(docs, embedding_model)

        # Perform similarity search
        results = faiss_index.similarity_search(query, k=10)

        st.subheader("📌 Top 10 Most Relevant Results")
        for i, doc in enumerate(results, 1):
            metadata = doc.metadata
            st.markdown(f"**{i}. PMID:** {metadata.get('PMID', 'N/A')}")
            st.markdown(f"**Title:** {metadata.get('Title', 'N/A')}")
            st.markdown(f"**Summary:** {doc.page_content}")
            st.markdown("---")

    except Exception as e:
        st.error(f"❌ Error: {str(e)}")

elif run_button:
    st.warning("⚠️ Please upload a file and enter a query.")
