from langchain_ollama import ChatOllama
from langchain_huggingface import HuggingFaceEmbeddings
from langchain_core.prompts import ChatPromptTemplate
from langchain.document_loaders import DataFrameLoader
import pandas as pd
import numpy as np
from langchain.vectorstores import FAISS
import time
from tqdm import tqdm
model=ChatOllama(model="qwen3:0.6b",temperature=0.1,reasoning=False,max_token=500,num_ctx=2058)
prompt = ChatPromptTemplate.from_messages([
    ("system", "You are a helpful assistant that summarizes medical abstracts.Just provide a concise summary that captures the main points of the abstract and procide summary dont include any additional information or context."),
    ("human", "Summarize the following abstract: {abstract}")
])
 
# Load Qwen model
model = ChatOllama(
    model="qwen3:0.6b",
    reasoning=False,
    temperature=0.1,
    max_tokens=500,
    num_ctx=2048
)
 
# Improved prompt for medical abstract summarization
prompt = ChatPromptTemplate.from_messages([
    ("system", "You are a helpful assistant that summarizes medical abstracts.Just provide a concise summary that captures the main points of the abstract and procide summary dont include any additional information or context."),
    ("human", "Summarize the following abstract: {abstract}")
])
 
# Load CSV
input_csv = r"C:\Users\shivam.mishra2\Downloads\Literature\Esomeprazole\esomeprazole.csv"
df = pd.read_csv(input_csv)
df.dropna(subset=['Abstract'], inplace=True)  # Drop rows with NaN abstracts
df['Abstract_Summary_Qwen'] = ''
 
# Summarize abstracts with tqdm progress bar
start = time.time()
for idx, row in tqdm(df.iterrows(), total=len(df), desc="Summarizing abstracts"):
    abstract = row.get('Abstract', '')
    if pd.notna(abstract) and abstract.strip() != '' and abstract != 'N/A':
        try:
            summary = model.invoke(prompt.invoke({"abstract": abstract}))
            df.at[idx, 'Abstract_Summary_Qwen'] = summary.content
        except Exception as e:
            df.at[idx, 'Abstract_Summary_Qwen'] = f"Error: {str(e)[:100]}"
    else:
        df.at[idx, 'Abstract_Summary_Qwen'] = 'No abstract available for summarization'
end = time.time()

 
embedding_model=HuggingFaceEmbeddings(model_name='NeuML/pubmedbert-base-embeddings')
model=ChatOllama(model="qwen3:0.6b",temperature=0.1,num_ctx=2000,reasoning=False,max_token=500)

df=df[["PMID","Abstract_Summary_Qwen","Title"]]
Dataloader=DataFrameLoader(df,page_content_column='Abstract_Summary_Qwen')
docs=Dataloader.load()
faiss_index = FAISS.from_documents(docs, embedding_model)
query='Use of esomeprazole for the not for children prolonged treatment  after i.v. induced prevention of re-bleeding of peptic ulcers'
results = faiss_index.similarity_search(query, k=10) 
results
