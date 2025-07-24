import pandas as pd
from langchain_core.prompts import ChatPromptTemplate
from langchain_community.chat_models import ChatOllama



# Ask the user for the indication
indication = input("Enter the therapeutic indication to evaluate against the abstract: ")

# Define the prompt template
template = ChatPromptTemplate.from_messages([
    ("system", 
     """You are a pharmacovigilance expert generating the PSUR sub-section titled “Characterisation of Benefits”.

Evaluate the abstract based on the following criteria:
1) Strength of evidence: comparator(s), effect size, statistical significance, consistency, methodology.
2) Clinical relevance of effect size.
3) Generalisability to full indicated population (e.g., any sub-population not showing benefit).
4) Dose-response characterisation.
5) Duration of effect.
6) Comparative efficacy.

Also:
- Verify whether the abstract supports the following therapeutic indication: "{Indication}"
- Identify any other indications mentioned in the abstract that differ from the one provided.

Respond in the following format:
Relevance: Relevant or Not Relevant
Indication Match: Yes or No
Other Indications: [List of other indications or 'None']
"""),
    ("user", "{Abstract}")
])

# Initialize the language model
llm = ChatOllama(model='gemma3:4b', temperature=0.1, num_ctx=1000)

# Create the pipeline
chain = template | llm

# Process each abstract
relevance_results = []
indication_results = []
other_indications_results = []

for abstract in df['Abstract']:
    result = chain.invoke({
        'Abstract': abstract,
        'Indication': indication
    }).content.strip()

    # Default values
    relevance = "Not Available"
    indication_match = "Not Available"
    other_indications = "NaN"

    # Parse the result
    for line in result.splitlines():
        if "Relevance:" in line:
            relevance = line.split(":", 1)[1].strip()
            if relevance.lower() == "irrelevant":
                relevance = "Not Relevant"
        elif "Indication Match:" in line:
            indication_match = line.split(":", 1)[1].strip()
        elif "Other Indications:" in line:
            other_indications = line.split(":", 1)[1].strip()
            if other_indications.lower() in ["none", "n/a", "not applicable"]:
                other_indications = "NaN"

    relevance_results.append(relevance)
    indication_results.append(indication_match)
    other_indications_results.append(other_indications)

# Add results to DataFrame
df['Relevance'] = relevance_results
df['Indication Match'] = indication_results
df['Other Indications'] = other_indications_results
df.to_csv("hello.csv")
# Output the result
print(df)
