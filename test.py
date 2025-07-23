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
        """Search PubM
        
