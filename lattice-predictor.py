# KB 2023
# hH - lattice parameter predictor v1
# PhD stuff

import streamlit as st
st.set_page_config(page_title='half-Heulser lattice predictor | Kaja Bilińska', page_icon='.png', layout='centered', initial_sidebar_state='auto')

# from sentence_transformers import SentenceTransformer
# from sklearn.metrics.pairwise import cosine_similarity # in requirements given as scikit-learn
# from transformers import AutoTokenizer, AutoModel
# from urllib.request import urlopen # no need to install it - standard library
# from bs4 import BeautifulSoup
import numpy as np
import pandas as pd

st.title('Lattice parameter predictor')
# st.subheader('In order to improve Your website...')

st.write('Hello there...')