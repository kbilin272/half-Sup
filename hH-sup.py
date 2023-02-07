# KB 2023
# hH - lattice parameter predictor v1
# PhD stuff

import streamlit as st
st.set_page_config(page_title='half-Sup | half-Heusler Support', page_icon='hS.jpeg', layout='centered', initial_sidebar_state='auto')

# from sentence_transformers import SentenceTransformer
# from sklearn.metrics.pairwise import cosine_similarity # in requirements given as scikit-learn
# from transformers import AutoTokenizer, AutoModel
# from urllib.request import urlopen # no need to install it - standard library
# from bs4 import BeautifulSoup
import numpy as np
import pandas as pd
import sklearn

def _read_vasprun():
    """
    Args:
        vasprun.xml / path
        
    Returns:
        input for calculate_kappa()
    """
    
    st.write('')
    st.write('Not ready yet...')
    st.write('')

def calculate_kappa(C11, C12, C44, V, n):
    """
    it shall be in the separated file, but will it streamlite handle it?
    """
    
    # csv with mass of the particular elements is required + some conditions to stop calculating
    
    # st.write('')
    
    # some in-between calculations for bulk and Shear
    B_V = ( C11 + 2*C12 ) / 3
    B_R = ( ( C11 + C12 ) * C11 - 2 * C12**2 ) / ( 3 * ( C11 - C12 ) )

    G_V = ( ( C11 - C12 ) + 3*C44 ) / 5
    G_R_1 = ( 18 * B_V ) / ( (C11 + C12) * C11 - 2 * C12**2 )
    G_R_2 = 6 / ( C11 - C12 )
    G_R_3 = 9 / C44
    G_R = 15 *  ( G_R_1 + G_R_2 + G_R_3 )**(-1)

    # bulk modulus
    B = np.round( ( B_V + B_R ) / 2, 2)

    # Shear modulus
    G = np.round( ( G_V + G_R ) / 2, 2)

    # Shear anisotropy factor
    A = np.round( ( 2 * C44 ) / ( C11 - C12 ), 2)

    # Poisson's ratio
    v = np.round( ( 3 * B - 2 * G ) / ( 6 * B + 2 * G ), 2)

    # Young's modulus
    Y = np.round( ( 9 * B * G ) / ( 3 * B + G ), 2)

    st.write(f'Bulk modulus, B = {B} GPa')
    st.write(f'Shear modulus, G = {G} ???')
    st.write(f'Shear anisotropy factor, A = {A} ???')
    st.write(f"Poisson's ratio, v = {v} ???")
    st.write(f"Young's modulus, Y = {Y} ???")
    
    # outpt data given as a latex table?
    
    st.balloons()

def main():
    """NLP App with Streamlit and TextBlob"""

    # st.title("half-Support")
    
    def _not_supported():

        title_templ = """
        <h1 style="color:#737373">Caution!</h1>
        </div>
        """

        st.markdown(title_templ,unsafe_allow_html=True)

        subheader_templ = """
        <div style="background-color:#737373;padding:8px;">
        <h3 style="color:#202020">Not supported yet...
        </h3>
        </div>
        """
        

        st.markdown(subheader_templ,unsafe_allow_html=True)
        st.write('')
            
    st.sidebar.image("half-Sup.jpeg", use_column_width=True)
    
    activity = ['About half-Support', "Lattice parameter prediction", "Slack formula", "Credits"]
    choice = st.sidebar.selectbox('',activity) # 'menu'

	# Text Analysis CHOICE
    if choice == 'About half-Support':
    
        _not_supported()

        st.write("")
        st.subheader("About half-Support")        
        st.write("")
        st.write("")

        # raw_text = st.text_area("Write something","Enter a Text in English...",height=250)

        # if st.button("Analyze"):
            # if len(raw_text) == 0:
            	# st.warning("Enter a Text...")
            # else:
            	# blob = TextBlob(raw_text)
            	# st.write("")
                
    elif choice == 'Lattice parameter prediction':
    
        _not_supported
    
    elif choice == "Slack formula":
    
        st.write("")
        st.subheader("Lattice thermal calculations")        
        st.write("")
        # st.write("")
        
        # st.write('In order to calculate lattice thermal conductivity - provide sufficient data.\n')
        # st.write('Decide on manual / supported version of file input. (on the left)')
        # st.write('')
        
        Slack_options = ['Manual input', "VASP input"]
        choice = st.sidebar.selectbox("Choose the input data manner:",Slack_options)
        
        if choice == Slack_options[0]:
        
            st.write("")
            st.write("Manual input required data below (mind units given):")
            st.write("")            
            
            coli1, coli2, coli3 = st.columns(3)
        
            with coli1:
                ion1 = st.text_input("ion 1")
                
            with coli2:
                ion2 = st.text_input("ion 2")
                
            with coli3:
                ion3 = st.text_input("ion 3")
            
            col4, col5 = st.columns(2)
        
            col1, col2, col3 = st.columns(3)
        
            with col1:
                C11 = st.text_input("C11 (GPa)")
                
            with col2:
                C12 = st.text_input("C12 (GPa)")
                
            with col3:
                C44 = st.text_input("C44 (GPa)")
            
            col4, col5 = st.columns(2)

            with col4:
                V = st.text_input('Volume of the cell ($\AA$³)')
            
            with col5:
                n = st.text_input('Number of atoms in the cell')
                
            def is_type(num, type='float'):
                try:
                    if type == 'float':
                        float(num)
                        return True
                    elif type == 'int':
                        int(num)
                        return True
                except ValueError:
                    return False
          

            if ion1 != '' and ion2 != '' and ion3 != '' and C11 != '' and C12 != '' and C44 != '' and V != '' and n != '':
                
                st.write('')
                st.write('')
                
                # map str on float or scream
                try:
                    C11 = float(C11)
                    C12 = float(C12)
                    C44 = float(C44)
                    V = float(V)
                    n = float(n)
                except:
                    st.write('')
                    st.write('Wrong type of data given. C11, C12, C44, V and n must be the numbers!')
                
                # if button is pushed
                if st.button('Calculate lattice thermal conductivity for given parameters', type='primary'):
                
                    st.write('')
                    st.write(f'Lattice thermal calculations for {ion1+ion2+ion3} basing on Slack equation in progress...')
                    st.write('')
                    
                    calculate_kappa(C11, C12, C44, V, n)
                    st.balloons()
                    
                # if not pushed
                else:
                    st.write('')
            

            
            # if some data is missing...
            elif C11 != '' or C12 != '' or C44 != '' or V != '' or n != '':
                st.write('At least one value is not given...')
            
        elif choice == Slack_options[1]:
        
            st.warning('For thermoelectic calculations following Slack formula,\
            vasprun.xml (cell volume & number of atoms in a cell) and OUTCAR  with elastic calculations performed are required.')

            # path to elats OUTCAR / vasprun:
            # "C:\Users\kajka\PhD\hH\HfRuTe\elast\OUTCAR"
            # "C:\Users\kajka\PhD\hH\HfRuTe\elast\vasprun.xml"

            with st.expander("INPUT for elastic VASP calculations"):
                st.write('')
                st.write('ENCUT = 500      # high plane-wave basis cutoff')
                st.write('NELM = 500       # high maximum number of SCF steps')
                st.write('ISMEAR = 0       # gaussian smearing for semiconductors')
                st.write('SIGMA = 0.01     # small for semiconductors')
                st.write('LMAXMIX = 4      # mixing for d-electron systems')
                st.write('LASPH = .TRUE.   # aspherical PAW for metaGGA')
                st.write('IBRION = 6   # calculations of elastic constants')
                st.write('ISIF = 3     # stress tensor etc.')

            col6, col7 = st.columns(2)

            with col6:
                f_vasprun = st.file_uploader("runvasp.xml file 👇", accept_multiple_files=False, type=['xml'])
            
            with col7:
                f_output = st.file_uploader("OUTPUT file 👇", accept_multiple_files=False)

            # uploaded_files = st.file_uploader("Choose a runvasp.xml file", accept_multiple_files=True)
            
            # here some condition is required - if runvasp.xml / OUTCAR are not complete
            # for uploaded_file in uploaded_files:
                # bytes_data = uploaded_file.read()
                # st.write("filename:", uploaded_file.name)
                # st.write(bytes_data)
            
            # implement method for gathering from this file:
            # ion1, ion2, ion3, volume, nr of atomns, C11, C12, C44 + lattice param?
                
                        
            _not_supported()
        
    elif choice == 'Credits':
    
        _not_supported()
        
        st.write("")
        st.write("")
        st.write("")
        
        st.subheader("Credits")

        st.write("")
        st.write("")
        
        st.image('git_qr.png', width = 250, caption = 'Git repository of half-Sup')

        st.markdown("""
        
        ##### By:
        + [Kaja Bilińska](mailto:k.bilinska@intibs)
        """)
        
        # add github QR code to this repo
        #         + **[Kaja Bilińska](https://www.youtube.com/channel/UCDn-FahQNJQOekLrOcR7-7Q)**
                ### NLP Simple Examples (App with Streamlit and TextBlob)

    

if __name__ == '__main__':
	main()
