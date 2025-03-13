##   PYTHON FILE HEADER #
##
##   File:         [frontend.py]
##
##   Author(s):    ['Pedro H.F Matias']
##   Site(s):      ['https://github.com/phfmatias']
##   Email(s):     ['phfmatias@discente.ufg.br']
##   Credits:      ['Copyright © 2024 LEEDMOL. All rights reserved.']
##   Date:         ['14.11.2024']
##   Version:      ['1.0.0']
##   Status:       ['Development']
##   Language:     ['Python']
##   Description:  ['Frontend of the application.']

### IMPORTS ###

import streamlit as st
import matplotlib.pyplot as plt
import shap
import numpy as np
from pandas import DataFrame
from streamlit_option_menu import option_menu
from streamlit_theme import st_theme
from streamlit_ketcher import st_ketcher
from warnings import filterwarnings
from .drawSmiles import *
import cairosvg
filterwarnings('ignore')

### SELF IMPORTS ###

from srcCode.texts import Texts
from srcCode.backend import Backend

### CODE ###

class Frontend(Backend):
    def __init__(self):

        super().__init__()

        if 'clicked' not in st.session_state:
            st.session_state.clicked = False

        if 'clickedSmiles' not in st.session_state:
            st.session_state.clickedSmiles = False

        if 'display_substructures' not in st.session_state:
            st.session_state.display_substructures = False

        texts = Texts()

        self.text1 = texts.text1()
        self.text2 = texts.text2()
        self.text3 = texts.text3()
        self.text4 = texts.text4()
        self.text5 = texts.text5()

        Frontend.main(self)

    def main(self):

        st.set_page_config(page_title='BTD_ML', page_icon='🔬', layout='centered', initial_sidebar_state='expanded')
        navbar = Frontend.NavSideBar(self)

        if navbar == 'Home':
            st.header('Integrating Machine Learning and SHAP Analysis to Advance the Rational Design of Benzothiadiazole Derivatives with Tailored Photophysical Properties')
            st.markdown('---')
            st.markdown(self.text1, unsafe_allow_html=True)

        elif navbar == 'Predictors':
            
            st.header('Predictors 📈')
            st.markdown('---')            
            self.mod = Frontend.PredictorsMenu(self)
            if self.mod == 'Random Forest':
                st.markdown('### {} 🎲🌲 selected to predict the properties.'.format(self.mod))
            elif self.mod == 'XGBoost':
                st.markdown('### {} ⚡🌲 selected to predict the properties.'.format(self.mod))
            elif self.mod == 'LightGBM':
                st.markdown('### {} 💡 selected to predict the properties.'.format(self.mod))

            st.markdown('---')
            st.markdown('### Insert SMILES 😀')
            Frontend.inputSmiles(self)       
            st.markdown('---')
            st.markdown('### Select Solvent 💧')
            Frontend.selectSolvent(self)
            st.markdown('---')
            if st.button('Predict'):
                Frontend.predict(self)
                Frontend.shap_analysis(self)

        elif navbar == 'About':
            st.header('About ℹ️')
            st.markdown('---')
            st.markdown(self.text2, unsafe_allow_html=True)
            st.markdown('---')
            st.header('Cite Us 📚')
            st.markdown(self.text3, unsafe_allow_html=True)
            st.code(self.text4)

        elif navbar == 'Data':
            st.header('Data 📊')
            st.markdown('---')
            st.write(self.text5)

            st.write(self.dataDL)

            with open('Data/no_missing_data.csv', 'r') as file:
                st.download_button('Download Data', file, 'no_missing_data.csv', 'csv')
                

    def NavSideBar(self):

        themeAsw = st_theme()

        with st.sidebar:
            if themeAsw['base'] == 'light':
                st.image('static/logoleed_cor.png')

                nav = option_menu("Menu", ["Home", "Predictors", "About", "Data"],
                                icons=["house", "cpu", "info-circle", "journal-bookmark-fill", "database-fill"],
                                menu_icon="grid", default_index=0, orientation="vertical",
                                styles={
                            "container": {"background-color": "#fafafa"},
                            "icon": {"color": "#ff6600", "font-size": "25px", "margin-left": "2px"}, 
                            "nav-link": {"font-size": "20px", "text-align": "left", "margin":"0px", "--hover-color": "#d7d7d7", "border-radius":"5px", "padding":"5px", "margin-bottom":"10px"},
                            "nav-link-selected": {"background-color": "#007bff", "text-transform": "none"},
                            })

            else:   
                st.image('static/logoleed.png')

                nav = option_menu("Menu", ["Home", "Predictors", "About", "Data"],
                                icons=["house", "cpu", "info-circle", "journal-bookmark-fill"],
                                menu_icon="grid", default_index=0, orientation="vertical",
                                styles={
                            "menu-title": {"color": "#ffffff"},
                            "container": {"background-color": "#0f0f0f"},
                            "icon": {"color": "#ff6600", "font-size": "25px", "margin-left": "2px"}, 
                            "nav-link": {"font-size": "20px", "text-align": "left", "margin":"0px", "--hover-color": "#d7d7d7", "border-radius":"5px", "padding":"5px", "color":"#ffffff", "margin-bottom":"10px"},
                            "nav-link-selected": {"background-color": "#007bff", "text-transform": "none", "color":"#FFFFFF"},
                            })

        return nav


    def PredictorsMenu(self):

        modelsNav = option_menu("Model", ["Random Forest", "XGBoost", "LightGBM"],
                    icons=["tree-fill", "lightning-charge-fill", "lightbulb-fill"], menu_icon="cpu-fill", default_index=0, orientation='vertical',
                    styles = {"icon": {"color": "#ff6600", "font-size": "25px", "margin-left": "2px"}, 
                              "nav-link-selected": {"background-color": "#007bff", "color":"#FFFFFF"}})
            
        return modelsNav
    
    def click_button(self):
        st.session_state.clicked = True

    def click_button_substructure(self):
        st.session_state.display_substructures = True

    def inputSmiles(self):
        
        st.button("Draw molecule instead...", on_click=self.click_button)

        if st.session_state.clicked:
            self.smiles = st_ketcher('C1C=CC2C(=NSN=2)C=1')

            if self.smiles:
                st.info('Current SMILES: {}'.format(self.smiles))         

                similarities, max_percentage, is_trustable = Backend._tanimoto(self, self.smiles)  
                self.max_percentage = max_percentage     

                st.markdown('#### Want to continue?')

                col1, col2 = st.columns(2)

                with col1:
                    if st.button('Yes'):
                        asw=True
                with col2:
                    if st.button('No'):
                        st.error('Please, draw again or try to insert a SMILES instead.')
                        st.session_state.clicked = False 
                        self.smiles = ''

                        if st.button('Click here to insert a SMILES'):
                            st.session_state.clicked = False 
                            
        else:
            self.smiles = st.text_input('Smiles Input', placeholder='Insert SMILES here')
            
            if self.smiles:
                self.smiles = st_ketcher(self.smiles)

                st.info('Current SMILES: {}'.format(self.smiles))

                similarities, max_percentage, is_trustable = Backend._tanimoto(self, self.smiles)  
                self.max_percentage = max_percentage     

                st.markdown('#### This is the molecule that you inserted. Do you want to continue?') 
                            
                col1, col2 = st.columns(2)

                with col1:
                    if st.button('Yes'):  
                        asw=True
                with col2:
                    if st.button('No'):
                        st.error('Please insert another SMILES, or try to draw the molecule instead.')
                        self.smiles = ''

    def selectSolvent(self):
        solvents = ['Dimethylsulfoxide (0.444)', 'Water (1.0)', 'Dichloromethane (0.309)', 'Acetonitrile (0.46)', 'Toluene (0.099)', 'Ethanol (0.654)', 'Chloroform (0.259)', 'Methanol (0.762)', 'Dimethylformamide (0.386)', 'Cyclochexane (0.006)', 'Hexane (0.009)', 'Tetrahidrofurane (0.605)', 'Methyl Cianide (0.46)', 'Acetate (0.355)', 'Methyl Phenil (0.099)', 'Isopropanol (0.546)', 'Dioxane (0.164)', 'mXylene (0.074)', 'Chlorobenzene (0.333)', 'Ethyl Acetate (0.228)', 'DiethylEter (0.117)', 'Octanol (0.537)', 'Nitrobenzila (0.333)', 'Benzene (0.111)', 'Dichloroethane (0.194)']

        solvent = st.selectbox('Select Solvent', solvents)

        self.etn = Backend._solventToEtn(self, solvent)
        
    def predict(self):
        st.markdown('### Predicting... 🧠⚙️')

        fp, bit = Backend._generateMF(self, self.smiles)
        fp = fp.reshape(1, -1)


        self.finger_print = fp
        self.bit = bit

        input = Backend._PrepareInput(self, fp, self.etn)

        if self.mod == 'Random Forest':
            my_em = self.rf_maxem.predict(input)
            my_abs = self.rf_maxabs.predict(input)

            my_shap_em = self.shap_rf_maxem.shap_values(input)
            my_shap_abs = self.shap_rf_maxabs.shap_values(input)

            my_explainer_em = self.shap_rf_maxem.expected_value
            my_explainer_abs = self.shap_rf_maxabs.expected_value

        elif self.mod == 'XGBoost':
            my_em = self.xgb_maxem.predict(input)
            my_abs = self.xgb_maxabs.predict(input)

            my_shap_em = self.shap_xgb_maxem.shap_values(input)
            my_shap_abs = self.shap_xgb_maxabs.shap_values(input)

            my_explainer_em = self.shap_xgb_maxem.expected_value
            my_explainer_abs = self.shap_xgb_maxabs.expected_value

        elif self.mod == 'LightGBM':
            my_em = self.lgbm_maxem.predict(input)
            my_abs = self.lgbm_maxabs.predict(input)

            my_shap_em = self.shap_lgbm_maxem.shap_values(input)
            my_shap_abs = self.shap_lgbm_maxabs.shap_values(input)

            my_explainer_em = self.shap_lgbm_maxem.expected_value
            my_explainer_abs = self.shap_lgbm_maxabs.expected_value

        self.my_shap_em = my_shap_em
        self.my_shap_abs = my_shap_abs
        self.my_explainer_em = my_explainer_em
        self.my_explainer_abs = my_explainer_abs

        self.mf_used_in_model = input

        col1,col2 = st.columns(2)

        with col1:
            st.markdown('#### Predicted max abs: {:.3f}nm'.format(my_abs[0]))
        with col2:
            st.markdown('#### Predicted max em: {:.3f}nm'.format(my_em[0]))

        st.markdown('<div style="text-align: center;"><h4>Predicted stokes shift: {:.3f} nm</div><h4>'.format(my_em[0] - my_abs[0]), unsafe_allow_html=True)

        if self.max_percentage >= 70:        
            st.markdown('<div style="text-align: center;"><h4>Similarity: <span style="color:green">{:.2f}%</span></h4>'.format(self.max_percentage), unsafe_allow_html=True)  

        elif self.max_percentage < 70 and self.max_percentage > 50:
            st.markdown('<div style="text-align: center;"><h4>Similarity: <span style="color:orange">{:.2f}%</span></h4>'.format(self.max_percentage), unsafe_allow_html=True)              
        
        else:
            st.markdown('<div style="text-align: center;"><h4>Similarity: <span style="color:red">{:.2f}%</span></h4></div>'.format(self.max_percentage), unsafe_allow_html=True)

        st.markdown('---')


    def shap_analysis(self):
        st.markdown('### SHAP Analysis 📊')

        col3, col4 = st.columns(2)

        st.markdown('#### Max Absorption')

        features_names = [
        f"{str(self.mf_used_in_model[0][i])[0]} {self.features_names[i]}" if self.features_names[i] != "etn" 
        else self.features_names[i]
        for i in range(len(self.features_names))
        ]

        shap_values_abs_plotter = shap.Explanation(self.my_shap_abs, self.my_explainer_abs, feature_names=features_names)
        
        fig,ax = plt.subplots()
        shap.waterfall_plot(shap_values_abs_plotter[0], max_display=10, show=False)
        st.pyplot(fig)

        top_indices_abs = np.argsort(-np.abs(shap_values_abs_plotter.values[0]))[:10]
        ten_features_abs = [self.features_names[i] for i in top_indices_abs]
        bits = [x.split('_')[1] for x in ten_features_abs if '_' in x]

        st.markdown('---')

        images_presence = []
        images_absence = []
        
        labels_presence = []
        labels_absence = []

        dl = 4
        
        for j in bits:
            bit_number = int(j)
            smiles = self.smiles
            mol = Chem.MolFromSmiles(smiles)
            addcords(mol)
            bi = self.bit
            fp = self.finger_print
            try:
                img = draw_substructure(mol, bi, bit_number, 400, detail_level=dl)
                images_presence.append(img)
                labels_presence.append('MF_{}'.format(bit_number))
            except KeyError:
                molNumber = self.df_drawn.loc[self.df_drawn['MF_{}'.format(bit_number)] == 1].index[0]
                smiles = self.data['smiles'].iloc[molNumber]
                mol = Chem.MolFromSmiles(smiles)
                addcords(mol)
                bi = {}
                fp = AllChem.GetMorganFingerprintAsBitVect(mol, 4, bitInfo=bi, nBits=2048)
                
                img = draw_substructure(mol, bi, bit_number, 400, detail_level=dl)
                images_absence.append(img)
                labels_absence.append('MF_{}'.format(bit_number))
        
        st.markdown('#### Presence of MF')
        
        max_columns = 5
        num_columns = min(len(images_presence), max_columns)

        cols_presence = st.columns(num_columns)

        for i, img in enumerate(images_presence):
            with cols_presence[i % num_columns]:  
                st.text(labels_presence[i])  
                st.image(img, use_container_width=True)

        st.markdown('---')

        st.markdown('#### Absence of MF')

        cols_absence = st.columns(num_columns)

        for i, img in enumerate(images_absence):
            with cols_absence[i % num_columns]:  
                st.text(labels_absence[i])  
                st.image(img, use_container_width=True)

        st.markdown('---')

        st.markdown('#### Max Emission')
        
        shap_values_em_plotter = shap.Explanation(self.my_shap_em, self.my_explainer_em, feature_names=features_names)
        fig,ax = plt.subplots()
        shap.waterfall_plot(shap_values_em_plotter[0], max_display=10, show=False)
        st.pyplot(fig)

        top_indices_abs = np.argsort(-np.abs(shap_values_em_plotter.values[0]))[:10]
        ten_features_abs = [shap_values_em_plotter.feature_names[i] for i in top_indices_abs]
        bits = [x.split('_')[1] for x in ten_features_abs if '_' in x]

        images_presence = []
        images_absence = []
        
        labels_presence = []
        labels_absence = []

        for j in bits:
            bit_number = int(j)
            smiles = self.smiles
            mol = Chem.MolFromSmiles(smiles)
            addcords(mol)
            bi = self.bit
            fp = self.finger_print
            
            
            try:
                img = draw_substructure(mol, bi, bit_number, 400, detail_level=dl)
                images_presence.append(img)
                labels_presence.append('MF_{}'.format(bit_number))
            except KeyError:
                molNumber = self.df_drawn.loc[self.df_drawn['MF_{}'.format(bit_number)] == 1].index[0]
                smiles = self.data['smiles'].iloc[molNumber]
                mol = Chem.MolFromSmiles(smiles)
                addcords(mol)
                bi = {}
                fp = AllChem.GetMorganFingerprintAsBitVect(mol, 4, bitInfo=bi, nBits=2048)
                
                img = draw_substructure(mol, bi, bit_number, 400, detail_level=dl)
                images_absence.append(img)
                labels_absence.append('MF_{}'.format(bit_number))
            
        st.markdown('#### Presence of MF')
        
        max_columns = 5
        num_columns = min(len(images_presence), max_columns)

        cols_presence = st.columns(num_columns)

        for i, img in enumerate(images_presence):
            with cols_presence[i % num_columns]:  
                st.text(labels_presence[i])  
                st.image(img, use_container_width=True)

        st.markdown('---')

        st.markdown('#### Absence of MF')

        cols_absence = st.columns(num_columns)

        for i, img in enumerate(images_absence):
            with cols_absence[i % num_columns]:  
                st.text(labels_absence[i])  
                st.image(img, use_container_width=True)

        st.markdown('---')

