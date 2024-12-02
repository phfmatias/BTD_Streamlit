##   PYTHON FILE HEADER #
##
##   File:         [texts.py]
##
##   Author(s):    ['Pedro H.F Matias']
##   Site(s):      ['https://github.com/phfmatias']
##   Email(s):     ['phfmatias@discente.ufg.br']
##   Credits:      ['Copyright © 2024 LEEDMOL. All rights reserved.']
##   Date:         ['14.11.2024']
##   Version:      ['1.0.0']
##   Status:       ['Development']
##   Language:     ['Python']
##   Description:  ['Here are the texts of the application.']

### IMPORTS ###

import lorem

### CODE ###

class Texts():
    def __init__(self) -> None:
        pass

    def text1(self):
        text1 = """
        <div>
            <h3>Welcome</h3>
            <p>
                Welcome to our web application for predicting the photophysical properties of benzothiazole-diazole (BTD) derivatives! 
                This tool is designed to aid researchers and scientists in exploring the optical behavior of BTD compounds, widely recognized 
                for their applications in bioimaging, optoelectronics, and material sciences.
            </p>
            <p>
                By leveraging state-of-the-art machine learning algorithms, the application offers reliable predictions of key properties, 
                such as absorption and emission wavelengths and Stokes shifts. Users can input the SMILE code for a specific BTD derivative 
                or draw the molecule directly on the interface. The model will then process the data and obtain instant results, enabling 
                rapid screening of candidate molecules for various applications.
            </p>
            <h4>Main features:</h4>
            <ul>
                <li>Drawing of BTD molecules directly on the interface</li>
                <li>Prediction of absorption and emission wavelengths</li>
                <li>Prediction of Stokes shifts</li>
            </ul>
            <h4>How to use:</h4>
            <ol>
                <li>Go to the <span style="background-color: #007bff; color: white; padding: 2px 5px; border-radius: 3px;">Predictor</span> tab</li>
                <li>Select the model you want to use</li>
                <li>Input the SMILE code or draw the molecule</li>
                <li>Select the solvent</li>
                <li>Click the <span style="background-color: #007bff; color: white; padding: 2px 5px; border-radius: 3px;">Predict</span> button</li>
            </ol>
        </div>
        """
        return text1

    def text2(self):
        text2 = '''
            <div style="text-align: justify;">
                Benzothiadiazole (BTD) derivatives show promise in advanced photophysical applications, but designing molecules with optimal desired properties remains challenging due to complex structure-property relationships. Existing computational methods have a high cost to predict precise photophysical characteristics. This study employs machine learning with Morgan fingerprints to forecast BTD derivatives' maximum absorption and emission wavelengths. Random forest achieved R2 values of 0.92 for absorption and 0.89 for emission, validated internally with 10-fold cross-validations and externally with experimental datasets. SHAP analysis revealed critical design insights, highlighting tertiary amine presence and solvent polarity as key drivers of red-shifted emissions. By developing a web-based predictive tool, we demonstrate machine learning's potential to accelerate molecular design, providing researchers a powerful approach to engineer BTD derivatives with enhanced photophysical properties.
            </div>
        '''
        return text2
    
    def text3(self):
        text3 = '''
        If you use this application, please cite the following paper:
        '''
        return text3
    
    def text4(self):
        text4 = '''
    @article{example2024,
    title={Example Title},
    author={Author, A. and Author, B.},
    journal={Journal of Example Research},
    volume={42},
    pages={123-456},
    year={2024},
    publisher={Example Publisher}
    }
        '''
        return text4

    
    def text5(self):
        text5 = '''
        The dataset used in this study is available below and can be downloade as a CSV file. Also, the DOI of the papers used in the dataset are available as well.        
    '''
        
        return text5