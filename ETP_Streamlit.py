"""
# 1. Importation des librairies necessaires

1. pandas for data processing
2. numpy for matrix algebra
3. matplotlib for graph plotting

"""


from cgitb import html
from multiprocessing.sharedctypes import Value
import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
from math import*
import math
import numpy as np
from base64 import b64decode
import base64
#from sklearn.datasets  import load_diabetes
#from keybert import KeyBERT
from PIL import Image
from io import BytesIO
#from pyxlsb import open_workbook as open_xlsb
from openpyxl import Workbook
import plotly.express as px



#------------------------------------------------------------#
# Les fonctions pour l'ETP

# U list en m/s
# Zv int en m
# return list

def U2(U,Zv):
    U2 = []
    for i in U:
        U2.append(i*(4.87/(np.log(67.8*Zv-5.42))))
    return U2

# Z altitude de la station
def Pa(Z):
    return 101.3 * (((293-0.0065*Z)/293)**5.26)

def Tmoy(Tmin, Tmax):
    Tm = []
    for i in range(0,len(Tmin)):
        Tm.append((Tmin[i] + Tmax[i])*0.5)
    return Tm

def delta(Tmoy):
    delt = []
    for i in Tmoy:
        delt.append(4098*(0.6108*np.exp(17.27*i/(i+237.3)))/((i+237.3)**2))
    return delt

#P Pression atmosphérique
def gamma(P):
    return 0.000665 * P

# T peut être Tmin ou Tmax

def eo(T):
    return 0.6108 * np.exp((17.27*T)/(T+237.3))

def es(Tmin,Tmax):
    e=[]
    for i in range(0,len(Tmin)):
        e.append((eo(Tmax[i])+eo(Tmin[i]))/2)
    return e

def ea(Tmin,Tmax,RHmin,RHmax):
    e = []
    for i in range(0, len(Tmin)):
        e.append((eo(Tmin[i])*(RHmax[i]/100) + eo(Tmax[i])*(RHmin[i]/100))*0.5)
    return e

# Il peut être déterminé pour chaque jour D du mois M
def J(D,M,bissextile = False):
    e=[]
    for i in range(0, len(D)):
        r = int(275*M[i]/9-30+D[i])-2
        if M[i]<3:
            r =r+2
        if (bissextile == True and M[i] > 2):
            r=r+1
        e.append(r)
    return e
    
def dr(J):
    e=[]
    for i in J:
        e.append(1 + 0.033 * np.cos(2 * np.pi * i / 365))
    return e


def sigma(J):
    e=[]
    for i in J:
        e.append(0.409 * np.sin(2 * np.pi * i / 365 - 1.39))
    return e

def phi(Degre,Minute, Hemisphere):
    if Hemisphere == 'N':
        degre_deci = Degre + Minute/60
    if Hemisphere == 'S':
        degre_deci = -Degre - Minute/60
    return (np.pi / 180) * degre_deci


def ws(phi, sigma):
    e=[]
    for i in sigma:
        e.append(np.arccos(- np.tan(phi) * np.tan(i)))
    return e


def Ra(dr,ws,phi,sigma):
    Gsc = 0.0820
    e=[]
    for i in range(0,len(sigma)):
        e.append((24*60/np.pi) * Gsc * dr[i] * (ws[i]* np.sin(phi) * np.sin(sigma[i]) + np.cos(phi) * np.cos(sigma[i]) * np.sin(ws[i])))
    return e

def Rso(Z, Ra):
    e=[]
    for i in range(0,len(Ra)):
        e.append((0.75 + 2 * 10**-5 * Z) * Ra[i])
    return e

def Rs(Insolation,Ra):
    e = []
    for i in range(0,len(Ra)):
        e.append((0.25+0.5*7.1/Insolation[i])*Ra[i])
    return e


def Rns(Rs, albedo = 0.23):
    e = []
    for i in Rs:
        e.append((1- albedo)*i)
    return e

def Rnl(Tmax,Tmin,Rs,Rso,ea):
    e = []
    for i in range(0,len(Tmin)):
        e.append(4.903 * 10**-9 * (((Tmax[i]+273.16)**4 + (Tmin[i]+273.16)**4)/2)*(0.34 - 0.14 * np.sqrt(ea[i])) * (1.35 * (Rs[i]/Rso[i]) - 0.35))
    return e

def Rn(Rns, Rnl):
    e = []
    for i in range(0,len(Rns)):
        e.append(Rns[i] - Rnl[i])
    return e

def Lambda(Tmean):
    e = []
    for i in range(0,len(Tmean)):
        e.append(2.501-(2.361*10**-3)*Tmean[i])
    return e

# Flux thermique du sol. Pour le pas journalier, G=0

def G():
    return 0

def RH(RHmin, RHmax):
    Tm = []
    for i in range(0,len(RHmin)):
        Tm.append((RHmin[i] + RHmax[i])*0.5)
    return Tm

# La formule de FAO Penman-Monteith (----) Journalier

def ETP_Penman_M_FAO(Tmax,Tmin,RHmax,RHmin,Rs,Vent10,Degre,Minute,Hemisphere,altitude,Jour,Mois,albedo=0.23):
    Z=altitude
    P=Pa(Z)
    g = gamma(P)
    U=Vent10
    Vent2=U2(U,10)
    Tmean=Tmoy(Tmin, Tmax)
    es1 = es(Tmin,Tmax)
    d = delta(Tmean)
    ea1= ea(Tmin,Tmax,RHmin,RHmax)
    Jo=J(Jour,Mois,bissextile = False)
    phi1=phi(Degre,Minute, Hemisphere)
    s=sigma(Jo)
    ws1=ws(phi1, s)
    dr1=dr(Jo)
    Ra1=Ra(dr1,ws1,phi1,s)
    Rso1=Rso(Z, Ra1)
    Rnl1=Rnl(Tmax,Tmin,Rs,Rso1,ea1)
    Rns1=Rns(Rs, albedo)
    Rn1=Rn(Rns1, Rnl1)
    G=0
    ETP=[]
    for i in range(0,len(Tmin)):
        ETP.append((0.408*d[i]*(Rn1[i]-G)+g*(900/(Tmean[i]+273))*Vent2[i]*(es1[i]-ea1[i]))/(d[i]+g*(1+0.34*Vent2[i])))
    return ETP

# La méthode de Romanenko (1961) Journalier

def ETP_Romanenko(Tmax,Tmin,RHmax,RHmin):
    Tmean=Tmoy(Tmin,Tmax)
    es1 = es(Tmin,Tmax)
    ea1= ea(Tmin,Tmax,RHmin,RHmax)
    ETP = []
    for i in range(0,len(Tmin)):
        ETP.append(4.5*(1+Tmean[i]/25)*(1-ea1[i]/es1[i]))
    return ETP

# La méthode de Hargreaves (1975)

def ETP_Hargreaves(Tmax,Tmin,Degre,Minute,Hemisphere,Jour,Mois):
    Tmean=Tmoy(Tmin,Tmax)
    la=Lambda(Tmean)
    Jo=J(Jour,Mois,bissextile = False)
    phi1=phi(Degre,Minute, Hemisphere)
    s=sigma(Jo)
    ws1=ws(phi1, s)
    dr1=dr(Jo)
    Ra1=Ra(dr1,ws1,phi1,s)
    ETP =[]
    for i in range(0,len(Tmin)):
        ETP.append(0.0023*(Ra1[i]/la[i])*np.sqrt(Tmax[i]-Tmin[i])*(Tmean[i] + 17.8))
    return ETP

# La méthode de 

def ETP_Jensen_Haise(Tmax,Tmin,Rs):
    Tmean=Tmoy(Tmin,Tmax)
    la=Lambda(Tmean)
    ETP =[]
    for i in range(0,len(Tmin)):
        ETP.append((Rs[i]/la[i])*(0.025*Tmean[i]+0.08))
    return ETP

def ETP_Oudin(Tmax,Tmin,Rs):
    Tmean=Tmoy(Tmin,Tmax)
    ETP =[]
    for i in range(0,len(Tmin)):
        ETP.append(Rs[i]*(Tmean[i]+5)/100)
    return ETP

def ETP_Turc(Tmin,Tmax,RHmin,RHmax,Rs):
    Tmean=Tmoy(Tmin,Tmax)
    RHmean = RH(RHmin, RHmax)
    ETP =[]
    for i in range(0,len(Tmin)):
        if RHmean[i] > 50 or RHmean[i] == 50:
            ETP.append(0.01333*(23.9001*Rs[i]+50)*(Tmean[i]/(Tmean[i]+15)))
        else:
            ETP.append(0.01333*(23.9001*Rs[i]+50)*(Tmean[i]/(Tmean[i]+15))*(1+(50-RHmean[i])/70))
    return ETP
    
def ETP_Dalton(Tmax,Tmin,RHmax,RHmin,Vent10):
    U=Vent10
    Vent2=U2(U,10)
    es1 = es(Tmin,Tmax)
    ea1= ea(Tmin,Tmax,RHmin,RHmax)
    ETP=[]
    for i in range(0,len(Tmin)):
        ETP.append((0.3648+0.07223*Vent2[i])*(es1[i]-ea1[i]))
    return ETP






#------------------------------------------------------------#





# Construction de l'application web


st.set_page_config(
    page_title="Calcul de l'évapotranspiration avec plusieurs méthodes", layout="wide"
    )

c30, c31, c32 = st.columns([50, 1, 3])
with c30:
    st.title("🌧️ Calcul de l'évapotranspiration avec plusieurs méthodes")
    #st.image("Evapotranspiration.png", width=100)
    st.header("")

with st.expander(" ℹ️  Informations", expanded=True):
    st.write(
        """
- Cette page permet le calcul de l'évapotranspiration à l'aide de plusieurs méthodes
- Pour ce faire, vous devez avoir un fichier excel au format csv contenant les données nécessaire en fonction de la formule utilisée.
- Vous pouvez vous référer au fichier exemple dans le menu.
            """

    )
    st.markdown("")

st.markdown("")
st.markdown("## **📇Jeu de données**")



# Bar de Menu pour importer les données

st.sidebar.header('Choisir un fichier Excel')
uploaded_file = st.sidebar.file_uploader('Choisir le fichier au format csv')
st.sidebar.markdown("""
[Exemple de fichier](https://drive.google.com/file/d/16AmTiI6sTxFg1Zs18JprwL1O8H8BUxSA/view?usp=sharing)

""")
st.sidebar.header('Méthode')

Methode = st.sidebar.selectbox("Choisir la méthode : ", ['Formule de Penman (FAO)', 'Formule de Hargreaves', 'Formule de Jensen-Haise', 'Formule de Romanenko', 'Formule de Oudin', 'Formule de Turc', 'Formule de Dalton'])


st.sidebar.header('Paramètres du model')

Degree = st.sidebar.number_input('Latitude de la station (degré)')
Minute = st.sidebar.number_input('Latitude de la station (minute)')
albedo = st.sidebar.number_input("Entrez l'albedo")
Altitude = st.sidebar.number_input('Altitude de la station (m)')

Hemisphere = st.sidebar.radio("Choisir l'Hemisphere", ['N', 'S'])

# Pour ajouter une image
if Methode == 'Formule de Penman (FAO)':
    st.markdown(f"""
    La {Methode} est la suivante : 
    """)
    st.latex(r'''
    ETP = \frac{0.408 (R_n - G) + \gamma \frac{900}{T_{moy} + 273} u_2 (e_s - e_a)}{\Delta + \gamma (1 + 0.34u_2) }
    ''')
    image = Image.open('Les Formules référence.png')
    st.image(image)

if Methode == 'Formule de Hargreaves':
    st.markdown(f"""
    La {Methode} est la suivante : 
    """)
    st.latex(r'''
    ETP = 0.0023(\frac{R_a}{\lambda}) \sqrt{T_{max} - T_{min}} (T_{moy} + 17.8)
    ''')
    image = Image.open('Les Formules référence.png')
    st.image(image)

if Methode == 'Formule de Jensen-Haise':
    st.markdown(f"""
    La {Methode} est la suivante : 
    """)
    st.latex(r'''
    ETP = \frac{R_s}{\lambda} (0.025 T_{moy} +0.08) 
    ''')
    image = Image.open('Les Formules référence.png')
    st.image(image)

if Methode == 'Formule de Romanenko':
    st.markdown(f"""
    La {Methode} est la suivante : 
    """)
    st.latex(r'''
    ETP = 4.5(1 + \frac{T_{moy}}{25})(1-\frac{e_a}{e_s})
    ''')
    image = Image.open('Les Formules référence.png')
    st.image(image)


if Methode == 'Formule de Oudin':
    st.markdown(f"""
    La {Methode} est la suivante : 
    """)
    st.latex(r'''
    ETP = R_s \frac{(T_{moy} + 5)}{100}    
    ''')
    image = Image.open('Les Formules référence.png')
    st.image(image)

if Methode == 'Formule de Turc':
    st.markdown(f"""
    La {Methode} est la suivante : 
    """)
    st.write("Si l'humidité est supérieure ou égale à 50 %")
    st.latex(r'''
    ETP = 0.01333 (23.9001 R_s + 50) (\frac{T_{moy}}{T_{moy} + 15})
    ''')
    st.write("Si l'humidité est inférieure à 50 %")
    st.latex(r'''
    ETP = 0.01333 (23.9001 R_s + 50) (\frac{T_{moy}}{T_{moy} + 15}) (1+ (50 - RH_{moy}))
    ''')  
    image = Image.open('Les Formules référence.png')
    st.image(image)

if Methode == 'Formule de Dalton':
    st.markdown(f"""
    La {Methode} est la suivante : 
    """)
    st.latex(r'''
    ETP = 0.3648 + 0.07223 * u_2 (e_s - e_a)
    ''')
    image = Image.open('Les Formules référence.png')
    st.image(image)

# Page Principale

#st.subheader('Jeu de données')

# Implémentation de l'ETP

#@st.cache


def plotting(Result, Methode):
    plt.figure(figsize=(19.2, 10.8))
    plt.plot(Result, label=f'{Methode}')
    plt.legend()
    plt.xlabel('Jour')
    plt.ylabel('ET (mm/day)')
    plt.savefig(f'{Methode}.png', dpi=200)
    return st.pyplot()

def plotting1(data,Result):
    date = list(range(len(Result)))
    d = {'Temps' : date, 'ETP (mm/Jour)' : Result}
    #fig = px.line(d, x='Date',y='ETP',text='ETP')
    st.line_chart(data = d, x='Temps', y='ETP (mm/Jour)')

def load_data():
    url = 'https://drive.google.com/file/d/16AmTiI6sTxFg1Zs18JprwL1O8H8BUxSA/view?usp=sharing'
    html = pd.read_html(url, header = 0)
    data = html[0]
    return data
def load_data2():
    data = pd.read_csv("Fichier_test.csv", sep=",")
    return data

def filedownloadexcel(data):
    writer = pd.ExcelWriter('update2.xlsx')
    data.to_excel(writer, index=False, header=True, encoding='utf-8')
    with open(writer, 'rb') as f:
        b64 = base64.b64encode(f.read()).decode()
        href = f'<a href="data:file/xls;base64,{b64}" download="resultat.xlsx">Download xlsx</a>'
    st.write(href, unsafe_allow_html=True)


def filedownloadcsv(data):
    csv = data.to_csv(index=False)
    b64 = base64.b64encode(csv.encode()).decode()
    href = f'<a href="data:file/csv;base64,{b64}" download="model_performance.csv">Download CSV'


def build_model(data, Methode):
    Tmax=data['Tmin']
    Tmin=data['Tmin']
    RHmax=data['Rhmax']
    RHmin=data['Rhmin']
    Rs=data['Rs']
    Vent10 = data['Vent10']
    Jour = data['Jour']
    Mois = data['Mois']
    Degre = Degree
    altitude = Altitude

    st.markdown("")
    
    if Methode == 'Formule de Penman (FAO)':
        Result = ETP_Penman_M_FAO(Tmax,Tmin,RHmax,RHmin,Rs,Vent10,Degre,Minute,Hemisphere,altitude,Jour,Mois,albedo=0.23)
    
    if Methode == 'Formule de Hargreaves':
        Result = ETP_Hargreaves(Tmax,Tmin,Degre,Minute,Hemisphere,Jour,Mois)
    
    if Methode == 'Formule de Jensen-Haise':
        Result = ETP_Jensen_Haise(Tmax,Tmin,Rs)
    
    if Methode == 'Formule de Romanenko':
        Result = ETP_Romanenko(Tmax,Tmin,RHmax,RHmin)

    if Methode == 'Formule de Oudin':
        Result = ETP_Oudin(Tmax,Tmin,Rs)

    if Methode == 'Formule de Turc':
        Result = ETP_Turc(Tmin,Tmax,RHmin,RHmax,Rs)

    if Methode == 'Formule de Dalton':
        Result = ETP_Dalton(Tmax,Tmin,RHmax,RHmin,Vent10)



    def _max_width_():
        max_width_str = f"max-width: 1400px;"
        st.markdown(
            f"""
        <style>
        .reportview-container .main .block-container{{
            {max_width_str}
        }}
        </style>
        """,
            unsafe_allow_html=True,
    )
    _max_width_()   

    with st.form(key="my_form"):

        st.markdown("## **⏬ Affichage des résultats**")


        if st.form_submit_button('Appuyer pour afficher le resultats'):
            d = {'Date' : data['Mois'], 'ETP' : Result}
            st.dataframe(d)
          
        st.markdown("## **🖨️ Telechargement vers un fichier excel**")

        if st.form_submit_button('Appuyer pour telecharger le fichier'):
            list1 =Result[col1].tolist()
            col1 = "ETP"
            donnee = pd.DataFrame({col1: list1})
            donnee.to_excel('fichier_result.xlsx', sheet_name='sheet1', index=False)
            st.markdown('Le fichier a été telecharger avec succès')
            st.dataFrame(Result)
        st.markdown("## **📈 Affichage du graphe**")

        if st.form_submit_button('Appuyer pour afficher la représentations graphique'):
            st.header("Graphe montrant la variation de l'Etp en fonction du temps")
            #plotting(Result, Methode)
            st.write("Appuiyez sur le bouton dans le coin en haut à droite sur le carte pour sauvegarder l'image")
            plotting1(data,Result)




#------------------------------#
if uploaded_file is not None:
    data = pd.read_excel(uploaded_file)
    if st.button('Appuyer pour afficher les données'):
        st.write(data)
    st.metric(label= "Albédo : ", value = f"{albedo}")
    #st.latex(f'''
    #Albédo : {albedo}
    #''')   
    st.metric(label= "Latitude : ", value = f"{Degree}° {Minute}'")
    #st.latex(fr'''
    #Latitude : {math.floor(Degree)}°{math.floor(Minute)}'
    #''')
    if Hemisphere == 'N':
        st.write('Hemisphère : Nord')
        #st.latex(r'''
        #Hemisphère : Nord
        #''')   
    if Hemisphere == 'S':
        st.write('Hemisphère : Sud')
        #st.latex(r'''
        #Hemisphère : Sud
        #''')   
    build_model(data, Methode)
else:
    st.info('En attente de chargement du fichier.')
    data = load_data2()
    if st.button('Appuyer pour utiliser des données par défaut'):
        st.markdown('Les donnees par défaut sont utilisés')
        st.write(data)
        # .head(5)
    st.metric(label= "Albédo : ", value = f"{albedo}")
    st.metric(label= "Latitude : ", value = f"{Degree}° {Minute}'")
    if Hemisphere == 'N':
        st.write('Hemisphère : Nord') 
    if Hemisphere == 'S':
        st.write('Hemisphère : Sud')

    build_model(data, Methode)
