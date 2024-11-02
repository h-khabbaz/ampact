from flask import Flask, request, render_template

app = Flask(__name__)

@app.route('/')
def index():
    return render_template('index.html', length=None)

@app.route('/calculate', methods=['POST'])
def calculate():
    input_string = request.form['input_string']
    f = open("input_peptide_sequences.txt", 'w')
    f.write(input_string)
    f.close()



    #Start Classifier#Start Classifier#Start Classifier#Start Classifier#Start Classifier#Start Classifier#Start Classifier#Start Classifier
    #from __future__ import print_function
    import datetime
    import ipc
    #import matplotlib
    #import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import csv
    import statistics
    #from propy import PyPro #Requires Installing Propy Package
    #import propy #Requires Installing Propy Package
    #from propy.PyPro import GetProDes #Requires Installing Propy Package
    import peptides
    import hmomentk

    import pickle


    import argparse
    import math
    import os
    import time
    labelDict={'negative':'non-AMP' , 'positive': 'AMP'}

    aa_charge = {'E': -1, 'D': -1, 'K': 1, 'R': 1}

    #
    b= open('input_peptide_sequences.txt').read().splitlines()
    def Merge(dict1, dict2):
        return(dict2.update(dict1))
    pp=[]

    #Create Dictionary of Sequences
    c={}
    x = range(1,len(b)+1)
    dic={}
    for i in b:
        dic[b.index(i)+1]=i
    #Calculate Properties
    aggDict = pickle.load(open('aggDict.pkl', 'rb'))#Requires aggDict File
    chargeDict = pickle.load(open('chargeDict.pkl', 'rb'))#Requires chargeDict File
    #isoElP=[]
    netcharge=[]
    agg=[]
    chargeDensities=[]
    mw=[]
    hm=[]
    hp=[]
    isop=[]

    dsss={"seq":b}
    for key, val in dic.items():
        agList=[]
        chargeList=[]
        mwList=[]
        for j in list(val):
            agList.append(aggDict[j])
            chargeList.append(chargeDict[j])
        agg.append(statistics.mean(agList))
        chargeDensities.append(sum(chargeList)/ipc.calculate_molecular_weight(i))
        mw.append(ipc.calculate_molecular_weight(val))
        netcharge.append(sum(chargeList))
        pepseq = peptides.Peptide(val)
        hp.append(-1*peptides.Peptide.hydrophobicity(pepseq,scale="KyteDoolittle"))
        hm.append(hmomentk.hydrop_moment(pepseq))
        isop.append(peptides.Peptide.isoelectric_point(pepseq,pKscale="Lehninger"))

    dsss["molecularWeight"]=mw
    dsss["netCharge"]=netcharge
    dsss["chargeDensity"]=chargeDensities
    dsss["aggregationPropensityInVivo"]=agg
    dsss["hydrophobicity"]=hp
    dsss["hydrophobicMoment"]=hm
    dsss["isoelectricPoint"]=isop
    df=pd.DataFrame(data=pp)
    p5=pd.DataFrame(data=dsss)
    p6=pd.concat([p5, df.reindex(p5.index)], axis=1)
    #Rescaling Data
    #Needs Final Features MinMaxes
    p=[]
    features=["molecularWeight","netCharge","chargeDensity","aggregationPropensityInVivo","hydrophobicity","hydrophobicMoment","isoelectricPoint"]
    d=[]
    newds={"id":list(dic.keys())}
    fMinMax=pd.read_csv("ampNonAMPminMaxesTable50.csv")
    for j in features:
        p=[]
        for i in p6.loc[:,j]:
            p.append((i-fMinMax.loc[0, j])/(fMinMax.loc[1, j]-fMinMax.loc[0, j]))
        newds[j]=p
    newdf=pd.DataFrame(data=newds)
    newdf
    rfModel = pickle.load(open('rf_AMP_classifier_final_model.sav', 'rb'))
    #svcModel = pickle.load(open('svc_AMP_classifier_final_model.sav', 'rb'))
    mr=newdf.iloc[:, 0].values
    print('-------AMP/non-AMP Classification Results-------')
    print('Random Forest Classifier')
    staphActCheckList=[]
    staphActCheckListNo=[]
    for i in mr:
        take=np.where(mr == i)
        pred=rfModel.predict(newdf.iloc[take[0],1:].values)
        #print(i,pred[0])
        print(i,labelDict[pred[0]],dsss['seq'][i-1])
        if labelDict[pred[0]]=='AMP':
            staphActCheckList.append(dsss['seq'][i-1])
            staphActCheckListNo.append(i)





    #End Classifier#End Classifier#End Classifier#End Classifier#End Classifier#End Classifier#End Classifier#End Classifier#End Classifier


    string_length = len(input_string)+5
    return render_template('index.html', label=labelDict[pred[0]], input_string=input_string)

if __name__ == '__main__':
    app.run(debug=True)