#! /usr/bin/python2

import os
import pickle
import pylab

from sympy import symbols
from sympy import sympify

import routes
import catSubType
import libWholeSet
import HPlibraryReader
import resultsReader



def setOfTypesUpdFilesMaker(Nmin,Nmax,kinParameters):
    """
    integer, integer, list of floats -> list of lists of CatSubType
    """
    route=routes.route
    for maxLength in range(Nmin,Nmax+1):
        for parameter in kinParameters:
            kP=parameter
            a=kP[0]
            alpha=kP[1]
            d=kP[2]
            Km=kP[3]
            eH=kP[4]
            kinParStr='_a'+str(a)+'-alpha'+str(alpha)+'-d'+str(d)+'-Km'+str(Km)+'-eH'+str(eH)
            setFile2preKinetics(maxLength)
            preKinetics2Workfile(maxLength,kP)
            setOfTypesUpd=setPickle2setWithConc(maxLength,kP)
    return None

def bestSeqsFinder(Nmin,Nmax,kinParameters):
    """
    
    """
    route=routes.route
    bestSeqsDict={}
    for maxLength in range(Nmin,Nmax+1):
        bestSeqsDict[maxLength]=[]
        for parameter in parameters:
            kP=parameter
            a=kP[0]
            alpha=kP[1]
            d=kP[2]
            Km=kP[3]
            eH=kP[4]
            kinParStr='_a'+str(a)+'-alpha'+str(alpha)+'-d'+str(d)+'-Km'+str(Km)+'-eH'+str(eH)
            setOfTypesUpd=pickle.load(open(route+'/results/'+str(maxLength)+'/'+"setOfTypes"+kinParStr+".p","rb"))
            DictOfTypes=lengthSorter(setOfTypesUpd)
            DictOfTypesSorted=DictOfTypesConcSorter(DictOfTypes)
            maxConcTypeDict=maxConcTypeFinder(DictOfTypesSorted)
            bestSeqsDict[maxLength].append(maxConcTypeDict)
            
    return bestSeqsDict


###debugging###
setOfTypesUpdFilesMaker(Nmin,Nmax,kinParameters)






