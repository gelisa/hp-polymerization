#! /usr/bin/python2

import nativeChain

def classesGenerator():
    classes=[]
    for first in [4,3,2]:
        for second in range(first,-1,-1):
            if not second==0:
                classes.append('H'*first+'P'+'H'*second)
            else:
                classes.append('H'*first+'PP')
    return classes   




#def decomposition2Patter((primaryH, Pnumber, secondaryH)):
    #return primaryH*'H'+Pnumber*'P'+secondaryH*'H'



def thisPatternSumGenerator(pattern, patternDict,catalysts,substrateConcentrations):
    patternSum=0
    #for every element in the list of working patterns for a given pattern
    for currPattern in patternDict[pattern]:
        #list of all applicable catalysts in terms of substrates
        catalystTmp=seqType.CatalystType(currPattern)
        if catalystTmp in catalysts.catalystList:
            elementNumber=catalysts.catalystList.index(catalystTmp)
            allCatsInSubsTerms=catalysts.catalystList[elementNumber].substrateTypeDecomposition
            #let's sum all the concentrations corresponding to a given pattern
            currSum=0
            for item in allCatsInSubsTerms:
                currSum+=item.quantity*substrateConcentrations[item.number-1]
            patternSum+=currSum
    thisPatternSum=patternSum
    return thisPatternSum

def patternSumGenerator(patternDict,catalysts,substrateConcentrations):
    patternSum={}
    for item in patternDict:
        patternSum[item]=thisPatternSumGenerator(item, patternDict,catalysts,substrateConcentrations)
    return patternSum


def answDictGenerator(answ, substrates):
    answDict={}
    for i in range(len(answ)):
        answDict[substrates.substrateList[i]]=answ[i]
    
    return answDict








