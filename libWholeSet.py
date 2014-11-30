#! /usr/bin/python2

import os
import pickle

from sympy import symbols
from sympy import sympify

import routes

def possiblePatternsGenerator():
    """
    None -> list of strings
    Returns a list of the patterns of catalysis, which seem possible from my point of view.
    """
    classes=[]
    for first in [4,3,2]:
        for second in range(first,-1,-1):
            if not second==0:
                classes.append('H'*first+'P'+'H'*second)
            else:
                classes.append('H'*first+'PP')
    return classes  

def catalystsPopulationsGenerator(setOfTypes):
    """
    list of CatSubType -> dictionary
    For every catalytic pattern returns a list of numbers of Types of catSubs which have this catalytic pattern
    numbers are given in the form of a la x135
    """
    patterns=possiblePatternsGenerator()
    catalystsPopulations={}
    
    for CSType in setOfTypes:
        if CSType.Type[2] in patterns:
            if CSType.Type[2] in catalystsPopulations.keys():
                catalystsPopulations[CSType.Type[2]].append('x'+str(CSType.indx))
            else:
                catalystsPopulations[CSType.Type[2]]=[]
                catalystsPopulations[CSType.Type[2]].append('x'+str(CSType.indx))
            
    
    return catalystsPopulations

def pattern2Decomposition(pattern):
    primaryH=0
    for i in range(len(pattern)):
        if not pattern[i]=='P':
            primaryH+=1
            Pnumber=pattern.count('P')
            secondaryH=pattern[i+1:].count('H')
        else:
            break
        
    return (primaryH, Pnumber, secondaryH)

def DictOfCatalystsForPatternGenerator(catalystsPopulations):
    """
    Say we have hHpp pattern and interested how H is catalyzed,
    Addition of this monomer is catalyzed not only by HHPP catalysts, but also by HHHPP,HHPH and some other.
    For this situation HHPP is a catUp or catDown (see letter.py).
    This function for each catUp or catDown returns a list of catalytic patterns which correspond to this catUp/Down.
    """
            
    #initiate a dictionary which will be {catUp/Down aka pattern: list of catalytic patterns which can catalyze catUp/Down}
    patternDict={}
    #what are possible patterns?
    patterns=possiblePatternsGenerator()
    #for every pattern in the list of possible patterns
    for item in patterns:
        #initiate a list
        patternDict[item]=[]
        #what kind of patern is it? how many primary H and secondary H does it have?
        decomp=pattern2Decomposition(item)
        #if the pattern has 2 Ps it means it doesn't have secondary H
        if decomp[1]==2:
            #for numbers in primaryH..4
            for first in range(decomp[0],5):
                #the pattern can be catalyzed by patterns with secondary Hs as well. 
                #They just need to have the same number of primary Hs.
                #So for every number 0..first
                for second in range(0,first+1):
                    #if this number isn't 0
                    if not second==0:
                        #current pattern is 'H'*first+'P'+'H'*second
                        currPatt='H'*first+'P'+'H'*second
                        #if the number of types which have this catalytic pattern isn't 0
                        if currPatt in catalystsPopulations.keys():
                            #add to the list
                            patternDict[item].append(currPatt)
                    #if second = 0
                    else:
                        #current pattern is 'H'*first+'PP'
                        currPatt='H'*first+'PP'
                        if currPatt in catalystsPopulations.keys():
                            #add to the list
                            patternDict[item].append(currPatt)
        else:
            first=decomp[0]
            for second in range(decomp[2],first+1):
                #current pattern is 'H'*first+'P'+'H'*second
                currPatt='H'*first+'P'+'H'*second
                if currPatt in catalystsPopulations.keys():
                    #add to the list
                    patternDict[item].append(currPatt)
            
    return patternDict


def popsOfCatsDictGenerator(patternDict,setOfTypes):
    popsOfCatsDict={}
    catPops=catalystsPopulationsGenerator(setOfTypes)
    for pattern in patternDict.keys():
        if not patternDict[pattern]==[]:
            #with each element of patternDict[pattern], 
            #so element is one of the catalytic patterns which work for a given catUp/catDown 
            #aka patern in patternDict.keys()
            listOfSums=map(lambda element: reduce(lambda x, y: x+'+'+y,catPops[element]),patternDict[pattern])
            pops=reduce(lambda x, y: x+'+'+y,listOfSums)
            popsOfCatsDict[pattern]=pops
    
    return popsOfCatsDict

def set2EquationsHOM(setOfTypes, a, alpha, d, Km, eH, popsOfCatsDict):
    setOfEquations=[CSType.type2FormulaHOM(a, alpha, d, Km, eH, popsOfCatsDict) for CSType in setOfTypes]
    
    return setOfEquations

def setOfVarsGenerator(setOfTypes):
    setOfVars=map(lambda x: 'x'+str(x.indx),setOfTypes)
    
    return setOfVars

def vars2Symbols(setOfVars):
    varsComma=reduce(lambda x, y: x+', '+y, setOfVars)
    MyInput=varsComma+'= symbols(\''+varsComma+'\')'
    #f = open('conversionFile.py', 'w')
    #f.write('#! /usr/bin/python2')
    #f.write(MyInput)
    #f.close()
    
    return MyInput

def equationsToHOMStrings(setOfEquations):
    HOMStringsTmp=map(lambda x: '('+x[0]+'-'+x[1]+')',setOfEquations)
    print('before sympify')
    HOMStrings=[sympify(x) for x in HOMStringsTmp]
    print('after sympify')
    #HOMStrings=[x.expand(basic=True) for x in HOMStrings]
    for i in range(len(HOMStrings)):
        HOMStrings[i]=HOMStrings[i].expand(basic=True)
        print(str(i))
    print('expanded')
    HOMStrings=[str(x)+';' for x in HOMStrings]
    print('; added')
    HOMStrings=[x.replace("**","^") for x in HOMStrings]
    print('done')
    
    return HOMStrings

def makeSymFile(HOMStrings,workfile):
    f = open(workfile, 'w')
    f.write('{\n')
    for aString in HOMStrings:
        f.write(aString+'\n')
    f.write('}\n')
    f.close()
        
    return None

def kinParamReader(filename):
    """
    string -> list of lists of floats
    """
    inFile = open(filename, 'r', 0)
    lines=inFile.readlines()
    lines=lines[1:]
    paramLines=[line.split(' ') for line in lines]
    paramLines=[[float(param) for param in line] for line in paramLines]
    
    return paramLines


###Unifying functions###

def setFile2preKinetics(maxLength):
    """
    integer -> pickle; returns None
    Fog a given set of kinetic parameters and for a given set of types of chains 
    creates a file for the HOM solver
    """
    #load set of Types of certain length to the memory
    try:
        popsOfCatsDict=pickle.load(open("popsOfCatsDict"+str(maxLength)+".p","rb"))
        setOfVars=pickle.load(open("setOfVars"+str(maxLength)+".p","rb"))
    except:
        setOfTypes=pickle.load(open("setOfTypes"+str(maxLength)+".p","rb"))
        print('setOfTypes done')
        catalystsPopulations=catalystsPopulationsGenerator(setOfTypes)
        print('catalystsPopulations done')
        patternDict=DictOfCatalystsForPatternGenerator(catalystsPopulations)
        print('patternDict done')
        popsOfCatsDict=popsOfCatsDictGenerator(patternDict,setOfTypes)
        print('popsOfCatsDict done')
        setOfVars=setOfVarsGenerator(setOfTypes)
        print('setOfVars done')

        pickle.dump(popsOfCatsDict,open("popsOfCatsDict"+str(maxLength)+".p","wb"))
        pickle.dump(setOfVars,open("setOfVars"+str(maxLength)+".p","wb"))
        print("popsOfCatsDict"+str(maxLength)+" pickled")
        print("setOfVars"+str(maxLength)+" pickled")
    else:
        print("popsOfCatsDict"+str(maxLength)+" already exists")
        print("setOfVars"+str(maxLength)+" already exists")
    
    return None

def preKinetics2Workfile(maxLength,kP):
    """
    number, list -> file; returns None
    """
    setOfTypes=pickle.load(open("setOfTypes"+str(maxLength)+".p","rb"))
    popsOfCatsDict=pickle.load(open("popsOfCatsDict"+str(maxLength)+".p","rb"))
    setOfVars=pickle.load(open("setOfVars"+str(maxLength)+".p","rb"))
    
    MyInput=vars2Symbols(setOfVars)
    print('MyInput done')
    exec MyInput
    print('execed')
    a=kP[0]
    alpha=kP[1]
    d=kP[2]
    Km=kP[3]
    eH=kP[4]
    kinParStr='_a'+str(a)+'-alpha'+str(alpha)+'-d'+str(d)+'-Km'+str(Km)+'-eH'+str(eH)
    filename=routes.routeHOM+'workfile'+str(maxLength)+kinParStr+'.sym'
    try:
        workfile=open(filename, 'r')
    except:
        setOfEquations=set2EquationsHOM(setOfTypes,a, alpha, d, Km, eH, popsOfCatsDict)
        print('setOfEquations done')
        HOMstrings=equationsToHOMStrings(setOfEquations)
        print('HOMstrings done')
        makeSymFile(HOMstrings,filename)
        print('workfile'+str(maxLength)+kinParStr+'.sym'+' written')
    else:
        workfile.close()
        print('workfile'+str(maxLength)+kinParStr+'.sym'+' already exists')
    
    return None

def setPickle2setWithConc(maxLength,kP):
    setOfTypes=pickle.load(open("setOfTypes"+str(maxLength)+".p","rb"))
    setOfVars=pickle.load(open("setOfVars"+str(maxLength)+".p","rb"))
    a=kP[0]
    alpha=kP[1]
    d=kP[2]
    Km=kP[3]
    eH=kP[4]
    kinParStr='_a'+str(a)+'-alpha'+str(alpha)+'-d'+str(d)+'-Km'+str(Km)+'-eH'+str(eH)
    #check if there is a file
    try:
        setOfTypesUpd=pickle.load(open(routes.route+'/results/'+str(maxLength)+'/'+'setOfTypes'+kinParStr+'.p',"rb"))
    #if there is no such file create it
    except:
        
        #preKinetics2Workfile(maxLength,kP)
        os.system('cd ../HOM/HOM4PS2 && ./hom4ps2 workfile'+str(maxLength)+kinParStr+'.sym')
        os.system('cd ../HOM/HOM4PS2 && cp data.roots data'+str(maxLength)+kinParStr+'.roots')
        
        filename='/data/research/06.origins of life/codes/HOM/HOM4PS2/data'+str(maxLength)+kinParStr+'.roots'
        nOfRoots=resultsReader.nOfRootsChecker(filename,len(setOfTypes))
        roots=resultsReader.goodRootsChecker(filename,nOfRoots,len(setOfTypes))
        cons=resultsReader.roots2concs(roots)
        varOrder=resultsReader.rootOrderFinder(filename)
        consSorted=resultsReader.rootSorter(cons,varOrder)
        
        setOfTypesUpd=[]
        for pair in zip(setOfTypes,consSorted):
            pair[0].concentration=pair[1]
            setOfTypesUpd.append(pair[0])
                    
        route=routes.route
        pickle.dump(consSorted,open("consSorted"+str(maxLength)+kinParStr+".p","wb"))
        pickle.dump(setOfTypesUpd,open(route+'results/'+str(maxLength)+'/'+'setOfTypes'+kinParStr+'.p','wb'))
        print('consSorted'+str(maxLength)+kinParStr+' dumped')
        print('setOfTypes'+kinParStr+' dumped')
    else:
        print('File '+'setOfTypes'+kinParStr+' already exists')
        #consSorted=[CSType.concentration for CSType in setOfTypesUpd]
        #pickle.dump(consSorted,open("consSorted"+str(maxLength)+kinParStr+".p","wb"))
        #print('consSorted'+str(maxLength)+kinParStr+' dumped')
            
    return setOfTypesUpd


###After HOM solver###

def lengthSorter(setOfTypes=0,filename=None):
    """
    list of CatSubType -> dictionary {integer: list of CatSubType}
    given a list of the chain types (setOfTypes) return dictionary, 
    where each key is the length of a chains and entry for this key is the list of types,which have this length
    """
    if setOfTypes==0:
        setOfTypes=setOfTypes=pickle.load(open(filename,"rb"))
    DictOfTypes={}
    DictOfTypes=dict((CSType.Type[0],[]) for CSType in setOfTypes)
    for CSType in setOfTypes:
        DictOfTypes[CSType.Type[0]].append(CSType) 
    
    return DictOfTypes

def DictOfTypesConcSorter(DictOfTypes):
    """
    dictionary {integer: list of CatSubType} -> dictionary {integer: list of CatSubType}
    given dictionary of setOfTypes sorted by lengths, 
    returns a dictionary: {length: list of CatSubType sorted by concentration, descending}
    """
    DictOfTypesSorted=DictOfTypes.copy()
    for key in DictOfTypesSorted.keys():
        val=DictOfTypesSorted[key]
        concentrations=[CSType.concentration for CSType in val]
        DictOfTypesSorted[key]=[CSType for (conc,CSType) in sorted(zip(concentrations,val), reverse=True)]
    
    return DictOfTypesSorted

def maxConcTypeFinder(DictOfTypesSorted):
    """
    dictionary {integer: list} -> dictionary {integer: CatSubType}
    from the dictionary of concentration-wise sorted Types 
    returns dictionary {length: Type with the biggest concentration}
    """
    maxConcTypeDict={}
    for key in DictOfTypesSorted.keys():
            maxConcTypeDict[key]=DictOfTypesSorted[key][0]
                
    return maxConcTypeDict

###Debugging###
#maxLength=11
#kP=[1,1,1,1,1.5]
#setFile2preKinetics(maxLength)
#preKinetics2Workfile(maxLength,kP)
#setOfTypesUpd=setPickle2setWithConc(maxLength,kP)
#filename=routes.route+'results/10/'+'setOfTypes_a1-alpha1-d1-Km1-eH0.9.p'
#DictOfTypes=lengthSorter(setOfTypesUpd)
#DictOfTypesSorted=DictofTypesConcSorter(DictOfTypes)
##for key in DictOfTypesSorted.keys():
    ##for val in DictOfTypesSorted[key]:
        ##print(val)
        ##print(val.concentration)
#maxConcTypeDict=maxConcTypeFinder(DictOfTypesSorted)
#print maxConcTypeDict










