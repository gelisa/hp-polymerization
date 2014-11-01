 #! /usr/bin/python2

import sys
import os
import nativeChain
import catSubType

def HPlibrary2sequencesList(path,subClasses=None,catClasses=None):
    """
    
    dictionary looks like 'URDDL...': [
                                    [HPHPH...., '<native energy>']
                                    .........
                                    ]
    dictEntry is the Value
    
    This function generates a list of sequences (type NativeChain) capable to fold into one native structure
    """
    def HPlibraryReader(path):
        """
        From the data files generates the dictionary, which later is used to gererate
        a list of NativeChain objects
        The files look the following way
        URULLD
            PHPPHPH 2
            HHPPHPH 2
        URRDLD
            HPHPPHP 2
        lines starting with U represent the conformation of the sequence and called 'vector'
        """
        chainDict={}
        for dirname, dirnames, filenames in os.walk(path):
            # print path to all subdirectories first.
            # for subdirname in dirnames:
                #    print os.path.join(dirname, subdirname)
        
                # From the file
            for filename in filenames:
                #full path
                seqFile =  os.path.join(dirname, filename)
                #print seqFile
                #read the file
                inFile = open(seqFile, 'r', 0)
                #for every line in the file
                for line in inFile:
                    #if line starts with the letter U
                    if line[0]=='U':
                        #then this is the vector string
                        #we create a dictionary entry with this vector being a key and entry an empty list
                        chainDict[line.rstrip('\n')]=[]
                        currentConf=line.rstrip('\n')
                    #if no, this is the sequence which has current conformation in the native state
                    else:
                        #in this case we add a list to an entry. 
                        #the first element of the list is hpstring, the second native energy
                        chainDict[currentConf].append(line.rstrip('\n').split(' ')[2:])
                            
        return chainDict
    
    #make temporary dictionary from which hpstrings and conformations can be extracted
    seqDict=HPlibraryReader(path)
    #this will be a list of all the sequences in the library
    nativeList=[]
    #this will be a list only of sequences which are catalysts and substrates simultaneously.
    CSList=[]
    count=0
    #for every key (which is vector 'URDDL...' describing the conformation)
    for key in seqDict:
        #assign a name to an entry, which is a list of lists of the following type [HPHPH...., '<native energy>']
        #several sequences can have the same conformation.
        dictEntry=seqDict[key]
        #then we take each list representing one sequence
        for i in range(len(dictEntry)):
            #assign a temporary name to it
            sequence=nativeChain.NativeChain(dictEntry[i][0])
            #save some obvious information
            sequence.nativeEnergy=int(dictEntry[i][1])
            sequence.vec=key
            sequence.coordinates=sequence.vec2coords()
            sequence.generatedNumber=count
            #check if this sequence is a substrate to one of the possible catalysts
            #some of the catalysts can be actually absent in the system\
            if not subClasses==None:
                sequence.ifSub=sequence.ifSubstrate(subClasses)[0]#!!!changed
                #by what type of catalyst it can be catalyzed?
                sequence.substratePattern=sequence.ifSubstrate(subClasses)[1]#!!!changed
            #if it's a catalyst
            if not catClasses==None:
                sequence.ifCat=sequence.ifCatalyst(catClasses)[0]
                #what kind of catalyst
                sequence.catalystPattern=sequence.ifCatalyst(catClasses)[1]
            #if it's a catalyst
                if sequence.ifCat:
                    #and if it's a substrate
                    if sequence.ifSub:
                        count+=1
                        #then assign a Type for it (see catSubType.py)
                        sequence.Type=(sequence.n,sequence.substratePattern,sequence.catalystPattern)
                        #assign them numbers
                        sequence.generatedNumber=count
                        #append to the list of substrates which are catalysts
                        CSList.append(sequence)
            #append to the list of 
            nativeList.append(sequence)
            
    return (nativeList,CSList)
 
def setOfTypesGenerator(CSList):
    setOfTypes=[]
    count=0
    for chain in CSList:
        typeOfThisChain=catSubType.CatSubType(chain.Type)
        if typeOfThisChain in setOfTypes:
            indx=setOfTypes.index(typeOfThisChain)
            setOfTypes[indx].numberOfSeqs+=1
            setOfTypes[indx].sequencesContent.append(chain.generatedNumber)
        else:
            count+=1
            typeOfThisChain.indx=count
            typeOfThisChain.sequencesContent.append(chain.generatedNumber)
            setOfTypes.append(typeOfThisChain)
        
    return setOfTypes


def libraryChoicer(route,libraryChoice):
    """
    string, string -> string
    given the route to the library and a type of a library returns name of the folder containing the library
    """
    if (libraryChoice=='s' or libraryChoice=='short'):
        seqChoice="HP_designing_test/"
    elif (libraryChoice=='m' or libraryChoice=='medium'):
        seqChoice="HP_designing_medm/"
    elif (libraryChoice=='l' or libraryChoice=='long'):
        seqChoice="HP_designing_long/"
    else:
        print('Please type s, m or l to choose a library of the HP sequences')
        
    return seqChoice

def maxLengthIdentifier(path):
    """
    string, string -> integer
    given a route to the directory containing the library and a choice of the directory returns the maximum length of the
    sequences in the library
    """
    for dirname, dirnames, filenames in os.walk(path):
        sLengths=[int(filename.rstrip('.txt').lstrip('HPn')) for filename in filenames]
        maxLength=max(sLengths)
    
    return maxLength


        
 
 
 
 
 
 
 
 
