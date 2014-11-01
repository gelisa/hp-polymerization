#! /usr/bin/python2

import pickle

answer=raw_input('Corrections in libraries?: y/n ')
if (answer=='n'):
    import catSubType
    import libWholeSet
elif (answer=='y'):
    import catSubType
    import libWholeSet
    import HPlibraryReader
    reload(catSubType)
    reload(libWholeSet)
else:
    print('??')
    import catSubType
    import libWholeSet
Nmer=raw_input('What is the length of the longest chain in the set? ')

setOfTypes=pickle.load(open("setOfTypes"+str(N)+".p","rb"))
print('setOfTypes done')
catalystsPopulations=libWholeSet.catalystsPopulationsGenerator(setOfTypes)
print('catalystsPopulations done')
patternDict=libWholeSet.DictOfCatalystsForPatternGenerator(catalystsPopulations)
print('patternDict done')
popsOfCatsDict=libWholeSet.popsOfCatsDictGenerator(patternDict,setOfTypes)
print('popsOfCatsDict done')
setOfVars=libWholeSet.setOfVarsGenerator(setOfTypes)
print('setOfVars done')

pickle.dump(popsOfCatsDict,open("popsOfCatsDict"+str(Nmer)+".p","wb"))
pickle.dump(setOfVars,open("setOfVars"+str(Nmer)+".p","wb"))
print("popsOfCatsDict"+str(Nmer)+" pickled")
print("setOfVars"+str(Nmer)+" pickled")

