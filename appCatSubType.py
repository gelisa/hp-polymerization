#! /usr/bin/python2

import pickle

answer=raw_input('Corrections in libraries?: y/n ')
if (answer=='n'):
    import HPlibraryReader
    import libWholeSet
elif (answer=='y'):
    import catSubType
    import libWholeSet
    dreload(catSubType)
    dreload(libWholeSet)
else:
    print('??')
    import HPlibraryReader

setOfTypes6=pickle.load(open("setOfTypes6.p","rb"))
catalystsPopulations=libWholeSet.catalystsPopulationsGenerator(setOfTypes6)
patternDict6=libWholeSet.DictOfCatalystsForPatternGenerator(catalystsPopulations)
popsOfCatsDict=libWholeSet.popsOfCatsDictGenerator(patternDict6,setOfTypes6)
result=[CSType.patternInType2FormulaHOM(1.0, 0.5, 1.0, 1.4, popsOfCatsDict) for CSType in setOfTypes6] 
formula=[CSType.type2FormulaHOM(1.0, 0.5, 1.0, 1.4, popsOfCatsDict) for CSType in setOfTypes6] 

