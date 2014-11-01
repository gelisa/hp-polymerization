#! /usr/bin/python2
import pickle

answer=raw_input('Corrections in libraries?: y/n ')
if (answer=='n'):
    import resultsReader
elif (answer=='y'):
    import resultsReader
    dreload(resultsReader)
else:
    print('??')
    import resultsReader
    
setOfTypes=pickle.load(open("setOfTypes10.p","rb"))
filename='/data/research/06.origins of life/codes/HOM/HOM4PS2/data.roots'
nOfRoots=resultsReader.nOfRootsChecker(filename,len(setOfTypes))
roots=resultsReader.goodRootsChecker(filename,nOfRoots,len(setOfTypes))
cons=resultsReader.roots2concs(roots)
varOrder=resultsReader.rootOrderFinder(filename)
consSorted=resultsReader.rootSorter(cons,varOrder)
print consSorted

