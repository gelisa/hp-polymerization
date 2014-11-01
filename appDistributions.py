#! /usr/bin/python2

import distributions
#import appResultsReader
import pylab
import pickle
#from appResultsReader import *
import routes

#totalCatCons=distributions.totalCatConsMaker(setOfTypes,consSorted)
#distr, nOfSeqsD=distributions.lengthDistrCatMaker(setOfTypes,totalCatCons)
#normCatDistr=distributions.normCatDistrMaker(setOfTypes,totalCatCons)
#simpDistr=distributions.simpDistrMaker(distr,a, alpha, d)
maxLength=10
kP=[1,1,1,1,1.5]
a=float(kP[0])
alpha=float(kP[1])
d=float(kP[2])
filename=routes.route+'results/10/'+'setOfTypes_a1-alpha1-d1-Km1-eH1.5.p'
setOfTypesUpd=pickle.load(open(filename,"rb"))
lenBestConcTuple=distributions.bestChainsDistrMaker(setOfTypesUpd)
simpDistr=distributions.simpDistrMaker(maxLength,a, alpha, d)

pylab.plot(lenBestConcTuple[0],lenBestConcTuple[1], label='catalyzed polymers')
pylab.plot(simpDistr[0],simpDistr[1], label='non-catalyzed polymers')
pylab.legend()
pylab.show()
