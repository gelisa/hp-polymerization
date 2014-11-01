#! /usr/bin/python2

import os
import pickle

import libExploration
import resultsReader
#import appWholeSet


filename='kineticParams.data'
Nmin=10
Nmax=10

parameters=libExploration.kinParamReader(filename)

libExploration.consentrationsFilesMaker(Nmin,Nmax,parameters)









