import sys
import os

from Config import *
from Chain import *
from Monty import *
from Replica import *
from Trajectory import *

from coordinateReader import *

allSeqs=coordinateReader('./sequences/conf/hp14')

#Finds the coordinates of non-'./sequences/conf/''./sequences/conf/'polar monomers
def hCoordFinder(sequence):
    hList=[]
    sequenceListed=list(sequence)
    #print sequenceListed
    for i in range(len(sequenceListed)):
        if sequenceListed[i]=='H':
            #print 1
            hList.append(allSeqs[sequence][i])
    return hList

#for a given coord tuple returns the list of contacting cells
def primaryNeighborsFinder(coordTuple):
    contactList=[]
    contactList.append((coordTuple[0],coordTuple[1]+1))
    contactList.append((coordTuple[0],coordTuple[1]-1))
    contactList.append((coordTuple[0]+1,coordTuple[1]))
    contactList.append((coordTuple[0]-1,coordTuple[1]))
    return contactList

def allNeighborsFinder(coordTuple):
    neighborsList=primaryNeighborsFinder(coordTuple)
    neighborsList.append((coordTuple[0]+1,coordTuple[1]+1))
    neighborsList.append((coordTuple[0]+1,coordTuple[1]-1))
    neighborsList.append((coordTuple[0]-1,coordTuple[1]+1))
    neighborsList.append((coordTuple[0]-1,coordTuple[1]-1))
    return neighborsList

#Finds coordinates of the free cells contacting hydrophobes            
def hContactFinder(sequence):
    hList=hCoordFinder(sequence)
    hContactList=[]
    for item in hList:
        contactList=primaryNeighborsFinder(item)
        for pair in contactList:
            if not pair in allSeqs[sequence]:
                if not pair in hContactList:
                    hContactList.append(pair)
    return hContactList
        
#Finds coordinates of the all cells, which are in primary contact with the given sequence
def primaryContactFinder(sequence):
    coordinates=allSeqs[sequence]
    primaryContacts=[]
    for item in coordinates:
        contactList=primaryNeighborsFinder(item)
        for pair in contactList:
            if not pair in allSeqs[sequence]:
                if not pair in primaryContacts:
                    primaryContacts.append(pair)
    return primaryContacts

def allContactFinder(sequence):
    coordinates=allSeqs[sequence]
    allContacts=[]
    for item in coordinates:
        contactList=allNeighborsFinder(item)
        for pair in contactList:
            if not pair in allSeqs[sequence]:
                if not pair in allContacts:
                    allContacts.append(pair)
    return allContacts
    
def stepUp(coordinate):
    return (coordinate[0],coordinate[1]+1)
def stepDown(coordinate):
    return (coordinate[0],coordinate[1]-1)
def stepRight(coordinate):
    return (coordinate[0]+1,coordinate[1])
def stepLeft(coordinate):
    return (coordinate[0]-1,coordinate[1])
def stepUpRight(coordinate):
    return (coordinate[0]+1,coordinate[1]+1)
def stepUpLeft(coordinate):
    return (coordinate[0]-1,coordinate[1]+1)
def stepDownRight(coordinate):
    return (coordinate[0]+1,coordinate[1]-1)
def stepDownLeft(coordinate):
    return (coordinate[0]-1,coordinate[1]-1)


    
    
#Returns the list of coordinates of the surface of a given sequence
def surfaceConnector(sequence):
    def allowedStepsFinder(coord):
        return [stepDown(coord),stepLeft(coord),stepUp(coord),stepRight(coord)]
    
    allContacts=allContactFinder(sequence)
    surfaceCoordinates=[]
    surfaceCoordinates.append(allContacts[0])
    currentStep=allContacts[0]
    nextStep=0
    count=0
    while count<=len(allContacts):
        count+=1
        allowedSteps=allowedStepsFinder(currentStep)
        for step in allowedSteps:
            if (step in allContacts) and (not step in surfaceCoordinates):
                surfaceCoordinates.append(step)
                currentStep=step
                #print currentStep
                break
            
    return surfaceCoordinates
    

#From the coordinates of the hydrophobes and all surface coordinates calculates surface sequence
def surfaceSequenceMaker(sequence):
    surfaceCoordinates=surfaceConnector(sequence)
    hContactList=hContactFinder(sequence)
    surfaceSequence=''
    for i in range(len(surfaceCoordinates)):
        if surfaceCoordinates[i] in hContactList:
            surfaceSequence=surfaceSequence+'H'
        else:
            surfaceSequence=surfaceSequence+'P'
    return surfaceSequence

#Makes a dictionary: {sequence: {surface sequence: surface coordinates}}
def surfaceDictMaker():
    allSurfaces={}
    for seq in allSeqs:
        surfaceCoordinates=surfaceConnector(seq)
        surfaceSequence=surfaceSequenceMaker(seq)
        innerDict={}
        innerDict[surfaceSequence]=surfaceCoordinates
        allSurfaces[seq]= innerDict
    return allSurfaces
    
def surfaceSmallDictMaker():
    allSurfacesNoCoord={}
    for seq in allSeqs:
        surfaceCoordinates=surfaceConnector(seq)
        surfaceSequence=surfaceSequenceMaker(seq)
        allSurfacesNoCoord[seq]= surfaceSequence
    return allSurfacesNoCoord


def catalyzersFinder(pattern):
    catalyzersDict={}
    allSurfacesNoCoord=surfaceSmallDictMaker()
    for sequence in allSurfacesNoCoord:
        if (allSurfacesNoCoord[sequence]+allSurfacesNoCoord[sequence]).find(pattern)>-1:
            catalyzersDict[sequence]=allSurfacesNoCoord[sequence]
    return catalyzersDict


            








