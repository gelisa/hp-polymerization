#! /usr/bin/python2

import sys
import coordinate

sys.path.append("/data/research/06.origins_of_life/codes/hp-model")

def HP2Bin(HPstring):
    """Convert a string of type 'HPHPHPPPHHP' to a string of 1s and 0s.""" 
    binseq =''
    for i in range(0,len(HPstring)):
            if HPstring[i] == 'H':
                binseq = binseq+'1'    
            else:
                binseq = binseq+'0'     
    return binseq


class NativeChain(object):
    """
    An object to represent the 2D HP lattice chain's sequence, conformation and free energy
    """ 
    def __init__(self,hpstring):
        self.hpstring = hpstring.upper()  # The HP sequence as a string
        self.n = len(self.hpstring)                 # the chain length
        self.vec = ''                      # an (n-1)-dimensional vector representation of the chain
        self.coordinates=[]
        self.nativeEnergy=None
        self.generatedNumber=None
        self.ifSub=None
        self.substratePattern=None
        self.ifCat=None
        self.catalystPattern=None
        self.Type=None
        self.hpbinary = self.hpstr2bin()                # The HP seq in binary rep (H=1 P=0)
        
        # The chain is represented by an (n-1)-dimensional vector, which is stored in a python list.
        # Each element represents the direction of successive links (bonds) along the chain length:
        #
        #     U     up
        #     R     right
        #     D     down
        #     L     left
        #
        # By convention, the first monomer (bead) in the chain is fixed at the origin on
        # a two-dimensional square lattice.

        
        
    def __repr__(self):
        if self.vec==[]:
            return "#%d, n=%d, %s: vector isn't determined" % (self.generatedNumber, self.n, self.hpstring)
        else:
            line1="#%d, n=%d, %s: %s" % (self.generatedNumber, self.n, self.hpstring, self.vec)
            line2="\n Sub: "+str(self.substratePattern)+" Cat: "+str(self.catalystPattern)
            return line1+line2
    
    def __str__(self):
        if self.vec==[]:
            return "#%d, n=%d, %s: vector isn't determined" % (self.generatedNumber, self.n, self.hpstring)
        else:
            line1="#%d, n=%d, %s: %s" % (self.generatedNumber, self.n, self.hpstring, self.vec)
            line2="\n Sub: "+str(self.substratePattern)+" Cat: "+str(self.catalystPattern)
            return line1+line2
          
###From here changed to Class attributes    
    #Comments
    def vec2coords(self):
        """Convert a list of chain vectors to a list of coordinates (duples).""" 
        tmp = [coordinate.Coordinate((0,0))]
        x = 0
        y = 0
        if self.vec==[]:
            print("no vector is given!")
        else:
            for i in range(0,len(self.vec)):
                if self.vec[i] == 'U':
                    y = y + 1
                if self.vec[i] == 'R':
                    x = x + 1
                if self.vec[i] == 'D':
                    y = y - 1
                if self.vec[i] == 'L':
                    x = x - 1
                tmp.append(coordinate.Coordinate((x,y)))
        return tmp



    ###
###Comments, debugging

    def allContactFinder(self):
        allContacts=[]
        for item in self.coordinates:
            contactList=item.allNeighborsFinder()
            #print contactList
            for pair in contactList:
                if not (pair in self.coordinates):
                    #print(str(pair)+' '+str(self.coordinates))
                    if not pair in allContacts:
                        allContacts.append(pair)
        return allContacts
    
    def hpstr2bin(self):
        """Convert a string of type 'HPHPHPPPHHP' to a list of 1s and 0s.""" 
        binseq =''
        for i in range(0,len(self.hpstring)):
             if self.hpstring[i] == 'H':
                binseq = binseq+'1'       
             else:
                binseq = binseq+'0' 
        return binseq
    
###Comments, debugging
    #Returns the list of coordinates of the surface of a given sequence
    
    def surfaceConnector(self):
        #First, looking for the all points contacting a given sequence
        allContacts=self.allContactFinder()
        #Here will be a list of surface coordinates in the order, in which one can walk the whole surface around with only allowedSteps
        surfaceCoordinates=[]
        #Starting from the first point in the list
        surfaceCoordinates.append(allContacts[0])
        currentStep=allContacts[0]
        nextStep=None
        count=0
        #While uncounted contacts present
        while count<=len(allContacts):
            #count counted contacts
            count+=1
            allowedSteps=currentStep.allowedStepsFinder()
            #Looking for all possible steps from the given point
            for step in allowedSteps:
                # if the step in the surface and not mentioned yet
                if (step in allContacts) and (not step in surfaceCoordinates):
                    #the first one is fine
                    surfaceCoordinates.append(step)
                    #moving to this coordinate
                    currentStep=step
                    #print currentStep
                    break
                
        return surfaceCoordinates

    #Finds the coordinates of non-polar monomers
    def hCoordFinder(self):
        hList=[]
        sequenceListed=list(self.hpstring)
        #Check every letter in the sequence
        for i in range(len(sequenceListed)):
            if sequenceListed[i]=='H':
                #append a coordinate of another H
                hList.append(self.coordinates[i])
        return hList

    #Finds coordinates of the free cells contacting hydrophobes            
    def hContactFinder(self):
        #hList is the list of all H's in the sequence
        hList=self.hCoordFinder()
        hContactList=[]
        #For every H in sequence
        for item in hList:
            #Find primary contacts of this H. This is a list
            HContacts=item.primaryNeighborsFinder()
            #For every pair in the contact list
            for pair in HContacts:
                #If this pair doesn't belong to the sequence
                if not pair in self.coordinates:
                    #and if it's not in the list of contacts of all H's
                    if not pair in hContactList:
                        #Then add it to the list
                        hContactList.append(pair)
        return hContactList
###Comments, debugging
    #From the coordinates of the hydrophobes and all surface coordinates calculates surface sequence
    def surfaceSequenceMaker(self):
        #Get the list of coordinates of the surface of a given sequence
        surfaceCoordinates=self.surfaceConnector()
        #Find H-contacting coordinates
        hContactList=self.hContactFinder()
        #Initialize the surface sequence
        surfaceSequence=''
        #for every coordinate in the surface coordinate list
        for i in range(len(surfaceCoordinates)):
            #if the coordinate belongs to H-contacting coordinates
            if surfaceCoordinates[i] in hContactList:
                #Then add 'H' to the sequence
                surfaceSequence=surfaceSequence+'H'
            # and if the coordinate doesn't contact 'H'
            else:
                #add 'P'
                surfaceSequence=surfaceSequence+'P'
                
        return surfaceSequence

    #Comments, debugging
    def ifCatalyst(self,catClassesHnP):#!!!Change all dependencies!!!
        surfaceSeq=self.surfaceSequenceMaker()
        for pattern in catClassesHnP:
            if (surfaceSeq+surfaceSeq).find(pattern)>-1:
                answer = (True, pattern)
                break
            else: 
                answer = (False,'')
        return answer
    
    #def ifCatalystBin(self,catClasses):#!!!Change all dependencies!!!
        #surfaceSeq=HP2Bin(self.surfaceSequenceMaker())
        #print surfaceSeq
        #for pattern in catClasses:
            #if (surfaceSeq+surfaceSeq).find(pattern)>-1:
                #answer = (True, pattern)
                #print answer
                #break
            #else: 
                #answer = (False,'')
        #return answer

    #def nativeToCatalyst(self):
        #if self.ifCatalyzer()[0]:
            #ct=Catalyzer(self)
            #ct.catalyzerClass=self.ifCatalyzer()[1]
        #return ct

    def ifSubstrate(self,subClasses):#TEST BUG?
        for pattern in subClasses:
            if (self.hpstring).find(pattern)>-1:
                answer = (True, pattern)
                break
            else: 
                answer = (False,'')
        return answer




#   







