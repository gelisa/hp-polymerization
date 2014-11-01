#! /usr/bin/python2

import letter

class CatSubType(object):
    """
    An object which represents a type of sequence which has certain kinetic properties:
    All sequences which have the same length, the same catalytic pattern and the same substrate pattern
    should behave the same according to our kintetic model.
    We treat them as identical to each other.
    Object CatSubType also has an iformation on how many sequences belong to a given type.
    """
    def __init__(self,Type,numberOfSeqs=1):
        '''
        '''
        #Type=(length,substrateType,catalystType)
        self.Type=Type
        #number of sequences belonging to this type
        self.numberOfSeqs=numberOfSeqs
        #generated number for the type
        self.indx=None
        #total concentration of all the sequences in this type
        self.concentration=None
        #what sequences belong to this type
        self.sequencesContent=[]
        
    def __str__(self):
        if self.numberOfSeqs==1:
            return '#'+str(self.indx)+' Type: '+str(self.Type)
        else:
            return '#'+str(self.indx)+' Type: '+str(self.Type)+',\t'+str(self.numberOfSeqs)+': '+str(self.sequencesContent)
    
    def __repr__(self):
        if self.numberOfSeqs==1:
            return '#'+str(self.indx)+' Type: '+str(self.Type)
        else:
            return '#'+str(self.indx)+' Type: '+str(self.Type)+',\t'+str(self.numberOfSeqs)+': '+str(self.sequencesContent)
    
    def __eq__(self,other):
        return self.Type==other.Type
    
    def type2lettersList(self):
        lettersList=[]
        for i in range(len(self.Type[1])):
            lett=letter.LetterInPattern(self.Type[1],i+1)
            lettersList.append(lett)
        
        return lettersList
    
    def patternInType2FormulaHOM(self,alpha, d, Km, eH, popsOfCatsDict):
        lettersList=self.type2lettersList()
        letFormulas=[lett.letterToFormulaWithCatsPopulations(alpha, d, Km, eH, popsOfCatsDict) for lett in lettersList]
        patternInTypeFormulaHOM=reduce(lambda x,y: (x[0]+'*'+y[0],(x[1]+'*'+y[1])),letFormulas)
        
        return patternInTypeFormulaHOM
    
    def type2FormulaHOM(self, a, alpha, d, Km, eH, popsOfCatsDict):
        """
        This function takes the Type and converts it to formula which shows how this Type is being polymerized.
        a, alpha, d, Km, eH  are kinetic parameters
        popsOfCatsDict are from libWholeSet.py
        """
        initialPol=str(a/(2.0*alpha))+'*'
        def simpleStep(alpha, d):
            return ('('+str(alpha)+')','('+str(2*alpha+d)+')')
        
        #formula for polymerization of the substrate pattern
        patternInTypeFormulaHOM=self.patternInType2FormulaHOM(alpha, d, Km, eH, popsOfCatsDict)
        
        #the number of monomers not in the pattern
        simpleStepsNumber=self.Type[0]-len(self.Type[1])
        #if the whole sequence is pattern
        if simpleStepsNumber==0:
            #type formula is pattern formula * a/2alpha
            typeFormulaHOM=(initialPol+patternInTypeFormulaHOM[0],patternInTypeFormulaHOM[1]+'*'+'x'+str(self.indx))
        elif simpleStepsNumber==1:
            leftPart=patternInTypeFormulaHOM[0]+'*'+str(simpleStep(alpha,d)[0])
            rightPart=patternInTypeFormulaHOM[1]+'*'+str(simpleStep(alpha,d)[1]+'*'+'x'+str(self.indx))
            typeFormulaHOM=(initialPol+leftPart, rightPart)
        else:
            leftPart=patternInTypeFormulaHOM[0]+'*'+str(simpleStep(alpha,d)[0])+'**'+str(simpleStepsNumber)
            rightPart=patternInTypeFormulaHOM[1]+'*'+str(simpleStep(alpha,d)[1])+'**'+str(simpleStepsNumber)+'*'+'x'+str(self.indx)
            typeFormulaHOM=(initialPol+leftPart, rightPart)
            
        return typeFormulaHOM
        
    
    
    
    