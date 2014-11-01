#! /usr/bin/python2

import letter
from sympy import *

class SequenceType(object):
    def __init__(self, length):
        """
        numberOfSeqs -- number of seqs which have the same pattern and the same length
        """
        self.length=length
        
class SubstrateType(SequenceType):
    def __init__(self, length,Spattern):
        self.length=length
        self.Spattern=Spattern
        self.SType=(self.length,self.Spattern)
        self.number=None
        #number of sequences which belong to this type
        self.quantity=None
        self.concentration=None
        
        #primaryH=0
        #for i in range(len(self.substrateType)):
            #if not self.substrateType[i]=='P':
                #primaryH+=1
                #Pnumber=self.substrateType.count('P')
                #secondaryH=self.substrateType[i+1:].count('H')
            #else:
                #break
        #self.substrateDecomposition=(primaryH, Pnumber, secondaryH)
        
    def __str__(self):
        return '#'+str(self.number)+'\t'+str(self.SType)+': '+str(self.quantity)
    
    def __repr__(self):
        return '#'+str(self.number)+'\t'+str(self.SType)+': '+str(self.quantity)
    
    def __eq__(self,other):
        return self.SType==other.SType
        
    def Type2Formula(self,alpha, d,a,Km,eH,Test=True, patternSum=None):
        def pattern2Formula(pattern,alpha, d,Test,patternSum):
            formula=1.0
            for i in range(1,len(pattern)+1):
                lett=letter.LetterInPattern(pattern,i)
                formula=formula*lett.letterToFormula(alpha, d,Km,eH,Test,patternSum)
            
            return formula
        
        def simpleStep(alpha, d,Test=True):
            if Test==True:
                return Symbol('alpha')/(2*Symbol('alpha')+Symbol('d'))
            else:
                return (alpha)/(2*alpha+d)
        
        pattern=self.Spattern
        formula=pattern2Formula(pattern,alpha, d,Test,patternSum)
        theRest=self.length-len(self.Spattern)
        formula=formula*((simpleStep(alpha, d,Test))**(theRest))*a/2.0/alpha
        
        return formula
        
class CatalystType(SequenceType):
    def __init__(self, Cpattern):
        self.CType=Cpattern
        self.substrateTypeDecomposition=[]
        self.quantity=None
            
        #primaryH=0
        #for i in range(len(self.substrateType)):
            #if not self.substrateType[i]=='P':
                #primaryH+=1
                #Pnumber=self.substrateTypeself.substrateType.count('P')
                #secondaryH=self.substrateType[i+1:].count('H')
            #else:
                #break
        #self.catalystDecomposition=(primaryH, Pnumber, secondaryH)
        
    def __str__(self):
        return str(self.CType)+': '+str(self.substrateTypeDecomposition)
    
    def __repr__(self):
        return str(self.CType)+': '+str(self.substrateTypeDecomposition)
    
    def __eq__(self,other):
        return self.CType==other.CType
  

    







