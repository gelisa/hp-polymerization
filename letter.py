#! /usr/bin/python2

#from sympy import *
import math
import libWholeSet

class LetterInPattern(object):
    def __init__(self,pattern,position):
        """
        position runs from 1.
        pattern is the sequence that specifies
        """
        self.pattern=pattern
        self.position=position
        self.letter=pattern[position-1]
        #what if the letter is the first letter in the pattern
        if self.position==1:
            self.catUp=None
            self.etaUp=None
            self.catDown=None
            self.etaDown=None
        #what if the letter is not the first
        else:
            #what if this letter is H
            if self.letter=='H':
                #what if it's a secondary H, meaning that there is a P in front of it
                if self.pattern[0:self.position].find('P')>-1:
                    #what if it's the last H
                    if self.position==len(self.pattern):
                        self.catUp=self.pattern
                        self.etaUp=self.pattern.count('H')
                        self.catDown=None
                        self.etaDown=None
                    #what if it's not the last letter
                    else:
                        self.catUp=self.pattern[0:position]
                        self.etaUp=self.pattern[0:position].count('H')
                        self.catDown=self.pattern[0:position+1]
                        self.etaDown=self.pattern[0:position+1].count('H')
                #what if there is no P before current H, it means it's primary H
                else:
                    #what if the next letter is P
                    if self.pattern[self.position]=='P':
                        self.catUp=self.pattern[0:position]+'PP'
                        self.etaUp=self.pattern[0:position].count('H')
                        self.catDown=None
                        self.etaDown=None
                    #what if the next letter is H
                    else:
                        self.catUp=self.pattern[0:position]+'PP'
                        self.etaUp=self.pattern[0:position].count('H')
                        self.catDown=self.pattern[0:position+1]+'PP'
                        self.etaDown=self.pattern[0:position+1].count('H')
            #what if the letter is P
            else:
                #what if it's the last letter or the next letter is P. 
                #It means that this P is in the PP tail
                if self.position==len(self.pattern) or self.pattern[self.position]=='P':
                    self.catUp=None
                    self.etaUp=None
                    self.catDown=None
                    self.etaDown=None
                #what if this P isn't in the PP tail
                else:
                    self.catUp=None
                    self.etaUp=None
                    self.catDown=self.pattern[0:position+1]
                    self.etaDown=self.pattern[0:position+1].count('H')
                    
    def __repr__(self):
        if self.position==1:
            line1 = self.letter+self.pattern[1:].lower()
        elif self.position==len(self.pattern):
            line1 = self.pattern[0:self.position-1].lower()+self.letter
        else:
            line1 = self.pattern[0:self.position-1].lower()+self.letter+self.pattern[self.position:].lower()
        line2 = '\n'+'catUp:'+str(self.catUp)+'\t catDown: '+str(self.catDown)
        line3 = '\n'+'etaUp:'+str(self.etaUp)+'\t etaDown: '+str(self.etaDown)
        return line1+line2+line3
    
    def __str__(self):
        if self.position==1:
            line1 = self.letter+self.pattern[1:].lower()
        elif self.position==len(self.pattern):
            line1 = self.pattern[0:self.position-1].lower()+self.letter
        else:
            line1 = self.pattern[0:self.position-1].lower()+self.letter+self.pattern[self.position:].lower()
        line2 = '\n'+'catUp:'+str(self.catUp)+'\t catDown: '+str(self.catDown)
        line3 = '\n'+'etaUp:'+str(self.etaUp)+'\t etaDown: '+str(self.etaDown)
        return line1+line2+line3
    
    def letterEta2RealEta(self, alpha, Km,eH):
        """
        Eta in formulas above is the amount of "active" monomers for the current letter.
        This function converts old Eta into real Eta, which has kinetic meaning.
        """
        #if etaDown is present
        if not self.etaDown==None:
            #convert it to formula
            realEtaDown=alpha/Km*math.exp(self.etaDown*eH)
        else:
            realEtaDown=None
        #etc
        if not self.etaUp==None:
            realEtaUp=alpha/Km*math.exp(self.etaUp*eH)
        else:
            realEtaUp=None
        return (realEtaUp,realEtaDown)
    
    def letterEta2HOMEta(self, alpha, Km,eH):
        """
        same as letterEta2RealEta, but returns a tuple of strings
        """
        if not self.etaDown==None:
            #convert it to formula
            realEtaDown=str(alpha/Km*math.exp(self.etaDown*eH))
        else:
            realEtaDown=None
        #etc
        if not self.etaUp==None:
            realEtaUp=str(alpha/Km*math.exp(self.etaUp*eH))
        else:
            realEtaUp=None
        return (realEtaUp,realEtaDown)
    
    def letterToFormula(self,alpha, d, Km, eH, catUp, catDown,  test=True):
        """
        This function return a tuple. 
        The first element is the nominator of the steady state formula for this letter.
        The second is the denominator.
        """
        
        def simpleStep(alpha, d):
            return (alpha,(2*alpha+d))
            
        def HHstep(alpha, d, etaUp, catUp, etaDown, catDown):
            return ((alpha+etaUp*catUp),(2*alpha+d+etaDown*catDown))
            
        def HPstep(alpha, d, etaUp, catUp):
            return ((alpha+etaUp*catUp),(2*alpha+d))
            
        def PHstep(alpha, d, etaDown, catDown):
            return ((alpha),(2*alpha+d+etaDown*catDown))
        
        if test==False:
            if self.catDown==None:
                if self.catUp==None:
                    formula=simpleStep(alpha, d)
                else:
                    formula=HPstep(alpha, d, self.letterEta2RealEta(alpha, Km,eH)[0], catUp)
            else:
                if self.catUp==None:
                    formula=PHstep(alpha, d, self.letterEta2RealEta(alpha, Km,eH)[1], catDown)
                else:
                    formula=HHstep(alpha, d, self.letterEta2RealEta(alpha, Km,eH)[0], catUp, self.letterEta2RealEta(alpha, Km,eH)[1], catDown)
                    
        else:
            if self.catDown==None:
                if self.catUp==None:
                    formula=simpleStep(Symbol('alpha'), Symbol('d'))
                else:
                    formula=HPstep(Symbol('alpha'), Symbol('d'), Symbol('etaUp'), Symbol('catUp'))
            else:
                if self.catUp==None:
                    formula=PHstep(Symbol('alpha'), Symbol('d'), Symbol('etaDown'), Symbol('self.catDown'))
                else:
                    formula=HHstep(Symbol('alpha'), Symbol('d'), Symbol('etaUp'), Symbol('catUp'), Symbol('etaDown'), Symbol('self.catDown'))
        return formula
    
    def letterToHOMFormula(self,alpha, d, Km, eH, catUp, catDown):
        """
        Same as previous, but returns a tuple of strings
        """
        def simpleStep(alpha, d):
            return ('('+str(alpha)+')','('+str(2*alpha+d)+')')
            
        def HHstep(alpha, d, etaUp, catUp, etaDown, catDown):
            return ('('+str(alpha)+'+'+str(etaUp)+'*('+str(catUp)+'))','('+str(2*alpha+d)+'+'+str(etaDown)+'*('+str(catDown)+'))')
            
        def HPstep(alpha, d, etaUp, catUp):
            return ('('+str(alpha)+'+'+str(etaUp)+'*('+str(catUp)+'))','('+str(2*alpha+d)+')')
            
        def PHstep(alpha, d, etaDown, catDown):
            return ('('+str(alpha)+')','('+str(2*alpha+d)+'+'+str(etaDown)+'*('+str(catDown)+'))')
        
        
        if self.catDown==None:
            if self.catUp==None:
                formula=simpleStep(alpha, d)
            else:
                formula=HPstep(alpha, d, self.letterEta2HOMEta(alpha, Km,eH)[0], catUp)
        else:
            if self.catUp==None:
                formula=PHstep(alpha, d, self.letterEta2HOMEta(alpha, Km,eH)[1], catDown)
            else:
                formula=HHstep(alpha, d, self.letterEta2HOMEta(alpha, Km,eH)[0], catUp, self.letterEta2HOMEta(alpha, Km,eH)[1], catDown)
                
        return formula

    def letterToFormulaWithCatsPopulations(self,alpha, d, Km, eH, cat2PopulationsDict):
        if self.catUp in cat2PopulationsDict.keys():
            catUp=cat2PopulationsDict[self.catUp]
        else:
            catUp='0'
        if self.catDown in cat2PopulationsDict.keys():
            catDown=cat2PopulationsDict[self.catDown]
        else:
            catDown='0'
        formulaWithPops=self.letterToHOMFormula(alpha, d, Km, eH, catUp, catDown)
        
        return formulaWithPops







