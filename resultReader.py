#! /usr/bin/python2
"""
File has the following structure (example for 2 vars, 2 roots):
answer for var1 root1
answer for var2 root1
 
 residue          =    ...
 condition number =    ...
 ---------------------------------------
 answer for var1 root2
answer for var2 root2
 
 residue          =    ...
 condition number =    ...
 ---------------------------------------
  The order of variables : 
 x1
 x2
 """
def nOfRootsChecker(filename,nOfVars):
    """
    string, number -> number
    calculates total number of root produced by the HOM solver
    """
    #open the file
    inFile = open(filename, 'r', 0)
    #assign number of lines variables
    lCount=0
    #for every line
    for line in inFile:
        #increase "number of lines" by 1
        lCount+=1
        #if encountered 'The order of variables'
        if line.find('The order of variables')>-1:
            #stop
            break
    #number of roots is the whole part of "number of lines" divided by "number of vars"+4
    nOfRoots=(lCount)/(nOfVars+4)
    
    return nOfRoots

def goodRootsChecker(filename,nOfRoots,nOfVars):
    """
    string,number,number -> list of lists of strings
    returns the root for which all vars are real and positive
    """
    inFile = open(filename, 'r', 0)
    lines=inFile.readlines()
    #initiale roots list
    roots=[]
    #for every number from 0 to the "number of roots"-1
    for i in range(nOfRoots):
        #root is a list of lines, which correspond to the answers for variables in the file
        #find the root
        root=lines[i*(nOfVars+4):(i+1)*(nOfVars+4)-4]
        #set quality of a root to "good"=1
        mark=1
        #every line in the root corresponds to an answer for a variable
        #for every variable
        for Line in root:
            #if the answer is negative
            if Line[0:3]=='( -':
                #the root is "bad"
                mark=0
                #stop checking this root
                break
            #if the line has an imaginary part
            elif not Line.find('0.0000000000000000E+00'):
                #thin is a bad root
                mark=0
                #stop checking it
                break
        #if after checking all the lines of root, the root is still good
        if mark==1:
            #write this root as a list of strings
            posRoot=[Line.lstrip('( ').split(' , ')[0] for Line in root]
            #append to the list of all roots
            roots.append(posRoot)
    
    return roots

def roots2concs(roots):
    """
    list of lists of strings -> list of floats
    """
    #initiate list of concentrations
    concs=[]
    
    #if there are more then 1 positive root for the solution
    if len(roots)>1:
        #I don't know what to do with it
        print('Warning! There are 2 or more positive concentrations!')
        #don't do anything
        return None
    #if there are no roots
    elif len(roots)==0:
        print('No roots!')
        #don't do anything
        return None
    #otherwise there is only one root
    else:
        #extract it
        Roots=roots[0]
        #for every variable in this root
        for root in Roots:
            #extract the base of the exponent
            numb=float(root.split('E')[0])
            #extract power of the exponent
            expt=root.split('E')[1]
            #convert base and exponent notation into a float-number. This is concentration.
            conc=numb*10.0**(int(expt))
            #append to the list of the concentrations
            concs.append(conc)
    
    return concs

def rootOrderFinder(filename):
    """
    string -> list
    returns a list of the indexes of variables in order in which HOM treats them
    """
    inFile = open(filename, 'r', 0)
    #for every line
    varOrder=[]
    for line in inFile:
        #if encountered 'The order of variables'
        if line.find('The order of variables')>-1:
            #stop
            break
    for line in inFile:
        if line[0:2]==' x':
            varOrder.append(int(line.lstrip(' x')))
        else:
            break
    
    return varOrder

def rootSorter(concs,varOrder):
    """
    list, list -> list
    returns concentrations of spices ordered, so that their indexes ascending
    """
    consSorted=[x for (y,x) in sorted(zip(varOrder,concs))]
        
    return consSorted





