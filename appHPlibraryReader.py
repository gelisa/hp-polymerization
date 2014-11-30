#! /usr/bin/python2

import pickle
import routes
import generators

route=routes.route

answer=raw_input('Corrections in libraries?: y/n ')
if (answer=='n'):
    import HPlibraryReader
elif (answer=='y'):
    import HPlibraryReader
    reload(HPlibraryReader)
else:
    print('??')
    import HPlibraryReader
#libraryChoice=libraries/#raw_input('Please choose library: short(s),medium(m),long(l) ')
seqChoice='libraries/'#HPlibraryReader.libraryChoicer(route,libraryChoice)
path=route+seqChoice

catClasses=generators.classesGenerator()
print catClasses
subClasses=generators.classesGenerator()
print subClasses
(nativeList,CSList)=HPlibraryReader.HPlibrary2sequencesList(path,subClasses,catClasses)
setOfTypes=HPlibraryReader.setOfTypesGenerator(CSList)
 
maxLength=HPlibraryReader.maxLengthIdentifier(path)
print('max length is '+str(maxLength))


pickle.dump(nativeList,open("nativeList"+str(maxLength)+".p","wb"))
print("nativeList"+str(maxLength)+".p created")
pickle.dump(CSList,open("CSList"+str(maxLength)+".p","wb"))
print("CSList"+str(maxLength)+".p created")
pickle.dump(setOfTypes,open("setOfTypes"+str(maxLength)+".p","wb"))
print("setOfTypes"+str(maxLength)+".p created")
