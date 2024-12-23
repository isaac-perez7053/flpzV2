#!/usr/bin/env python3
import numpy as np
import subprocess
import xml.etree.ElementTree as ET
import sys
np.set_printoptions(precision=10)
#takes as input an smodes input file and irrep name and creates the distortions needed to do a symmetry adapted modes calculation
#we expect a POSCAR file in this directory to get the header from


#old version set this to your acell, as if you had normalized rprim. Ignore for now
#precLatParam=np.array([7.3477206390E+00,  7.3477206390E+00,  7.3477206390E+00]) #this is in Bohr
#######################################################################################################


#all the atom labels and their atomic masses
atomList=[('1',   'H',    '1.008'), ('2',   'He'   ,'4.002'),('3',  'Li',   '6.94'),('4', 'Be', ' 9.012'),('5', 'B', ' 10.81'),('6', 'C', ' 12.011'), ('7', 'N', ' 14.007'),('8', 'O', ' 15.999'),('9', 'F', ' 18.998'),('10','Ne', ' 20.180'),('11','Na', ' 22.990'),('12','Mg', ' 24.305'),('13','Al', ' 26.982'),('14','Si', ' 28.085'),('15','P', ' 30.974'),('16','S', ' 32.06'),('17','Cl', ' 35.45'),('18','Ar', ' 39.948'),('19','K', ' 39.098'),('20','Ca', ' 40.078'),('21','Sc', ' 44.956'),('22','Ti', ' 47.867'),('23','V', ' 50.942'),('24','Cr', ' 51.996'),('25','Mn', ' 54.938'),('26','Fe', ' 55.845'),('27','Co', ' 58.933'),('28','Ni', ' 58.693'),('29','Cu', ' 63.546'),('30','Zn', ' 65.38'),('31','Ga', ' 69.723'),('32','Ge', ' 72.630'),('33','As', ' 74.922'),('34','Se', ' 78.971'),('35','Br', ' 79.904'),('36','Kr', ' 83.798'),('37','Rb', ' 85.468'),('38','Sr', ' 87.62'),('39','Y', ' 88.906'),('40','Zr', ' 91.224'),('41','Nb', ' 92.906'),('42','Mo', ' 95.95'),('43','Tc', ' 98'),('44','Ru', '101.07'),('45','Rh', '102.91'),('46','Pd', '106.42'),('47','Ag', '107.87'),('48','Cd', '112.41'),('49','In', '114.82'),('50','Sn', '118.71'),('51','Sb', '121.76'),('52','Te', '127.60'),('53','I ', '126.90'),('54','Xe', '131.29'),('55','Cs', '132.91'),('56','Ba', '137.33'),('57','La', '138.91'),('58','Ce', '140.12'),('59','Pr', '140.91'),('60','Nd', '144.24'),('61','Pm', '145'),('62','Sm', '150.36'),('63','Eu', '151.96'),('64','Gd', '157.25'),('65','Tb', '158.93'),('66','Dy', '162.50'),('67','Ho', '164.93'),('68','Er', '167.26'),('69','Tm', '168.93'),('70','Yb', '173.05'),('71','Lu', '174.97'),('72','Hf', '178.49'),('73','Ta', '180.95'),('74','W', '183.84'),('75','Re', '186.21'),('76','Os', '190.23'),('77','Ir', '192.22'),('78','Pt', '195.08'),('79','Au', '196.97'),('80','Hg', '200.59'),('81','Tl', '204.38'),('82','Pb', '207.2'),('83','Bi', '208.98'),('84','Po', '209'),('85','At', '210'),('86','Rn', '222'),('87','Fr', '223'),('88','Ra', '226'),('89','Ac', '227'),('90','Th', '232.04'),('91','Pa', '231.04'),('92','U ', '238.03'),('93','Np', '237'),('94','Pu', '244'),('95','Am', '243'),('96','Cm', '247'),('97','Bk', '247'),('98','Cf', '251'),('99','Es', '252'),('100','Fm', '257'),('101','Md', '258'),('102','No', '259'),('103','Lr', '266'),('104','Rf', '267'),('105','Db', '268'),('106','Sg', '269'),('107','Bh', '270'),('108','Hs', '277')]

dispMag=0.001 #In angstroms!!
symmPrec=1e-5 #this is the best smodes can do
smodesFile=sys.argv[1]
targetIrrep=sys.argv[2]

#get acell from smodes file
s=open(smodesFile)
sLines=s.readlines()
s.close()
sWords=sLines[0].split()
precLatParam = [float(sWords[0]), float(sWords[1]), float(sWords[2])]
print("!!!!!")
print("REMEMBER TO SET ACELL IN 1st LINE OF SMODES FILE, right now is:")
print(precLatParam)
print("!!!!!")


smodes_path = '/home/iperez/isobyu/smodes'  # Replace with the actual path you found
args = f'{smodes_path} <{smodesFile}'
proc = subprocess.Popen(args, shell=True, stdout=subprocess.PIPE)
output = proc.stdout.read()
parsedOutput=output.decode('ascii')

targetOutput=[]
startTarget=999
endTarget=0
#endPositionis=0
outList=parsedOutput.split("\n")
for l in range(len(outList)):
    thisLine=outList[l].split()
    #print([l,thisLine])
    if len(thisLine)>1:
        if (    thisLine[0]=="Irrep") and (thisLine[1]==targetIrrep):
            startTarget=l
        
    if (len(thisLine)>0) and (startTarget<999):
        if thisLine[0]=="***********************************************":
            endTarget=l
            break
targetOutput=outList[startTarget:endTarget]
#check if it's raman or IR active, because we'll need to remove that line
if (targetOutput[3].split()[0]=='These'):
    print('Includes translational modes')
    del targetOutput[3] 

if (targetOutput[3].split()[0]=='IR'):
    print('IR Active')
    del targetOutput[3] 
if (targetOutput[3].split()[0]=='Raman'):
    print('Raman Active')
    del targetOutput[3] 



degeneracy=int(targetOutput[1].split()[-1])
print("Degeneracy: ", degeneracy)

numModes_withoutDegen=int(targetOutput[2].split()[-1])
print("Number of Modes: ", numModes_withoutDegen)

numModes=int(numModes_withoutDegen/degeneracy)

print("(Meaning ", numModes,"modes to find)")

v1=[float(i) for i in targetOutput[4].split()]
v2=[float(i) for i in targetOutput[5].split()]
v3=[float(i) for i in targetOutput[6].split()]
shapeCell=np.matrix([v1,v2,v3])
atom_list=[]
atom_pos_raw=[]

for l in range(8,len(targetOutput)):
    thisLine=targetOutput[l].split()
    if (targetOutput[l]=="Symmetry modes:"):
        break
    if(len(thisLine)>=4):
        atom_list.append(thisLine[1])
#atom_pos.append([float(thisLine[2])*precLatParam[0],float(thisLine[3])*precLatParam[1],float(thisLine[4])*precLatParam[2]])
        atom_pos_raw.append([float(thisLine[2]),float(thisLine[3]),float(thisLine[4])])

numAtoms=len(atom_list)
print("Number of Atoms: ",numAtoms)


#find number of atom types and count them
typeList=[]
for a in range(numAtoms):
    if (atom_list[a] in typeList) == False:
        typeList.append(atom_list[a])

typeCount=np.zeros((len(typeList)))
for a in range(numAtoms):
    for type in range(len(typeList)):
        if (atom_list[a] == typeList[type]):
            typeCount[type]=typeCount[type]+1
typeString=' '
for type in range(len(typeList)):
    typeString=typeString+str(int(typeCount[type]))+' '+typeList[type]+'\t '
print('We use '+typeString)
#clean the cell and convert it, searching for irrational numbers and even fractions
cleanList=[float(1.0/3.0), float(2.0/3.0)]

for i in range(1,10):
    cleanList.append(np.sqrt(3)/float(i))
    cleanList.append(2*np.sqrt(3)/float(i))
    cleanList.append(3*np.sqrt(3)/float(i))
    cleanList.append(4*np.sqrt(3)/float(i))
    cleanList.append(5*np.sqrt(3)/float(i))
    cleanList.append(np.sqrt(2)/float(i))
    cleanList.append(2*np.sqrt(2)/float(i))
    cleanList.append(3*np.sqrt(2)/float(i))
    cleanList.append(4*np.sqrt(2)/float(i))
    cleanList.append(5*np.sqrt(2)/float(i))
    cleanList.append(float(i)/6.0)
    cleanList.append(float(i)/8.0)

for n in range(3):
    for i in range(3):
        for c in range(len(cleanList)):
            if (np.absolute(np.absolute(shapeCell[n,i])-np.absolute(cleanList[c])))<symmPrec:
                shapeCell[n,i]=np.sign(shapeCell[n,i])*cleanList[c] 
Cell=np.multiply(shapeCell,np.matrix([precLatParam,precLatParam,precLatParam]))
#Cell=np.vstack((shapeCell[0,:]*precLatParam[0],shapeCell[1,:]*precLatParam[1],shapeCell[2,:]*precLatParam[2]))
print("Cell:")
print(Cell)

#also clean the atom positions

atom_pos=atom_pos_raw
for n in range(numAtoms):
    for i in range(3):
        for c in range(len(cleanList)):
            if (np.absolute(np.absolute(atom_pos[n][i])-np.absolute(cleanList[c])))<symmPrec:
#                print(str(n)+" "+str(i)+" "+str(atom_pos[n][i])+" to: ")
                atom_pos[n][i]=np.sign(atom_pos[n][i])*cleanList[c] 
#                print(str(atom_pos[n][i]))
#convert to cartesian
#print(precLatParam)
for n in range(numAtoms):
    for i in range(3):
        atom_pos[n][i]=atom_pos[n][i]*precLatParam[i]
posMatCart=np.matrix(atom_pos)
#print(posMatCart)

#fine all the maxses
atomicNumList=[]
atomicMassList=[]
for i in range(len(typeList)):
        foundAtom=0
        for k in range(len(atomList)):
            if atomList[k][1]==typeList[i]:
                atomicNumList.append(atomList[k][0])
                atomicMassList.append(atomList[k][2])
                foundAtom=foundAtom+1
                break
        if foundAtom==0:
            print("Could not recognize atom "+str(typeList[i])+", you'll have to manually enter mass and atomic number")
            atomicNumList.append(0)
            atomicMassList.append(0)


#write header file and displacements
headerName="headerFile_"+str(targetIrrep)+".dat"
h=open(headerName,'w')
h.write("Irrep: "+str(targetIrrep)+"\n")
h.close()
h=open(headerName,'a')
h.write("NumSAM: "+str(numModes)+"\n")
h.write("NumAtomTypes: "+str(len(typeList))+"\n")
h.write("NumAtoms: "+str(numAtoms)+"\n")
h.write("DispMag(angstrom): "+str(dispMag)+"\n")
for i in range(len(typeList)):
      h.write(f"{typeList[i]} {int(typeCount[i])} {atomicMassList[i]}\n")

distMat=np.zeros((numAtoms,3,numModes+1))
modeInd=1
modeTypeList=[" "]*(numModes+1)

startLine=8+numAtoms+3
if(targetOutput[startLine-2].split()[0]=="The"): #"The following include N translational modes"
    startLine=startLine+1 #we want to skip this line if it exists

for l in range(startLine,len(targetOutput)):
    thisLine=targetOutput[l].split()
    if (targetOutput[l]=='------------------------------------------'):
        modeInd=modeInd+1
    else:    
        atom=int(thisLine[0])-1
        modeTypeList[modeInd]=thisLine[1]
        disp1=float(thisLine[2])
        disp2=float(thisLine[3])
        disp3=float(thisLine[4])
        distMat[atom,0,modeInd]=disp1        
        distMat[atom,1,modeInd]=disp2        
        distMat[atom,2,modeInd]=disp3        

#now normalize the SAMs:
for m in range(1,numModes+1):
    distMat[:,:,m]=distMat[:,:,m]/np.sqrt(np.sum(np.multiply(distMat[:,:,m],distMat[:,:,m])))

#now orthogonalize them
orthMat=np.zeros((numAtoms,3,numModes+1))
for m in range(1,numModes+1):
    SAM=distMat[:,:,m]
    for n in range(m+1,numModes+1):
        #subtract off projection for every other mode we haven't orthogonalized yet:
        SAM=SAM-distMat[:,:,n]*np.sum(np.multiply(distMat[:,:,n],SAM))
    #re-normalize
    orthMat[:,:,m]=SAM/np.sqrt(np.sum(np.multiply(SAM,SAM)))
distMat=orthMat

#abinit needs cross(R1,R2)*R3 to be positive
crossDot=np.dot(np.cross(Cell[0,:],Cell[1,:]),np.transpose(Cell[2,:]))
print("cross(R1,R2)*R3 is: ", crossDot, ", must be positive or Abinit will freak out")
for m in range(numModes+1):
    filename="dist_"+str(targetIrrep)+"_"+str(m)
    thisDispCart=posMatCart+1.88973*dispMag*distMat[:,:,m] #the 1.88973 converts from angstrom to bohr
    thisDispFrac=np.matmul(thisDispCart,np.linalg.inv(Cell))
    f=open(filename,'w')
#    f.write("Automatically generated with SMODES for irrep "+str(targetIrrep)+"\n")
    f.close()
    f=open(filename,'a')
    f.write("natom "+str(numAtoms)+"\n")
    f.write("ntypat "+str(len(typeCount))+"\n")
    f.write("typat ")
    for i in range(len(typeCount)):
        for j in range(int(typeCount[i])):
            f.write(str(int(i)+1)+'  ')
    f.write("\n")
    f.write("znucl ")
    for i in range(len(atomicNumList)):
        f.write(str(atomicNumList[i])+" ")

    f.write("\n")
    f.write("xred \n")
    #we're also going to use this loop to write displacements to the header file
    if m>0:
        h.write("SAM_"+str(m)+": "+modeTypeList[m]+"\n")
    for i in range(numAtoms):
        for j in range(3):
            f.write('{:.10f}'.format(thisDispFrac[i,j])+" ")
            if m>0:
                h.write('{:.10f}'.format(distMat[i,j,m])+" ")
        f.write("\n")
        if m>0:
            h.write("\n")
    f.write("acell 3*1.0\nrprim\n") 
    for i in range(3):
        for j in range(3):
            f.write('{:.10f}'.format(Cell[i,j])+" ")
        f.write("\n")

    f.close()
h.close()


