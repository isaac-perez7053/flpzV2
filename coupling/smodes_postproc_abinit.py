#!/usr/bin/env python3
import numpy as np
import xml.etree.ElementTree as ET
import sys
np.set_printoptions(precision=10)
#takes as input an smodes input file and irrep name and creates the distortions needed to do a symmetry adapted modes calculation
#we expect a POSCAR file in this directory to get the header from
targetIrrep=str(sys.argv[1])

headerFile="SMODES_"+targetIrrep+"/headerFile_"+targetIrrep+".dat"
h=open(headerFile)
headerLines=h.readlines()
h.close()
Irrep=headerLines[0].split()[1]
NumSAM=int(headerLines[1].split()[1])
NumAtomTypes=int(headerLines[2].split()[1])
NumAtoms=int(headerLines[3].split()[1])
DispMag=float(headerLines[4].split()[1])

print("Irrep: "+Irrep)
print("NumSAM: "+str(NumSAM))
print("NumAtomTypes: "+str(NumAtomTypes))
print("NumAtoms: "+str(NumAtoms))
print("DispMag: "+str(DispMag))

typeList=[]
typeCount=[]
massList=[]
for lineNum in range(5,5+NumAtomTypes):
	words = headerLines[lineNum].split()
	typeList.append(words[0])
	typeCount.append(int(words[1]))
	if words[2]=="MASS":
		print("You forgot to set the mass in "+headerFile+" (in AMU)! We can't go any further.")
		print("Quitting...")
		quit()
	else:
		massList.append(float(words[2]))
	
typeText="Unit cell consists of "
for i in range(len(massList)):
	typeText=typeText+str(typeCount[i])+" "+str(typeList[i])+" atoms of mass "+str(massList[i])+", "

print(typeText[:-2])

#we're just going to leave the first matrix blank so we don't have to reindex the files
SAMmat=np.zeros((NumAtoms,3,NumSAM+1))
SAMatomLabel=[]
modeIndex=0
atomInd=0
for lineNum in range(5+NumAtomTypes,len(headerLines)):
	words = headerLines[lineNum].split()
	if words[0][0:3]=="SAM":
		modeIndex=modeIndex+1
		atomInd=0
		SAMatomLabel.append(words[1])
	else:
		SAMmat[atomInd,0,modeIndex]=float(words[0])
		SAMmat[atomInd,1,modeIndex]=float(words[1])
		SAMmat[atomInd,2,modeIndex]=float(words[2])
		atomInd=atomInd+1

#now get the forces from the runs

#its one bigger because we have the initial forces we're going to subtract off
forceMat_raw=np.zeros((NumAtoms,3,NumSAM+1))
for SAM in range(NumSAM+1):
	thisOUTCAR="SMODES_"+targetIrrep+"/dist_"+str(SAM)+"/dist_"+str(SAM)+".abo"
	f=open(thisOUTCAR)
	OUTCARLines=f.readlines()
	f.close()	
	lineStart=0
	atomInd=0
	for lineNum in range(len(OUTCARLines)):
		words = OUTCARLines[lineNum].split()
		if len(words)>=1:
			if (words[0]=="cartesian") and (words[1]=="forces" ) and (words[2]=="(eV/Angstrom)" ):
				lineStart=lineNum+1
				break
	for lineNum in range(lineStart,lineStart+NumAtoms):
		words = OUTCARLines[lineNum].split()
		forceMat_raw[atomInd,0,SAM]=float(words[1])
		forceMat_raw[atomInd,1,SAM]=float(words[2])
		forceMat_raw[atomInd,2,SAM]=float(words[3])
		atomInd=atomInd+1

#subtract off initial forces and now fix the indexing to go from zero again
forceList=np.zeros((NumAtoms,3,NumSAM))
for SAM in range(NumSAM):
	for i in range(NumAtoms):
		for j in range(3):
			forceList[i,j,SAM]=forceMat_raw[i,j,SAM+1]-forceMat_raw[i,j,0]			
#build force matrix
forceMat=np.zeros((NumSAM,NumSAM))
SAMmat=SAMmat[:,:,1:]
for f in range(NumSAM):
	for s in range(NumSAM):
		forceVal=np.multiply(forceList[:,:,f],SAMmat[:,:,s])
		forceMat[f,s]=forceVal.sum()
#build mass matrix
#for m in range(NumSAM):
#	thisMass=0
#	for n in range(len(typeList)):
#		if SAMatomLabel[m]==typeList[n]:
#			thisMass=massList[n]
#	if thisMass==0:
#		print("Problem with building mass matrix. Quitting...")
#		quit()
#	else:
#		M[m,m]=np.sqrt(thisMass)

#make a vector thats the list of masses in order of the typelist
massVec=np.zeros((NumSAM))
for m in range(NumSAM):
	thisMass=0
	for n in range(len(typeList)):
		if SAMatomLabel[m]==typeList[n]:
			thisMass=massList[n]
			massVec[m]=thisMass
	if thisMass==0:
		print("Problem with building mass matrix. Quitting...")
		quit()

MM=np.zeros((NumSAM,NumSAM))
for m in range(NumSAM):
	for n in range(NumSAM):
		MM[m,n]=np.sqrt(massVec[m])*np.sqrt(massVec[n])
#MM=np.matmul(M,npM)
#divide by disp mag to get fc matrix

FC_mat=-forceMat/DispMag

#symmetrize

FC_mat=(FC_mat+np.transpose(FC_mat))/2.0
Dyn_mat=np.divide(FC_mat,MM)
#solve for the eigensystem, evals are still in SAM basis	
FCevals, FCevecs_SAM = np.linalg.eig(FC_mat)
Dynevals, Dynevecs_SAM = np.linalg.eig(Dyn_mat)
#convert Dyn evals to frequency in THz
eV_to_J=1.602177E-19
ang_to_m=1.0E-10
AMU_to_kg=1.66053E-27
c= 2.9979458E10 #speed of light
#c= 3E10 #speed of light
#make it so imaginary frequencies look negative
Freq_THz=np.multiply(np.sign(Dynevals),np.sqrt(np.absolute(Dynevals)*eV_to_J/(ang_to_m**2*AMU_to_kg))*1.0E-12)
FC_eval=np.multiply(np.sign(FCevals),np.sqrt(np.absolute(FCevals)))

#sort them greatest to least 
idx_Dyn = np.flip(Freq_THz.argsort()[::-1] ) 
Freq_THz = Freq_THz[idx_Dyn]/(2*np.pi) #convert from 2piTHz
Dynevecs_SAM = Dynevecs_SAM[:,idx_Dyn]

#Freq_THz = Freq_THz/(2*np.pi) #convert from 2piTHz (did this above)
#convert to cm^-1 
Freq_cm=Freq_THz*1.0E12/(c)


idx_FC = np.flip(FC_eval.argsort()[::-1] ) 
FC_eval = FC_eval[idx_FC]
FCevecs_SAM = FCevecs_SAM[:,idx_FC]

#print(Dyn_mat)
#print(Dynevecs_SAM[:,0]) #this one is the eigenvector, with second index corresponding to eig index
#print(Dynevecs_SAM[0,:])

#now convert evecs to real space basis
Dynevecs=np.zeros((NumAtoms,3,NumSAM))
Fcevecs=np.zeros((NumAtoms,3,NumSAM))
for evec in range(NumSAM):
	realDynEvec=np.zeros((NumAtoms,3))
	realFCEvec=np.zeros((NumAtoms,3))
	for s in range(NumSAM):
		#realDynEvec=realDynEvec+Dynevecs_SAM[evec,s]*SAMmat[:,:,s]#wrong bc you used 1st index as eig ind	
		#realFCEvec=realFCEvec+FCevecs_SAM[evec,s]*SAMmat[:,:,s]	#wrong bc you used 1st index as eig ind
		realDynEvec=realDynEvec+Dynevecs_SAM[s,evec]*SAMmat[:,:,s]	
		realFCEvec=realFCEvec+FCevecs_SAM[s,evec]*SAMmat[:,:,s]	
	Dynevecs[:,:,evec]=realDynEvec
	Fcevecs[:,:,evec]=realFCEvec

#make mass matrix for defining reduced mass and phonon displacement eigenvectors
MassCol=np.zeros((NumAtoms,3))
atomind=0
for atype in range(NumAtomTypes):
	for j in range(typeCount[atype]):
		MassCol[atomind,:]=np.sqrt(massList[atype])
		atomind=atomind+1

#define phonon displacement eigenvectors, u
#these are solutions to the generalized eigenvalue problem with
#omega^2M^(1/2)u=M^(-1/2)Fu, where F is the FC matrix and
#M^(1/2) is a matrix with sqrt mass along the diagonal. 
#They are equal to M^(-1/2)e, where
#e is a phonon eigenvector of the dynamical matrix. These
#are the things you'd freeze in to find phonon frequencies
#from displacement curvature (with help from reduced mass):
#frequency=sqrt(deriv2/(redmass*AMU_to_kg)); %in 2piTHz
#frequency=frequency/(2*pi*3e10); %in cm^-1

PhonDispEigs=np.zeros((NumAtoms,3,NumSAM))
redmassvec=np.zeros((NumAtoms,1)) #the reduced mass associated with the phonon mode
for mode in range(NumSAM):
	PhonDispEigs[:,:,mode]=np.divide(Dynevecs[:,:,mode],MassCol)
	magSquared=np.sum(np.sum(np.multiply(PhonDispEigs[:,:,mode],PhonDispEigs[:,:,mode])))
	redmassvec[mode]=1.0/magSquared
	#normalize the new eigenvector
	PhonDispEigs[:,:,mode]=PhonDispEigs[:,:,mode]/np.sqrt(magSquared)

#now write these to files and put them in the original directory

DynFreqsFileName="SMODES_"+targetIrrep+"/DynFreqs.dat"
DynevecFileName="SMODES_"+targetIrrep+"/DynEvecs.dat"
FCevalFileName="SMODES_"+targetIrrep+"/FCEvals.dat"
FCevecFileName="SMODES_"+targetIrrep+"/FCEvecs.dat"
PhonDispFileName="SMODES_"+targetIrrep+"/PhonDispVecs.dat"
RedMassFileName="SMODES_"+targetIrrep+"/RedMass.dat"

#write the Dynmat files
Dval=open(DynFreqsFileName,'w')
Dvec=open(DynevecFileName,'w')
Dval.write("THz \t cm^-1 \n")
#Dv.write("Irrep: "+str(targetIrrep)+"\n")
for mode in range(NumSAM):
	valstring="%.2f \t %.2f \n" % (Freq_THz[mode], Freq_cm[mode])
	Dval.write(valstring)
	for atom in range(NumAtoms):
		for j in range(3):
			Dvec.write("%.5f \t"%(Dynevecs[atom,j,mode]))
#			print("Dynevecs["+str(atom)+","+str(j)+","+str(mode)+"] = "+str(Dynevecs[atom,j,mode]))
	Dvec.write("\n")
Dval.close()
Dvec.close()

#write the FC files
FCval=open(FCevalFileName,'w')
FCvec=open(FCevecFileName,'w')
FCval.write("eV/A^2 \n")
#Dv.write("Irrep: "+str(targetIrrep)+"\n")
for mode in range(NumSAM):
	valstring="%.2f \n" % (FC_eval[mode])
	FCval.write(valstring)
	for atom in range(NumAtoms):
		for j in range(3):
			FCvec.write("%.5f \t"%(Fcevecs[atom,j,mode]))
	
	FCvec.write("\n")

FCval.close()
FCvec.close()

#write the phonon dist amplitude files and reduced masses
Phonvec=open(PhonDispFileName,'w')
Redmass=open(RedMassFileName,'w')
Redmass.write("AMU \n")
for mode in range(NumSAM):
	valstring="%.4f \n" % (redmassvec[mode])
	Redmass.write(valstring)

	for atom in range(NumAtoms):
		for j in range(3):
			Phonvec.write("%.5f \t"%(PhonDispEigs[atom,j,mode]))
	
	Phonvec.write("\n")

Phonvec.close()
Redmass.close()



