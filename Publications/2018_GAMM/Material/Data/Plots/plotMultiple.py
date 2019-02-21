#############################################
# Calculation of Energy release rate
# ----------------------
# by
#
# Dr.-Ing. Christian Willberg
# christian.willberg@dlr.de
#
# ----------------------
#
# Execution: python lowPass.py
# 
##############################################
import numpy as np
import matplotlib.pyplot as plt
import csv

def readData(name='noName.csv', nodesName='N', displacementName = "avg(Displacement (1))", forceName = "avg(Force (1))"):

    infile = open(name, "r")
    csvdata = csv.reader(infile, delimiter=",")
    data = []
    for row in csvdata:
        data.append(row)
    checkHeader = data[0]
    count = 0
    nnID = -1
    d1AVGID = 0
    f1AVGID = 0
    for check in checkHeader:
        if check == nodesName:
            nnID = count
        if check == displacementName:
            d1AVGID = count
        if check == forceName:
            f1AVGID = count
        count += 1
    
    disp = []
    force = []
    # if nNodes not given, because global values are used the value is set to one 
    
    nNodes = 1

    for dat in data[1:]:
        if nnID!=-1:
            nNodes = int(dat[nnID])
        disp.append(float(dat[d1AVGID]))
        force.append(-float(dat[f1AVGID]))
    return disp, force, nNodes
def linearFunctions(disp, force, n = 1):
    if n > len(force)-1:
        n=len(force)-1
    if n<0:
        n=0
    
    linY = []
    for i in range(0,len(force)):
        
        if i > n:
            linY.append(force[i])
        else:
            linY.append(force[n]/disp[n]*disp[i])
    
    return linY
def integration(disp, force, linY):
    
    IntVal = 0.0
    for i in range(0,len(force)):
        IntVal += force[i]-linY[i]
    IntVal = IntVal * disp[1]    
        
    
    return IntVal    

#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
# INPUT DEFINITION
#*******************************************************************************
#*******************************************************************************
# DEFINE INPUT AND OUTPUT FILE NAME
input = 'h2d.csv'
output = 'output'
# DEFINE UNIT DEFINITION
unitDisp = "m"
unitForce = "N"
# DEFINE WIDTH OF THE DCB PROBE
width = 0.003

# DEFINE THE TWO CRACK LENGHTS AT SOLUTION STEP N1
#n1 = 400
#l1 = 16.*0.00033
n1= 1500
l1 = 15.*0.00033
# AND SOLUTION STEP N2
n2 = 2600
l2 = 56*0.00033

########################################################################
#DEFINE SOLUTION ROW NAMES TO FIND THE DATA EXPECTED
# IN THIS CASE THE AVERAGE FORCE, THE NUMBER OF NODES SELECTED AND 
# THE DISPLACEMENT y DIRECTION IS USED.
########################################################################
nodesName='N'
displacementName = "External_Displacement:1"
forceName = "Reaction_Force:1"
########################################################################
# CALCULATION
########################################################################
outputEnergy = output + "EnergyValues"
# Input data file name; row name
input = 'dx00005hor001005.csv'
disp1, force1, nNodes = readData(name=input, nodesName=nodesName, displacementName = displacementName, forceName = forceName)
input = 'dx00005hor001505.csv'
disp2, force2, nNodes = readData(name=input, nodesName=nodesName, displacementName = displacementName, forceName = forceName)
input = 'dx00005hor002005.csv'
disp3, force3, nNodes = readData(name=input, nodesName=nodesName, displacementName = displacementName, forceName = forceName)
input = 'dx00005hor002505.csv'
disp4, force4, nNodes = readData(name=input, nodesName=nodesName, displacementName = displacementName, forceName = forceName)

input = 'dx000033hor0006633.csv'
disp5, force5, nNodes = readData(name=input, nodesName=nodesName, displacementName = displacementName, forceName = forceName)
input = 'dx000033hor0009933.csv'
disp6, force6, nNodes = readData(name=input, nodesName=nodesName, displacementName = displacementName, forceName = forceName)
input = 'dx000033hor0013233.csv'
disp7, force7, nNodes = readData(name=input, nodesName=nodesName, displacementName = displacementName, forceName = forceName)
input = 'dx000033hor00167.csv'
disp8, force8, nNodes = readData(name=input, nodesName=nodesName, displacementName = displacementName, forceName = forceName)
input = 'dx000025hor0005025.csv'
disp9, force9, nNodes = readData(name=input, nodesName=nodesName, displacementName = displacementName, forceName = forceName)
input = 'dx000025hor0007525.csv'
disp10, force10, nNodes = readData(name=input, nodesName=nodesName, displacementName = displacementName, forceName = forceName)
input = 'dx000025hor0010025.csv'
disp11, force11, nNodes = readData(name=input, nodesName=nodesName, displacementName = displacementName, forceName = forceName)
input = 'dx000025hor0012525.csv'
disp12, force12, nNodes = readData(name=input, nodesName=nodesName, displacementName = displacementName, forceName = forceName)
input = 'dx0000125hor0000251.csv'
disp13, force13, nNodes = readData(name=input, nodesName=nodesName, displacementName = displacementName, forceName = forceName)
input = 'dx0000125hor00050125.csv'
disp14, force14, nNodes = readData(name=input, nodesName=nodesName, displacementName = displacementName, forceName = forceName)

#*******************************************************************************
#*******************************************************************************

input = 'plotCrackBranchCritStretch.csv'
disp15, force15, nNodes = readData(name=input, nodesName=nodesName, displacementName = displacementName, forceName = forceName)

input = 'plotCrackBranchCritEnergy075.csv'
disp16, force16, nNodes = readData(name=input, nodesName=nodesName, displacementName = displacementName, forceName = forceName)
input = 'plotCrackBranchCritEnergy081.csv'
disp17, force17, nNodes = readData(name=input, nodesName=nodesName, displacementName = displacementName, forceName = forceName)
input = 'plotCrackBranchCritEnergy084.csv'
disp18, force18, nNodes = readData(name=input, nodesName=nodesName, displacementName = displacementName, forceName = forceName)

#################
# PLOT
#################

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.tick_params(axis='both', which='major', labelsize=40, label2On=False,color='black')
ax.plot(disp1, force1, 'k-',linewidth=2, label='dx=0.5mm; h=1.005mm')
ax.plot(disp2, force2, 'r-',linewidth=2, label='dx=0.5mm; h=1.505mm')
ax.plot(disp3, force3, 'k--',linewidth=2, label='dx=0.5mm; h=2.005mm')
ax.plot(disp4, force4, 'r--',linewidth=2, label='dx=0.5mm; h=2.405mm')
ax.legend(fontsize = 30, loc = 4)
plt.grid(color='k', linestyle='--')
plt.xlabel('Displacement [m]', fontsize=40)
plt.ylabel('summedForce [N]', fontsize=40)
#labels = [item.get_text() for item in ax.get_xticklabels()]
from matplotlib.ticker import FormatStrFormatter 
ax.xaxis.set_major_formatter(FormatStrFormatter('%1.1E'))

plt.show()
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.tick_params(axis='both', which='major', labelsize=40, label2On=False,color='black')
ax.plot(disp5, force5, 'k-',linewidth=2, label='dx=0.33mm; h=0.663mm')
ax.plot(disp6, force6, 'r-',linewidth=2, label='dx=0.33mm; h=0.993mm')
ax.plot(disp7, force7, 'k--',linewidth=2, label='dx=0.33mm; h=1.323mm')
ax.plot(disp8, force8, 'k--',linewidth=2, label='dx=0.33mm; h=1.1653mm')
ax.legend(fontsize = 30, loc = 4)
plt.grid(color='k', linestyle='--')
plt.xlabel('Displacement [m]', fontsize=40)
plt.ylabel('summedForce [N]', fontsize=40)
#labels = [item.get_text() for item in ax.get_xticklabels()]

ax.xaxis.set_major_formatter(FormatStrFormatter('%1.1E'))
plt.show()
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.tick_params(axis='both', which='major', labelsize=40, label2On=False,color='black')
ax.plot(disp9, force9,   'k-',linewidth=2, label='dx=0.25mm; h=0.5025mm')
ax.plot(disp10, force10, 'r-',linewidth=2, label='dx=0.25mm; h=0.7525mm')
ax.plot(disp11, force11, 'b-',linewidth=2, label='dx=0.25mm; h=1.0025mm')
ax.plot(disp12, force12, 'k--',linewidth=2, label='dx=0.25mm; h=1.2525mm')
ax.plot(disp13, force13, 'r--',linewidth=2, label='dx=0.125mm; h=0.251mm')
ax.plot(disp14, force14, 'b--',linewidth=2, label='dx=0.125mm; h=0.501mm')
ax.legend(fontsize = 30, loc = 4)
plt.grid(color='k', linestyle='--')
plt.xlabel('Displacement [m]', fontsize=40)
plt.ylabel('summedForce [N]', fontsize=40)
#labels = [item.get_text() for item in ax.get_xticklabels()]

ax.xaxis.set_major_formatter(FormatStrFormatter('%1.1E'))
#newLabels = ['%E'.format(label) for label in labels]
plt.show()
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.tick_params(axis='both', which='major', labelsize=40, label2On=False,color='black')
ax.plot(disp3, force3,   'k-',linewidth=2, label='dx=0.5mm; h=2.005mm')
ax.plot(disp4, force4,   'r-',linewidth=2, label='dx=0.5mm; h=2.405mm')
ax.plot(disp6, force6,   'b-',linewidth=2, label='dx=0.33mm; h=0.993mm')
ax.plot(disp7, force7,   'm-',linewidth=2, label='dx=0.33mm; h=1.323mm')
ax.plot(disp11, force11, 'k--',linewidth=2, label='dx=0.25mm; h=1.0025mm')
ax.plot(disp12, force12, 'r--',linewidth=2, label='dx=0.25mm; h=1.2525mm')
ax.plot(disp14, force14, 'b--',linewidth=2, label='dx=0.125mm; h=0.501mm')
ax.legend(fontsize = 30, loc = 4)
plt.grid(color='k', linestyle='--')
plt.xlabel('Displacement [m]', fontsize=40)
plt.ylabel('summedForce [N]', fontsize=40)
#labels = [item.get_text() for item in ax.get_xticklabels()]

ax.xaxis.set_major_formatter(FormatStrFormatter('%1.1E'))
plt.show()

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.tick_params(axis='both', which='major', labelsize=40, label2On=False,color='black')
ax.plot(disp11, force11, 'k--',linewidth=2, label='dx=0.25mm; h=1.0025mm')
ax.plot(disp12, force12, 'r--',linewidth=2, label='dx=0.25mm; h=1.2525mm')
ax.plot(disp14, force14, 'b--',linewidth=2, label='dx=0.125mm; h=0.501mm')
#ax.set_xticklabels(newLabels)
ax.legend(fontsize = 30, loc = 4)
plt.grid(color='k', linestyle='--')
plt.xlabel('Displacement [m]', fontsize=40)
plt.ylabel('summedForce [N]', fontsize=40)
#labels = [item.get_text() for item in ax.get_xticklabels()]

ax.xaxis.set_major_formatter(FormatStrFormatter('%1.1E'))

plt.show()

print("plot branching")
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.tick_params(axis='both', which='major', labelsize=40, label2On=False,color='black')
ax.plot(disp15, force15, 'k-',linewidth=2, label='Critical Stretch sc = 0.000433593')
ax.plot(disp16, force16, 'k--',linewidth=2, label='G0 = 0.75 N/m')
ax.plot(disp17, force17, 'r--',linewidth=2, label='G0 = 0.81 N/m')
ax.plot(disp18, force18, 'r-',linewidth=2, label='G0 = 0.84 N/m')
#ax.set_xticklabels(newLabels)
ax.legend(fontsize = 30, loc = 1)
plt.grid(color='k', linestyle='--')
plt.xlabel('Displacement [m]', fontsize=40)
plt.ylabel('summedForce [N]', fontsize=40)
#labels = [item.get_text() for item in ax.get_xticklabels()]

ax.xaxis.set_major_formatter(FormatStrFormatter('%1.1E'))

plt.show()
