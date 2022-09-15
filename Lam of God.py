import numpy as np
from colorama import Fore, Back, Style
from math import log10, floor
import itertools

def round_sig(x, sig=3):
    return round(x, sig-int(floor(log10(abs(x))))-1)
def biggestgroup(input):
    grouped = ([(list(g)) for k, g in itertools.groupby(input)])
    list_len = [len(i) for i in grouped]
    # print(grouped)
    Frequency = max(list_len)
    Value_index = list_len.index(Frequency)
    # print(Value_index)
    Value = grouped[Value_index][0]
    return Value,Frequency
def all_equal(iterable):
    g = itertools.groupby(iterable)
    return next(g, True) and not next(g, False)

layup = np.array([0,45,90,-45])  # Orientation of the fibres for each layer/ply of the laminate
plythickness = 0.25
E11 = 125
E22 = 8
G12 = 5
v12 = 0.3



plycount = layup.size
result = biggestgroup(layup)
max_angle = max(abs(np.diff(layup%180)))

split = np.array_split(layup,2)
if plycount%2==1:
    lastElementIndex = len(split[0]) - 1
    split[0] = split[0][:lastElementIndex]
flip = np.flip(split[1])
symmetry = np.sum(abs(split[0]-flip))

if (symmetry == 0 and plycount%2==0):
    anglehop =np.diff(split[0]%180)
    QI = all_equal(anglehop)
else:
    anglehop = np.diff(layup % 180)
    QI1 = all_equal(anglehop)
    if(180/plycount-biggestgroup(anglehop)[0]==0 and QI1==True):
        QI_Result=True



E1 = np.full((layup.shape[0], 1), 125)[:, 0]  # Axial stiffness of each ply [GPa]
# E1 = np.array([45.6,45.6,45.6])                 # (Use this instead if layer properties are not the same for all layers)

E2 = np.full((layup.shape[0], 1), 8)[:, 0]  # transverse stiffness of each ply [GPa]
# E2 = np.array([16.2,16.2.16,2])                 # (Use this instead if layer properties are not the same for all layers)

v12 = np.full((layup.shape[0], 1), 0.3)[:, 0]  # Poisson's ratio
# v12 = np.array([0.278,0.278,0.278])             # (Use this instead if layer properties are not the same for all layers)
v21 = (v12 * E2) / E1  # (Since the compliance matrix is symmetric, v12/E1=v21/E2)

G12 = np.full((layup.shape[0], 1), 5)[:, 0]  # Shear modulus [GPa]
# G12 = np.array([5.83,5.83,5.83])                # (Use this instead if layer properties are not the same for all layers)

h = np.full((plycount),plythickness)
h_mod = np.full((plycount+1),plythickness)
for i in range(plycount+1):
    h_mod[i] = (i*plythickness)-((plycount/2)*plythickness)

# Initiating Q as a list that can contain the stiffness matrices transformed into the global coordinates for each layer
Q=[]
# Looping through each layer of the laminate
for i in range(layup.shape[0]):
    theta = layup[i]*np.pi/180 # Current ply angle changed into radians
    # assigning the current material properties to temporary variables (could also be used directly):
    E1l=E1[i]
    E2l=E2[i]
    v12l=v12[i]
    v21l=v21[i]
    G12l=G12[i]
    # print('Ply no. ' + str(i+1) + ': theta='+ str(layup[i]) +', E1l=' + str(E1l), ' E2l=' + str(E2l))
    # Establishing current local stiffness matrix (of the ply), Ql:
    Ql = 1/(1-v12l*v21l)*np.array([[E1l,v21l*E1l,0],\
                                   [v12l*E2l,E2l,0],\
                                   [0,0,G12l*(1-v12l*v21l)]])


    # Transformation matrix:
    T=np.array([[np.cos(theta)**2,np.sin(theta)**2,-2*np.sin(theta)*np.cos(theta)], \
                [np.sin(theta)**2,np.cos(theta)**2,2*np.sin(theta)*np.cos(theta)],\
                [np.sin(theta)*np.cos(theta),-np.sin(theta)*np.cos(theta),np.cos(theta)**2-np.sin(theta)**2]])
    # Adding the current stiffness matrix in the global coordinate system to the Q-list variable:
    Q.append(np.dot(np.dot(T,Ql),np.transpose(T)))


A=np.zeros((3,3))
RealA =np.zeros((3,3))
RealB =np.zeros((3,3))
RealD =np.zeros((3,3))
for i in range(layup.shape[0]):
    A=A+Q[i]*h[i]
for i in range(plycount):
    RealA = RealA.round(2) +Q[i]*(h_mod[i+1]-h_mod[i])
    RealB = RealB.round(2)+ (Q[i] * (h_mod[i + 1]**2 - h_mod[i]**2))*0.5
    RealD = RealD.round(2) + (Q[i] * (h_mod[i + 1] ** 3 - h_mod[i] ** 3))*0.33333333
RealA = RealA.round(2)
RealB = RealB.round(2)
RealD = RealD.round(2)

Qstar = A/sum(h)
Sstar = np.linalg.inv(Qstar)
Ex = round_sig(1/Sstar[0,0])
Ey = round_sig(1/Sstar[1,1])
Gxy = round_sig(1/Sstar[2,2])
vxy=-round_sig(Sstar[1,0]/Sstar[0,0])
print(Fore.BLUE+"LAYUP"+Style.RESET_ALL)
print(layup)
print('\n'+Fore.BLUE+"LAMINATE MATERIAL PROPERTIES"+Style.RESET_ALL)
print('Ex = '+ str(Ex)+ '        (Primary Dir)\nEy = ' + str(Ey) + '        (Secondary Dir)\nGxy = '+ str(Gxy) + '       (Shear Modulus)\nvxy = '+str(vxy)+'      (poissons ratio)\n')
top=np.concatenate((RealA,RealB),0)
bot=np.concatenate((RealB,RealD),0)
ABD=np.concatenate((top,bot),1)
print(Fore.BLUE+"ABD MATRIX"+Style.RESET_ALL)
print(ABD)
print("\n")

print(Fore.BLUE+"REPORT"+Style.RESET_ALL)
# print(Fore.GREEN + "VERIFIED: Laminate is Balanced"+Style.RESET_ALL)
if symmetry==0:
    print(Fore.GREEN + "VERIFIED: Layup is Symmetric"+Style.RESET_ALL)
else:
    print(Fore.YELLOW + "WARNING: Layup is Asymmetric" + Style.RESET_ALL)
# print(Fore.GREEN + "VERIFIED: Laminate is Balanced-Symmetric"+Style.RESET_ALL)
# print(Fore.GREEN + "VERIFIED: Laminate is QI"+Style.RESET_ALL)

if result[1]>=3:
    print(Fore.YELLOW + "WARNING: Layup has "+str(result[1])+" plies at "+str(result[0])+" degrees next to each other"+Style.RESET_ALL)
if max_angle>45:
    print(Fore.YELLOW + "WARNING: Large angular difference between plies ("+str(max_angle)+" degrees)"+Style.RESET_ALL)
if (RealA[0,2]or RealA[1,2])!=0:
    print(Fore.RED + "ACHTUNG: In-plane shear present (A16,A26)"+Style.RESET_ALL)
if (RealB[0,0]or RealB[1,1])!=0:
    print(Fore.RED + "ACHTUNG: Tension-bending coupling present (B11,B22)"+Style.RESET_ALL)
if (RealB[0,1]or RealB[1,0])!=0:
    print(Fore.RED + "ACHTUNG: Out-of-plane tension-bending coupling present (B12)"+Style.RESET_ALL)
if (RealB[0,2]or RealB[1,2]or RealB[2,0]or RealB[2,1])!=0:
    print(Fore.RED + "ACHTUNG: Tension-twist coupling present (B16, B26)"+Style.RESET_ALL)
if (RealB[2,2])!=0:
    print(Fore.RED + "ACHTUNG: Shear-twist coupling present (B66)"+Style.RESET_ALL)
if (RealD[0,2]or RealD[1,2])!=0:
    print(Fore.RED + "ACHTUNG: Out-of plane shear present (D16,D26)"+Style.RESET_ALL)
