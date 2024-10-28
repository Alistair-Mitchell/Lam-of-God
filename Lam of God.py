import numpy as np
from colorama import Fore, Back, Style
from math import log10, floor
import itertools

np.set_printoptions(precision=3,suppress=True)

def round_sig(x, sig=3):
    return round(x, sig-int(floor(log10(abs(x))))-1)
def engineering_notation(x):
    if x == 0:
        return "0.000"
    
    exponent = int(np.floor(np.log10(abs(x)) / 3) * 3)
    coefficient = x / (10 ** exponent)
    
    return f"{coefficient:.3f}e{exponent}"
def print_engineering_matrix(matrix, use_table_format=False):
    # Apply engineering notation to the matrix
    formatted_matrix = np.vectorize(engineering_notation)(matrix)
    # Find the maximum width for consistent spacing
    max_width = max(len(entry) for row in formatted_matrix for entry in row)
    if use_table_format:
        # Create top border
        print("+" + "-" * (max_width + 2) * matrix.shape[1] + "+")
    # Print each row with consistent spacing
    for row in formatted_matrix:
        if use_table_format:
            row_str = " | ".join(f"{entry:>{max_width}}" for entry in row)
            print(f"| {row_str} |")
        else:
            row_str = "\t".join(f"{entry:>{max_width}}" for entry in row)
            print(row_str)
    if use_table_format:
        # Create bottom border after the last row
        print("+" + "-" * (max_width + 2) * matrix.shape[1] + "+")
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

def printCR(colour,string):
    colourdict = {
      "g": Fore.GREEN,
      "y": Fore.YELLOW,
      "r": Fore.RED,
      "b": Fore.BLUE
    }
    print(colourdict[colour],string,Style.RESET_ALL)
    
print("\nLam of God V1.4, AListair Mitchell (AJM) 2022")
layup = np.array([])  # Orientation of the fibres for each layer/ply of the laminate
plycount = layup.size

if plycount == 0:
    print("\nWhere a layup isn't specified, please enter a layup seperated by /")
    arr = input("Enter layup : ")  # takes the whole line of n numbers
    l = list(map(int, arr.split('/')))  # split those numbers with space( becomes ['2','3','6','6','5']) and then map every element into int (becomes [2,3,6,6,5])
    layup = np.array(l)
plycount = layup.size
#Please use SI units if you want anything useful out
plythickness = 0.005
YM11 = 20010000
YM22 = 1301000
Shear = 1001000
Poiss = 0.3

#Check for largest grouping of plies
grouping = biggestgroup(layup)

#Check for balanced laminate
balance=np.zeros((plycount))
for i in range(plycount):
    if layup[i]!=layup[i]%90:
        balance[i]=-(layup[i]%90)
    else:
        balance[i] = layup[i]
balance_check = np.sum(balance)

#Check for symetry
split = np.array_split(layup,2)
if plycount%2==1:
    lastElementIndex = len(split[0]) - 1
    split[0] = split[0][:lastElementIndex]
flip = np.flip(split[1])
symmetry = np.sum(abs(split[0]-flip))

#Check for quasi isotropy
if (symmetry == 0 and plycount%2==0 and plycount>2):
    max_angle = max(abs(np.diff(layup % 180)))
    anglehop =np.diff(split[0]%180)
    # QI = all_equal(anglehop)
    # print(QI)
    QI1 = all_equal(anglehop)
    if (180 / (plycount/2) - biggestgroup(anglehop)[0] == 0 and QI1 == True and (plycount/2)>=3):
        QI_Result = True
    else:
        QI_Result = False
elif plycount == 1:
    max_angle = 0
    QI_Result = False
else:
    max_angle = max(abs(np.diff(layup % 180)))
    anglehop = np.diff(layup % 180)
    QI1 = all_equal(anglehop)
    if(180/plycount-biggestgroup(anglehop)[0]==0 and QI1==True and plycount>=3):
        QI_Result=True
    else:
        QI_Result = False

# E1 = np.array([45.6,45.6,45.6])                 # (Use this instead if layer properties are not the same for all layers)
E1 = np.full((layup.shape[0], 1), YM11)[:, 0]  # Axial stiffness of each ply [GPa]
E2 = np.full((layup.shape[0], 1), YM22)[:, 0]  # transverse stiffness of each ply [GPa]
v12 = np.full((layup.shape[0], 1), Poiss)[:, 0]  # Poisson's ratio
v21 = (v12 * E2) / E1  # (Since the compliance matrix is symmetric, v12/E1=v21/E2)
G12 = np.full((layup.shape[0], 1), Shear)[:, 0]  # Shear modulus [GPa]

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
printCR('b',"\nLAYUP")
print(layup)
printCR('b',"\nLAMINATE MATERIAL PROPERTIES")
print('Ex = '+ engineering_notation(Ex)+ '\t\t(Primary Dir)\n'
      'Ey = ' + engineering_notation(Ey) + '\t\t(Secondary Dir)\n'
      'Gxy = '+ engineering_notation(Gxy) + '\t\t(Shear Modulus)\n'
      'vxy = '+str(vxy)+'\t\t(poissons ratio)\n')
top=np.concatenate((RealA,RealB),0)
bot=np.concatenate((RealB,RealD),0)
ABD=np.concatenate((top,bot),1)
printCR('b',"ABD MATRIX")
print_engineering_matrix(ABD,use_table_format=False)


#This is all reporting based on the ABD matrix, helpful to highlight how a laminate may not be ideal
printCR('b',"\nREPORT")
if (symmetry==0 and balance_check==0):
    printCR('g', "VERIFIED: Laminate is Balanced-Symmetric")
else:
    if symmetry==0:
        printCR('g',"VERIFIED: Layup is Symmetric")
    if balance_check==0:
        printCR('g',"VERIFIED: Laminate is Balanced")

if QI_Result==True:
    printCR('g',"VERIFIED: Laminate is quasi-isotropic")
if symmetry != 0:
    printCR('y',"WARNING: Layup is Asymmetric")
if result[1]>=3:
    printCR('y',"WARNING: Layup has "+
            str(grouping[1])+" plies at "+
            str(grouping[0])+" degrees next to each other")
if max_angle>45:
    printCR('y',"WARNING: Large angular difference between plies ("+
          str(max_angle)+" degrees)")
if (RealA[0,2]or RealA[1,2])!=0:
    printCR('r',"ACHTUNG: In-plane shear present (A16,A26)")
if (RealB[0,0]or RealB[1,1])!=0:
    printCR('r',"ACHTUNG: Tension-bending coupling present (B11,B22)")
if (RealB[0,1]or RealB[1,0])!=0:
    printCR('r',"ACHTUNG: Out-of-plane tension-bending coupling present (B12)")
if (RealB[0,2]or RealB[1,2]or RealB[2,0]or RealB[2,1])!=0:
    printCR('r',"ACHTUNG: Tension-twist coupling present (B16, B26)")
if (RealB[2,2])!=0:
    printCR('r',"ACHTUNG: Shear-twist coupling present (B66)")
if (RealD[0,2]or RealD[1,2])!=0:
    printCR('r',"ACHTUNG: Out-of plane shear present (D16,D26)")
