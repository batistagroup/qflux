import numpy as np
from qiskit import QuantumRegister, QuantumCircuit

from .params import I,X,Y,Z

##=======the Walsh operator scheme===========
#1. functions for Walsh representation
#(a). the binary and qubit state representation of a integer
#(b). the function to calculate walsh coefficient
#this is binary representation of a integer j!
#j=[j1,j2,j3] for example: 3=[1,1,0]
def binary(j,nbits):
  arr = np.zeros(nbits,dtype=int)
  res = j
  for i in range(nbits):
    if(res>1):
        arr[i] = res%2
        res=res//2
    else:
        arr[i] = res
        break
  return arr

#input: array of 0 or 1
#output: max index of 1
#e.g.: j=4 binary is [0,0,1,0], max_bit is 2
def max_bit(arr):
  iflag = 0
  for i in range(len(arr)-1,-1,-1):
    if(arr[i]==1):
        iflag = i
        break
  return iflag

#j-th gray code:
#first do the right shift of the binary representation
#then do the bit wise XOR
def gray(j,nbits):
  arr_n1 = binary(j//2,nbits) #
  arr_j = binary(j,nbits) #
  arr_res = np.zeros_like(arr_n1)
  for i in range(nbits):
    if(arr_n1[i]==arr_j[i]):
        arr_res[i] = 0
    else:
        arr_res[i] = 1
  return arr_res

# return the decimal of arr_j
# [1,1,0] is j0=1,j1=1, which is 3
def decimal(arr_j):
  res = 0
  nbits = len(arr_j)
  for i in range(nbits):
    res+= 2**i*arr_j[i]
  return res

#2. The Walsh coefficient and operator
#input an fk array
#(that's the discretize array in qubit statevector) |qn,...,q1>
#e.g. f3 is the coefficient of |0,1,1>, which is q0=1, q1=1
#this just match the definition of function binary
#output the coefficient in the walsh function representation
def walsh_coef(fk,nbits):
  nlen = 2**nbits
  if(nlen!=len(fk)): print('walsh_coef: err in array length')

  a_arr = np.zeros(nlen)
  for i in range(nlen):
    arr_j = binary(i,nbits)
    for k in range(nlen):
      fac = (-1)**(np.dot(arr_j,binary(k,nbits)))
      a_arr[i] += fk[k]*fac
    a_arr[i] /= nlen
  return a_arr

#walsh operator, \Prod \otimes Z^(j)
#j=[1,1,0] return IZZ (q0,q1 is Z, q2 is I, in qiskit convention order)
def walsh_oprQ(arr_j,nbits):
  res = 1.0
  for i in range(nbits-1,-1,-1):
    if(arr_j[i]==0):
      res = np.kron(res,I)
    elif(arr_j[i]==1):
      res = np.kron(res,Z)
    else:
      print('err')
  return res


#3. Walsh circuit list and Gray code optimize for diagonal unitaries
#  functions to generate the quantum circuit of the diagnoal unitaries in the Walsh Operator representation
#input an walsh coefficient (arra) of f_k array,
#output the list that contain all the parameters to construct the circuit for U=diag(e^(i f_k))
def cirq_list_walsh(arra,Nqb,epsilon):
  if(len(arra)!=2**Nqb): print('err')

  U=[]

  for j in range(1,2**Nqb):

    #encoding: Gray Code or Binary Code
    #arr_j = binary(j,Nqb)
    arr_j = gray(j,Nqb)
    j_code = decimal(gray(j,Nqb))

    if(abs(arra[j_code])<epsilon): continue

    max_index = max_bit(arr_j) #the R gate will located at the max index

    for ibit in range(max_index):
      if(arr_j[ibit]==1):
        U.append(('C',ibit,max_index))

    U.append(('R',max_index,-2.0*arra[j_code]))

    for ibit in range(max_index-1,-1,-1):
      if(arr_j[ibit]==1):
        U.append(('C',ibit,max_index))
  return U

#generate circuit from list of turple U
def cirq_from_U(U,Nqb):
  nlen = len(U)

  #the quantum circuit
  qr = QuantumRegister(Nqb)
  #cr = ClassicalRegister(1)
  qc = QuantumCircuit(qr)

  for i in range(nlen):
    gate_str = U[i]
    if(gate_str[0]=='C'):
      qc.cx(gate_str[1],gate_str[2])
    elif(gate_str[0]=='R'):
      qc.rz(gate_str[2],gate_str[1])
    else:
      print('err in reading circuit')
  return qc


#Optimize the quantum circuit
def optimize(U):
  U0,iflag0 = scanI(U)
  U0,iflag1 = scanRule2(U0)

  #some possible exchange and check
  i = 0
  while(i<len(U0)-1):

    if(U0[i+1][0] == 'C' and U0[i][0] =='C'):

      ictrl = U0[i][1]; itarg = U0[i][2]
      if(ictrl == U0[i+1][1] or itarg == U0[i+1][2]):
        Utmp = U0.copy()
        Utmp[i] = U0[i+1]
        Utmp[i+1] = U0[i]
        Utmp,iflag0 = scanI(Utmp)
        Utmp,iflag1 = scanRule2(Utmp)

        #if(iflag0 or iflag1): U0 = Utmp
        U0 = Utmp

    i += 1
  return U0

#scan if there are identical adjacent CNOT gate
def scanI(U):
  nlen = len(U)
  newU = []
  i=0
  iflag = 0
  while(i<nlen-1):
    if(U[i]==U[i+1] and U[i][0] =='C'):
      i += 2
      iflag = 1
    else:
      newU.append(U[i])
      i += 1
  #end
  if(i==nlen-1): newU.append(U[nlen-1])
  return newU,iflag

#scan if there have 3 gate can reduced to:
#the target of one CNOT is the control of another.
#CijCjk = CjkCikCij
def scanRule2(U):
  nlen = len(U)
  newU = []
  i=0
  iflag = 0
  while(i<nlen-2):

    if(U[i][0] == U[i+1][0] == U[i+2][0] =='C'):
      sj = U[i][1]
      sk = U[i][2]
      si = U[i+1][1]
      if(sj == U[i+2][2] and si == U[i+2][1] and sk == U[i+1][2]):
        i += 3
        newU.append(('C',si,sj))
        newU.append(('C',sj,sk))
        print('Rule2')
        iflag = 1
      else:
        newU.append(U[i])
        i += 1
    else:
      newU.append(U[i])
      i += 1
  #end
  if(i==nlen-2):
    newU.append(U[nlen-2])
    newU.append(U[nlen-1])
  return newU,iflag
##=======end the Walsh operator scheme===========

##=======test functions for Walsh operator scheme=========
#extented matrix of the quantum gate
#note that this is the reverse order with qiskit convention
def ext_cx(ictrl,itarg,Nqb):
  if(ictrl>=Nqb or itarg>=Nqb): print('err')
  Proj00 = np.array([[1.0,0],[0,0]])
  Proj11 = np.array([[0.0,0],[0,1]])

  res00 = 1
  res11 = 1
  for i in range(Nqb-1,-1,-1):
  #for i in range(Nqb):
    if(i!=ictrl and i!=itarg):
      res00 = np.kron(res00,I)
      res11 = np.kron(res11,I)
    if(i==ictrl):
      res00 = np.kron(res00,Proj00)
      res11 = np.kron(res11,Proj11)
    if(i==itarg):
      res00 = np.kron(res00,I)
      res11 = np.kron(res11,X)
  return res00+res11

def ext_Rz(phi,iqb,Nqb):
  if(iqb>=Nqb): print('err')
  mat_rz = np.array([[np.exp(-1j*phi/2),0],[0,np.exp(1j*phi/2)]])
  res = 1
  for i in range(Nqb-1,-1,-1):
  #for i in range(Nqb):
    if(i!=iqb):
      res = np.kron(res,I)
    if(i==iqb):
      res = np.kron(res,mat_rz)
  return res

def ext_Z(iqb,Nqb):
  if(iqb>=Nqb): print('err')
  res = 1
  for i in range(Nqb-1,-1,-1):
  #for i in range(Nqb):
    if(i!=iqb):
      res = np.kron(res,I)
    if(i==iqb):
      res = np.kron(res,Z)
  return res

#generate matrix of quantum circuit from list of turple U
#convient for verify the result
def cirqmat_from_U(U,Nqb):
  nlen = len(U)

  Umat = []

  for i in range(nlen):
    gate_str = U[i]
    if(gate_str[0]=='C'):
      Umat.append(ext_cx(gate_str[1],gate_str[2],Nqb))
    elif(gate_str[0]=='R'):
      Umat.append(ext_Rz(gate_str[2],gate_str[1],Nqb))
    else:
      print('err in reading circuit')

  #connect all the gates
  res = Umat[-1]
  for i in range(nlen-2,-1,-1):
    res = res@Umat[i]
  return res

#input an walsh coefficient (arra) of f_k array, output the circuit matrix for U=diag(e^(i f_k))
def cirqmat_walsh(arra,Nqb,epsilon):
  if(len(arra)!=2**Nqb): print('err')

  U = []

  for j in range(1,2**Nqb):
    if(abs(arra[j])<epsilon): continue
    arr_j = binary(j,Nqb)
    max_index = max_bit(arr_j) #the R gate will located at the max index

    for ibit in range(max_index):
      if(arr_j[ibit]==1):
        U.append(ext_cx(ibit,max_index,Nqb))

    U.append(ext_Rz(-2*arra[j],max_index,Nqb))

    for ibit in range(max_index-1,-1,-1):
      if(arr_j[ibit]==1):
        U.append(ext_cx(ibit,max_index,Nqb))

  #connect all the gates
  res = U[-1]
  for i in range(len(U)-2,-1,-1):
    res = res@U[i]
  return res
#============================