import numpy as np
import functools
import sys

from mpsqd.utils import MPS
from mpsqd.utils.split import split_rq,split_qr

from mpsqd.tdvp import ttfunc as ttf


"""
Modified from Mpsqd (https://github.com/qiangshi-group/MPSQD)

The MIT License (MIT)
Copyright (c) 2020 mpsqd team, Qiang Shi group

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

By Weizhong Guan; Peng Bao; Jiawei Peng; Zhenggang Lan and Qiang Shi 
J. Chem. Phys. 161, 122501 (2024)
Changes:
    add rk4slices as argument (for update_type='rk4')
"""
def tdvp1site(rin,pall,dt,update_type='krylov',mmax=30,nsteps=1,rk4slices = 10):
  argdict = {'dt':dt,'mmax':mmax,'nsteps':nsteps,'update_type':update_type, 'rk4slices':rk4slices}
  try:
    for key in ['nsteps','mmax']:
      argdict[key] = int(argdict[key])
    for key in ['dt']:
      argdict[key] = float(argdict[key])
    for key in ['dt','nsteps','mmax']:
      if (argdict[key] <= 0):
        print("Wrong value in tdvp: ",key," = ",argdict[key])
        sys.exit()
    key = 'update_type'
    argdict[key] = argdict[key].lower()
    if (not argdict[key] in ['krylov','rk4']):
      print("Wrong value in tdvp: update_type = ",update_type)
      sys.exit()
    if(rin.__class__.__name__ != "MPS"):
      print("Wrong value in tdvp: improper MPS")
      sys.exit()
    if(pall.__class__.__name__ != "MPO"):
      print("Wrong value in tdvp: improper MPO")
      sys.exit()
  except (ValueError):
    print("Wrong value in tdvp: ",key," = ",argdict[key])
    sys.exit()
  for istep in range(nsteps):
    rin = _tdvp1(rin,pall,argdict['dt'],argdict['update_type'],argdict['mmax'],argdict['rk4slices'])
  return rin


"""
Modified from Mpsqd (https://github.com/qiangshi-group/MPSQD)

The MIT License (MIT)
Copyright (c) 2020 mpsqd team, Qiang Shi group

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

By Weizhong Guan; Peng Bao; Jiawei Peng; Zhenggang Lan and Qiang Shi 
J. Chem. Phys. 161, 122501 (2024)
Changes:
    add rk4slices as argument (for update_type='rk4')
"""
def _tdvp1(rin,pall,dt,update_type='krylov',mmax=30, rk4slices = 10):

  dt2 = 0.5*dt
  if(update_type=='krylov'):
    update_v = functools.partial(ttf.expmv,mmax=mmax,dt=dt2)
  else:
    update_v = functools.partial(ttf.update_rk4,dt=dt2,rk4slices=rk4slices)

  nlen = len(rin.nodes)
  r1 = MPS(nlen,rin.nb)
  r2 = MPS(nlen,rin.nb)
  r3 = MPS(nlen,rin.nb)

  phia = []

  # phia[nlevel+2]
  vtmp = np.ones((1,1,1),dtype=np.complex128)
  phia.append(vtmp)

#======================================================
# ortho from right
  r, q = split_rq(rin.nodes[nlen-1])
  
  r1.nodes.append(q)
  u1 = r
# phia[nlevel+1]
  phia.append(ttf.phia_next(phia[0],pall.nodes[nlen-1],q,q,1))
#-------------------------------------------------------
# intermediate terms
  for i in range(nlen-2,0,-1):

    rtmp = np.tensordot(rin.nodes[i], u1, axes=((2),(0)))

    r, q = split_rq(rtmp)
    r1.nodes.append(q)
    u1 = r

    phia.append(ttf.phia_next(phia[nlen-1-i],pall.nodes[i],q,q,1))

# the left matrix

  rtmp = np.tensordot(rin.nodes[0], u1, axes=((2),(0)))
  r1.nodes.append(rtmp)

  r1.nodes.reverse()

# add phia[0] and reverse to normal
  vtmp = np.ones((1,1,1),dtype=np.complex128)
  phia.append(vtmp)
  phia.reverse()

#===============================================
###### the first part of KSL, from left to right
  for i in range(nlen-1):
    phi1 = phia[i]
    phi2 = phia[i+1]

    if (i == 0):
      ksol = r1.nodes[i].copy()
    else:
      ksol = r2.nodes[i].copy()


    #for islice in range(pa.rk4slice):
    ksol = update_v(yy=ksol, phi1=phi1, phi2=phi2, mat1=pall.nodes[i])
    q, r = split_qr(ksol)


    if (i == 0):
      r2.nodes.append(q)
    else:
      r2.nodes[i] = q


    phi1 = ttf.phia_next(phi1,pall.nodes[i],q,q,0)
    phia[i+1] = phi1

# ?? need to copy
    ssol = r
    ssol = update_v(yy=ssol, phi1=phi1, phi2=phi2)

    rtmp1 = np.tensordot(ssol,r1.nodes[i+1],axes=((1),(0)))
    r2.nodes.append(rtmp1)

#--------------------------------------------------------
### right most part

  phi1 = phia[nlen-1]
  phi2 = phia[nlen]
  ksol = r2.nodes[nlen-1].copy()
  ksol = update_v(yy=ksol, phi1=phi1, phi2=phi2, mat1=pall.nodes[nlen-1])

  r2.nodes[nlen-1] = ksol
#===================================================================
###### the second part of KSL, from right to left
  for i in range(nlen-1,0,-1):

    phi1 = phia[i]
    phi2 = phia[i+1]

    if (i==nlen-1):
      ksol = r2.nodes[i].copy()
    else:
      ksol = r3.nodes[nlen-1-i].copy()

    #for islice in range(pa.rk4slice):
    ksol = update_v(yy=ksol, phi1=phi1, phi2=phi2, mat1=pall.nodes[i])
    r, q = split_rq(ksol)

    if (i==nlen-1):
      r3.nodes.append(q)
    else:
      r3.nodes[nlen-1-i] = q

    phi2 = ttf.phia_next(phi2,pall.nodes[i],q,q,1)
    phia[i] = phi2

    ssol = r

    #for islice in range(pa.rk4slice):
    ssol = update_v(yy=ssol, phi1=phi1, phi2=phi2)

    rtmp1 = np.tensordot(r2.nodes[i-1],ssol,axes=((2),(0)))
    r3.nodes.append(rtmp1)

  r3.nodes.reverse()
  ### the left most matrix
  phi1 = phia[0]
  phi2 = phia[1]
  ksol = r3.nodes[0].copy()
  ksol = update_v(yy=ksol, phi1=phi1, phi2=phi2, mat1=pall.nodes[0])

  r3.nodes[0] = ksol

  return r3