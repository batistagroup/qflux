'''
MIT License

Copyright (c) 2024 Saurabh Shivpuje

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
'''

import numpy as np
from itertools import combinations, product
from scipy.linalg import expm, kron

# Ansatz Class and Operators

class Ansatz_class:
    def __init__(self, nqbit, u0, relrcut, pool, theta=[], ansatz=[]):
        self.theta = np.array(theta)
        self.A = ansatz
        self.state = u0.copy()
        self.ref = u0.copy()
        self.relrcut = relrcut
        self.nqbit = nqbit
        self.pool = pool

class AnsatzOperatorBase:
    pass

class PauliOperator_class(AnsatzOperatorBase):
    def __init__(self, mat, tag, nqbit):
        self.mat = mat
        self.tag = tag
        self.nqbit = nqbit

def single_clause(ops, q_ind, weight, num_qubit):
    si = np.eye(2)
    res = weight * np.eye(1)
    for i in range(1, num_qubit + 1):
        if i in q_ind:
            op2 = eval(f"{ops[q_ind.index(i)]}")
            res = kron(res, op2)
        else:
            res = kron(res, si)
    return res

def PauliOperator(ops, idx, w, nqbit):
    sortedRep = sorted(zip(ops, idx), key=lambda x: x[1])
    tag = ''.join([f"{o}{i}" for o, i in sortedRep])
    mat = single_clause(ops, idx, w, nqbit)
    return PauliOperator_class(mat, tag, nqbit)

# Defining pool of operators/gates
sx = np.array([[0, 1], [1, 0]])
sy = np.array([[0, -1j], [1j, 0]])
sz = np.array([[1, 0], [0, -1]])
# S = np.array([[1, 0], [0, 1j]])


def build_pool(nqbit):
    pauliStr = ["sx", "sz", "sy"]
    res = []
    for order in range(1, nqbit+1):
        for idx in combinations(range(1, nqbit + 1), order):
            for op in product(pauliStr, repeat=order):
                res.append(PauliOperator(op, list(idx), 1, nqbit))
    return res

def Ansatz(u0, relrcut, theta=[], ansatz=[]):
    nqbit = int(np.log2(len(u0)))
    pool_qubit =  nqbit
    # u0 = np.outer(u0, u0).flatten()
    pool = build_pool(pool_qubit)
    return Ansatz_class(nqbit, u0, relrcut, pool, theta, ansatz)
