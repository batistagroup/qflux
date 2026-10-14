import qflux

from qflux.closed_systems import DynamicsCS
print("Checking Closed Systems...")
print(dir(DynamicsCS))

from qflux.variational_methods import qmad
print("Checking Variational Methods...")
print(dir(qmad))

print("Checking qmad modules...")

from qflux.variational_methods.qmad import effh, ansatz, solver
print("Methods from effh:")
print(dir(effh))

print("Methods from ansatz:")
print(dir(ansatz))

print("Methods from solver:")
print(dir(solver))

from qflux.error_mitigation.rbm import RBM
print("Checking RBM ...")
print(dir(RBM))

from qflux.error_mitigation import mitigation
print("Checking Second Module in Error Mitigation...")
print(dir(mitigation))

print("Testing citation output...")
qflux.cite()
