import qflux

def test_imports():
    from qflux.closed_systems import DynamicsCS
    from qflux.variational_methods import qmad
    from qflux.variational_methods.qmad import effh, ansatz, solver
    from qflux.error_mitigation.rbm import RBM
    from qflux.error_mitigation import mitigation

    assert qflux is not None
    assert DynamicsCS is not None
    assert qmad is not None
    assert effh is not None
    assert ansatz is not None
    assert solver is not None
    assert RBM is not None
    assert mitigation is not None