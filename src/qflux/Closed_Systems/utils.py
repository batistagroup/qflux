# Utilities
# Plotting utilities, conversion factors, etc.

# Conversion Factors:
def convert_au_to_eV(input_val):
    au2ev = 27.21138602
    return(au2ev * input_val)


def convert_eV_to_au(input_val):
    au2ev = 27.21138602
    ev2au = 1/au2ev
    return(ev2au * input_val)


def convert_bohr_to_au(input_val):
    bohr2au = 0.52917721092
    return(bohr2au * input_val)


def convert_au_to_bohr(input_val):
    bohr2au = 0.52917721092
    au2bohr = 1/bohr2au
    return(input_val * au2bohr)


def convert_fs_to_au(input_val):
    fs2au = 41.3414
    return(fs2au * input_val)


def convert_au_to_fs(input_val):
    fs2au = 41.3414
    au2fs = 1/fs2au
    return(au2fs * input_val)


def get_proton_mass():
    proton_mass = 1836.15267343 # proton-electron mass ratio
    return(proton_mass)


def execute(self, QCircuit, backend=None, shots=None):
    '''
        Function to replace the now-deprecated Qiskit
        `QuantumCircuit.execute()` method.

        Input:
          - `QCircuit`: qiskit.QuantumCircuit object
          - `Backend`: qiskit.Backend instance
          - `shots`: int specifying the number of shots
          - `real_backend`: bool specifying whether the provided backend is
            a real device (True) or not (False)
    '''
    if shots:
        n_shots = shots
    else:
        n_shots = 1024 # Use the qiskit default if not specified
    backend_type = type(backend)
    sv_type = qiskit_aer.backends.statevector_simulator.StatevectorSimulator
    if backend_type == sv_type:
        real_backend = False
    else:
        real_backend = True

    if real_backend:
        QCircuit.measure_all()
        qc = transpile(QCircuit, backend=backend)
        sampler = Sampler(backend)
        job = sampler.run([qc], shots=n_shots)
    else:
        # Transpile circuit with statevector backend
        tmp_circuit = transpile(QCircuit, backend)
        # Run the transpiled circuit
        job = backend.run(tmp_circuit, n_shots=shots)
    return(job)
