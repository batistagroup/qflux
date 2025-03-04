# Utility Functions/Patches:
from qiskit import transpile


def execute(QCircuit, backend=None, shots=None):
    '''
        Function to replace the now-deprecated Qiskit
        `QuantumCircuit.execute()` method.

        Input:
          - `QCircuit`: qiskit.QuantumCircuit object
          - `Backend`: qiskit.Backend instance
          - `shots`: int specifying the number of shots
    '''
    # Transpile circuit with statevector backend
    tmp_circuit = transpile(QCircuit, backend)
    # Run the transpiled circuit
    if shots:
        job = backend.run(tmp_circuit, n_shots=shots)
    else:
        job = backend.run(tmp_circuit)
    return(job)
