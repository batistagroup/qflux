import numpy as np
from . import params as pa
from typing import List, Tuple

def output_superoper_array(
    time: List[float],
    s_array: np.ndarray,
    prefix: str
    ) -> None:
    """
    Write the time-dependent superoperator array to disk.

    Args:
        time (List[float]): List of time values.
        s_array (np.ndarray): Superoperator array (shape (T, N^2, N^2)).
        prefix (str): Output file name prefix.
    """
    Nlen = len(time)
    for j in range(pa.DOF_E_SQ):
        a, b = divmod(j, pa.DOF_E)
        for k in range(pa.DOF_E_SQ):
            c, d = divmod(k, pa.DOF_E)

            filename = f"{prefix}{a}{b}{c}{d}{pa.PARAM_STR}.txt"
            with open(filename, "w") as f:
                for i in range(Nlen):
                    real_part = s_array[i, j, k].real
                    imag_part = s_array[i, j, k].imag
                    f.write(f"{time[i]}\t{real_part}\t{imag_part}\n")


def read_superoper_array(
    Nlen: int,
    prefix: str
    ) -> Tuple[np.ndarray, np.ndarray]:
    """
    Read the time-dependent superoperator array from disk.

    Args:
        Nlen (int): Number of time steps.
        prefix (str): Input file name prefix.

    Returns:
        Tuple[np.ndarray, np.ndarray]: 
            - time (array of shape (Nlen,))
            - S_array (array of shape (Nlen, N^2, N^2))
    """
    S_array = np.zeros((Nlen, pa.DOF_E_SQ, pa.DOF_E_SQ), dtype=np.complex128)
    time = np.zeros(Nlen, dtype=np.float64)

    for j in range(pa.DOF_E_SQ):
        a, b = divmod(j, pa.DOF_E)
        for k in range(pa.DOF_E_SQ):
            c, d = divmod(k, pa.DOF_E)

            filename = f"{prefix}{a}{b}{c}{d}{pa.PARAM_STR}.txt"
            data = np.loadtxt(filename)
            time_read, real_part, imag_part = np.hsplit(data, 3)

            for i in range(Nlen):
                time[i] = time_read[i]
                S_array[i, j, k] = real_part[i] + 1j * imag_part[i]

    return time, S_array


def output_operator_array(
    time: List[float],
    sigma: np.ndarray,
    prefix: str
    ) -> None:
    """
    Write the time-dependent vectorized operator array to disk.

    Args:
        time (List[float]): List of time values.
        sigma (np.ndarray): Operator array (shape (T, N^2)).
        prefix (str): Output file name prefix.
    """
    for j in range(pa.DOF_E_SQ):
        a, b = divmod(j, pa.DOF_E)
        filename = f"{prefix}{a}{b}{pa.PARAM_STR}.txt"
        with open(filename, "w") as f:
            for i in range(len(time)):
                real_part = sigma[i, j].real
                imag_part = sigma[i, j].imag
                f.write(f"{time[i]}\t{real_part}\t{imag_part}\n")


def read_operator_array(
    Nlen: int,
    prefix: str
    ) -> Tuple[np.ndarray, np.ndarray]:
    """
    Read the time-dependent vectorized operator array from disk.

    Args:
        Nlen (int): Number of time steps.
        prefix (str): Input file name prefix.

    Returns:
        Tuple[np.ndarray, np.ndarray]: 
            - time (array of shape (Nlen,))
            - O_array (array of shape (Nlen, N^2))
    """
    O_array = np.zeros((Nlen, pa.DOF_E_SQ), dtype=np.complex128)
    time = np.zeros(Nlen, dtype=np.float64)

    for j in range(pa.DOF_E_SQ):
        a, b = divmod(j, pa.DOF_E)
        filename = f"{prefix}{a}{b}{pa.PARAM_STR}.txt"
        data = np.loadtxt(filename)
        time_read, real_part, imag_part = np.hsplit(data, 3)

        for i in range(Nlen):
            time[i] = time_read[i]
            O_array[i, j] = real_part[i] + 1j * imag_part[i]

    return time, O_array
