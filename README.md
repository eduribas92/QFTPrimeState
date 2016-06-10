# Quantum Fourier Transform of a Prime State
Simulation of a quantum algorithm to compute the Quantum Fourier Transform (QFT) of a Prime State. The result of the QFT gives information about Prime Numbers: number of primes, biases and ditribution.

Final Bachelor Project, Faculty of Physics, University of Barcelona, Spain

## Usage

On a Python Shell:

  `import QFTPrimeState as QFTPS`
  
Call the function:

  `QFTPS.QFT_Simulation(1000, 10, True)`

where:
* `N`: only the primes less than *N* will be considered
* `qubits`: if `True`, instead of *N*, it will consider the first greater power of *2*
* `num_peaks`: sets the number of highests peaks whose numerical value will be given
