# xas-convolution
## Code for computing the convolution of the single particle XAS with the XPS

This will contain code written in C++ that is designed to perform the convolution of the spectral function (XPS) and the XAS, both as functions of frequency. 

To use this program, one must first write a parameter file (see sample_parameters.txt for more details). 

The program is expecting that the data files be in the format of one line per data point with fields written in columns delimited by spaces.

The XAS file should have 2 columns. The first is the frequency and the second is the value at that frequency.

The SPEC file should have 3 columns. The first is the frequency, the second is the spectral density function, and the third column is the spectral cumulative function.

## Usage
The program simply accepts a parameterfile with the name of the input and output files. From there it will automatically detect the size of the arrays for the data sets. It will also automatically perform peak shifting so that the XAS spectrum has it's first peak occuring at zero. It expects the spectral function to already have its main peak at zero, with the satellites to the right (higher energies).

The output file will have 5 columns, the first is frequency, the second is the convolution, the third is the interpolated XAS, the fourth column is the interpolated spectral density. The fifth column is the normalization as a function of frequency.

It will also print out useful diagnostic information as it runs, so you may want to dump the standard output into a log file.

## Operation
The operation that this file performs is a modified convolution. Let the XAS be mu(w) and the spectral function be A(w). Then this will compute (and output)

MU(w) = 1/(integral_{-infintiy}^{w-Efermi} A(w')dw' ) integral_{-infinity}^{w-Efermi} A(w')mu(w-w')dw'

where Efermi is the Fermi energy of the XAS.
