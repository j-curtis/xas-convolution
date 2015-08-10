# xas-convolution
## Code for computing the convolution of the single particle XAS with the XPS

This will contain code written in C++ that is designed to perform the convolution of the spectral function (XPS) and the XAS, both as functions of frequency. 

To use this program, one must first write a parameter file (see sample_parameters.txt for more details). 

The program is expecting that the data files be in the format of one linea per data point with columns delimited by spaces.

The XAS file should have 2 columns. The first is the frequency and the second is the value at that frequency.

The XPS file should have 3 columns. The first is the frequency, the second is the real part of the XPS, and the third column is the imaginary part (the spectral function).

The program xas-convolve is expecting that the spectral function have it's **main peak centered at zero with satellite peaks located at higher energies**. 

If this is not the case, run the python shell script flip-sign.py on the data first. 

The usage of flip-sign.py is "./flip-sign.py inputxps.dat outputxps.dat"

Once this is done, run the program on the parameter file with the XPSFILE flag set to the name of the output xps file.
