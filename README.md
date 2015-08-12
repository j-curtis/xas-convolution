# xas-convolution
## Code for computing the convolution of the single particle XAS with the XPS

This will contain code written in C++ that is designed to perform the convolution of the spectral function (XPS) and the XAS, both as functions of frequency. 

To use this program, one must first write a parameter file (see sample_parameters.txt for more details). 

The program is expecting that the data files be in the format of one line per data point with columns delimited by spaces.

The XAS file should have 2 columns. The first is the frequency and the second is the value at that frequency.

The XPS file should have 3 columns. The first is the frequency, the second is the real part of the XPS, and the third column is the imaginary part (the spectral function).

## Usage
The program simply accepts a parameterfile with the name of the input and output files. From there it will automatically detect the size of the arrays for the data sets. It will also automatically perform peak shifting so that:
+ The XAS spectrum has it's first peak occuring at zero
+ The XPS spectrum has it's first (main) peak occuring at zero
+ THe XPS spectrum has it's satellites to the right (positive frequencies) of the main (first) peak

Once it computes the convolution, it will find the location of the new first peak and shift this back to zero. After that, it will print to the output file the (shifted) frequency in the first column, the convolution in the second column, the interpolated XPS in the third column, and the interpolated XAS in the fourth column. It is delimited by spaces. 
