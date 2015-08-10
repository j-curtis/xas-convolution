# xas-convolution
Code for computing the convolution of the single particle XAS with the XPS

This will contain code written in C++ that is designed to perform the convolution of the spectral function (XPS) and the XAS, both as functions of frequency. 

It will scan a parameter file for the parameters of the convolution as well as the names of the data files containing the XAS and XPS. 

It will then read in these files, with the XAS being a two column data file of frequency and absorbtion, delimited by a space.

The XPS is a three column file with the frequency, then the real and imaginary parts as the second and third columns, also delimted by a space. We are interested in the third column. 


It will then print the ouput to the specified output file with the first column the frequency, the second column the convolution, the third column the XPS for that frequency, and the fourth (last) column being the XAS at that frequency.

The code implements a linear interpolation of the input data onto a single frequency array and then convolutes over this array.
Currently, the way the array is generated is by using user input for the size (number of samples) and then making the largest possible array out of the user given bounds, and the bounds extracted from the frequencies of each input data set. 

Currently, there is no algorithm implemented to determine a good range, so it is somewhat arbitrary and must be user tuned to get a good plot range. 

It also does not implement trapezoidal rule and may also be prone to Gibbs phenomenon near large jumps in the spectrum. 
