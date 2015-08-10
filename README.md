# xas-convolution
Code for computing the convolution of the single particle XAS with the XPS

This will contain code written in C++ that is designed to perform the convolution of the spectral function (XPS) and the XAS, both as functions of frequency. 

It will scan a parameter file for the parameters of the convolution as well as the names of the data files containing the XAS and XPS. 

It will then read in these files, with the XAS being a two column data file of frequency and absorbtion, delimited by a space.

The XPS is a three column file with the frequency, then the real and imaginary parts as the second and third columns, also delimted by a space. We are interested in the third column. 
