//This is the main file for the convolution script
//Jonathan Curtis
//08/09/15

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <string>
#include <cmath>
#include <complex>
#include <cstdlib>
using namespace std;

//These are methods for parsing the data file
//They will evaluate label=="flag" and if it does, will assign variable to value (appropriately typed).
//Each one parses for a different data type. 
void parseStr(string label, string value, string flag, string & variable){
  //This parses for a string type variable
  if(label == flag){
    variable = value;
  }
}
void parseDbl(string label, string value, string flag, double & variable){
  //This parses for a double type variable
  if(label == flag){
    variable = atof(value.c_str() );  //We convert the value to a double
  }
}
void parseInt(string label, string value, string flag, int & variable){
  //This parses for a double type variable
  if(label == flag){
    variable = atoi(value.c_str() );  //We convert the value to a double
  }
}

//This function will linearly interpolate between the two values 
//NOTE: this is not safe against division by zero resulting from x2 = x1
double linearInterpolate(double x, double x1, double f1, double x2, double f2){
  //Given x, x1, x2, f1, f2, this will interpolate f(x) as:
  //f(x) = f1 if x<=x1
  //f(x) = f1 + (f2-f1)*(x-x1)/(x2-x1) if x1<x<x2
  //f(x) = f2 if x>=x2
  if(x<=x1){
    return f1;
  }
  else if(x>=x2){
    return f2;
  }
  else{
    return f1+(f2-f1)*(x-x1)/(x2-x1);  //Note this is not safe against x2=x1==> division by zero error
  }
}

//This function will snap a given point onto a grid 
//It will return the largest array index that has a value less than the passed value
//If the element is out of bounds, it will return 0 or arrysize-1 depending on which end of the array it crosses
//It assumes that the array is uniformly spaced
int snapToGrid(double x, int size, double * arry){
  //First we deal with the case that x<values[0]
  if(x<arry[0]){
    return 0;
  }
  //Next, if it is in the array range
  else if(x<arry[size-1]){
    return floor( (x-arry[0])/(arry[1]-arry[0]) ); //We assume the array is uniform and has spacing values[1]-values[0]
  }
  else{
    return size-1;  //It is larger than the largest entry so we return the largest index
  }
}

//This function generates an array of N uniformly spaced values from min to max
//It returns a pointer to the first element of the array
//The returned array must be freed at the end 
double * genGrid(int Num, double min, double max){
  double * val = new double[Num];
  double delta = (max-min)/double(Num-1);
  
  for(int i = 0; i < Num; i++){
    val[i] = min + double(i)*delta;
  }
  
  return val;
}

//This function will locate the first local max in an array using a simple "walkthrough" technique
//It will return the index of the first max. If there is no local max, it will return 0
//To filter out roundoff error and small fluctuations, we require it to be larger than MULTIPLE neighbors
int locateFirstMax(int size, double * values){
  for(int i = 2; i < size-2; i++){
    //We step through and check to see if an element is greater than all of its neighbors
    //To minimize floating point erros we will multiply each element by 10^12 and then compare the results cast as integers
    //We do this by multiplying by 10000 and then casting as an int and then comparing the integers
    bool nneighbor_check = ( ( int(1e12*values[i]) > int(1e12*values[i-1]) ) && ( int(1e12*values[i]) > int(1e12*values[i+1]) ) );
    bool nnneighbor_check = ( ( int(1e12*values[i]) > int(1e12*values[i-2]) ) && ( int(1e12*values[i]) > int(1e12*values[i+2]) ) );

    if(nneighbor_check && nnneighbor_check){
      return i;
    } 
  }
  return 0;
}


//Now the main method
int main(int argc, char * argv []){
  //We check to make sure they have the right number of arguments
  if(argc!=2){
    cout<<"Incorrect usage. Usage is: 'program_name' INPUTFILE "<<endl;
    return 0;
  }
  
  //Now we define the variables we will use 
  int num_w_steps;//The number of points we compute the convolution for
  int num_xps_steps;  //The number of points in the xps array
  int num_xas_steps;  //The number of points in the xas array
  
  //These are the arrays for the input and output data
  double * xps_freqs; //The frequencies for the input xps 
  double * xps;   //The input xps values for each frequency
  double * xas_freqs; //The frequencies for the input xas
  double * xas;   //The input xas values for each frequency
  double * out_freqs; //The output frequencies
  double * out;   //The output values at each frequency

  string XPSFILE; //The file where the XPS is input from
  string XASFILE; //The file where the XAS is input from
  string OUTFILE; //The file where the output is written to
  
  //Now we parse for these values
  ifstream infile(argv[1]);
  string infileline;
  while(getline(infile,infileline)){
    //We parse for tokens of label, values pairs
    string label, value;
    istringstream ss(infileline); //This is what actually breaks the line into tokens
    ss>>label>>value;
    //We use the variable name as a label 
    parseStr(label,value,"XPSFILE",XPSFILE);
    parseStr(label,value,"XASFILE",XASFILE);
    parseStr(label,value,"OUTFILE",OUTFILE);
    
    //parseInt(label,value,"num_w_steps",num_w_steps); 
    //parseDbl(label,value,"min_w",min_w);
    //parseDbl(label,value,"max_w",max_w);   
    parseInt(label,value,"num_xps_steps",num_xps_steps);
    parseInt(label,value,"num_xas_steps",num_xas_steps);
  }
  infile.close();
  
  //Allocation of data arrays
  xps_freqs = new double[num_xps_steps];
  xps = new double[num_xps_steps];

  xas_freqs = new double[num_xas_steps];
  xas = new double[num_xas_steps];

  //Open up XPS and XAS files and read the data in
  ifstream xpsfile(XPSFILE.c_str());
  string xpsline;
  int xpscount = 0; //To increment through each line of the array that we are filling
  while(getline(xpsfile,xpsline)){
    //We now parse for 3 tokens, each a double 
    string token1, token2, token3;
    istringstream ss(xpsline);  //Extracts each token
    ss>>token1>>token2>>token3;
    //We only want the first token (frequency), and the third token (xps)
    xps_freqs[xpscount] = atof(token1.c_str()); //We must also convert it to a double 
    xps[xpscount] = atof(token3.c_str()); //We must also convert this to a double 
    xpscount++; //Increment the counter
  }
  xpsfile.close();  //Close the file
  
  ifstream xasfile(XASFILE.c_str());
  string xasline;
  int xascount = 0; //To increment through each line of the array that we are filling
  while(getline(xasfile,xasline)){
    //We now parse for 2 tokens, each a double 
    string token1, token2;
    istringstream ss(xasline);  //Extracts each token
    ss>>token1>>token2;
    xas_freqs[xascount] = atof(token1.c_str()); //We must also convert it to a double 
    xas[xascount] = atof(token2.c_str()); //We must also convert this to a double 
    xascount++; //Increment the counter
  }
  xasfile.close();  //Close the file

  //Now we need to flip the xps and center it and center the xas as well
  //First we locate the first peak of the xas
  int xas_fp = locateFirstMax(num_xas_steps,xas);
  double xas_fpw = xas_freqs[xas_fp];
  cout<<xas_fpw<<endl;
  //We subtract the peak frequency from each frequency
  for(int i = 0; i < num_xas_steps; i++){
    xas_freqs[i] -= xas_fpw;
  }

  //Now we reflect the xps about zero (so that the main peak is also the first peak)
  //To do this, we first reflect the frequencies by multipliying by -1
  //Then we re-order both the xps and xps_freqs array so that they read in the opposite order as they do now 
  for(int i = 0; i < num_xps_steps; i++){
    xps_freqs[i] *= -1.0;
  }
  //Now we reorder
  for(int i = 0; i < int(num_xps_steps/2.0); i++){
    double tmp_w = xps_freqs[i]; //Used to store the frequency of the element we are swapping out
    double tmp_x = xps[i]; //Used to store the xps of the element we are swapping out

    xps_freqs[i] = xps_freqs[num_xps_steps-1-i];  
    xps[i] = xps[num_xps_steps-1-i];  

    //Finally, we replace the last values with the temp values
    xps_freqs[num_xps_steps-1-i] = tmp_w;
    xps[num_xps_steps-1-i] = tmp_x;
  }

  //Now we locate the max and subtract it off
  int xps_fp = locateFirstMax(num_xps_steps,xps);
  double xps_fpw = xps_freqs[xps_fp];
  cout<<xps_fpw<<endl;
  //We subtract off using a loop
  for(int i = 0; i < num_xps_steps; i++){
    xps_freqs[i] -= xps_fpw;
  }
  //Now both functions should have their first peak at zero
  
  //Next we must linearly interpolate them onto the same frequency arrays 
  //We interpolate the XAS onto the XPS grid
  //First we allocate the output arrays

  num_w_steps = num_xps_steps;  //We use the same size array, but this new variable will make it more clear when we are looping

  out_freqs = new double[num_w_steps];
  out = new double[num_w_steps];

  for(int i = 0; i < num_w_steps; i++){
    out_freqs[i] = xps_freqs[i];
  }

  //Now we snap the xas onto this grid 
  //We use linear interpolation
  double * xas_snapped = new double[num_w_steps];
 
  for(int i = 0; i < num_w_steps; i++){
    //First we compute the snapped to grid index for the xas and xps 
    double w = out_freqs[i];  //The value of the frequency we are trying to interpolate to
    int xas_index = snapToGrid(w,num_xas_steps,xas_freqs);  //The corresponding snapped index on the xas grid
    //Can range from 0 to num_xas_steps-1

    //Repeat for the xas
    if(xas_index != num_xas_steps-1){
      xas_snapped[i] = linearInterpolate(w,xas_freqs[xas_index],xas[xas_index],xas_freqs[xas_index+1],xas[xas_index+1]);
    }
    else{
      xas_snapped[i] = xas[num_xas_steps-1];
    }

  }

  //Now all the data arrays have been filled, we may perform the convolution
  //We use a double loop to implement the integral
  //We do NOT use trapezoidal rule
  //Now we use the formula
  //out(w) = int_{wmin}^{wmax} du' xps(u)xas(w-u)
  //Because we have already interpolated xas and xps onto the same grid, we simply implement discrete convolution of 
  //out[i] = sum_j delta xps[j]*xas[i-j]

  double delta_w = out_freqs[1] - out_freqs[0];
  for(int i = 0; i < num_w_steps; i++){
    out[i] = 0.0; //We initialize to zero
    for(int j = 0; j < num_w_steps; j++){
      //We only add if i - j > 0. Otherwise we approximate xas[- |i-j| ] = xas[0]
      if((i-j)<0){
        out[i] += delta_w * xps[j] * xas_snapped[0];
      }
      else{
        out[i] += delta_w * xps[j] * xas_snapped[i-j];
      }
    }
  }

  //Now we find the location of the first peak and we subtract off the frequency
  int peak_index = locateFirstMax(num_w_steps,out);
 
  double peak_freq = out_freqs[peak_index]; //This is where the peak occurs
  //We use this loop to output the results to the output file, to save on the number of loops
  //We also use this loop to subtract off the first peaks location
  ofstream outfile(OUTFILE.c_str());
  for(int i = 0; i < num_w_steps; i++){
    //We also print the snapped xas and xps
    outfile<<out_freqs[i]-peak_freq<<" "<<out[i]<<" "<<xps[i]<<" "<<xas_snapped[i]<<endl;
  }
  //And we close the output file
  outfile.close();
  
  //Delete allocated memory
  delete [] xps_freqs;
  delete [] xps;
  delete [] xas_freqs;
  delete [] xas;
  delete [] out_freqs;
  delete [] out;
  delete [] xas_snapped;
  
  return 0;
}
