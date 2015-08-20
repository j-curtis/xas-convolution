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
  int num_spec_steps;  //The number of points in the xps array
  int num_xas_steps;  //The number of points in the xas array
  
  //These are the arrays for the input and output data
  double * spec_freqs; //The frequencies for the input spectral function 
  double * spec_den;   //The input spectral function density values for each frequency
  double * spec_cum;
  double * spec_den_snapped;  //The spectral density snapped to grid
  double * spec_cum_snapped;  //The spectral cumulative snapped to grid   
  
  double * xas_freqs; //The frequencies for the input xas
  double * xas;   //The input xas values for each frequency
  double * xas_snapped; //This will be the xas interpolated onto the output grid size
  
  double * out_freqs; //The output frequencies
  double * out;   //The output values at each frequency

  string SPECFILE; //The file where the spectral data is input from
  string XASFILE; //The file where the XAS is input from
  string OUTFILE; //The file where the output is written to
  
  cout<<"Parsing "<<argv[1]<<endl;

  //Now we parse for these values
  ifstream infile(argv[1]);
  string infileline;
  while(getline(infile,infileline)){
    //We parse for tokens of label, values pairs
    string label, value;
    istringstream ss(infileline); //This is what actually breaks the line into tokens
    ss>>label>>value;
    //We use the variable name as a label 
    parseStr(label,value,"SPECFILE",SPECFILE);
    parseStr(label,value,"XASFILE",XASFILE);
    parseStr(label,value,"OUTFILE",OUTFILE);
  }
  infile.close();

  cout<<"Done"<<endl;
  cout<<"Determining size of files"<<endl;

  //Automatic detection of array sizes
  //This increments a counter for every line in the file, which will yield the array size for us
  ifstream specsize(SPECFILE.c_str());
  string line;
  int speclinecounter = 0;
  while(getline(specsize,line)){
    speclinecounter++;
  }
  specsize.close();

  num_spec_steps = speclinecounter;

  //Same for xas
  //We use the same line string to save on memory and variable number
  ifstream xassize(XASFILE.c_str());
  int xaslinecounter = 0;
  while(getline(xassize,line)){
    xaslinecounter++;
  }
  xassize.close();

  num_xas_steps = xaslinecounter;

  cout<<"Done"<<endl;
  cout<<"SPECFILE has "<<num_spec_steps<<" entries"<<endl;
  cout<<"XASFILE has "<<num_xas_steps<<" entries"<<endl;
  
  //Allocation of data arrays
  spec_freqs = new double[num_spec_steps];
  spec_den = new double[num_spec_steps];
  spec_cum = new double[num_spec_steps];

  xas_freqs = new double[num_xas_steps];
  xas = new double[num_xas_steps];

  cout<<"Reading in from "<<SPECFILE<<endl;

  //Open up SPEC and XAS files and read the data in
  ifstream specfile(SPECFILE.c_str());
  string specline;
  int speccount = 0; //To increment through each line of the array that we are filling
  while(getline(specfile,specline)){
    //We now parse for 3 tokens, each a double 
    string token1, token2, token3;
    istringstream ss(specline);  //Extracts each token
    ss>>token1>>token2>>token3;
    //First token is frequency, second is density, third is cumulative
    spec_freqs[speccount] = atof(token1.c_str());
    spec_den[speccount] = atof(token2.c_str());
    spec_cum[speccount] = atof(token3.c_str());
    
    speccount++; //Increment the counter
  }
  specfile.close();  //Close the file

  cout<<"Done"<<endl;
  cout<<"Reading in from "<<XASFILE<<endl;
  
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

  cout<<"Done"<<endl;
  cout<<"Interpolating"<<endl;

  //Now we need to center the xas 
  int xas_fp = locateFirstMax(num_xas_steps,xas);
  double xas_fpw = xas_freqs[xas_fp];
  //We subtract the peak frequency from each frequency
  for(int i = 0; i < num_xas_steps; i++){
    xas_freqs[i] -= xas_fpw;
  }
  //Now both functions should have their first peak at zero

  //Next we must linearly interpolate them onto the same frequency arrays 
  //We use the relation that, if one has arrays 
  // h[n] = f[h]*g[n] 
  //With sizes Nh, Nf, and Ng, respectively, then 
  //Nh = Nf+Ng-1
  //In this case, we will be using Nf=Ng=num_spec_steps
  num_w_steps = 2*num_spec_steps-1;  //We use the same size array, but this new variable will make it more clear when we are looping

  out_freqs = new double[num_w_steps];
  out = new double[num_w_steps];

  //We generate the output frequencies by starting from spec_freqs[0] and using delta_w_spec to go 2*num_spec_steps-1
  double delta_w = spec_freqs[1] - spec_freqs[0];

  for(int i = 0; i < num_w_steps; i++){
    out_freqs[i] = spec_freqs[0] + double(i)*delta_w;
    out[i] = 0.0; //We also initialize the ouput array values to zero, just in case they aren't already.
  }


  //To implement this algorithm we snap the XAS and SPEC onto the output grid
  //For the spectral function this is just zero padding the end of the array (so it has the correct size)
  //We also padd the cumulative function with 1s where the spectral density is padded with zeros
  //For this XAS this is interpolating and then padding with the last value
  //We use linear interpolation when appropriate
  spec_den_snapped = new double[num_w_steps];
  spec_cum_snapped = new double[num_w_steps];
  xas_snapped = new double[num_w_steps];

  for(int i = 0; i < num_w_steps; i++){
    //First we compute the snapped to grid index for the xas
    double w = out_freqs[i];  //The value of the frequency we are trying to interpolate to
    int xas_index = snapToGrid(w,num_xas_steps,xas_freqs);  //The corresponding snapped index on the xas grid
    //Can range from 0 to num_xas_steps-1
    //Interpolate the XAS onto the value w
    if(xas_index != num_xas_steps-1){
      xas_snapped[i] = linearInterpolate(w,xas_freqs[xas_index],xas[xas_index],xas_freqs[xas_index+1],xas[xas_index+1]);
    }
    else{
      xas_snapped[i] = xas[num_xas_steps-1];
    }

    //Copy the XPS if it is in the original range, otherwise we just zero pad
    if(i<num_spec_steps){
      spec_den_snapped[i] = spec_den[i];
      spec_cum_snapped[i] = spec_cum[i];
    } 
    else{
      spec_den_snapped[i] = 0.0;
      spec_cum_snapped[i] = 1.0;
    }

  }

  cout<<"Done"<<endl;
  cout<<"Convolving"<<endl;

  //Now all the data arrays have been filled, we may perform the convolution
  //We will use trapezoidal rule
  //First we will need to determine the Fermi-energy
  //Then we will need to compute the function
  //MU(w) = 1/(1-spec_cum(Efermi-w))int_{Efermi-w}^infty A(w-u)mu(u)du
  //Frist we must compute the fermi energy
  //We will define the Fermi energy to be the first energy where the xas reaches a certain % of the first peak height
  double fermi_peak_percent = .1;  //The percent of the first peak height we use to determine the Fermi energy

  double xas_fph = xas[xas_fp]; //We don't need to use the snapped value because this won't be very much different from it

  double fermi_value = xas_fph * fermi_peak_percent;  //This is the value that we will compare against to determine the fermi enrgy

  int iFermi = 0; //The index of the fermi energy
  double wFermi = 0;  //The fermi energy
  double hFermi = 0;  //The height of the XAS at the Fermi energy

  for(int i = 1; i < num_w_steps-1; i++){
    if(xas_snapped[i] >= fermi_value){
      //We have reached the peak, we save the values and break the loop
      iFermi = i;
      wFermi = out_freqs[iFermi];
      hFermi = xas_snapped[iFermi];
      break;
    }
  }

  cout<<"Fermi index "<<iFermi<<endl;
  cout<<"Fermi energy "<<wFermi<<" eV"<<endl;
  cout<<"XAS(Fermi energy) "<<hFermi<<endl;

  cout<<"Fermi threshold percent "<<fermi_peak_percent<<endl;
  cout<<"Fermi threshold value "<<fermi_value<<endl;

  cout<<"Done"<<endl;
  cout<<"Printing"<<endl;

  ofstream outfile(OUTFILE.c_str());
  for(int i = 0; i < num_w_steps; i++){
    //We also print the snapped xas and xps
    outfile<<out_freqs[i]<<" "<<out[i]<<" "<<xas_snapped[i]<<" "<<spec_den_snapped[i]<<" "<<spec_cum_snapped[i]<<endl;
  }
  //And we close the output file
  outfile.close();
  
  cout<<"Finished."<<endl;

  //Delete allocated memory
  delete [] spec_freqs;
  delete [] spec_den;
  delete [] spec_den_snapped;
  delete [] spec_cum;
  delete [] spec_cum_snapped;

  delete [] xas_freqs;
  delete [] xas;
  delete [] xas_snapped;

  delete [] out_freqs;
  delete [] out;
  
  return 0;
}
