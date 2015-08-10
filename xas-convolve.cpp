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
  //f(x) = f1 + f2*(x-x1)/(x2-x1) if x1<x<x2
  //f(x) = f2 if x>=x2
  if(x<=x1){
    return f1;
  }
  else if(x>=x2){
    return f2;
  }
  else{
    return f1+f2*(x-x1)/(x2-x1);  //Note this is not safe against x2=x1==> division by zero error
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

//Now the main method
int main(int argc, char * argv []){
  //We check to make sure they have the right number of arguments
  if(argc!=2){
    cout<<"Incorrect usage. Usage is: 'program_name' INPUTFILE "<<endl;
    return 0;
  }
  
  //Now we define the variables we will use 
  double min_w = 0.0; //The minimum frequency we compute the convolution for
  double max_w = 0.0; //The maximum frequency we compute the convolution for
  int num_w_steps = 3;//The number of points we compute the convolution for
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
  string inline;
  while(getline(infile,inline)){
    //We parse for tokens of label, values pairs
    string label, value;
    istringstream ss(inline); //This is what actually breaks the line into tokens
    ss>>label>>value;
    //We use the variable name as a label 
    parseStr(label,value,"XPSFILE",XPSFILE);
    parseStr(label,value,"XASFILE",XASFILE);
    parseStr(label,value,"OUTFILE",OUTFILE);
    
    parseInt(label,value,"num_w_steps",num_w_steps);
    parseDbl(label,value,"min_w",min_w);
    parseDbl(label,value,"max_w",max_w);
    
    parseInt(label,value,"num_xps_steps",num_xps_steps);
    parseInt(label,value,"num_xas_steps",num_xas_steps);
  }
  infile.close();
  
  //Now we do value checking 
  if(num_w_steps <= 2){
    cout<<"Error: insufficient number of output frequency steps (<3)."<<endl;
    return 0;
  }
  if(num_xps_steps <= 2){
    cout<<"Error: insufficient number of input xps steps (<3)."<<endl;
    return 0;
  }
  if(num_xas_steps <= 2){
    cout<<"Error: insufficient number of input xas steps (<3)."<<endl;
    return 0;
  }
  if(max_w<=min_w){
    cout<<"Error: max_w <= min_w is not permissible."<<endl;
    return 0;
  }
  
  //Allocation
  xps_freqs = new double[num_xps_steps];
  xps = new double[num_xps_steps];
  xas_freqs = new double[num_xas_steps];
  xas = new double[num_xas_steps];
  out_freqs = genGrid(num_w_steps,min_w,max_w); //Generate a uniform list of frequencies from scratch
  out = new double[num_w_steps];  //This will be generated from the convolution
  
  //Open up XPS and XAS files and read the data in
  ifstream xpsfile(XPSFILE);
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
    xpscounter++; //Increment the counter
  }
  xpsfile.close();  //Close the file
  
  ifstream xasfile(XASFILE);
  string xasline;
  int xascount = 0; //To increment through each line of the array that we are filling
  while(getline(xasfile,xasline)){
    //We now parse for 2 tokens, each a double 
    string token1, token2;
    istringstream ss(xasline);  //Extracts each token
    ss>>token1>>token2;
    xas_freqs[xascount] = atof(token1.c_str()); //We must also convert it to a double 
    xas[xascount] = atof(token2.c_str()); //We must also convert this to a double 
    xascounter++; //Increment the counter
  }
  xasfile.close();  //Close the file
  
  //Now all the data arrays have been filled, we may perform the convolution
  //We use a double loop to implement the integral
  //We do NOT use trapezoidal rule
  //We precompute delta_w for convinience
  //We also use this loop to output the results to the output file, to save on the number of loops
  ofstream outfile(OUTFILE);
  
  double delta_w = out_freqs[1]-out_freqs[0];
  
  for(int i = 0; i < num_w_steps; i++){
    out[i] = 0.0; //We initialize to zero
    //Now we use the formula
    //out(w) = int_{wmin}^{wmax} dw' xps(w')xas(w-w')
    //We use linear interpolation to comute the integrands at each integration point 
    for(int j = 0; j < num_w_steps; j++){
      //We compute the snapped to grid value of j for the xas and xps grids 
      //We also need i snapped to the xas grid (to compute xas(w-w') )
      double w = out_freqs[i];  //The value of w
      double u - out_freqs[j];  //The value of w'
      
      int j_xps = snapToGrid(u,num_xps_steps,xps_freqs); //The value of j snapped to the xps grid
      //Thus, j_xps is between 0 and num_xps_steps - 1
      int i_xas = snapToGrid(w-u,num_xas_steps,xas_freqs); //The value of i-j snapped to the xas grid
      //Thus, this takes values between 0 and num_xas_steps
      
      //Now we compute the linearly extrapolated xas and xps for w and u
      //We have to make sure we don't go out of array bounds when computing the right endpoint for the linear extrapolations
      double xps_term = 0.0;
      double xas_term = 0.0;
      
      if(j_xps < num_xps_steps-2){
        xps_term = linearInterpolate(u,xps_freqs[j_xps],xps_freqs[j_xps+1],xps[j_xps],xps[j_xps+1]);
      } 
      else{
        xps_term = linearInterpolate(u,xps_freqs[j_xps],xps_freqs[j_xps],xps[j_xps],xps[j_xps]);
        //Note that is j_xps is always returned to be, at most num_xps_steps - 1, then when x = x2 = x1 it will return f2
      }
      
      //Now we must do the same for the xas term
      if(i_xas < num_xas_steps - 2){
        xas_term = linearInterpolate(w-u,xas_freqs[i_xas],xas_freqs[i_xas+1],xas[i_xas],xas[i_xas+1]);
      }
      else{
        xas_term = linearInterpolate(w-u,xas_freqs[i_xas],xas_freqs[i_xas],xas[i_xas],xas[i_xas]);
      }
      
      out[i] += delta_w *xas_term*xps_term; //Add to the integral
    }
    
    //And we write this term to a file
    outfile<<out_freq[i]<<" "<<out[i]<<endl;
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
  
  return 0;
}



