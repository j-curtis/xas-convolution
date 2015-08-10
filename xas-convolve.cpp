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
  //f(x) = f1 + f2*(x-x1)/(x2-x1) if x1<x<=x2
  //f(x) = f2 if x>x2
  if(x<=x1){
    return f1;
  }
  else if(x>x2){
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
}



