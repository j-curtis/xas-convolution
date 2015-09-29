//This is an updated and improved XAS convolution calculator
//Jonathan Curtis
//09/16/15

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <string>
#include <cmath>
#include <complex>
#include <cstdlib>
using namespace std;

//PARSING METHODS 
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

//This function counts the number of lines in a data file 
int lineCount(string FILENAME){
	//Open an ifstream
	ifstream file(FILENAME.c_str());
	//Start a counter
	int counter = 0;
	string line;	//To hold the value of each line

	while(getline(file,line)){
		counter++;
	}
	file.close();

	return counter;
}

//MAX FINDING FUNCTIONS
//This function will locate the first local max in an array using a simple "walkthrough" technique
//It will return the index of the first max. If there is no local max, it will return 0
//To filter out roundoff error and small fluctuations, we require it to be larger than MULTIPLE neighbors
int locateFirstMax(int size, double * values){
  for(int i = 2; i < size-2; i++){
    //We step through and check to see if an element is greater than all of its neighbors
    //To minimize floating point erros we will multiply each element by 10^12 and then compare the results cast as longs (too avoid overflow for type int)
    //We do this by multiplying by a large integer and then casting as an int and then comparing the integers
    bool nneighbor_check = ( ( long(1e12*values[i]) > long(1e12*values[i-1]) ) && ( long(1e12*values[i]) > long(1e12*values[i+1]) ) );
    bool nnneighbor_check = ( ( long(1e12*values[i]) > long(1e12*values[i-2]) ) && ( long(1e12*values[i]) > long(1e12*values[i+2]) ) );

    if(nneighbor_check && nnneighbor_check){
      return i;
    } 
  }
  return 0;
}

//This function will locate the global max of an array using a simple loop through elements
int locateGlobalMax(int size, double * values){
    int index = 0;
    double val = values[index];
    
    //If we find a bigger value, we save it
    for(int i = 0; i < size; i++){
        if(values[i]>val){
            index = i;
            val = values[i];
        }
    }

    return index;
}

//INTERPOLATION FUNCTIONS
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


int main(int argc, char * argv[]){
	//We check to make sure they have the right number of arguments
 	if(argc!=2){
    	cout<<"Incorrect usage. Usage is: 'program_name' INPUTFILE "<<endl;
    	return 0;
  }
 
  string SPECFILE; //The file where the spectral data is input from
  string XASFILE; //The file where the XAS is input from
  string OUTFILE; //The file where the output is written to

	double fermi_ratio = 0.1;	//The ratio of the max peak height we define as being the Fermi energy (default of .1)

  //We read in the data 
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
    parseDbl(label,value,"fermi_ratio",fermi_ratio);
  }
  infile.close();

  cout<<"Done"<<endl;
  cout<<"Determining size of files"<<endl;

  int num_spec_steps = lineCount(SPECFILE);
  int num_xas_steps = lineCount(XASFILE);

  cout<<"   SPECFILE has "<<num_spec_steps<<" entries"<<endl;
  cout<<"   XASFILE has "<<num_xas_steps<<" entries"<<endl;
  cout<<"Done"<<endl;

  //Now we read in the data
  //First we must allocate arrays 
  double * spec_w = new double[num_spec_steps];
  double * spec = new double[num_spec_steps];
  double * xas_w = new double[num_xas_steps];
  double * xas = new double[num_xas_steps];

  //We read in the data 
	cout<<"Reading in from "<<SPECFILE<<endl;

  //Open up SPEC and XAS files and read the data in
  ifstream specfile(SPECFILE.c_str());
  string specline;
  int speccount = 0; //To increment through each line of the array that we are filling
  while(getline(specfile,specline)){
  //We now parse for 2 tokens, each a double 
    string token1, token2, token3;
    istringstream ss(specline);  //Extracts each token
    ss>>token1>>token2>>token3;
    //First token is frequency, second is density, third is cumulative (we ignore)
    spec_w[speccount] = atof(token1.c_str());
    spec[speccount] = atof(token2.c_str());
    
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
    xas_w[xascount] = atof(token1.c_str()); //We must also convert it to a double 
    xas[xascount] = atof(token2.c_str()); //We must also convert this to a double 
  	xascount++; //Increment the counter
  }
  xasfile.close();  //Close the file

  cout<<"Done"<<endl;
  cout<<"Interpolating"<<endl;

  //We will define the output grid to be the XPS grid 
  int num_out_steps = 2*num_spec_steps-1;	//We use the formula for convolution of Nf1*f2 = Nf1 + Nf2 - 1
  double * out_w = new double[num_out_steps];	
  double * out = new double[num_out_steps];
  double * xas_interp = new double[num_out_steps];	//This will be the xas interpolated onto the output grid 
  double * spec_interp = new double[num_out_steps]; //Same for spectral function 

  //We compute the spacing and generate a new grid of size 2*num_spec_steps - 1
  double delta_w = spec_w[1] - spec_w[0];

  //We build our frequencies from spec_w[0] + i*delta_w 
  for(int i = 0; i < num_out_steps; i++){
    out_w[i] = spec_w[0] + double(i)*delta_w;
  }

  //Now we interpolate the xas onto the spectral function grid 
  //The XAS is not, in general, centered. We do this by centering the highest peak
  int xas_gp = locateGlobalMax(num_xas_steps,xas);
  //We subtract the frequency this is from each frequency in the xas frequencies
  for( int i = 0; i < num_xas_steps; i++){
  	xas_w[i] -= xas_w[xas_gp];
  }

  //Now that this is centered, we interpolate it onto the output grid 
  for(int i = 0; i < num_out_steps; i++){
  	//First we go through the frequencies out_w and compute which xas_w they correspond to
  	double w = out_w[i];
  	int xas_index = snapToGrid(w,num_xas_steps,xas_w);
    int spec_index = snapToGrid(w,num_spec_steps,spec_w);

  	//Now we get the value of the xas and spec at this value, unless it is passed the last one, then we set it to the last value
  	if(xas_index != num_xas_steps-1){
  		xas_interp[i] = linearInterpolate(w,xas_w[xas_index],xas[xas_index],xas_w[xas_index+1],xas[xas_index+1]);
  	}
  	else{
  		xas_interp[i] = xas[num_xas_steps-1];
		}
    if(i < num_spec_steps){
      spec_interp[i] = spec[i];
    }
    else{
      spec_interp[i] = 0.0;
    }
 	}

  //Now we have the data on the same grid, we may compute the integrals 
  //First we must compute the fermi energy 
  //We calculate the global max of the XAS
  xas_gp = locateGlobalMax(num_out_steps,xas_interp);	//This variable is being repurposed
  double xas_qpw = out_w[xas_gp];	//The location of the global max
  double xas_gpv = xas_interp[xas_gp];	//The value of the global max

  //Now the threshold value
  double fermi_threshold = xas_gpv*fermi_ratio;	//The ratio of the max times the value of the max is the fermi value we look for
  //Now we loop until we reach this value
  int fermi_index = 0;
  
  //We look for the first value that passes the fermi threshold. Then we break the loop and return these values
  for(int i = 0; i < num_out_steps; i++){
    if(xas_interp[i] > fermi_threshold){
  		fermi_index = i;
      break;
  	}
  }

  //We also save the frequency and value of the XAS at the fermi level
  double fermi_w = out_w[fermi_index];
  double fermi_xas_value = xas_interp[fermi_index];

  //Print out the results 
  cout<<"XAS global max index = "<<xas_gp<<endl;
  cout<<"XAS global max w = "<<xas_qpw<<" eV"<<endl;
  cout<<"XAS global max value = "<<xas_gpv<<endl;
  cout<<"fermi_ratio = "<<fermi_ratio<<endl;
  cout<<"Fermi XAS threshold = "<<fermi_threshold<<endl;
  cout<<"Fermi XAS value = "<<fermi_xas_value<<endl;
  cout<<"Fermi index = "<<fermi_index<<endl;
  cout<<"Fermi Energy = "<<fermi_w<<" eV"<<endl;

  //Now we compute the integral (yikes!)
  //The integral is given by
  //out(w) = integral_{-infty}^{w-wFermi}spec(w')XAS(w-w')dw'
  //		   ------------------------------------------------
  //		   integral_{-infty}^{w-wFermi}spec(w')dw'

  //We print the results to the output file 

  //We calculate both the numerator and denominator in the same loop
  double * numerator = new double[num_out_steps];
  double * denominator = new double[num_out_steps];

  //We begin the calculation of the integrals
  for(int i = 0; i < num_out_steps; i++){
    //First we calculate the index that corresponds to w-fermi_w
    double w = out_w[i];  //The frequency we calculate for 

    //We integrate up to this energy 
    int shifted_w_index = snapToGrid(w-fermi_w,num_out_steps,out_w);

    //We use trapezoidal rule
    for(int j = 1; j < shifted_w_index; j++){
      numerator[i] += 0.5 * delta_w * (spec_interp[j-1] * xas_interp[i-j+1] + spec_interp[j] * xas_interp[i-j]);
      denominator[i] += 0.5 * delta_w * ( spec_interp[j-1] + spec_interp[j] );
    }
    
    //And we calculate the integral as numerator/denominator
    out[i] = numerator[i]/denominator[i];
  }

  //We now compute the maximum and its index and value 
  int out_gp = locateGlobalMax(num_out_steps,out);
  double out_gpw = out_w[out_gp];
  double out_gpv = out[out_gp];

  cout<<"out_global_max_index = "<<out_gp<<endl;
  cout<<"out_global_max_w = "<<out_gpw<<" eV"<<endl;
  cout<<"out_global_value = "<<out_gpv<<endl;

  //We also print the total integral of the XAS before and after
  double before_integral = 0.0;
  double after_integral = 0.0;

  for(int i = 1; i < num_out_steps; i++){
    before_integral += 0.5 * delta_w * (xas_interp[i-1]+ xas_interp[i]);
    after_integral += 0.5 * delta_w * (out[i-1] + out[i]);
  }

  cout<<"Integral of single particle XAS = "<<before_integral<<endl;
  cout<<"Integral of many particle XAS = "<<after_integral<<endl;

  ofstream outfile(OUTFILE.c_str());
  for(int i = 0; i < num_out_steps; i++){
  	outfile<<out_w[i]<<" "<<out[i]<<" "<<xas_interp[i]<<" "<<spec_interp[i]<<endl;
  }
  outfile.close();

  delete [] xas;
  delete [] xas_w;
  delete [] spec;
  delete [] spec_w;
  delete [] out_w;
  delete [] out;
  delete [] xas_interp;
  delete [] spec_interp;

  return 0;
}



