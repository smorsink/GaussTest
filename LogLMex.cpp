/***************************************************************************************/
/*                                   LogLike1d.cpp

One Dimensional Hastings-Metropolis Algorithm

*/
/***************************************************************************************/

#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <exception>
#include <vector>
#include <string>
#include <string.h>
#include "nrutil.h"
#include "bayes.h"
#include "Struct.h"
#include "io64.h"
#include "mex.h"

// MAIN

void mexFunction ( int numOutputs, mxArray *theOutput[], int numInputs, const mxArray *theInput[]){

    double *curveOut, *chiOut;

    // Setting up the output parameters
    int dimSize[2];
    dimSize[0] = 1;
    dimSize[1] = 1;

    // output needs to be an array so set it up as a [1,1] array holding one value
    theOutput[0] = mxCreateNumericArray(2, dimSize, mxDOUBLE_CLASS, mxREAL);

    chiOut = mxGetPr(theOutput[0]);

  /*********************************************/
  /* VARIABLE DECLARATIONS AND INITIALIZATIONS */
  /*********************************************/
    
   // double *probmean = dvector(1,NN);
    //double *probsig = dvector(1,NN);
    //double *meanvals = dvector(1,NN);
    //double *sigvals = dvector(1,NN);
  //std::ifstream in;      // output stream; printing information to the output file
  //char in_file[256] = "Nothing!";

  
  unsigned int nsteps(100);


  double
    xlo(-5.0), xhi(5.0);
  int
    numbins(100);

  double
    meanlo(-4.0), meanhi(4.0), siglo(0.0), sighi(6.0);
  int NN(4000);

	double xmean, xsigma;


 
	nsteps = mxGetScalar(theInput[0]); // number of data bins
	numbins = mxGetScalar(theInput[1]); // number of bins for histograms
	xmean = mxGetScalar(theInput[2]); // mean value of x for gaussian	
	xsigma = mxGetScalar(theInput[3]); // standard deviation of x for gaussian	

	double *data = dvector(1,nsteps); 

    	data = mxGetPr(theInput[4]); // array of double

  
 
    long int *histogram = ivector(1,numbins);
    double *probability = dvector(1,numbins);
    double *prob1 = dvector(1,numbins);
    double *xvals = dvector(1,numbins);
    double dx;

    dx = (xhi - xlo)/(1.0*numbins);
    for (unsigned int i(1); i<=numbins; i++){
      xvals[i] = xlo + (0.5 + i - 1)*dx;
      histogram[i] = 0;
    }

    int bin;

    double dmean, dsig;
    dmean = (meanhi-meanlo)/(1.0*NN);
    dsig = (sighi-siglo)/(1.0*NN);
    double *probmean = dvector(1,NN);
    double *probsig = dvector(1,NN);
    double *meanvals = dvector(1,NN);
    double *sigvals = dvector(1,NN);

    for (unsigned int i(1); i<=NN; i++){
      meanvals[i] = meanlo + (0.5+i-1)*dmean;
      sigvals[i] = siglo + (0.5+i-1)*dsig;
      probmean[i]=0.0;
      probsig[i] = 0.0;
    }


    //in.open(in_file);

    double x1(0.0);
    //double xave(0.0), yave(0.0);
  
    for ( unsigned int i(1); i <= nsteps; i++){

	x1 = data[i];

     //in >> x1;     
     // Increment the correct bin of the histogram

     if (x1 < xlo || x1 > xhi)
       std::cout << "x1 = " << x1 << " is out of bounds! Increase them!" << std::endl;
     else{ // Within the bounds
       bin = (x1-xlo)/dx + 1;
       histogram[bin] += 1;
     }
    } // End of Data Histogram Loop

    // Integrate the histogram
    double totprob=0.0;
    for (unsigned int i(1); i<=numbins; i++){
      totprob += histogram[i]*dx;
    }
   // Normalize the histogram
   for (unsigned int i(1); i<=numbins; i++){
     probability[i] = histogram[i]/(1.0*totprob);
     //cumul += probability[i];
   }

 

    double ll1=0.0;

 
 
 

    // Compute histogram1 with numbin bins, mean1, and sig1
    GaussProbDist(xvals,prob1,xmean,xsigma,numbins);
    ll1 = ProbCompare(probability,prob1,numbins);

    std::cout << "GaussTest1d: "
		  << " mean = " << xmean
		  << " sig  = " << xsigma
		  << " ll1 = " << ll1
		  << std::endl;

    /*****************************************************************/
    /* DUMPING DATA INTO MATLAB OUTPUT                               */
    /* For filewriting, need to use fopen, fprintf, and fclose (etc) */
    /*****************************************************************/


 	
    chiOut[0] = -ll1; // saved this to theOutput[0] at top just after declarations





     //dimSize[1] = 3; // number of columns (1: time, 2: flux in 1st energy band, 3: flux in 2nd energy band)
    //std::cout << "Output 2." << std::endl;

    //dimSize[0] = (int)numbins; // number of rows ( = numbins)
	dimSize[0] = 1;

    theOutput[1] = mxCreateNumericArray(2, dimSize, mxDOUBLE_CLASS, mxREAL); // formatting/setting up the output
  

	// Free Memory!

	free_ivector(histogram,1,numbins);
	free_dvector(probability,1,numbins);
	free_dvector(prob1,1,numbins);
	free_dvector(xvals,1,numbins);

	free_dvector(probmean,1,NN);
	free_dvector(probsig,1,NN);
	free_dvector(meanvals,1,NN);
	free_dvector(sigvals,1,NN);

 } 


