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

// MAIN
int main ( int argc, char** argv ) try {  // argc, number of cmd line args; 
                                          // argv, actual character strings of cmd line args

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
    
  std::ofstream out;      // output stream; printing information to the output file
  std::ofstream histout;      // output stream; printing information to the output file
  char out_file[256] = "output.txt";
  char out_sig[256] = "output.txt";

  std::ifstream in;      // output stream; printing information to the output file
  char in_file[256] = "Nothing!";

  
  unsigned int nsteps(100);


  double
    xlo(-5.0), xhi(5.0);
  int
    numbins(100);

  double
    meanlo(-4.0), meanhi(4.0), siglo(0.0), sighi(6.0);
  int NN(4000);

	double xmean, xsigma;



  /*********************************************************/
  /* READING IN PARAMETERS FROM THE COMMAND LINE ARGUMENTS */
  /*********************************************************/
    
    for ( int i(1); i < argc; i++ ) {
        if ( argv[i][0] == '-' ) {  // the '-' flag lets the computer know that we're giving it information from the cmd line
            switch ( argv[i][1] ) {

	    case 'i':  // Name of output file
	                sscanf(argv[i+1], "%s", in_file);
	                break;
	      
	    case 'n': // Number of steps
	      sscanf(argv[i+1], "%u",&nsteps);
	      break;

	    case 'b': // Number of xbins
	      sscanf(argv[i+1], "%u",&numbins);
	      break;

	    case 'm': // mean
	      sscanf(argv[i+1], "%lf",&xmean);
	      break;

	    case 's': // sigma
	      sscanf(argv[i+1], "%lf",&xsigma);
	      break;


	                
            } // end switch	
        } // end if
    } // end for

    //    int **histogram = imatrix(1,numradius,1,nummass);
 
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


    in.open(in_file);

    double x1(0.0);
    //double xave(0.0), yave(0.0);
  
    for ( unsigned int i(0); i < nsteps; i++){

     in >> x1;     
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


    

    return 0;


    //free_ivector(histogram,1,numbins);
    free_dvector(xvals,1,numbins);

 } 

catch(std::exception& e) {
       std::cerr << "\nERROR: Exception thrown. " << std::endl
	             << e.what() << std::endl;
       return -1;
}
