/***************************************************************************************/
/*                                   Metrop.cpp

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

  /*********************************************/
  /* VARIABLE DECLARATIONS AND INITIALIZATIONS */
  /*********************************************/
    
  std::ofstream out;      // output stream; printing information to the output file
  std::ofstream histout;      // output stream; printing information to the output file
  char out_file[256] = "output.txt";
  std::ifstream in;      // output stream; printing information to the output file
  char in_file[256] = "Nothing!";

  
  unsigned int nsteps(100);


  double
    xlo(-5.0), xhi(5.0);
  int
    numbins(100);

  struct Prob1d xx;


  /*********************************************************/
  /* READING IN PARAMETERS FROM THE COMMAND LINE ARGUMENTS */
  /*********************************************************/
    
    for ( int i(1); i < argc; i++ ) {
        if ( argv[i][0] == '-' ) {  // the '-' flag lets the computer know that we're giving it information from the cmd line
            switch ( argv[i][1] ) {

	    case 'i':  // Name of output file
	                sscanf(argv[i+1], "%s", in_file);
	                break;
	      
	    case 'o':  // Name of output file
	                sscanf(argv[i+1], "%s", out_file);
	                break;

	    case 'n': // Number of steps
	      sscanf(argv[i+1], "%u",&nsteps);
	      break;

	    case 'b': // Number of xbins
	      sscanf(argv[i+1], "%u",&numbins);
	      break;
	      
	                
            } // end switch	
        } // end if
    } // end for

    //    int **histogram = imatrix(1,numradius,1,nummass);
 
    long int *histogram = ivector(1,numbins);
    double *probability = dvector(1,numbins);
    double *prob1 = dvector(1,numbins);
    double *prob2 = dvector(1,numbins);

    double *xvals = dvector(1,numbins);

    double dx;

    dx = (xhi - xlo)/(1.0*numbins);
    for (unsigned int i(1); i<=numbins; i++){
      xvals[i] = xlo + (0.5 + i - 1)*dx;
      histogram[i] = 0;
    }

    int bin;

    // Set up x-bins for storing the histogram


    in.open(in_file);

    double x1(0.0), x2(0.0);
    double xave(0.0), yave(0.0);
  
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
     probability[i] = histogram[i]/(1.0*totprob) * dx;
     //cumul += probability[i];
   }

    srand(time(NULL));   // should only be called once

    double mean1(0.0), mean2(0.0);
    double sig1(1.0),sig2;
    double ll1=0.0, ll2=0.0;

    double var(2.0);
    double mean(2.0);
    double sig;
    sig = sqrt(var);

    double meanave(0.0), sigave(0.0);

    double qmux, qmuy, qvar, qsig;
    double r;
    double yes(0);

    for (unsigned int i(0); i < 10; i++){
     
      meanave += mean1;
      sigave += sig1;

      // Compute histogram1 with numbin bins, mean1, and sig1
      GaussProbDist(xvals,prob1,mean1,sig1,numbins);

      // Compute ll1 = log(likelihood) of first proposal by comparing histogram1 with data histogram

      // Create Proposal with mean=x1; var=1
      qmux = mean1;
      qmuy = sig1;
      qvar = 1.0;
      qsig = sqrt(qvar);

      // Draw a new value for mean and sigma from the Proposal Distribution
      mean2 = NormalDev(qmux,qsig);
      sig2 = NormalDev(qmuy,qsig);

      // Compute histogram1 with numbin bins, mean1, and sig1
      // Compute ll22 = loglikelihood of second proposal by comparing histogram2 with data histogram
    
     r = Rand1();
     // Do we accept the new value?

      if ( ll2 - ll1 > log(r)){ // Accept the new step!

	mean1 = mean2;
	sig1 = sig2;
	ll1 = ll2;
	yes += 1;

	std::cout << "i = " << i
		  << " mean = " << mean1
		  << " sig  = " << sig1
		  << " ll1 = " << ll1
		  << std::endl;


      }
      else{ // Reject the new step!

      }



    }




   //std::cout << "total prob = " << totprob << std::endl;
   double xmed(xlo);
   double xl1(xlo), xl2(xlo), xl3(xlo);
   double xr1(xlo), xr2(xlo), xr3(xlo);
   double cumul(0.0);

   // Normalize the histogram
   for (unsigned int i(1); i<=numbins; i++){
     probability[i] = histogram[i]/(1.0*totprob) * dx;
     cumul += probability[i];
     /*std::cout << "i=" << i 
	       << " x[i]= " << xvals[i]
	       << " prob[i]=  " << probability[i]
	       << " cumul= " << cumul
	       << std::endl;*/
   }

      std::cout << "N steps = " << nsteps << "\t"
		<< " dx = " << dx
		<< " Integral of Prob = " << cumul 
		<< " xave = " << xave
	     << " Median x = " << xmed << std::endl;

      xx = ProbContours1D(xvals,probability,xave,numbins);
      std::cout << "Probability distribution  for " << nsteps << " values of x" << std::endl;

      std::cout << "2 sigma x = " << xx.xl2 << std::endl;
      std::cout << "1 sigma x = " << xx.xl1 << std::endl;
      std::cout << "Median  x = " << xx.median << std::endl;
      std::cout << "1 sigma x = " << xx.xr1 << std::endl;
      std::cout << "2 sigma x = " << xx.xr2 << std::endl;


      histout.open("gplotx.txt");
      histout << "unset arrow" << std::endl;
      histout << "set arrow 1 from " << xx.xl3 << ",0 to " << xx.xl3 << "," << xx.pl3 << " nohead lt 5" << std::endl;
      histout << "set arrow from " << xx.xl2 << ",0 to " << xx.xl2 << "," << xx.pl2 << " nohead lt 4" << std::endl;
      histout << "set arrow from " << xx.xl1 << ",0 to " << xx.xl1 << "," << xx.pl1 << " nohead lt 2" << std::endl;
      histout << "set arrow from " << xx.median << ",0 to " << xx.median << "," << xx.pmed << " nohead lt 1 lw 3 " << std::endl;
      histout << "set arrow from " << xx.xr1 << ",0 to " << xx.xr1 << "," << xx.pr1 << " nohead lt 2" << std::endl;
      histout << "set arrow from " << xx.xr2 << ",0 to " << xx.xr2 << "," << xx.pr2 << " nohead lt 4" << std::endl;
      histout << "set arrow from " << xx.xr3 << ",0 to " << xx.xr3 << "," << xx.pr3 << " nohead lt 5" << std::endl;
      histout << "set title  \"1D Gaussian Test with " << nsteps << " values \"" << std::endl;
      histout << "set xlabel \"x\"" << std::endl;
      histout << "set ylabel \"Probability\"" << std::endl;
      histout << "set key left" << std::endl;
      histout << "plot \"" << out_file << "\" with lines lt -1 notitle, 1/0 t \"2 sigma x = "
	      << xx.xl2 << "\" lt 4, 1/0 t \"1 sigma x = "
	      << xx.xl1 << "\" lt 2, 1/0 t \"Median  x = "
	      << xx.median << "\" lt 1 lw 3, 1/0 t \"1 sigma x = "
	      << xx.xr1 << "\" lt 2, 1/0 t \"2 sigma x = "
	      << xx.xr2 << "\" lt 4, \"" << out_file << "\" using 1:3 with lines lt 1 notitle"
	      << std::endl;     
      histout.close();

      
      mean = 0.0;
      sig = 1.0;
      
   // Output histogram
   out.open(out_file);
   for (unsigned int i(1); i<=numbins; i++){

     out << xvals[i] << " " << probability[i] << " " 
	 << Gaussian(xvals[i],mean,sig)*dx << " "
	 << i << std::endl;
   }


    return 0;


    //free_ivector(histogram,1,numbins);
    free_dvector(xvals,1,numbins);

 } 

catch(std::exception& e) {
       std::cerr << "\nERROR: Exception thrown. " << std::endl
	             << e.what() << std::endl;
       return -1;
}
