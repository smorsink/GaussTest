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



  int chainlength(1000);
  double fudge(1.0);


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

	    case 'O':  // Name of output file
	                sscanf(argv[i+1], "%s", out_sig);
	                break;


	    case 'n': // Number of steps
	      sscanf(argv[i+1], "%u",&nsteps);
	      break;

	    case 'b': // Number of xbins
	      sscanf(argv[i+1], "%u",&numbins);
	      break;
	      
	    case 'c': // Chainlength
	      sscanf(argv[i+1], "%u",&chainlength);
	      break;

	    case 'f': // Fudge Factor
	      sscanf(argv[i+1], "%lf",&fudge);
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

    srand(time(NULL));   // should only be called once

    double mean1(0.5), mean2(0.0);
    double sig1(1.5),sig2;
    double ll1=0.0, ll2=0.0;

    double var(2.0);
    //double mean(2.0);
    double sig;
    sig = sqrt(var);

    double meanave(0.0), sigave(0.0);

    double qmux, qmuy, qvar, qsig;
    double r;
    double yes(0),acceptance(0);

    // Compute histogram1 with numbin bins, mean1, and sig1
    GaussProbDist(xvals,prob1,mean1,sig1,numbins);
    ll1 = ProbCompare(probability,prob1,numbins);

    std::cout << "First Choice"
		  << " mean = " << mean1
		  << " sig  = " << sig1
		  << " ll1 = " << ll1
		  << std::endl;


    out.open("trace.txt");

    for (unsigned int i(0); i < chainlength; i++){ // M-T loop
     
      meanave += mean1;
      sigave += sig1;

      // Create Proposal with mean=x1; var=1
      qmux = mean1;
      qmuy = sig1;
      qvar = 1.0;
      qsig = sqrt(qvar);

      // Draw a new value for mean and sigma from the Proposal Distribution
      mean2 = NormalDev(qmux,qsig*fudge);
      sig2 = NormalDev(qmuy,qsig*fudge);


      //sig2 = 1.0; // This is only for debugging purposes

      // Compute histogram1 with numbin bins, mean1, and sig1
      // Compute ll22 = loglikelihood of second proposal by comparing histogram2 with data histogram

      if (sig2 >= 0.0){
	GaussProbDist(xvals,prob2,mean2,sig2,numbins);
	ll2 = ProbCompare(probability,prob2,numbins);
      }
      else
	ll2 = -100000.0;


      /*   std::cout << std::endl << "New Proposal i = " << i
		<< " mean2 = " << mean2
		<< " sig2  = " << sig2
		<< " ll2 = " << ll2
		<< std::endl;*/

     r = Rand1();
     // Do we accept the new value?

     //std::cout << "log(r) = " << log(r) 
     //	       << " ll2 - ll1 = " << ll2 - ll1 
     //	       << std::endl;

      if ( ll2 - ll1 > log(r)){ // Accept the new step!

	mean1 = mean2;
	sig1 = sig2;
	ll1 = ll2;
	yes += 1;
	acceptance = yes/(1.0*i);

	if (mean1 < meanlo || mean1 > meanhi)
	  std::cout << "mean1 = " << mean1 << " is out of bounds! Increase them!" << std::endl;
	else{ // Within the bounds
	  bin = (mean1-meanlo)/dmean + 1;
	  probmean[bin] += 1;
	}
	if (sig1 < siglo || sig1 > sighi)
	  std::cout << "sig1 = " << sig1 << " is out of bounds! Increase them!" << std::endl;
	else{ // Within the bounds
	  bin = (sig1-siglo)/dsig + 1;
	  probsig[bin] += 1;
	}
	
	/*	std::cout << "i = " << i
		  << " mean = " << mean1
		  << " sig  = " << sig1
		  << " ll2 = " << ll2
		  << " yes = " << yes
		  << std::endl;*/

	if (i%1000==0)
	  out << i << "\t"
	      << mean1 << "\t"
	      << sig1 << "\t"
	      << ll1 << "\t"
	      << acceptance << std::endl;
      }
      else{ // Reject the new step!
      }

    } //end M-T Loop

    out.close();

    std::cout << "Average Value of the Mean = " << meanave/(chainlength*1.0) << std::endl;
    std::cout << "Average Value of the StdDev = " << sigave/(chainlength*1.0) << std::endl;
      
    // Normalize the probability distributions for mean and sigma

    double totmean(0.0), totsig(0.0);
    for (unsigned int i(1);i<=NN;i++){
      totmean += probmean[i];
      totsig += probsig[i];
    }
    for (unsigned int i(1);i<=NN;i++){
      probmean[i] /= totmean;
      probsig[i] /= totsig;
    }

   // Output histogram
   out.open(out_file);
   for (unsigned int i(1); i<=NN; i++){
     out << meanvals[i] << " " << probmean[i] << " " 
	 << std::endl;
   }
   out.close();
   out.open(out_sig);
   for (unsigned int i(1); i<=NN; i++){
     out << sigvals[i] << " " << probsig[i] << " " 
	 << std::endl;
   }
   out.close();

   
    struct Prob1d Mean, Sig;
    Mean = ProbContours1D(meanvals,probmean,meanave,NN);
      std::cout << "Probability distribution  for the mean value of x. " << nsteps << " values of x" << std::endl;

      std::cout << "2 sigma xmean = " << Mean.xl2 << std::endl;
      std::cout << "1 sigma xmean = " << Mean.xl1 << std::endl;
      std::cout << "Median  xmean = " << Mean.median << std::endl;
      std::cout << "1 sigma xmean = " << Mean.xr1 << std::endl;
      std::cout << "2 sigma xmean = " << Mean.xr2 << std::endl;

      histout.open("gplotmean.txt");
      histout << "unset arrow" << std::endl;
      histout << "set arrow 1 from " << Mean.xl3 << ",0 to " << Mean.xl3 << "," << Mean.pl3 << " nohead lt 5" << std::endl;
      histout << "set arrow from " << Mean.xl2 << ",0 to " << Mean.xl2 << "," << Mean.pl2 << " nohead lt 4" << std::endl;
      histout << "set arrow from " << Mean.xl1 << ",0 to " << Mean.xl1 << "," << Mean.pl1 << " nohead lt 2" << std::endl;
      histout << "set arrow from " << Mean.median << ",0 to " << Mean.median << "," << Mean.pmed << " nohead lt 1 lw 3 " << std::endl;
      histout << "set arrow from " << Mean.xr1 << ",0 to " << Mean.xr1 << "," << Mean.pr1 << " nohead lt 2" << std::endl;
      histout << "set arrow from " << Mean.xr2 << ",0 to " << Mean.xr2 << "," << Mean.pr2 << " nohead lt 4" << std::endl;
      histout << "set arrow from " << Mean.xr3 << ",0 to " << Mean.xr3 << "," << Mean.pr3 << " nohead lt 5" << std::endl;
      histout << "set title  \"1D Gaussian Test with " << nsteps << " values of x \"" << std::endl;
      histout << "set xlabel \"Mean value of x\"" << std::endl;
      histout << "set ylabel \"Probability\"" << std::endl;
      histout << "set key right" << std::endl;
      histout << "plot \"" << out_file << "\" with lines lt -1 notitle, 1/0 t \"2 sigma mean x = "
	      << Mean.xl2 << "\" lt 4, 1/0 t \"1 sigma mean x = "
	      << Mean.xl1 << "\" lt 2, 1/0 t \"Median  mean x = "
	      << Mean.median << "\" lt 1 lw 3, 1/0 t \"1 sigma mean x = "
	      << Mean.xr1 << "\" lt 2, 1/0 t \"2 sigma mean x = "
	      << Mean.xr2 << "\" lt 4 " 
	      << std::endl;     
      histout.close();

      Sig = ProbContours1D(sigvals,probsig,sigave,NN);
      std::cout << "Probability distribution  for the mean value of x. " << nsteps << " values of x" << std::endl;

      std::cout << "2 sigma xsig = " << Sig.xl2 << std::endl;
      std::cout << "1 sigma xsig = " << Sig.xl1 << std::endl;
      std::cout << "Median  xsig = " << Sig.median << std::endl;
      std::cout << "1 sigma xsig = " << Sig.xr1 << std::endl;
      std::cout << "2 sigma xsig = " << Sig.xr2 << std::endl;

      histout.open("gplotsig.txt");
      histout << "unset arrow" << std::endl;
      histout << "set arrow 1 from " << Sig.xl3 << ",0 to " << Sig.xl3 << "," << Sig.pl3 << " nohead lt 5" << std::endl;
      histout << "set arrow from " << Sig.xl2 << ",0 to " << Sig.xl2 << "," << Sig.pl2 << " nohead lt 4" << std::endl;
      histout << "set arrow from " << Sig.xl1 << ",0 to " << Sig.xl1 << "," << Sig.pl1 << " nohead lt 2" << std::endl;
      histout << "set arrow from " << Sig.median << ",0 to " << Sig.median << "," << Sig.pmed << " nohead lt 1 lw 3 " << std::endl;
      histout << "set arrow from " << Sig.xr1 << ",0 to " << Sig.xr1 << "," << Sig.pr1 << " nohead lt 2" << std::endl;
      histout << "set arrow from " << Sig.xr2 << ",0 to " << Sig.xr2 << "," << Sig.pr2 << " nohead lt 4" << std::endl;
      histout << "set arrow from " << Sig.xr3 << ",0 to " << Sig.xr3 << "," << Sig.pr3 << " nohead lt 5" << std::endl;
      histout << "set title  \"1D Gaussian Test with " << nsteps << " values of x \"" << std::endl;
      histout << "set xlabel \"Standard Deviation of x\"" << std::endl;
      histout << "set ylabel \"Probability\"" << std::endl;
      histout << "set key right" << std::endl;
      histout << "plot \"" << out_sig << "\" with lines lt -1 notitle, 1/0 t \"2 sigma x sigma = "
	      << Sig.xl2 << "\" lt 4, 1/0 t \"1 sigma x sigma = "
	      << Sig.xl1 << "\" lt 2, 1/0 t \"Median  x sigma= "
	      << Sig.median << "\" lt 1 lw 3, 1/0 t \"1 sigma x sigma= "
	      << Sig.xr1 << "\" lt 2, 1/0 t \"2 sigma x sigma= "
	      << Sig.xr2 << "\" lt 4 " 
	      << std::endl;     
      histout.close();


      


    return 0;


    //free_ivector(histogram,1,numbins);
    free_dvector(xvals,1,numbins);

 } 

catch(std::exception& e) {
       std::cerr << "\nERROR: Exception thrown. " << std::endl
	             << e.what() << std::endl;
       return -1;
}
