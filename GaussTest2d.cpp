/***************************************************************************************/
/*                                   GaussTest2d.cpp

Two Dimensional Gaussian TEst

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
  char out_filex[256] = "output.txt";
  char out_filey[256] = "output.txt";
  std::ifstream in;      // output stream; printing information to the output file
  char in_file[256] = "Nothing!";

  unsigned int nsteps(100), NN(10000);
  unsigned int chainlength(10);

  double
    xlo(-23.0), xhi(23.0);
  double
    ylo(-23.0), yhi(23.0);
  int
    numbins(500);

  struct Prob1d xx, yy;
  struct Prob2d contours;

  double fudge(1.0);

  /*********************************************************/
  /* READING IN PARAMETERS FROM THE COMMAND LINE ARGUMENTS */
  /*********************************************************/
    
    for ( int i(1); i < argc; i++ ) {
        if ( argv[i][0] == '-' ) {  // the '-' flag lets the computer know that we're giving it information from the cmd line
            switch ( argv[i][1] ) {

	    case 'i':  // Name of input file
	                sscanf(argv[i+1], "%s", in_file);
	                break;
	      
	    case 'o':  // Name of output file
	                sscanf(argv[i+1], "%s", out_filex);
	                break;

	    case 'O':  // Name of output file
	      sscanf(argv[i+1], "%s", out_filey);
	      break;

	    case 'n': // Number of steps
	      sscanf(argv[i+1], "%u",&nsteps);
	      break;

	    case 'b': // Number of xbins
	      sscanf(argv[i+1], "%u",&numbins);
	      break;

	    case 'p': // Number of probability bins
	      sscanf(argv[i+1], "%u",&NN);
	      break;

	    case 'f': // Fudge Factor
	      sscanf(argv[i+1], "%lf",&fudge);
	      break;

	    case 'c': // Chainlength
	      sscanf(argv[i+1], "%u",&chainlength);
	      break;
	      
            } // end switch	
        } // end if
    } // end for

 
    long int **histogram = imatrix(1,numbins,1,numbins);
    double **probability = dmatrix(1,numbins,1,numbins);
    double *xprob = dvector(1,numbins);
    double *yprob = dvector(1,numbins);
    double *xvals = dvector(1,numbins);
    double *yvals = dvector(1,numbins);

    double dx, dy;

    dx = (xhi - xlo)/(1.0*numbins);
    dy = (yhi - ylo)/(1.0*numbins);

    for (unsigned int i(1); i<=numbins; i++){
      xvals[i] = xlo + (0.5 + i - 1)*dx;
      yvals[i] = ylo + (0.5 + i - 1)*dy;
    }
    std::cout << "dx = " << dx << " dy = " << dy << std::endl;

    for (unsigned int i(1);i<=numbins; i++)
      for (unsigned int j(1); j<=numbins; j++)
	histogram[i][j] = 0;
    
    int xbin, ybin;

    int bin;

    //double dmean, dsig;
    //double meanhi(10.0), meanlo(-10.0), sighi(50.0),siglo(0.0);

    double *dmean = dvector(1,2);
    double *dsig = dvector(1,2);
    double *meanhi = dvector(1,2);
    double *meanlo = dvector(1,2);
    double *sighi = dvector(1,2);
    double *siglo = dvector(1,2);

    // Initialize 1e2
    /* meanhi[1] = 0.1;
    meanlo[1] = -0.1;
    meanhi[2] = 5.0;
    meanlo[2] = -5.0;

    sighi[1] = 2.0;
    siglo[1] = 0.0;
    sighi[2] = 20.0;
    siglo[2] = 0.0;*/

    // Initialize 1e4
    meanhi[1] = 0.1;
    meanlo[1] = -0.1;
    meanhi[2] = 4.0;
    meanlo[2] = -4.0;

    sighi[1] = 1.0;
    siglo[1] = 0.0;
    sighi[2] = 15.0;
    siglo[2] = 5.0;


    

    for (unsigned int j(1); j<=2; j++){
      dmean[j] = (meanhi[j]-meanlo[j])/(1.0*NN);
      std::cout << "delta(mean["<< j << "]=" << dmean[j] << std::endl;
      dsig[j] = (sighi[j]-siglo[j])/(1.0*NN);
      std::cout << "delta(sig["<< j << "]=" << dsig[j] << std::endl;
    }
    
    double **probmean = dmatrix(1,2,1,NN);
    double **probsig = dmatrix(1,2,1,NN);
    double **meanvals = dmatrix(1,2,1,NN);
    double **sigvals = dmatrix(1,2,1,NN);

    // j=1 is "-" solution
    // j=2 is "+" solution
    for (unsigned int i(1); i<=NN; i++){
      for (unsigned int j(1);j<=2;j++){
      meanvals[j][i] = meanlo[j] + (0.5+i-1)*dmean[j];
      sigvals[j][i] = siglo[j] + (0.5+i-1)*dsig[j];
      probmean[j][i]=0.0;
      probsig[j][i] = 0.0;
      }
    }

    // Set up x-bins for storing the histogram

  double x1(0.0);
  double y1(0.0);
  double xave(0.0), yave(0.0);

  // Read in the Data file
  in.open(in_file);

   for ( unsigned int i(0); i < nsteps; i++){

     in >> x1;     
     in >> y1;
          
     xave += x1;
     yave += y1;

     // Increment the correct bin of the histogram

     if (x1 < xlo || x1 > xhi)
       std::cout << "x1 = " << x1 << " is out of bounds! Increase them!" << std::endl;
     else{ // Within the x bounds
       if (y1 < ylo || y1 > yhi)
	 std::cout << "y1 = " << y1 << " is out of bounds! Increase them!" << std::endl;
       else{ // Within the y bounds

	 xbin = (x1-xlo)/dx + 1;
	 ybin = (y1-ylo)/dy + 1;
	 
	 histogram[xbin][ybin] += 1;

       }
     }

   } // End of binning loop

   xave /= (1.0*nsteps);
   yave /= (1.0*nsteps);

   std::cout << "Average x = " << xave << std::endl;
   std::cout << "Average y = " << yave << std::endl;


   // Integrate the histogram
   double totprob=0.0;
   for (unsigned int i(1); i<=numbins; i++){
     for (unsigned int j(1); j<=numbins; j++){
       totprob += histogram[i][j];
     }
   }
   std::cout << "totprob = " << totprob << std::endl;
   // Normalize the probability
   for (unsigned int i(1); i<=numbins; i++){
     for (unsigned int j(1); j<=numbins; j++){
       probability[i][j] = 4.0*histogram[i][j]/(1.0*totprob*dx);
       /*std::cout << "hist[i][j] = " << histogram[i][j]
		 << " data[i][j] = " << probability[i][j]
		 << std::endl;*/
     }
   }
   // Probability is normalized so integral over x and y yields 1.0
   // totprob = 0.0;
   //Integrate over y to get the xprobability 
   for (unsigned int i(1); i<=numbins; i++){
     xprob[i] = 0.0;
     for (unsigned int j(1); j<=numbins; j++){
       xprob[i] += probability[i][j]*dx*dy;
     }
     //totprob += xprob[i]*dx;
   }
   for (unsigned int j(1); j<=numbins; j++){
     yprob[j] = 0.0;
     for (unsigned int i(1); i<=numbins; i++){
       yprob[j] += probability[i][j]*dx*dy;
     }
     //totprob += yprob[j]*dy;
   }


   // Create a proposal and evaluate log(likelyhood)
   
   double *mean1 = dvector(1,2);
   double *sig1 = dvector(1,2);
   double *mean2 = dvector(1,2);
   double *sig2 = dvector(1,2);
   double *meanave = dvector(1,2);
   double *sigave = dvector(1,2);

   double **prob1 = dmatrix(1,numbins,1,numbins);
   double **prob2 = dmatrix(1,numbins,1,numbins);

   double ll1,ll2;
   
   //Initialize
   mean1[1] = -0.00035;
   //mean1[1] = 0.001;
   mean1[2] = -0.1591;
   sig1[1] = 0.0987;
   sig1[2] = 10.025;

   // sig1[1] = 0.1;

   
   GaussProbDist2d(xvals, yvals, prob1,mean1,sig1,numbins);
   ll1 = ProbCompare2d(probability,prob1,numbins);

   out.precision(7);
    out.open("test.txt");
   for (unsigned int i(1);i<=numbins;i++){
     for (unsigned int j(1);j<=numbins;j++){

       out << xvals[i] << "\t"
	   << yvals[j] << "\t"
	   << probability[i][j] << "\t"
	   << prob1[i][j] <<"\t"
	   << histogram[i][j]
	   << std::endl;

     }
   }
   out.close();
   
   std::cout << "First Choice"
	      << " mean- = " << mean1[1] << " mean+ = " << mean1[2]
		  << " sig-  = " << sig1[1]
     		  << " sig+  = " << sig1[2]
		  << " ll1 = " << ll1
		  << std::endl;

    out.open("trace.txt");

   
    srand(time(NULL));   // should only be called once
    
    double *qsigmean = dvector(1,2);
    double *qsigsig = dvector(1,2);

    qsigmean[1] = 0.001;
    qsigmean[2] = 0.01;
    qsigsig[1] = 0.001;
    qsigsig[2] = 0.05;
    
    double r, yes(0), acceptance;

    for (unsigned int i(0); i < chainlength; i++){ // M-T loop
   
      // Draw a new value for mean and sigma from the Proposal Distribution
      for (unsigned int j(1);j<=2; j++){
	//Compute running totals for averages
	meanave[j] += mean1[j];
	sigave[j] += sig1[j];

	//mean2[j] = mean1[j];
	//sig2[j] = sig1[j];

	// Create a new proposal
	mean2[j] = NormalDev(mean1[j],qsigmean[j]*fudge);
	sig2[j] = NormalDev(sig1[j],qsigsig[j]*fudge);
	
      }

      //mean2[1] = NormalDev(mean1[1],qsigmean[1]*fudge);

      if ( sig2[1] > 0 && sig2[2] > 0){
	GaussProbDist2d(xvals, yvals, prob2,mean2,sig2,numbins);
	ll2 = ProbCompare2d(probability,prob2,numbins);
      }
      else
	ll2 = -10000000.0;

      for (unsigned int j(1);j<=2; j++){
	if (mean2[j] < meanlo[j] || mean2[j] > meanhi[j]){
	  std::cout << "mean2 = " << mean2[j] << " is out of bounds! Increase them!" << std::endl;
	  ll2 = -1000000.0;
	}
	if (sig2[j] < siglo[j] || sig2[j] > sighi[j]){
	  std::cout << "sig2 = " << sig2[j] << " is out of bounds! Increase them!" << std::endl;
	  ll2 = -1000000.0;
	}
      }
      
      /* std::cout << "i = " << i
	<< " Second Choice"
	      << " mean- = " << mean2[1] << " mean+ = " << mean2[2]
		  << " sig-  = " << sig2[1]
     		  << " sig+  = " << sig2[2]
		  << " ll2 = " << ll2
		  << std::endl;*/

      r = Rand1();

      if ( ll2 - ll1 > log(r)){ // Accept the new step!

	/*	std::cout << "Accept the new step: ll2 - ll1 = " << ll2 - ll1
		<< " >  log(r) = " << log(r) << std::endl;*/
	
	for (unsigned int j(1);j<=2;j++){
	  mean1[j] = mean2[j];
	  sig1[j] = sig2[j];
	}
	ll1 = ll2;
	yes += 1;
	acceptance = yes/(1.0*i);
	//para1 = para2;

	for (unsigned int j(1);j<=2;j++){
	  bin = (mean1[j]-meanlo[j])/dmean[j] + 1;
	  //std::cout << " bin mean = " << bin << std:endl;
	  probmean[j][bin] += 1;
	  bin = (sig1[j]-siglo[j])/dsig[j] + 1;
	  //std::cout << " bin mean = " << bin << std:endl;
	  probsig[j][bin] += 1;
	}
	if (i%10==0)
	  out << i << "\t"
	      << mean1[1] << "\t"
	      << sig1[1] << "\t"
	      << mean1[2] << "\t"
	      << sig1[2] << "\t"
	      << ll1 << "\t"
	      << acceptance << std::endl;
      }
      else{ // Reject the new step!

	/*	std::cout << "Reject the new step: ll2 - ll1 = " << ll2 - ll1
		<< " <  log(r) = " << log(r) << std::endl;*/

	
      }      
    } // End of M_T loop

      std::cout << "Average Value of the Mean-  = " << meanave[1]/(chainlength*1.0) << std::endl;
    std::cout << "Average Value of the StdDev- = " << sigave[1]/(chainlength*1.0) << std::endl;
     std::cout << "Average Value of the Mean+  = " << meanave[2]/(chainlength*1.0) << std::endl;
    std::cout << "Average Value of the StdDev+ = " << sigave[2]/(chainlength*1.0) << std::endl;
   

   out.open(out_filex);
   for (unsigned int i(1); i<=numbins; i++){
     out << xvals[i] << " "
	 << xprob[i] << std::endl;
   }
   out.close();

   out.open(out_filey);
   for (unsigned int i(1); i<=numbins; i++){
     out << yvals[i] << " "
	 << yprob[i] << std::endl;
   }
   out.close();

 

   xx = ProbContours1D(xvals,xprob,xave,numbins);
   yy = ProbContours1D(yvals,yprob,yave,numbins);

      histout.open("gplotx.txt");
      histout << "unset arrow" << std::endl;
      histout << "set arrow 1 from " << xx.xl3 << ",0 to " << xx.xl3 << "," << xx.pl3 << " nohead lt 5" << std::endl;
      histout << "set arrow from " << xx.xl2 << ",0 to " << xx.xl2 << "," << xx.pl2 << " nohead lt 4" << std::endl;
      histout << "set arrow from " << xx.xl1 << ",0 to " << xx.xl1 << "," << xx.pl1 << " nohead lt 2" << std::endl;
      histout << "set arrow from " << xx.median << ",0 to " << xx.median << "," << xx.pmed << " nohead lt 1 lw 3 " << std::endl;
      histout << "set arrow from " << xx.xr1 << ",0 to " << xx.xr1 << "," << xx.pr1 << " nohead lt 2" << std::endl;
      histout << "set arrow from " << xx.xr2 << ",0 to " << xx.xr2 << "," << xx.pr2 << " nohead lt 4" << std::endl;
      histout << "set arrow from " << xx.xr3 << ",0 to " << xx.xr3 << "," << xx.pr3 << " nohead lt 5" << std::endl;
      histout << "set title  \"2D Gaussian Test with " << nsteps << " values \"" << std::endl;
      histout << "set xlabel \"x\"" << std::endl;
      histout << "set ylabel \"Probability\"" << std::endl;
      histout << "set key left" << std::endl;
      histout << "plot \"" << out_filex << "\" with lines lt -1 notitle, 1/0 t \"2 sigma x = "
	      << xx.xl2 << "\" lt 4, 1/0 t \"1 sigma x = "
	      << xx.xl1 << "\" lt 2, 1/0 t \"Median  x = "
	      << xx.median << "\" lt 1 lw 3, 1/0 t \"1 sigma x = "
	      << xx.xr1 << "\" lt 2, 1/0 t \"2 sigma x = "
	      << xx.xr2 << "\" lt 4, \"" << out_filex << "\" using 1:3 with lines lt 1 notitle"
	      << std::endl;     
      histout.close();
   
      histout.open("gploty.txt");
      histout << "unset arrow" << std::endl;
      histout << "set arrow from " << yy.xl3 << ",0 to " << yy.xl3 << "," << yy.pl3 << " nohead lt 5" << std::endl;
      histout << "set arrow from " << yy.xl2 << ",0 to " << yy.xl2 << "," << yy.pl2 << " nohead lt 4" << std::endl;
      histout << "set arrow from " << yy.xl1 << ",0 to " << yy.xl1 << "," << yy.pl1 << " nohead lt 2" << std::endl;
      histout << "set arrow from " << yy.median << ",0 to " << yy.median << "," << yy.pmed << " nohead lt 1 lw 3 " << std::endl;
      histout << "set arrow from " << yy.xr1 << ",0 to " << yy.xr1 << "," << yy.pr1 << " nohead lt 2" << std::endl;
      histout << "set arrow from " << yy.xr2 << ",0 to " << yy.xr2 << "," << yy.pr2 << " nohead lt 4" << std::endl;
      histout << "set arrow from " << yy.xr3 << ",0 to " << yy.xr3 << "," << yy.pr3 << " nohead lt 5" << std::endl;
      histout << "set title  \"2D Gaussian Test with " << nsteps << " values \"" << std::endl;
      histout << "set xlabel \"y\"" << std::endl;
      histout << "set ylabel \"Probability\"" << std::endl;
      histout << "set key left" << std::endl;
      histout << "plot \"" << out_filex << "\" with lines lt -1 notitle, 1/0 t \"2 sigma y = "
	      << yy.xl2 << "\" lt 4, 1/0 t \"1 sigma y = "
	      << yy.xl1 << "\" lt 2, 1/0 t \"Median  y = "
	      << yy.median << "\" lt 1 lw 3, 1/0 t \"1 sigma y = "
	      << yy.xr1 << "\" lt 2, 1/0 t \"2 sigma y = "
	      << yy.xr2 << "\" lt 4, \"" << out_filex << "\" using 1:3 with lines lt 1 notitle"
	      << std::endl;     
      histout.close();

      //contours = ProbContours2D(xvals,yvals,probability,numbins,numbins);
      

   out.open("output.txt");
   for (unsigned int i(1); i<=numbins; i++){
     for (unsigned int j(1); j<=numbins; j++){

       out << xvals[i] << " " 
	   << yvals[j] << " " 
	   << probability[i][j] << " "
	   << std::endl;
     
     }
   }
   out.close();

 
      histout.open("gplot2d.txt");
      histout << "set xlabel \"x\"" << std::endl;
      histout << "set ylabel \"y\"" << std::endl;
      histout << "set dgrid3d 80,80" << std::endl;
      histout << "set contours base" << std::endl;
      histout << "set cntrparam levels discrete " << contours.sigma1 << "," << contours.sigma2 << std::endl;
      histout << "splot \"output.txt\" with lines notitle"<< std::endl;
      histout << "unset surface; set view map; replot" << std::endl;
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
