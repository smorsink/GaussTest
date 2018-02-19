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

  unsigned int nsteps(100);


  double
    xlo(-22.0), xhi(22.0);
  double
    ylo(-22.0), yhi(22.0);
  int
    numbins(500);

  struct Prob1d xx, yy;
  struct Prob2d contours;


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
      /* std::cout << " bin = " << i
		<< " bin centre = " << xvals[i] 
		<< " left edge = " << xvals[i] - 0.5*dx
		<< std::endl;*/

      //histogram[i] = 0;

    }
    for (unsigned int i(1);i<=numbins; i++)
      for (unsigned int j(1); j<=numbins; j++)
	histogram[i][j] = 0;


    int xbin, ybin;

    // Set up x-bins for storing the histogram

  double x1(0.0), x2(0.0);
  double y1(0.0), y2(0.0);

  double xave(0.0), yave(0.0);

  // Create Initial Step 
  x1=-2.0;
  y1=-2.0;

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

   // Normalize the probability
   for (unsigned int i(1); i<=numbins; i++){
     for (unsigned int j(1); j<=numbins; j++){
       probability[i][j] = histogram[i][j]/(1.0*totprob*dx*dy);
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
