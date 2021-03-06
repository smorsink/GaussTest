double Rand1();

double NormalDev(double mu, double sig);

double Gaussian(double x, double mu, double sig);

double Gaussian2d(double x, double y, double *mean, double **var);

// Compute the Probability Contours for a 1D distribution
struct Prob1d ProbContours1D(double *xvals, double *probability, double average, unsigned int numbins );

struct Prob2d ProbContours2D( double *rvals, double *mvals, double **probability, unsigned int numradius, unsigned int nummass);

void GaussProbDist(double *xvals, double *prob, double mean, double sig, unsigned int numbins);

double ProbCompare(double *prob1, double *prob2, int numbins);
double ProbCompare2d(double **prob1, double **prob2, int numbins);


void GaussProbDist2d(double *xvals, double *yvals, double **prob, double *mean, double *sig, unsigned int numbins);
