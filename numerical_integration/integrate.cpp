//  file: integrate.cpp               
//                                                                     
//  Programmer:  Dylan Wells (wells.1629@osu.edu)
//
//  Revision history:
//      2024-2-18 original version, for 5810 Computational Physics
//
//  Notes:  
//   
//
//************************************************************************
/*Analysis (b)
For milne and simpson, there seem to be 3 distinct regions in the log-log plot.
The first region begins when there are not enough points to accurately estimate the integral.
The lines are kind of jumpy, but the fitting lines show a slope of -3 for milne and -2 for simpson.
This indicates a loose initial power-law relation.
After this, we see a more consistent downward line, with a slope of -6 for milne and -4 for simpson.
These match the error terms of h^4 for simspon and h^6 for milne as discussed in the notes.
The region after this takes over when the step size is in the region of round-off errors and possibly even machine precision.
There isn't a set slope, but the estimating vary rapidly around the 10^-14 relative error.
I suspect that further increases in N could actually have an increase in the relative error as
these errors dominate further, but I wasn't able to run with more points.

The optimum number of points on the graph is between log(3.5) and log(4) or between 
~3,000 and 10,000.
Following the method in the notes, we want the approximation error to be
equal to the round off error for the best result. The approximation error
for milne is proportional to h^6, and the round off error goes as machine precision / h,
Equating these formulas, we can get that h should equal machine precision ^ 1/7
Then, we can find N by plugging this in for h. Approximating these values in desmos,
I get N = 1871.04307087. So, maybe I missed a coefficient for the approximation error
instead of just h^6?
*/


// include files
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
using namespace std;

#include "integ_routines.h"	// prototypes for integration routines

double test_integrand (double x);
double gsl_integrand (double x, void *params);
double relative_error(double result, double answer);

//************************************************************************

int main () {
    const int max_intervals = 50001;
    const double pi = 4.0*atan(1.0);
    const double lower = pi / 10.0;
    const double upper = 10.0;
    

    const double answer = 3.4422149121631; // From Mathematica
    double result = 0.0;

    ofstream simpson_out("simpson.dat");
    simpson_out << "# N " << "result" << endl;

    // Simpson needs 2i+1 intervals
    for (int i = 3; i <= max_intervals; i += 2){
        double h = (upper - lower) / (i - 1);
        simpson_out << setw(4) << log10(double(i));
        result = simpsons_rule(i, lower, upper, &test_integrand);
        simpson_out << " " << scientific << setprecision(15) << result << " " 
                    << log10(relative_error(result, answer))
                    << " " << log10(h) << endl;

    }

    simpson_out.close();

    ofstream milne_out("milne.dat");

    milne_out << "# N " << "result" << endl;
    // milne needs 4i+1 intervals
    for (int i = 5; i <= max_intervals; i += 4){
        double h = (upper - lower) / (i - 1);
        milne_out << setw(4) << log10(double(i));
        result = milne_rule(i, lower, upper, &test_integrand);
        milne_out << " " << scientific << setprecision(15) << result << " " 
                    << log10(relative_error(result, answer))
                    << " " << log10(h) << endl;
    }

    milne_out.close();

    ofstream gsl_out("gsl.dat");

    gsl_out << "# N " << "result" << endl;
    // gsl needs at least 10 intervals or it errors
    for (int i = 10; i <= max_intervals; i += 2){
        double h = (upper - lower) / (i - 1);
        gsl_out << setw(4) << log10(double(i));
        result = gsl_integration(i, lower, upper, &gsl_integrand);
        gsl_out << " " << scientific << setprecision(15) << result << " " 
                    << log10(relative_error(result, answer))
                    << " " << log10(h) << endl;
    }

    gsl_out.close();
    
    return (0);
}

// Function to integrate
double test_integrand(double x){
    // (sin(10x) * ln(x)) / (x * tan(x))
    return (sin(10.0 * x) * log(x)) / (x * tan(x));
}

// Function to integrate using gsl
double gsl_integrand(double x, void *){
    return test_integrand(x);
}

double relative_error(double result, double answer){
    return fabs((result - answer) / answer);
}
