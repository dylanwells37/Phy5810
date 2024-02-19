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

// include files
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
using namespace std;

#include "integ_routines.h"	// prototypes for integration routines

double test_integrand (double x);

double gsl_integrand (double x, void *params);

const double ME = 2.7182818284590452354E0;	// Euler's number 
//************************************************************************

int main () {
    const int max_intervals = 501;
    const double lower = 0.0;
    const double upper = 1.0;

    const double answer = 1. - (1. / ME);
    double result = 0.0;

    ofstream simpson_out("simpson.dat");
    simpson_out << "# N " << "result" << endl;

    for (int i = 3; i <= max_intervals; i += 2){
        result = simpsons_rule(i, lower, upper, &test_integrand);
        simpson_out << i << " " << result << " " << endl;
    }

    simpson_out.close();

    ofstream milne_out("milne.dat");

    milne_out << "# N " << "result" << endl;
    for (int i = 5; i <= max_intervals; i += 4){
        result = milne_rule(i, lower, upper, &test_integrand);
        milne_out << i << " " << result << endl;
    }

    milne_out.close();

    ofstream gsl_out("gsl.dat");

    gsl_out << "# N " << "result" << endl;
    for (int i = 3; i <= max_intervals; i += 2){
        result = gsl_integration(i, lower, upper, &gsl_integrand);
        gsl_out << i << " " << result << endl;
    }

    gsl_out.close();
    
    cout << answer << endl;
    return (0);
}

// Function to integrate
double test_integrand(double x){
    return (exp (-x));
}

// Function to integrate using gsl
double gsl_integrand(double x, void *){
    return test_integrand(x);
}

