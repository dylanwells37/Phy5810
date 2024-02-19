//  file: integ_routines.cpp
//
//  Integration routines for Simpson's, Milne, and GSL integration           
//                                                                     
//  Programmer:  Dylan Wells (wells.1629@osu.edu)
//
//  Revision history:
//      2024-2-18 original version, for 5810 Computational Physics
//
// 
//************************************************************************

// include files
#include <cmath>
#include "integ_routines.h"   // integration routine prototypes 
#include <gsl/gsl_integration.h>

double simpsons_rule (int num_points, double x_min, 
                      double x_max, double (*integrand) (double x) ) {
    /* Calcualte the simpson's rule estimation for a given integrand,
    num_points, and range */

    double h = (x_max - x_min) / (num_points - 1); // step size
    double sum = 0;

    // sum the even points
    for (int n = 2; n < num_points; n += 2){
        double x = x_min + h * (n - 1);
        sum += (4.0 / 3.0) * h * integrand(x);
    }
    // sum the odd points
    for (int n = 3; n < num_points; n += 2){
        double x = x_min + h * (n - 1);
        sum += (2.0 / 3.0) * h * integrand(x);
    }
    // add in the endpoint contributions
    sum += (h / 3.0) * (integrand(x_min) + integrand(x_max));

    return sum;
}

double milne_rule (int num_points, double x_min, 
                   double x_max, double (*integrand) (double x) ) {
    /* Calcualte the milne rule estimation for a given integrand,
    num_points, and range */
    
    double h = (x_max - x_min) / (num_points - 1); // step size
    double sum = 0;

    // add even contribution (always 64h/45 * f(x))
    for (int n = 2; n < num_points; n += 2){
        double x = x_min + h * (n - 1);
        sum += (64.0 / 45.0) * h * integrand(x);
    }
    // Add odd contributions
    for (int n = 3; n < num_points; n += 4){
        double x = x_min + h * (n - 1);
        sum += (24.0 / 45.0) * h * integrand(x);
    }
    // add odd contributions counted twice
    for (int n = 5; n < num_points; n+=4){
        double x = x_min + h * (n - 1);
        sum += (28.0 / 45.0) * h * integrand(x);
    }

    // add endpoint contribution
    sum += (14.0 / 45.0) * h * (integrand(x_min) + integrand(x_max));

    return sum;
}

double gsl_integration (int num_points, double x_min, 
                        double x_max, double (*gsl_integrand) (double x, void *) ){
    /* Calcualte the gsl integration estimation for a given integrand,
    num_points, and range */

    gsl_integration_workspace * w = gsl_integration_workspace_alloc (num_points);
    double result, error;

    gsl_function F;

    F.function = gsl_integrand;

    gsl_integration_qags (&F, x_min, x_max, 0, 1e-7, num_points, w, &result, &error);

    gsl_integration_workspace_free (w);

    return result;
}