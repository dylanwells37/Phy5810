//  file: integ_routines.h
// 
//  Header file for integ_routines.cpp
//
//
//  Programmer:  Dylan Wells (wells.1629@osu.edu)
//
//  Revision History:
//    2024-2-18 --- original version
//
//************************************************************************

// Function Prototypes

extern double simpsons_rule (int num_points, double x_min,
                            double x_max, double (*integrand) (double x));

extern double milne_rule (int num_points, double x_min,
                            double x_max, double (*integrand) (double x));

extern double gsl_integration (int num_points, double x_min,
                            double x_max, double (*integrand) (double x));
