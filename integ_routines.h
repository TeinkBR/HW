//  file: integ_routines.h
// 
//  Header file for integ_routines.cpp
//
//
//  Programmer:  Dick Furnstahl  furnstahl.1@osu.edu
//
//  Revision History:
//    05-Jan-2004 --- original version, based on C version
//    08-Jan-2005 --- function to be integrated now passed, changed names
//    09-Jan-2011 --- changed function names
//
//  To do:
//
//************************************************************************

//  begin: function prototypes 
 
extern float trapezoid_rule ( int num_pts, float x_min, float x_max, 
                         float (*integrand) (float x) );    // trapezoid rule 
extern float simpsons_rule ( int num_pts, float x_min, float x_max, 
                       float (*integrand) (float x) );    // Simpson's rule 
extern float gauss_quadrature( int num_pts, float x_min, float x_max, 
                       float (*integrand) (float x) );    // Gauss' rule 
   
extern void gauss(int npts, int job, double a, double b, 
                  double x[], double w[]);              // from gauss.cpp 

//  end: function prototypes 
