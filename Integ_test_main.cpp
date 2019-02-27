#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
using namespace std;

#include "integ_routines.h"	// prototypes for integration routines

float my_integrand (float x);

const double ME = 2.7182818284590452354E0;	// Euler's number

//************************************************************************

int
main ()
{
  // set up the integration specifiction
  const int max_intervals = 501;	// maximum number of intervals
  const float lower = 0.0;	// lower limit of integration
  const float upper = 1.0;	// upper limit of integration

  const double answer = 1. - 1. / ME;	// the "exact" answer for the test
  float result = 0.;  // approximate answer

  // open the output file stream
  ofstream integ_out ("integ.dat");	// save data in integ.dat
  integ_out << "#  N   Milne     Simpsons      GSL  " << endl;
  integ_out << "#-----------------------------------------" << endl;

  // Simpson's rule requires an odd number of intervals
  for (int i = 3; i <= max_intervals; i += 2)
  {
    integ_out << setw(4) << i;

    result = milne_rule (i, lower, upper, &my_integrand);
    integ_out << "  " << scientific << fabs (result - answer);

    result = simpsons_rule (i, lower, upper, &my_integrand);
    integ_out << "  " << scientific << fabs (result - answer);

    result = gsl_function (i, lower, upper, &my_integrand);
    integ_out << "  " << scientific << fabs (result - answer);

    integ_out << endl;
  }

  cout << "data stored in integ.dat\n";
  integ_out.close ();

  return (0);
}

//************************************************************************

// the function we want to integrate
float
my_integrand (float x)
{
  return (exp (-x));
}
