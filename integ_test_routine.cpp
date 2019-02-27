#include <cmath>
#include <iomanip>
#include "integ_routines.h"


//-----------------milne_rule--------------------------------------------------
float milne_rule ( int num_pts, float x_min, float x_max,
                      float (*integrand) (float x) )
{
   float interval = ((x_max - x_min)/float(num_pts - 1));  // called h in notes
   float sum=  0.;  // initialize integration sum to zero

   for (int n=2; n<num_pts; n+=4)                // loop for 2+4 points
   {
     float x = x_min + interval * float(n-1);
     sum += (32./45.)*interval * integrand(x);
   }

  for (int n=3; n<num_pts; n+=4)                // loop for 3+4 points
   {
     float x = x_min + interval * float(n-1);
     sum += (12./45.)*interval * integrand(x);
   }

   sum +=  (interval*14./45.) * (integrand(x_min) + integrand(x_max));
   for (int n=4; n<num_pts; n+=4)
   {
     float x = x_min + interval * float(n-1);
     sum += (12./45.)*interval * integrand(x);
    }
   for (int n=5; n<num_pts; n+=4)
   {
      float x = x_min + interval * float(n-1);
      sum += (32./45.)*interval * integrand(x);

   }

   // add in the endpoint contributions
   sum +=  (interval/45.) * (integrand(x_min) + integrand(x_max));

   return (sum);
}
//************************simpsons_rule***************************************

// Integration using Simpson's rule
float simpsons_rule ( int num_pts, float x_min, float x_max,
                      float (*integrand) (float x) )
{
   float interval = ((x_max - x_min)/float(num_pts - 1));  // called h in notes
   float sum=  0.;  // initialize integration sum to zero

   for (int n=2; n<num_pts; n+=2)                // loop for odd points
   {
     float x = x_min + interval * float(n-1);
     sum += (4./3.)*interval * integrand(x);
   }
   for (int n=3; n<num_pts; n+=2)                // loop for even points
   {
     float x = x_min + interval * float(n-1);
     sum += (2./3.)*interval * integrand(x);
   }
   // add in the endpoint contributions
   sum +=  (interval/3.) * (integrand(x_min) + integrand(x_max));

   return (sum);
}
//************************using_gsl*******************************************
{
  gsl_integration_workspace *work_ptr
    = gsl_integration_workspace_alloc (1000);

  double lower_limit = 0;	/* lower limit a */
  double upper_limit = 1;	/* upper limit b */
  double abs_error = 1.0e-8;	/* to avoid round-off problems */
  double rel_error = 1.0e-8;	/* the result will usually be much better */
  double result;		/* the result from the integration */
  double error;			/* the estimated error from the integration */

  double alpha = 1.0;		// parameter in integrand
  double expected = -4.0;	// exact answer

  gsl_function My_function;
  void *params_ptr = &alpha;

  My_function.function = &my_integrand;
  My_function.params = params_ptr;

  gsl_integration_qags (&My_function, lower_limit, upper_limit,
			abs_error, rel_error, 1000, work_ptr, &result,
			&error);

  cout.setf (ios::fixed, ios::floatfield);	// output in fixed format
  cout.precision (18);		// 18 digits in doubles

  int width = 20;  // setw width for output
  cout << "result          = " << setw(width) << result << endl;
  cout << "exact result    = " << setw(width) << expected << endl;
  cout << "estimated error = " << setw(width) << error << endl;
  cout << "actual error    = " << setw(width) << result - expected << endl;
  cout << "intervals =  " << work_ptr->size << endl;

  return 0;
}
