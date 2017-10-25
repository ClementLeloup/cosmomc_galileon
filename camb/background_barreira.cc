// Modified by Clement Leloup 26-04-2016 :
// Adding the initCond and ageOfUniverse functions
// Commenting previous main function and writing new one that just calculates initial conditions in x and H from initial conditions in rho and Omega_m
// Modified by Clement Leloup 13-04-2016 :
// Adding the "coord" parameter to calcHubbleGalileon, to take into account whether we use a or z to do the calculation. In case one use a to calculate, the calculation is done backward from a=1.
// Modification of the main function. Now, the behaviour is adapted whether the user prefers a or z to calculate the observables.

// Modified by Clement Leloup 07-04-2016 :
// Adding the contribution of orad to the evolution of Lander system when !useacoord
// Writing of the main function to test and compare with former results
// Commenting of every function except CalcValOmC2C3C4C5CGC0 and calcHubbleGalileon

// From modified CosFitter for Galileon Cosmology file lumdist.cc by Jeremy Neveu


#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_errno.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include "fstream"
#include <algorithm>
#include <cmath>
#include "string.h"

using namespace std;


/*! 
   \ingroup lumdist
  Not in lumdist class because the GSL does not accept pointers to member
  functions.
  \param[in] z Redshift
  \param[in] y Current value of integral : y[0] is h(z), y[1] is dpi/dz and y[2] is pi
  \param[out] f The value of the integral term : f[0] si dh/dz, f[1] is d^2pi/dz^2 and f[2] is dpi/dz
  \param[in] params Array of parameters 
     (\f$\Omega_m, c_2, c_3, c_4, c_5, c_G, c_0, orad, useacoord\f$)

  Tested with Mathematica with c0 and y0 = 0.

 */

/*
inline int calcPertOmC2C3C4C5CGC0(double z, const double y[3], 
			       double f[3], void* params) {
  double k1,k2,k3,k4,k5,k6;
  const double *in_params = static_cast<const double*>(params);
  //double om = in_params[0];
  double c2 = in_params[1];
  double c3 = in_params[2];
  double c4 = in_params[3];
  double c5 = in_params[4];
  double cG = in_params[5];
  double c0 = in_params[6];
  bool useacoord = (int)in_params[8];
  double zpo;
  if (useacoord==0) {
    zpo=(1.0+z);
  } else {
    // We consider zpo as -a :
    zpo=-z;
  }
  int status = 0;

  // system from lna in z : 
  // y[0](lna) -> y[0](z)
  // y[1](lna) -> -(1+z)*y[1](z)
  // f[0](lna) -> -(1+z)*f[0](z)
  // f[1](lna) -> (1+z)^2*f[1](z)+(1+z)*y[1](z)

  double prod = -zpo*y[0]*y[1];
  double prod2 = prod*prod;
  double prod3 = prod*prod2;
  double prod4 = prod*prod3;
  double h = y[0];
  double h2 = h*h;
  double h3 = h*h2;
  double dh = -zpo*f[0];
  double x = -zpo*y[1];
  double dx = zpo*zpo*f[1]+zpo*y[1];
  if (useacoord==1) dx += -2*zpo*y[1];

  k1 = -6*c4*h*prod2*(dh*x + h*dx + prod/3.0) -2*c0 + c5*h2*prod3*(12*h*dx + 15*dh*x + 3*prod) + 2.0*cG*(dh*prod+h2*dx+h*prod);

  k2 = -c2/2.0 + 6*c3*h*prod - 27*c4*h2*prod2 + 30*c5*h3*prod3 + 3.0*cG*h2;

  k3 = -(1.0 - 2.0*c0*y[2]) - c4/2.0*prod4 - 3*c5*h*prod4*(h*dx + dh*x) + cG*prod2;

  k4 = -2.0*(1.0 - 2.0*c0*y[2]) + 3.0*c4*prod4 - 6*c5*h*prod*prod4 - 2.0*cG*prod2;

  k5 = 2*c3*prod2 - 12*c4*h*prod3 -2.0*c0 + 15.0*c5*h2*prod4 + 4.0*cG*h*prod;

  k6 = c2/2.0 - 2.0*c3*( h2*dx + h*dh*x + 2*h*prod ) + c4*( 12*h3*prod*dx + 18*h*prod2*dh + 13*h2*prod2 ) - c5*( 18*h2*h2*prod2*dx + 30*h2*prod3*dh + 12*h3*prod3 ) -cG*( 2.0*h*dh + 3.0*h2 );

  double noghost = k2 + 1.5*k5*k5/k4;
  double cs2 = (4*k1*k4*k5 - 2*k3*k5*k5 - 2*k4*k4*k6)/(k4*(2*k4*k2 + 3*k5*k5));
  Geff = 4*( k3*k6 - k1*k1 ) / ( k5*(k4*k1-k5*k3) - k4*(k4*k6 - k5*k1) );
  // Here we can't have a condition noghost>0 because with some set of
  // parameters we have noghost=1e-16, condition passes anyway
  // and then rk4 fails and gives NaN values.
  if ( noghost>-1e-8 ) status = 6;
  if ( cs2<0 ) status = 7;

  // Check instabilities for tensorial modes (cf Felice & Tsujikawa)
  double noghostt = 0.5 - 0.75*c4*prod4 + 1.5*c5*prod4*prod*h + 0.5*cG*prod2 -c0*y[2];
  if ( noghostt < 0 )  status = 2;
  double ct2 = (0.5 + 0.25*c4*prod4 + 1.5*c5*prod4*h*(h*dx+dh*x) - 0.5*cG*prod2 -c0*y[2]) / (0.5 - 0.75*c4*prod4 + 1.5*c5*prod4*prod*h + 0.5*cG*prod2 -c0*y[2]);
  if ( ct2 < 0) status = 3;
  return status;

}
*/


/*! 
   \ingroup lumdist
  Not in lumdist class because the GSL does not accept pointers to member
  functions.
  \param[in] z Redshift
  \param[in] y Current value of integral : y[0] is h(z), y[1] is dpi/dz and y[2] is pi
  \param[out] f The value of the integral term : f[0] si dh/dz, f[1] is d^2pi/dz^2 and f[2] is dpi/dz
  \param[in] params Array of parameters 
     (\f$\Omega_m, c_2, c_3, c_4, c_5, c_G, c_0, orad, useacoord\f$)

  Tested with Mathematica with c0 and y0 = 0.

 */
int calcValOmC2C3C4C5CGC0(double z, const double y[3], 
			       double f[3], void* params) {

  double alpha,gamma,beta,sigma,lambda,omega,OmegaP,OmegaM,OmTest;
  const double *in_params = static_cast<const double*>(params);
  double om = in_params[0];
  double c2 = in_params[1];
  double c3 = in_params[2];
  double c4 = in_params[3];
  double c5 = in_params[4];
  double cG = in_params[5];
  double c0 = in_params[6];
  double orad = in_params[7];
  bool useacoord = (int)in_params[8];
  double *OmegaPiOverOrad = NULL;
  if (useacoord) OmegaPiOverOrad = const_cast<double*>(&in_params[9]);
  //if (useacoord) std::cout << "Using a instead of z" << std::endl; 
  double zpo;
  double a=1;
  if (!useacoord) {
    zpo=(1.0+z);
  } else {
    // We consider z as a :
    a = z;
    zpo=-a;
  }
  int status = 0;

  // system of Linder from lna to z :
  // y[0](lna) -> y[0](z)
  // y[1](lna) -> -zpo*y[1](z)
  // f[0](lna) -> -zpo*f[0](z)
  // f[1](lna) -> zpo*zpo*f[1](z)+zpo*y[1](z)

  // system of Linder from lna to a :
  // y[0](lna) -> y[0](z)
  // y[1](lna) -> a*y[1](z)
  // f[0](lna) -> a*f[0](z)
  // f[1](lna) -> a*a*f[1](z)+a*y[1](z)
  // OmegaP(lna) = 3*h2*rhoP(lna)
  // So :

  double prod = -zpo*y[0]*y[1]; // prod=h(z)x(z)
  double prod2 = prod*prod;
  double prod3 = prod*prod2;
  double prod4 = prod*prod3;
  double prod5 = prod*prod4;
  double h = y[0];
  double h2 = h*h;
  double h3 = h*h2;
  double h4 = h2*h2;

  OmegaP = ( 6.0*c0*(h*prod+y[2]*h2) + 0.5*c2*prod2 - 6.0*c3*h*prod3 + 22.5*c4*h2*prod4 - 21.0*c5*h3*prod5 - 9*cG*h2*prod2  )/(3.0*h2) ;
  if (!useacoord) {
    OmegaM = 1 - OmegaP - orad*pow(zpo, 4)/h2;
    OmTest = om*zpo*zpo*zpo/h2; 
  } else {
    OmegaM = 1 - OmegaP - orad/(pow(a,4)*h2);
    OmTest = om/(pow(a,3)*h2);
  }
  //if ( abs(OmegaM - OmTest)>1e-5  ) {
    //std::cout<<z<<" "<<om<<" "<<c2<<" "<<c3<<" "<<c4<<" "<<c5<<" "<<cG<<" "<<c0<<" "<<OmTest<<" "<<OmegaM<<" "<<OmTest<<" "<<OmegaP<<" "<<orad<<" "<<a<<" "<<h2<<" "<<abs(OmegaM - OmTest)<<std::endl;
    //return 4;
  //}
  //if ( OmegaP<0 ) 
  //return 5;
  // if ( useacoord && OmegaP>0.1 && abs(1/z-1000.0)<1e-8  ){
  //   std::cout<<"Warning : OmegaP is more than 10% of Omega_tot at z=1000 : OmegaP="<<OmegaP<<"\n\t with om="<<om<<" c2="<<c2<<" c3="<<c3<<" c4="<<c4<<" c5="<<c5<<" cG="<<cG<<std::endl;
  //   std::cout<<"\tBeware of theoretical validity of some computations!"<<std::endl;
  // }

  // The equations : 
  alpha = -3*c3*h*prod2 + 15*c4*h2*prod3 + c0*h + c2/6.0*prod - 17.5*c5*h3*prod4 - 3.0*cG*h2*prod;
  gamma = 2*c0*h2 - c3*h2*prod2 + c2/3.0*h*prod + 2.5*c5*h4*prod4 - 2.0*cG*h3*prod;
  beta = -2*c3*h3*prod + c2/6.0*h2 + 9*c4*h4*prod2 - 10*c5*h4*h*prod3 - cG*h4;
  sigma = 2.0*( 1.0 - 2*c0*y[2] )*h - 2.0*c0*prod + 2.0*c3*prod3 - 15.0*c4*h*prod4 + 21.0*c5*h2*prod5 + 6.0*cG*h*prod2;
  lambda = 3.0*( 1 - 2*c0*y[2] )*h2 - 2.0*c0*h*prod - 2.0*c3*h*prod3 + c2/2.0*prod2 + 7.5*c4*h2*prod4 - 9.0*c5*h3*prod5 - cG*h2*prod2;
  omega = -2*c0*h2 + 2*c3*h2*prod2 - 12*c4*h3*prod3 + 15*c5*h4*prod4 + 4.0*cG*h3*prod;
  if (useacoord) 
  {
    lambda += orad/pow(a,4);
  }
  else {
    lambda += orad*pow(zpo,4);
  }

  double denom = sigma*beta - alpha*omega;
  f[0] = (omega*gamma - lambda*beta) / (-zpo*denom);
  f[1] = (alpha*lambda - sigma*gamma) / (zpo*zpo*denom) ;
  f[2] = y[1];///(-zpo);
  if (useacoord) f[1] += -2.0*y[1]/a;

  // Compute and test perturbations :
  //status = calcPertOmC2C3C4C5CGC0(z, y, f, params) ; 
  //if (status != 0) return status;
  if (a==0.001 && useacoord) {
    *OmegaPiOverOrad = OmegaP / ( orad/(pow(a,4)*h2) );
  }
  if (useacoord) double dhda = f[0];
  return GSL_SUCCESS;
}



/*
  Evaluates the age of the universe in galileon cosmology with trapezoidal method.
  Be careful, in here, we only work with a coordinate !
  
  \param[in] cond \f$cond\f$, initial conditions array (H, xi, pi, age)
  \param[in] om \f$\Omega_m\f$
  \param[in] orad \f$\Omega_r\f$
  \param[in] c2 \f$c_2\f$, parameter of galileon model
  \param[in] c3 \f$c_3\f$, parameter of galileon model
  \param[in] c4 \f$c_4\f$, parameter of galileon model
  \param[in] c5 \f$c_5\f$, parameter of galileon model
  \param[in] cG \f$c_G\f$, parameter of galileon model
  \param[in] c0 \f$c_0\f$, parameter of galileon model
  \param[out] table of a coordinates where h(a) is evaluated
  \param[out] age of universe for given cosmology
  \returns p0_status :
                       \status = 0 : ok
		       \status = 1 : don't give a de Sitter solution (respect to equations A7 and A8 oof Gannoudji Sami)
		       \status = 2 : tensorial no ghosr condition not satisfied
		       \status = 3 : tensorial laplace condition not satisfied
		       \status = 4 : 00 Einstein equation not fulfilled
		       \status = 5 : rhoP<0
		       \status = 6 : scalar noghost condition not satified
		       \status = 7 : scalar imaginary speed of sound
*/


int ageOfUniverse(std::vector<double> &hubble, std::vector<double> &x, double &age, std::vector<double> &acoord, double cond[4], double om, double orad, double c2, double c3, double c4, double c5, double cG, double c0 ) {

  const double Gyr = 3.1556926e16; //Conversion factor from Gyr to sec
  const double Mpc = 3.085678e19; //Conversion factor from Mpc to km
  const double h0 = 71; // Present-day Hubble constant (WMAP7)

  double params[10];

  params[0] = om; params[1] = c2; params[2] = c3;
  params[3] = c4; params[4] = c5; params[5] = cG;
  params[6] = c0; params[7] = orad; params[8]= 1; params[9] = 0.0;
  int p0_status = 0;
  //p0_status = calcInitialValueForPi0(om,c2,c3,c4,c5,cG);	

  std::cout << cond[1] << endl;

  // Vector y = ( dh/dz, dx/dz, dy/dz )
  double y[3] = {cond[0], cond[1], cond[2]};  //Inital value of integral
  hubble[0] = cond[0];
  x[0] = cond[1];

  //std::cout<<"calcIntHubbleGalileon"<<std::endl;
  //std::cout<<"....Used memory = "<<ReadLinuxMem()<<" Mo"<<std::endl; 
  if ( p0_status > 1 ) return p0_status;

  const gsl_odeiv_step_type * T = gsl_odeiv_step_rkf45;
  gsl_odeiv_step * s = gsl_odeiv_step_alloc (T, 3);
  gsl_odeiv_control * c = gsl_odeiv_control_standard_new(1e-14, 1e-14, 1, 1);
  gsl_odeiv_evolve * e = gsl_odeiv_evolve_alloc (3);

  gsl_odeiv_system sys;
  sys.function = calcValOmC2C3C4C5CGC0;
  sys.dimension = 3;
  sys.params = &params;

  double a = acoord[1];
  int nstep = acoord.size();
  double h = 1e-6; //Initial step guess
  double acurrtarg;
  double prec_age = 1/(acoord[1]*cond[0]);
  int st;

  age = cond[3];
  //age = 0;

  for (int i = 2; i < nstep; i++) {
    acurrtarg = acoord[i];
    while (a < acurrtarg) {
      st = gsl_odeiv_evolve_apply(e, c, s, &sys, &a, acurrtarg, &h, y);
      double OmegaP = (0.5*c2*pow(y[0], 2)*pow(acoord[i]*y[1], 2) - 6*c3*pow(y[0], 4)*pow(acoord[i]*y[1], 3) + 22.5*c4*pow(y[0], 6)*pow(acoord[i]*y[1], 4) - 21*c5*pow(y[0], 8)*pow(acoord[i]*y[1], 5) - 9*cG*pow(y[0], 4)*pow(acoord[i]*y[1], 2))/(3.0*pow(y[0], 2)) ;
      double OmegaM = 1 - OmegaP - orad/(pow(a,4)*pow(y[0], 2));
      double OmTest = om/(pow(a,3)*pow(y[0], 2));
      if ( abs(OmegaM - OmTest)>1e-4  ) {
	std::cout<<om<<" "<<c2<<" "<<c3<<" "<<c4<<" "<<c5<<" "<<cG<<" "<<c0<<" "<<OmTest<<" "<<OmegaM<<" "<<OmTest<<" "<<OmegaP<<" "<<orad<<" "<<acoord[i]<<" "<<pow(y[0], 2)<<" "<<abs(OmegaM - OmTest)<<std::endl;
	st = 4;
      }
      if (st != 0) {
	gsl_odeiv_evolve_free (e);
	gsl_odeiv_control_free (c);
	gsl_odeiv_step_free (s);
	return st;
      }
    }
    if ( std::isnan(abs(y[0])) || std::isnan(abs(y[1])) || std::isnan(abs(y[2])) ) {
      gsl_odeiv_evolve_free (e);
      gsl_odeiv_control_free (c);
      gsl_odeiv_step_free (s);
      //std::cout<<"\nFailure with om="<<om<<" c2="<<c2<<" c3="<<c3<<" c4="<<c4<<" c5="<<c5<<" cG="<<cG<<" c0="<<c0<<" at z="<<zcurrtarg<<" h="<<y[0]<<" x="<<y[1]<<" y="<<y[2]<<std::endl;
      return 8;
    }
    age += 0.5*(1.0/(acoord[i]*y[0])+prec_age)*(acoord[i]-acoord[i-1]);
    prec_age = 1/(acoord[i]*y[0]);
    hubble[i-1] = y[0];
    x[i-1] = y[1];
  }

  gsl_odeiv_evolve_free (e);
  gsl_odeiv_control_free (c);
  gsl_odeiv_step_free (s);


  // Convert the age of the Universe in Gyr
  age *= Mpc/Gyr/h0;


  return p0_status;

}



/*
  Determine the initial condition in x from the initial conditions in rho and H at z=10^6
  cf Barreira et al. 2013 : Linear perturbations in Galileon gravity models
  
  \param[in] H \f$H\f$, { value of H one step before the initial step, initial value of H}
  \param[in] ratio_rho \f$ratio_\rho\f$, initial value of rho_p/rho_m
  \param[in] age \f$Age of Universe\f$, age of the universe in Gyr
  \param[in] om \f$\Omega_m\f$
  \param[in] orad \f$\Omega_r\f$
  \param[in] c2 \f$c_2\f$, parameter of galileon model
  \param[in] c3 \f$c_3\f$, parameter of galileon model
  \param[in] c4 \f$c_4\f$, parameter of galileon model
  \param[in] c5 \f$c_5\f$, parameter of galileon model
  \param[in] cG \f$c_G\f$, parameter of galileon model
  \param[in] c0 \f$c_0\f$, parameter of galileon model
  \param[out] initial value of x
*/

struct AgeofU
{
  int i; // Position in the vector
  double a; // Age of Universe

  AgeofU(int pos, double u_age) : i(pos), a(u_age) {}
};

struct less_than_key
{
  // To sort by the difference to the input age of universe, while keeping the initial position
  inline bool operator() (const AgeofU& struct1, const AgeofU& struct2)
  {
      return (struct1.a < struct2.a);
  }
};

int initCond(std::vector<double> &hubble, std::vector<double> &x, double &xi, std::vector<double> &acoord, double H[2], double ratio_rho, double age, double om, double orad, double c2, double c3, double c4, double c5, double cG, double c0)
{

  const double Gyr = 3.1556926e16; //Convertion factor from Gyr to sec
  const double Mpc = 3.085678e19; //Convertion factor from Mpc to km
  int status = 0;
  double coeff[6]; //array of coefficients of the polynomial
  size_t n = sizeof(coeff)/sizeof(double);
  double cond[4]; // initial conditions

  // Structure to store the ages of the Universe for all the roots of the polynomial
  //struct ageofU {
    //int i;
    //double a;
    //static bool less( const ageofU &a1, const ageofU &a2) {return abs(a1.a) < abs(a2.a);}
  //};

  AgeofU ini_vec_age = AgeofU(0, 1000);
  //ini_vec_age.i = 0; ini_vec_age.a = 1234569999;
  std::vector<AgeofU> aou(5, ini_vec_age);

  // coeff[0] = -3*ratio_rho*om/pow(acoord[1], 3); coeff[1] = 0; coeff[2] = c2/2*pow(H[1], 2)+9*cG*pow(H[1], 4);
  // coeff[3] = 6*c3*pow(H[1], 4); coeff[4] = 45/2*c4*pow(H[1], 6); coeff[5] = 21*c5*pow(H[1], 8);
  coeff[0] = -3*ratio_rho*om/pow(acoord[1], 3); coeff[1] = 0; coeff[2] = 0.5*c2*pow(H[1], 2)-9*cG*pow(H[1], 4);
  coeff[3] = -6*c3*pow(H[1], 4); coeff[4] = 22.5*c4*pow(H[1], 6); coeff[5] = -21*c5*pow(H[1], 8);

  gsl_poly_complex_workspace* w = gsl_poly_complex_workspace_alloc(n);
  double z[2*(n-1)]; // Array of complex numbers ordered like {Re[0],Im[0],Re[1],Im[1],...}

  status = gsl_poly_complex_solve(coeff, n, w, z); // Solve coeff[0]+coeff[1]*x+...+coeff[5]*x^5=0 and store sol in z

  cond[0] = H[1]; cond[2] = 0; cond[3] = Mpc/Gyr*(1/H[1]-pow(acoord[1], 2)/2*(1/(acoord[0]*H[0])-1/(acoord[1]*H[1]))/(acoord[0]-acoord[1])); // Initial value of the age of the Universe, extrapolating linearly the integral

  std::cout << "Initial value : " << cond[3] << endl;

  for(int j = 0; j<n-1; j++)
  {
    std::cout << "Root number " << j << " is " << z[2*j] << " + " << z[2*j+1] << "i" << endl;
    if(z[2*j+1]==0)
    {
      //cond[1] = z[2*j];
      cond[1] = z[2*j]/acoord[1];
      //cond[1] = 3.35e-14/acoord[1];
      aou[j].i = j;
      int status2 = ageOfUniverse(hubble, x, aou[j].a, acoord, cond, om, orad, c2, c3, c4, c5, cG, c0);
      aou[j].a -= age;
      if(status2 != 0) return status2;
    }
  }

  for(int j = 0; j<n-1; j++)
  {
    std::cout << "Root number " << aou[j].i << " gives age of universe equal to " << aou[j].a + age << endl; 
  }

  sort(aou.begin(), aou.end(), less_than_key());

  for(int j = 0; j<n-2; j++)
  {
    aou[j].a += age;
  }

  xi = z[2*(aou[0].i)];

  std::cout << "size of aout : " << aou.size() << endl;
  std::cout << "Initial condition for x is " << xi << ", it is the root number " << aou[0].i << " and it gives an age of the universe equal to " << aou[0].a << endl;

  return status;

} 



/*!
  Evaluates the h(z) in galileon cosmology
  
  \param[in] om \f$\Omega_m\f$
  \param[in] orad \f$\Omega_r\f$
  \param[in] c2 \f$c_2\f$, parameter of galileon model
  \param[in] c3 \f$c_3\f$, parameter of galileon model
  \param[in] c4 \f$c_4\f$, parameter of galileon model
  \param[in] c5 \f$c_5\f$, parameter of galileon model
  \param[in] cG \f$c_G\f$, parameter of galileon model
  \param[in] c0 \f$c_0\f$, parameter of galileon model
  \param[in] coord, parameter telling whether to use a or z coordinate
  \param[out] table of z coordinates where h(z) and x(z) are evaluated
  \param[out] table of values of h(z)
  \param[out] table of values of x(z)
  \returns p0_status :
                       \status = 0 : ok
		       \status = 1 : don't give a de Sitter solution (respect to equations A7 and A8 oof Gannoudji Sami)
		       \status = 2 : tensorial no ghost condition not satisfied
		       \status = 3 : tensorial laplace condition not satisfied
		       \status = 4 : 00 Einstein equation not fulfilled
		       \status = 5 : rhoP<0
		       \status = 6 : noghost condition not satisfied
		       \status = 7 : imaginary speed of sound
*/
int calcHubbleGalileon(std::vector<double>& hubble, std::vector<double>& x, std::vector<double>& zcoord, double om, double orad, double c2, double c3, double c4, double c5, double cG, double c0, double coord ) {

  double params[10];
  params[0] = om; params[1] = c2; params[2] = c3;
  params[3] = c4; params[4] = c5; params[5] = cG;
  params[6] = c0; params[7] = orad; params[8] = coord; params[9] = 0; 

  bool useacoord = (int)params[8];

  int p0_status = 0;
  //p0_status = calcInitialValueForPi0(om,c2,c3,c4,c5,cG);	
  if ( p0_status > 1 ) return p0_status;

  // Warning : in z x(0)/pi'0=-1
  //double y[3] = { 1.0, 1.0, 0.0 };  //Inital value of integral



  const gsl_odeiv_step_type * T = gsl_odeiv_step_rkf45;
  gsl_odeiv_step * s = gsl_odeiv_step_alloc(T, 3);
  //gsl_odeiv_control * c = gsl_odeiv_control_y_new (prec_lum, prec_lum_rel); // keep absolute error under first parameter, and relative eror under second parameter wrt y
  //gsl_odeiv_control * c = gsl_odeiv_control_y_new(1, 1e-10);
  gsl_odeiv_control * c = gsl_odeiv_control_standard_new(1e-14, 1e-14, 1, 1);
  gsl_odeiv_evolve * e = gsl_odeiv_evolve_alloc(3);
  gsl_odeiv_system sys;
  sys.function = calcValOmC2C3C4C5CGC0;
  sys.dimension = 3;
  sys.params = &params;
 
  if(useacoord)
  {
    double y[3] = { 1.0, 1.0, 0.0 };  //Inital value of integral
    double z = 1;
    int nstep = min(hubble.size(),x.size());
    double h = -1e-6; //Initial step guess
    double zcurrtarg;
    int st;
    for (int i = 0; i < nstep; ++i) 
    {
      zcurrtarg = zcoord[i];
      while (z > zcurrtarg) 
      {
	st = gsl_odeiv_evolve_apply(e, c, s, &sys, &z, zcurrtarg, &h, y);
	if (st != 0) 
	{
	  gsl_odeiv_evolve_free(e);
	  gsl_odeiv_control_free(c);
	  gsl_odeiv_step_free(s);
	  return st;
	}
      }
      if ( isnan(abs(y[0])) || isnan(abs(y[1]))  || isnan(abs(y[2])) ) 
      {
	gsl_odeiv_evolve_free(e);
	gsl_odeiv_control_free(c);
	gsl_odeiv_step_free(s);
	std::cout<<"\nFailure with om="<<om<<" c2="<<c2<<" c3="<<c3<<" c4="<<c4<<" c5="<<c5<<" cG="<<cG<<" c0="<<c0<<" at z="<<zcurrtarg<<" h="<<y[0]<<" x="<<y[1]<<" y="<<y[2]<<std::endl;
	// Careful, can't throw exceptions when linking to fortran with ifort
	//throw CosFitterExcept("lumdist","calcValOmC2C3C4C5CGC0", "NaN value!",1);
	return 8;
      }
      hubble[i] = y[0];
      x[i] = y[1];
    }
  } else 
  {
    double y[3] = { 1.0, -1.0, 0.0 };  //Inital value of integral
    double z = 0;
    int nstep = min(hubble.size(),x.size());
    double h = 1e-6; //Initial step guess
    double zcurrtarg;
    int st;
    for (int i = 0; i < nstep; ++i) 
    {
      zcurrtarg = zcoord[i];
      while (z < zcurrtarg) 
      {
	st = gsl_odeiv_evolve_apply(e, c, s, &sys, &z, zcurrtarg, &h, y);
	if (st != 0) 
	{
	  gsl_odeiv_evolve_free(e);
	  gsl_odeiv_control_free(c);
	  gsl_odeiv_step_free(s);
	  return st;
	}
      }
      if ( isnan(abs(y[0])) || isnan(abs(y[1]))  || isnan(abs(y[2])) ) 
      {
	gsl_odeiv_evolve_free(e);
	gsl_odeiv_control_free(c);
	gsl_odeiv_step_free(s);
	std::cout<<"\nFailure with om="<<om<<" c2="<<c2<<" c3="<<c3<<" c4="<<c4<<" c5="<<c5<<" cG="<<cG<<" c0="<<c0<<" at z="<<zcurrtarg<<" h="<<y[0]<<" x="<<y[1]<<" y="<<y[2]<<std::endl;
	// Careful, can't throw exceptions when linking to fortran with ifort
	//throw CosFitterExcept("lumdist","calcValOmC2C3C4C5CGC0", "NaN value!",1);
	return 8;
      }
      hubble[i] = y[0];
      x[i] = y[1];
    }    
  }

  gsl_odeiv_evolve_free(e);
  gsl_odeiv_control_free(c);
  gsl_odeiv_step_free(s);

  return p0_status;

}





/*!
  Evaluates the int(1/h(z)) in galileon cosmology with trapezoidal method.
  
  \param[in] om \f$\Omega_m\f$
  \param[in] c2 \f$c_2\f$, parameter of galileon model
  \param[in] c3 \f$c_3\f$, parameter of galileon model
  \param[in] c4 \f$c_4\f$, parameter of galileon model
  \param[in] c5 \f$c_5\f$, parameter of galileon model
  \param[in] cG \f$c_G\f$, parameter of galileon model
  \param[in] c0 \f$c_0\f$, parameter of galileon model
  \param[out] table of z coordinates where h(z) is evaluated
  \param[out] table of h(z) 
  \returns p0_status :
                       \status = 0 : ok
		       \status = 1 : don't give a de Sitter solution (respect to equations A7 and A8 oof Gannoudji Sami)
		       \status = 2 : tensorial no ghosr condition not satisfied
		       \status = 3 : tensorial laplace condition not satisfied
		       \status = 4 : 00 Einstein equation not fulfilled
		       \status = 5 : rhoP<0
		       \status = 6 : scalar noghost condition not satified
		       \status = 7 : scalar imaginary speed of sound
*/

/*
inline int lumdist::calcIntHubbleGalileon(std::vector<double>& intval, const std::vector<double>& zcoord, double om, double c2, double c3, double c4, double c5, double cG, double c0 ) const {

  double params[10];

  params[0] = om; params[1] = c2; params[2] = c3;
  params[3] = c4; params[4] = c5; params[5] = cG;
  params[6] = c0; params[7] = 0; params[8]= 0; params[9] = 0.0;
  int p0_status = 0;
  //p0_status = calcInitialValueForPi0(om,c2,c3,c4,c5,cG);	

  // Warning : in z coordinate x(0)/pi'0=-1
  // Vector y = ( dh/dz, dx/dz, dy/dz )
  double y[3] = { 1.0 , -1.0 , 0.0 };  //Inital value of integral

  //std::cout<<"calcIntHubbleGalileon"<<std::endl;
  //std::cout<<"....Used memory = "<<ReadLinuxMem()<<" Mo"<<std::endl; 
  if ( p0_status > 1 ) return p0_status;

  const gsl_odeiv_step_type * T = gsl_odeiv_step_rkf45;
  gsl_odeiv_step * s = gsl_odeiv_step_alloc (T, 3);
  gsl_odeiv_control * c = gsl_odeiv_control_y_new (prec_lum, prec_lum_rel);
  gsl_odeiv_evolve * e = gsl_odeiv_evolve_alloc (3);

  gsl_odeiv_system sys;
  sys.function = calcValOmC2C3C4C5CGC0;
  sys.dimension = 3;
  sys.params = &params;

  double z = 0;
  int nstep = intval.size();
  double h = 1e-6; //Initial step guess
  double zcurrtarg;
  double prec_invhubble = 1/y[0];
  int st;
  intval[0] = 0;
  for (int i = 1; i < nstep; ++i) {
    zcurrtarg = zcoord[i];
    while (z < zcurrtarg) {
      st = gsl_odeiv_evolve_apply(e, c, s, &sys, &z, zcurrtarg, &h, y);
      if (st != 0) {
	gsl_odeiv_evolve_free (e);
	gsl_odeiv_control_free (c);
	gsl_odeiv_step_free (s);
	return st;
      }
    }
    if ( std::isnan(abs(y[0])) || std::isnan(abs(y[1])) || std::isnan(abs(y[2])) ) {
      gsl_odeiv_evolve_free (e);
      gsl_odeiv_control_free (c);
      gsl_odeiv_step_free (s);
      //std::cout<<"\nFailure with om="<<om<<" c2="<<c2<<" c3="<<c3<<" c4="<<c4<<" c5="<<c5<<" cG="<<cG<<" c0="<<c0<<" at z="<<zcurrtarg<<" h="<<y[0]<<" x="<<y[1]<<" y="<<y[2]<<std::endl;
      return 8;
    }
    intval[i] = intval[i-1] + 0.5*(1.0/y[0]+prec_invhubble)*(zcoord[i]-zcoord[i-1]);
    prec_invhubble = 1/y[0];
  }

  gsl_odeiv_evolve_free (e);
  gsl_odeiv_control_free (c);
  gsl_odeiv_step_free (s);

  return p0_status;

}
*/


//int main()
//{

  //double coord = 0;  

  // File to store h and x as a function of z
  //ofstream f;
  //f.open("5params_all_2015_combined_noH0prior.dat", ios::out );
  
  /*
  // Set of parameters (Best fit 2015 BAO+Planck+Lya from Jeremy)
  double h = 0.762;
  double OmegaM0 = 0.261175;
  double OmegaR0 = 2.469*pow(10,-5)*(1+0.2271*3.04)/(pow(h,2));
  double om = OmegaM0;
  double orad = OmegaR0;
  double c2 = -5.25803;
  double c3 = -1.13137;
  double cG = 0;
  double c4 = 1.0/27*(30*(om+orad-1)-9*c2+24*c3-6*cG);
  double c5 = 1.0/3*(4*(om+orad-1)-c2+2*c3-2*cG);
  double c0 = 0;
  char* name = "tracker_test.dat";
  */

  /*
  // Best fit Unc 2015 All+JLA+noH0prior from Jeremy
  double h = 0.736;
  double OmegaM0 = 0.275;
  double OmegaR0 = 2.469*pow(10,-5)*(1+0.2271*3.04)/(pow(h,2));
  double om = OmegaM0;
  double orad = OmegaR0;
  double c2 = -4.145;
  double c3 = -1.545;
  double cG = 0;
  double c4 = -0.776;
  double c5 = (om+orad-1+c2/6-2*c3+7.5*c4-3*cG)/7;
  double c0 = 0;
  */
  

  /*
  // Best fit + cG +1 sigma 2015 All+JLA+noH0prior from Jeremy
  double h = 0.727;
  double OmegaM0 = 0.280;
  double OmegaR0 = 2.469*pow(10,-5)*(1+0.2271*3.04)/(pow(h,2));
  double om = OmegaM0;
  double orad = OmegaR0;
  double c2 = -3.434;
  double c3 = -1.062;
  double cG = 0.235;
  double c4 = -0.610;
  double c5 = (om+orad-1+c2/6-2*c3+7.5*c4-3*cG)/7;
  double c0 = 0;
  */


  /*
  if(coord)
  {
    f.open(name, ios::out );

    int n_steps = 10000000;
    double amax = 1;
    double amin = 0.000908;
    //double scale_factor[n_steps];
    //scale_factor[0] = amax;
    //for(int i=1; i<sizeof(scale_factor)/sizeof(double); i++) {scale_factor[i]=scale_factor[i-1]+(amin-amax)/n_steps;}
    //const std::vector<double> acoord(scale_factor, scale_factor+sizeof(scale_factor)/sizeof(double));
    std::vector<double> acoord;
    acoord.push_back(amax);
    for(int i=1; i<n_steps; i++) {acoord.push_back(acoord[0]+i*(amin-amax)/n_steps);}
    std::vector<double> hubble(acoord.size(), 999999);
    std::vector<double> x(acoord.size(), 999999);
    std::vector<double> OmegaM(acoord.size(), 999999);
    std::vector<double> OmegaR(acoord.size(), 999999);
    std::vector<double> OmegaP(acoord.size(), 999999);

    int status = calcHubbleGalileon(hubble, x, acoord, om, orad, c2, c3, c4, c5, cG, c0, coord);

    std::cout << "OmegaM0 = " << om << std::endl;
    std::cout << "OmegaR0 = " << orad << std::endl;
    std::cout << "c0 = " << c0 << std::endl;
    std::cout << "c2 = " << c2 << std::endl;
    std::cout << "c3 = " << c3 << std::endl;
    std::cout << "c4 = " << c4 << std::endl;
    std::cout << "c5 = " << c5 << std::endl;
    std::cout << "cG = " << cG << std::endl;
    std::cout << "\n\nstatus = " << status << "\n" << std::endl;

    if(status!=0) return 0;

    for(int i=0; i<hubble.size(); i++)
    {
      OmegaM[i]=om/(pow(hubble[i], 2)*pow(acoord[i], 3));
      OmegaR[i]=orad/(pow(hubble[i], 2)*pow(acoord[i], 4));
      OmegaP[i] = (c2/2.0*pow(acoord[i]*hubble[i]*x[i], 2) - 6.0*c3*hubble[i]*pow(acoord[i]*hubble[i]*x[i], 3) + 22.5*c4*pow(hubble[i], 2)*pow(acoord[i]*hubble[i]*x[i], 4) - 21.0*c5*pow(hubble[i], 3)*pow(acoord[i]*hubble[i]*x[i], 5) - 9*cG*pow(hubble[i], 2)*pow(acoord[i]*hubble[i]*x[i], 2) )/(3.0*pow(hubble[i], 2));

      f << 1/acoord[i]-1 << " ; " << hubble[i] << " ; " << x[i] << " ; " << OmegaP[i] << " ; " << OmegaM[i] << " ; " << OmegaR[i] << std::endl;
    }

  } else
  {
    f.open(name, ios::out );

    int n_steps = 1000000;
    double zmin = 0;
    double zmax = 1100;
    std::vector<double> zcoord;
    zcoord.push_back(zmin);
    for(int i=1; i<n_steps; i++) {zcoord.push_back(zcoord[0]+i*(zmax-zmin)/n_steps);}
    std::vector<double> hubble(zcoord.size(), 999999);
    std::vector<double> x(zcoord.size(), 999999);
    std::vector<double> OmegaM(zcoord.size(), 999999);
    std::vector<double> OmegaR(zcoord.size(), 999999);
    std::vector<double> OmegaP(zcoord.size(), 999999);

    int status = calcHubbleGalileon(hubble, x, zcoord, om, orad, c2, c3, c4, c5, cG, c0, coord);

    std::cout << "OmegaM0 = " << om << std::endl;
    std::cout << "OmegaR0 = " << orad << std::endl;
    std::cout << "c0 = " << c0 << std::endl;
    std::cout << "c2 = " << c2 << std::endl;
    std::cout << "c3 = " << c3 << std::endl;
    std::cout << "c4 = " << c4 << std::endl;
    std::cout << "c5 = " << c5 << std::endl;
    std::cout << "cG = " << cG << std::endl;
    std::cout << "\n\nstatus = " << status << "\n" << std::endl;

    if(status!=0) return 0;

    for(int i=0; i<hubble.size(); i++)
    {
      OmegaM[i]=om*pow(1+zcoord[i], 3)/pow(hubble[i], 2);
      OmegaR[i]=orad*pow(1+zcoord[i], 4)/pow(hubble[i], 2);
      OmegaP[i] = (c2/2.0*pow(-(1+zcoord[i])*hubble[i]*x[i], 2) - 6.0*c3*hubble[i]*pow(-(1+zcoord[i])*hubble[i]*x[i], 3) + 22.5*c4*pow(hubble[i], 2)*pow(-(1+zcoord[i])*hubble[i]*x[i], 4) - 21.0*c5*pow(hubble[i], 3)*pow(-(1+zcoord[i])*hubble[i]*x[i], 5) - 9*cG*pow(hubble[i], 2)*pow(-(1+zcoord[i])*hubble[i]*x[i], 2) )/(3.0*pow(hubble[i], 2));

      f << std::setprecision(12) << zcoord[i] << " ; " << hubble[i] << " ; " << x[i] << " ; " << OmegaP[i] << " ; " << OmegaM[i] << " ; " << OmegaR[i] << std::endl;
    }

  }
  */


  //return 0;

//}


int main()
{
  char* name = "barreira_background_sample/barreira_galileon_1_1e-4_test.dat";
  ofstream f;

  // Galileon parameters from Barreira
  // Galileon 1
  double dotphi = 3.25e-14;
  double c3 = -12.8;
  double c4 = -1.7;
  double c5 = -1.0;
  double cG = 0;
  double c2 = -27.00;
  double age = 13.978;
  double ratio_rho = 1e-4;
  //double c2 = -27.49;
  //double age = 14.317;
  //double ratio_rho = 1e-6;    
  //double c2 = -27.56;
  //double age = 14.366;
  //double ratio_rho = 1e-7;
  //double c2 = -27.58;
  //double age = 14.374;
  //double ratio_rho = 1e-8;
  //double c2 = -27.59;
  //double age = 14.375;
  double c0 = 0;
  double orad = 8e-5;
  double om = 0.265;

  double n_steps = 1000000;
  double apremin = 1e-7;
  double amin = 9.99999e-7;
  double amax = 1;
  std::vector<double> acoord;
  acoord.push_back(apremin);
  acoord.push_back(amin);
  for(int i=1; i<n_steps+1; i++) {acoord.push_back(acoord[1]+i*(amax-amin)/n_steps);}

  double xi = 0;
  double H[] = {sqrt(om/pow(apremin, 3)*(1+ratio_rho)+orad/pow(apremin, 4)), sqrt(om/pow(amin, 3)*(1+ratio_rho)+orad/pow(amin, 4))};
  //double H[] = {sqrt(om/pow(apremin, 3)+orad/pow(apremin, 4)), sqrt(om/pow(amin, 3)+orad/pow(amin, 4))};
  //double ratio_rho = pow(acoord[1], 3)/(3*om)*(0.5*c2*pow(dotphi*H[1], 2)-6*c3*H[1]*pow(dotphi*H[1], 3)+22.5*c4*pow(H[1], 2)*pow(dotphi*H[1], 4)-21*c5*pow(H[1], 3)*pow(dotphi*H[1], 5)-9*cG*pow(H[1], 2)*pow(dotphi*H[1], 2));

  std::vector<double> hubble(acoord.size()-1, 999999);
  std::vector<double> x(acoord.size()-1, 999999);

  int status = initCond(hubble, x, xi, acoord, H, ratio_rho, age, om, orad, c2, c3, c4, c5, cG, c0);

  std::cout << std::setprecision(16) << "Hi = " << H[1] << "xi = " << xi << endl;

  std::cout << "Status : " << status << endl;

  if(status != 0) return 0;

  f.open(name, ios::out );

  std::cout << "rho_phi/rho_m = " << ratio_rho << std::endl;
  std::cout << "c0 = " << c0 << std::endl;
  std::cout << "c2 = " << c2 << std::endl;
  std::cout << "c3 = " << c3 << std::endl;
  std::cout << "c4 = " << c4 << std::endl;
  std::cout << "c5 = " << c5 << std::endl;
  std::cout << "cG = " << cG << std::endl;
  std::cout << "\n\nstatus = " << status << "\n" << std::endl;

  for(int i=0; i<hubble.size(); i++)
  {
    double alpha = c2/6.0*hubble[i]*acoord[i+1]*x[i]-3*c3*pow(hubble[i], 3)*pow(acoord[i+1]*x[i], 2) + 15*c4*pow(hubble[i], 5)*pow(acoord[i+1]*x[i], 3) - 17.5*c5*pow(hubble[i], 7)*pow(acoord[i+1]*x[i], 4) - 3*cG*pow(hubble[i], 3)*acoord[i+1]*x[i];
    double gamma = c2/3.0*pow(hubble[i], 2)*acoord[i+1]*x[i]-c3*pow(hubble[i], 4)*pow(acoord[i+1]*x[i], 2) + 2.5*c5*pow(hubble[i], 8)*pow(acoord[i+1]*x[i], 4) - 2*cG*pow(hubble[i], 4)*acoord[i+1]*x[i];
    double beta = c2/6.0*pow(hubble[i], 2) -2*c3*pow(hubble[i], 4)*acoord[i+1]*x[i] + 9*c4*pow(hubble[i], 6)*pow(acoord[i+1]*x[i], 2) - 10*c5*pow(hubble[i], 8)*pow(acoord[i+1]*x[i], 3) - cG*pow(hubble[i], 4);
    double sigma = 2*hubble[i] + 2*c3*pow(hubble[i], 3)*pow(acoord[i+1]*x[i], 3) - 15*c4*pow(hubble[i], 5)*pow(acoord[i+1]*x[i], 4) + 21*c5*pow(hubble[i], 7)*pow(acoord[i+1]*x[i], 5) + 6*cG*pow(hubble[i], 3)*pow(acoord[i+1]*x[i], 2);
    double lambda = 3*pow(hubble[i], 2) + orad/pow(acoord[i+1], 4) + 0.5*c2*pow(hubble[i], 2)*pow(acoord[i+1]*x[i], 2) - 2*c3*pow(hubble[i], 4)*pow(acoord[i+1]*x[i], 3) + 7.5*c4*pow(hubble[i], 6)*pow(acoord[i+1]*x[i], 4) - 9*c5*pow(hubble[i], 8)*pow(acoord[i+1]*x[i], 5) - cG*pow(hubble[i], 4)*pow(acoord[i+1]*x[i], 2);
    double omega = 2*c3*pow(hubble[i], 4)*pow(acoord[i+1]*x[i], 2) - 12*c4*pow(hubble[i], 6)*pow(acoord[i+1]*x[i], 3) + 15*c5*pow(hubble[i], 8)*pow(acoord[i+1]*x[i], 4) + 4*cG*pow(hubble[i], 4)*acoord[i+1]*x[i];

    double x_prime = -acoord[i+1]*x[i]+(alpha*lambda-sigma*gamma)/(sigma*beta-alpha*omega);
    double h_prime = (omega*gamma-lambda*beta)/(sigma*beta-alpha*omega);

    double rho = 0.5*c2*pow(hubble[i], 2)*pow(acoord[i+1]*x[i], 2) - 6*c3*pow(hubble[i], 4)*pow(acoord[i+1]*x[i], 3) + 22.5*c4*pow(hubble[i], 6)*pow(acoord[i+1]*x[i], 4) - 21*c5*pow(hubble[i], 8)*pow(acoord[i+1]*x[i], 5) - 9*cG*pow(hubble[i], 4)*pow(acoord[i+1]*x[i], 2);
    double p = 0.5*c2*pow(hubble[i], 2)*pow(acoord[i+1]*x[i], 2) + 2*c3*pow(hubble[i], 3)*pow(acoord[i+1]*x[i], 2)*(h_prime*acoord[i+1]*x[i]+x_prime*hubble[i]) - c4*(4.5*pow(hubble[i], 6)*pow(acoord[i+1]*x[i], 4) + 12*pow(hubble[i], 6)*pow(acoord[i+1]*x[i], 3)*x_prime + 15*pow(hubble[i], 5)*pow(acoord[i+1]*x[i], 4)*h_prime) + 3*c5*pow(hubble[i], 7)*pow(acoord[i+1]*x[i], 4)*(5*hubble[i]*x_prime+7*h_prime*acoord[i+1]*x[i]+2*hubble[i]*acoord[i+1]*x[i]) + cG*(6*pow(hubble[i], 3)*pow(acoord[i+1]*x[i], 2)*h_prime + 4*pow(hubble[i], 4)*acoord[i+1]*x[i]*x_prime + 3*pow(hubble[i], 4)*pow(acoord[i+1]*x[i], 2));

    double w = p/rho;

    double hubble_LCDM = sqrt(om/pow(acoord[i+1], 3)+orad/pow(acoord[i+1], 4)+(1-om-orad));
    f << std::setprecision(12) << acoord[i+1] << " ; " << hubble[i] << " ; " << hubble_LCDM << " ; " << w << " ; " << p << " ; " << rho << " ; " << x[i] << std::endl;
  }

  return 0;


}


