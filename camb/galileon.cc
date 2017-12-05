// Modified by Clement Leloup 04-10-2017 :
// Corrected all functions calculating perturbations
// Added global variables to store current value of h,x,dh/dlna,dx/dlna
// Added calcPertOmC2C3C4C5CGC0 function to check if theoretical constraints are verified
// Removed the JNEVEUz way of solving background
// Added massive neutrinos
// Modified by Clement Leloup 09-12-2016 :
// Added the interpolation global variables and the interpolation in arrays_
// Modified by Clement Leloup 06-12-2016 :
// Added functions grhogal and gpresgal
// Changed name of file to "galileon.cc"
// Added all functions calculating perturbations part
// Modified by Clement Leloup 04-12-2016 :
// Created header and moved everything needed there
// Modified by Clement Leloup 05-10-2016 :
// Changed name of dtauda to arrays
// Added the function handxofa
// Added three global vectors intvar, hubble and x and modified everything to take it into account
// Moved common parameters to global
// Modified by Clement Leloup 29-09-2016 :
// Made the code a little more C-like (removing strings in favor of char*)
// Modified main and dtauda
// Changed type of readParamsFile type to array of char*
// Added the dotphi parameter
// Cleaned a little
// Modified by Clement Leloup 27-09-2016 :
// Changed type of dtauda to double array (here double pointer)
// Modified by Clement Leloup 22-07-2016 :
// Added the function dtauda to be linked with CAMB
// Cleaned a little
// Added the function readParamsFile
// Modified main
// Modified by Clement Leloup 23-05-2016 :
// Added the analytic tracker solution in main function
// Added the function calcHubbleTracker
// Modified by Clement Leloup 13-04-2016 :
// Added the "coord" parameter to calcHubbleGalileon, to take into account whether we use a or z to do the calculation. In case one use a to calculate, the calculation is done backward from a=1.
// Modified the main function. Now, the behaviour is adapted whether the user prefers a or z to calculate the observables.
// Modified by Clement Leloup 07-04-2016 :
// Added the contribution of orad to the evolution of Lander system when !useacoord
// Wrote of the main function to test and compare with former results
// Commented of every function except CalcValOmC2C3C4C5CGC0 and calcHubbleGalileon

// From modified CosFitter for Galileon Cosmology file lumdist.cc by Jeremy Neveu

#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>

#include <vector>

#include <iostream>
//#include <iomanip>
//#include "fstream"
//#include <algorithm>
//#include <math.h>
//#include "string.h"
//#include <stdio.h>
//#include <stdlib.h>

// using namespace std;

// Global vectors of the integration variable (a or z), h and x
std::vector<double> intvar;
std::vector<double> hubble;
std::vector<double> x;

// Background at current a point
double acurr = 0;
double hcurr = 0;
double xcurr = 0;
double dx = 0;
double dh = 0;

// Interpolation tools
gsl_interp_accel* acc;
gsl_spline* spline_h;
gsl_spline* spline_x;

// Global galileon parameters
double om = 0;
double orad = 0;
double h0 = 0;
double c2 = 0;
double c3 = 0;
double c4 = 0;
double c5 = 0;
double cG = 0;
double c0 = 0; // not implemented with c0 yet, but in case it is one day ...

// External functions from fortran
extern"C" void massivenu_mp_nu_rho_(double* am, double* rhonu);
extern"C" void massivenu_mp_nu_background_(double* am, double* rhonu, double* pnu);
extern"C" void massivenu_mp_nurhopres_(double* am, double* rhonu, double* pnu);

// Return absolute max of vector
double maxVec(std::vector<double> vec){

  double max = 1e-16;

  if(vec.size() > 0){
    for(int i = 0; i < vec.size(); i++){
      if(fabs(vec.at(i)) > max) max = fabs(vec.at(i));
    }
  }

  return max;

}

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
inline int calcPertOmC2C3C4C5CGC0(double a, const double y[3]) {

  fflush(stdout);

  double k1,k2,k3,k4,k5,k6;
  double alpha,gamma,beta,sigma,lambda,omega,denom;
  double dhdlna,dxdlna;

  // system from lna in z : 
  // y[0](lna) -> y[0](z)
  // y[1](lna) -> -(1+z)*y[1](z)
  // f[0](lna) -> -(1+z)*f[0](z)
  // f[1](lna) -> (1+z)^2*f[1](z)+(1+z)*y[1](z)

  double xgal = a*y[1];
  double prod = xgal*y[0];
  double prod2 = prod*prod;
  double prod3 = prod*prod2;
  double prod4 = prod*prod3;
  double prod5 = prod*prod4;
  double h = y[0];
  double h2 = h*h;
  double h3 = h*h2;
  double h4 = h2*h2;

  // The equations : 
  alpha = c0*h + c2/6.0*prod - 3*c3*h*prod2 + 15*c4*h2*prod3 - 17.5*c5*h3*prod4 - 3.0*cG*h2*prod;
  gamma = 2*c0*h2 + c2/3.0*h*prod - c3*h2*prod2 + 2.5*c5*h4*prod4 - 2.0*cG*h3*prod;
  beta = -2*c3*h3*prod + c2/6.0*h2 + 9*c4*h4*prod2 - 10*c5*h4*h*prod3 - cG*h4;
  sigma = 2.0*( 1.0 - 2*c0*y[2] )*h - 2.0*c0*prod + 2.0*c3*prod3 - 15.0*c4*h*prod4 + 21.0*c5*h2*prod5 + 6.0*cG*h*prod2;
  lambda = 3.0*( 1 - 2*c0*y[2] )*h2 + orad/(a*a*a*a) - 2.0*c0*h*prod - 2.0*c3*h*prod3 + c2/2.0*prod2 + 7.5*c4*h2*prod4 - 9.0*c5*h3*prod5 - cG*h2*prod2;
  omega = -2*c0*h2 + 2*c3*h2*prod2 - 12*c4*h3*prod3 + 15*c5*h4*prod4 + 4.0*cG*h3*prod;
  denom = sigma*beta - alpha*omega;

  dhdlna = (omega*gamma - lambda*beta)/denom;
  dxdlna = (alpha*lambda - sigma*gamma)/denom - xgal;

  k1 = -6*c4*h*prod2*(dhdlna*xgal + h*dxdlna + prod/3.0) -2*c0 + c5*h2*prod3*(12*h*dxdlna + 15*dhdlna*xgal + 3*prod) + 2.0*cG*(dhdlna*prod+h2*dxdlna+h*prod);
  k2 = -0.5*c2 + 6*c3*h*prod - 27*c4*h2*prod2 + 30*c5*h3*prod3 + 3.0*cG*h2;
  k3 = -(1.0 - 2.0*c0*y[2]) - 0.5*c4*prod4 - 3*c5*h*prod4*(h*dxdlna + dhdlna*xgal) + cG*prod2;
  k4 = -2.0*(1.0 - 2.0*c0*y[2]) + 3.0*c4*prod4 - 6*c5*h*prod5 - 2.0*cG*prod2;
  k5 = 2*c3*prod2 - 12*c4*h*prod3 -2.0*c0 + 15.0*c5*h2*prod4 + 4.0*cG*h*prod;
  k6 = 0.5*c2 - 2.0*c3*( h2*dxdlna + dhdlna*prod + 2*h*prod ) + c4*( 12*h3*prod*dxdlna + 18*h*prod2*dhdlna + 13*h2*prod2 ) - c5*( 18*h2*h2*prod2*dxdlna + 30*h2*prod3*dhdlna + 12*h3*prod3 ) - cG*( 2.0*h*dhdlna + 3.0*h2 );

  double noghost = k2 + 1.5*k5*k5/k4;
  double cs2 = (4*k1*k4*k5 - 2*k3*k5*k5 - 2*k4*k4*k6)/(k4*(2*k4*k2 + 3*k5*k5));
  // Here we can't have a condition noghost>0 because with some set of
  // parameters we have noghost=1e-16, condition passes anyway
  // and then rk4 fails and gives NaN values.

  // From Felice & Tsujikawa
  double noghostt = 0.5 - 0.75*c4*prod4 + 1.5*c5*prod5*h + 0.5*cG*prod2 - c0*y[2];
  double ct2 = (0.5 + 0.25*c4*prod4 + 1.5*c5*prod4*h*(h*dxdlna+dhdlna*xgal) - 0.5*cG*prod2 -c0*y[2]) / (0.5 - 0.75*c4*prod4 + 1.5*c5*prod5*h + 0.5*cG*prod2 - c0*y[2]);

  // // To test theoretical constraints
  // printf("a = %.12f \t Qs = %.12f \t cs2 = %.12f \t Qt = %.12f \t cT2 = %f \t OmegaPi = %.12f\n om = %.12f \t orad = %.12f \t c2 = %.12f \t c3 = %.12f \t c4 = %.12f \t cG = %.12f\n", a, noghost, cs2, noghostt, ct2, OmegaP, om, orad, c2, c3, c4, cG);

  if(noghost>-1e-8){
    fprintf(stderr, "Error : scalar no ghost constraint not verified : a = %.12f \t Qs = %.12f\n", a, noghost);
    return 1;
  }
  if(cs2<0){
    fprintf(stderr, "Error : complex sound speed : a = %.12f \t cs2 = %.12f\n", a, cs2);
    return 2;
  }
  if(noghostt < 0){
    fprintf(stderr, "Error : tensorial no ghost constraint not verified : a = %.12f \t Qt = %.12f\n", a, noghostt);
    return 3;
  }
  if(ct2 < 0){
    fprintf(stderr, "Error : Laplace stability condition not verified : a = %f \t cT2 = %f\n", a, ct2);
    return 4;
  }

  return 0;

}
/*
  \param[in] z Redshift
  \param[in] y Current value of integral : y[0] is h(z), y[1] is dpi/dz and y[2] is pi
  \param[in] params Array of parameters (\f$\Omega_m, c_2, c_3, c_4, c_5, c_G, c_0, orad, useacoord\f$)
  \param[out] f The value of the integral term : f[0] is dh/dz, f[1] is d^2pi/dz^2 and f[2] is dpi/dz
*/
// int calcValOmC2C3C4C5CGC0(double a, const double y[3], double f[3], void* params){
int calcValOmC2C3C4C5CGC0(double a, const double y[3], double f[3], void* params){

  double alpha,gamma,beta,sigma,lambda,omega,OmegaP,OmegaM,OmTest;
  const double *in_params = static_cast<const double*>(params);
  // bool store_derivs = (int)in_params[0];
  double grhormass[5];
  for(int i = 0; i<5; i++){ grhormass[i] = in_params[i]; }
  double nu_masses[5];
  for(int i = 0; i<5; i++){ nu_masses[i] = in_params[5+i]; }
  int nu_mass_eigenstates = in_params[10];

  double prod = a*y[0]*y[1]; // prod=h*x
  double prod2 = prod*prod;
  double prod3 = prod*prod2;
  double prod4 = prod*prod3;
  double prod5 = prod*prod4;
  double h = y[0];
  double h2 = h*h;
  double h3 = h*h2;
  double h4 = h2*h2;
  double a2 = a*a;

  // The equations : 
  alpha = c2/6.0*prod - 3*c3*h*prod2 + 15*c4*h2*prod3 - 17.5*c5*h3*prod4 - 3.0*cG*h2*prod;
  gamma = c2/3.0*h*prod - c3*h2*prod2 + 2.5*c5*h4*prod4 - 2.0*cG*h3*prod;
  beta = c2/6.0*h2 - 2*c3*h3*prod + 9*c4*h4*prod2 - 10*c5*h4*h*prod3 - cG*h4;
  sigma = 2.0*h + 2.0*c3*prod3 - 15.0*c4*h*prod4 + 21.0*c5*h2*prod5 + 6.0*cG*h*prod2;
  lambda = 3.0*h2 + orad/(a2*a2) + c2/2.0*prod2 - 2.0*c3*h*prod3 + 7.5*c4*h2*prod4 - 9.0*c5*h3*prod5 - cG*h2*prod2;

  // Contribution from massive neutrinos
  if(nu_mass_eigenstates>0){
    for(int i=0; i<nu_mass_eigenstates; i++){
      double rhonu = 0;
      double pnu = 0;
      double am = a*nu_masses[i];
      massivenu_mp_nu_background_(&am, &rhonu, &pnu);
      lambda += grhormass[i]*pnu/(a2*a2*h0*h0);
      // printf("i : %d \t a : %.16f \t am : %.16f \t rhonu : %.16f \t lambda_numass : %.16f\n", i, a, am, rhonu, grhormass[i]*pnu/(a2*a2*h0*h0));
    }
  }

  omega = 2*c3*h2*prod2 - 12*c4*h3*prod3 + 15*c5*h4*prod4 + 4.0*cG*h3*prod;

  double denom = sigma*beta - alpha*omega;
  f[0] = (omega*gamma - lambda*beta) / (a*denom);
  f[1] = (alpha*lambda - sigma*gamma) / (a2*denom) - 2.0*y[1]/a;
  f[2] = y[1];

  return GSL_SUCCESS;
}


/*!
  Evaluates the h(z) in galileon cosmology

  \param[in] om \f$\Omega_m\f$
  \param[in] orad \f$\Omega_{rad}\f$
  \param[in] c2 \f$c_2\f$, parameter of galileon model
  \param[in] c3 \f$c_3\f$, parameter of galileon model
  \param[in] c4 \f$c_4\f$, parameter of galileon model
  \param[in] c5 \f$c_5\f$, parameter of galileon model
  \param[in] cG \f$c_G\f$, parameter of galileon model
  \param[in] c0 \f$c_0\f$, parameter of galileon model
  \param[in] coord, whether to use a or z coordinate
  \param[out] table of a or z coordinates where h is evaluated
  \param[out] table of h
  \param[out] table of x
  \returns p0_status :
                       \status = 0 : ok
		       \status = 1 : don't give a de Sitter solution (respect to equations A7 and A8 oof Gannoudji Sami)
		       \status = 2 : tensorial no ghost condition not satisfied
		       \status = 3 : tensorial laplace condition not satisfied
		       \status = 4 : 00 Einstein equation not fulfilled
		       \status = 5 : rhoP<0
		       \status = 6 : noghost condition not satisfied
		       \status = 7 : imaginary speed of sound
		       \status = 8 : failed to integrate
		       \status = 9 : one of the global vectors is empty
*/
int calcHubbleGalileon(double* grhormass, double* nu_masses, int* nu_mass_eigenstates){

  fflush(stdout);

  if(intvar.size() == 0 || hubble.size() == 0 || x.size() == 0){
    printf("One of the global arrays is empty\n");
    fflush(stdout);
    return 9;
  }

  int testPert;
  double h2, x1, x2, x3, a3;
  double OmegaP, OmegaM, OmTest;
  double params[11]; // parameters for massive neutrinos
  for(int i = 0; i<5; i++){ params[i] = grhormass[i]; }
  for(int i = 0; i<5; i++){ params[5+i] = nu_masses[i]; }
  params[10] = (*nu_mass_eigenstates);

  const gsl_odeiv_step_type * T = gsl_odeiv_step_rkf45;
  gsl_odeiv_step * s = gsl_odeiv_step_alloc(T, 3);
  gsl_odeiv_control * c = gsl_odeiv_control_y_new(1e-18, 1e-16);
  gsl_odeiv_evolve * e = gsl_odeiv_evolve_alloc(3);
  gsl_odeiv_system sys;
  sys.function = calcValOmC2C3C4C5CGC0;
  sys.dimension = 3;
  sys.params = &params;

  double y[3] = { 1.0, 1.0, 0.0 };  //Inital value of integral
  double a = 1;
  // int nstep = std::min(hubble.size(),x.size());
  int nstep = intvar.size();
  double h = -1e-6; //Initial step guess
  double acurrtarg;
  int st;
  for(int i = 0; i < nstep; ++i){
    acurrtarg = intvar[i];
    while(a > acurrtarg){
      st = gsl_odeiv_evolve_apply(e, c, s, &sys, &a, acurrtarg, &h, y);
      }
   
    if(isnan(fabs(y[0])) || isnan(fabs(y[1]))  || isnan(fabs(y[2]))){
      gsl_odeiv_evolve_free(e);
      gsl_odeiv_control_free(c);
      gsl_odeiv_step_free(s);
      fprintf(stderr, "\nFailure with om = %f, c2 = %f, c3 = %f, c4 = %f, c5 = %f, cG = %f and c0 = %f at a = %f, h = %f, x = %f and y = %f\n", om, c2, c3, c4, c5, cG, c0, acurrtarg, y[0], y[1], y[2]);
      return 8;
    }

    testPert = calcPertOmC2C3C4C5CGC0(intvar[i], y);
    if(testPert != 0){
      return 7;
    }
    
    // Variables to optimize computation
    h2 = y[0]*y[0];
    // x1 = intvar[i]*y[1];
    x1 = a*y[1];
    x2 = x1*x1;
    x3 = x2*x1;
    // a3 = intvar[i]*intvar[i]*intvar[i];
    a3 = a*a*a;

    // Friedmann equation for background
    OmegaP = (0.5*c2*x2 - 6.0*c3*h2*x3 + 22.5*c4*h2*h2*x2*x2 - 21.0*c5*h2*h2*h2*x3*x2 - 9.0*cG*h2*x2)/3.0;
    // OmegaM = 1 - OmegaP - orad/(a3*intvar[i]*h2);
    OmegaM = 1 - OmegaP - orad/(a3*a*h2);

    // Contribution from massive neutrinos
    if((*nu_mass_eigenstates)>0){
      for(int j=0; j<(*nu_mass_eigenstates); j++){
    	double rhonu = 0;
	double pnu = 0;
    	double am = intvar[i]*nu_masses[j];
    	massivenu_mp_nu_background_(&am, &rhonu, &pnu);
    	OmegaM -= grhormass[j]*rhonu/(a3*a*3.*h2*pow(h0, 2));
      }
    }

    OmTest = om/(a3*h2);

    // //For test
    // printf("Test on OmegaM : %f %f %f %f %f %f %f %f %f %f %f %f %f %f %.12f\n", om, orad, c2, c3, c4, c5, cG, c0, OmTest, OmegaM, orad/(a3*a*h2), OmegaP, a3*h2, intvar[i], fabs((OmegaM - OmTest)/OmTest));
    // // printf("lambda : %.12f\tusual : %.12f\tfrom massive neutrinos : %.12f\n", lambda, 3.0*h2 + orad/(a3*intvar[i]) - 2.0*c0*h2*x1 - 2.0*c3*h2*h2*x3 + c2/2.0*h2*x2 + 7.5*c4*h2*h2*h2*x2*x2 - 9.0*c5*h2*h2*h2*h2*x3*x2 - cG*h2*h2*x2, grhormass[0]*nuPres(intvar[i]*nu_masses[0])/(a3*intvar[i]*h0*h0));

    if(OmegaP<0){
    //   // fprintf(stderr, "Negative galileon energy density : a = %.12f \t %.12f\n", intvar[i], OmegaP);
      fprintf(stderr, "Negative galileon energy density : a = %.12f \t %.12f\n", a, OmegaP);
      return 5;
    }

    if ( fabs((OmegaM - OmTest)/OmTest)>1e-4 ) {
      // fprintf(stderr, "Integration error : %f %f %f %f %f %f %f %f %f %f %f %f %.16f %.16f %f\n", om, orad, c2, c3, c4, c5, cG, c0, OmTest, OmegaM, OmegaP, OmegaM - 1 + OmegaP + orad/(a3*intvar[i]*h2), intvar[i], a, fabs((OmegaM - OmTest)/OmTest));
      fprintf(stderr, "Integration error : %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", om, orad, c2, c3, c4, c5, cG, c0, OmTest, OmegaM, OmegaP, OmegaM - 1 + OmegaP + orad/(a3*intvar[i]*h2), a, fabs((OmegaM - OmTest)/OmTest));
      return 4;
    }

    hubble[i] = y[0];
    x[i] = y[1];
  }

  gsl_odeiv_evolve_free(e);
  gsl_odeiv_control_free(c);
  gsl_odeiv_step_free(s);

  return 0;

}


/*!
  Calculate analytically the h(z) and x(z) tracker solution
  
  \param[in] om \f$\Omega_m\f$
  \param[in] orad \f$\Omega_{rad}\f$
  \param[in] c2 \f$c_2\f$, parameter of galileon model
  \param[in] c3 \f$c_3\f$, parameter of galileon model
  \param[in] c4 \f$c_4\f$, parameter of galileon model
  \param[in] c5 \f$c_5\f$, parameter of galileon model
  \param[in] cG \f$c_G\f$, parameter of galileon model
  \param[in] c0 \f$c_0\f$, parameter of galileon model
  \param[in] coord, whether to use a or z coordinate
  \param[out] table of a or z coordinates where h is evaluated
  \param[out] table of h 
  \param[out] table of x
*/
void calcHubbleTracker(){

  fflush(stdout);
 
  if(intvar.size() == 0 || hubble.size() == 0 || x.size() == 0){
    printf("One of the global arrays is empty\n");
    return ;
  }

  double op = (c2/2.0 - 6.0*c3 + 22.5*c4 - 21.0*c5 - 9*cG )/3.0; //Omega_phi

  int nstep = std::min(hubble.size(),x.size());
  for(int i = 0; i < nstep; ++i){
   double a2 = intvar[i]*intvar[i];
    hubble[i] = sqrt(0.5*(om/(a2*intvar[i])+orad/(a2*a2)-3*cG)+sqrt(op/9+3*cG+0.25*pow(3*cG-om/(a2*intvar[i])-orad/(a2*a2),2))); // Analytical solution for H
    x[i] = intvar[i]/(hubble[i]*hubble[i]); // since for tracker, h^2*x = 1
  }

}


// Function that calculates dtau/da = 1/(a^2*H) (tau being confromal time)
// The function is the one linked to CAMB in order to calculate the background contribution
extern "C" int arrays_(double* omegar, double* omegam, double* H0in, double* c2in, double* c3in, double* c4in, double* cGin, double* grhormass, double* nu_masses, int* nu_mass_eigenstates){

  fflush(stdout);

  // Clear three global vectors
  intvar.clear();
  hubble.clear();
  x.clear();

  om = (*omegam);
  orad = (*omegar);
  h0 = (*H0in);
  c2 = (*c2in);
  c3 = (*c3in);
  c4 = (*c4in);
  cG = (*cGin);
  c5 = (-1. + om + orad + c2/6. - 2*c3 + 7.5*c4 - 3*cG)/7.;

  // Add massive neutrinos
  double rhonu = 0;
  double grhom = 3*pow(h0, 2); // critical density
  if((*nu_mass_eigenstates)>0){
    for(int i = 0; i<(*nu_mass_eigenstates); i++){
      massivenu_mp_nu_rho_(&(nu_masses[i]), &rhonu);
      c5 += rhonu*grhormass[i]/(7.*grhom);
      // printf("Omeganu of mass eigenstate %d = %.16f\n", i, rhonu*grhormass[i]*0.726*0.726*94.07/grhom);
    }
  }

  // The status of the integration, 0 if everything ok
  int status = 0;

  // Where to put initial condition
  // double amin = 9.99999e-7;
  double amin = 1e-10;
  double amax = 1.;

  // Fill the vector of a with a geometric sequence
  double nb = 30000; // number of a points using q0
  // double nb = 300000; // number of a points using q0
  double q = pow(amin/amax, 1./nb);
  for(int i = 0; i<=nb; i++){
    intvar.push_back(amax*pow(q, i));
  }

  printf("Number of points : %i\n", intvar.size());
  fflush(stdout);

  printf("OmegaM0 = %.16f\nOmegaR0 = %.16f\nc0 = %.16f\nc2 = %.16f\nc3 = %.16f\nc4 = %.16f\nc5 = %.16f\ncG = %.16f\nh0 = %.16f km/s/Mpc\n", om, orad, c0, c2, c3, c4, c5, cG, h0*2.99792458e8/1000);
  fflush(stdout);

  hubble.resize(intvar.size(), 999999);
  x.resize(intvar.size(), 999999);

  // printf("Tracker criterion : %.12f\n", fabs(c2-6*c3+18*c4-15*c5-6*cG));

  // Integrate and fill hubble and x both when tracker and not tracker
  if(fabs(c2-6*c3+18*c4-15*c5-6*cG)>1e-8)
    {
      printf("nu_mass_eigenstates : %d\n", (*nu_mass_eigenstates));
      status = calcHubbleGalileon(grhormass, nu_masses, nu_mass_eigenstates);
    }
  else {
    calcHubbleTracker();
  }

  printf("status = %i\n", status);
  fflush(stdout);

  // if(status != 0){
  //   hubble.clear();
  //   x.clear();
  // }

  // FILE* f = fopen(outfile, "w");
  // for(int i = 0; i<intvar.size()-1; i++){
  //   double alpha = c2/6*hubble[i]*intvar[i+1]*x[i]-3*c3*pow(hubble[i], 3)*pow(intvar[i+1]*x[i], 2) + 15*c4*pow(hubble[i], 5)*pow(intvar[i+1]*x[i], 3) - 17.5*c5*pow(hubble[i], 7)*pow(intvar[i+1]*x[i], 4) - 3*cG*pow(hubble[i], 3)*intvar[i+1]*x[i];
  //   double gamma = c2/3*pow(hubble[i], 2)*intvar[i+1]*x[i]-c3*pow(hubble[i], 4)*pow(intvar[i+1]*x[i], 2) + 2.5*c5*pow(hubble[i], 8)*pow(intvar[i+1]*x[i], 4) - 2*cG*pow(hubble[i], 4)*intvar[i+1]*x[i];
  //   double beta = c2/6*pow(hubble[i], 2) -2*c3*pow(hubble[i], 4)*intvar[i+1]*x[i] + 9*c4*pow(hubble[i], 6)*pow(intvar[i+1]*x[i], 2) - 10*c5*pow(hubble[i], 8)*pow(intvar[i+1]*x[i], 3) - cG*pow(hubble[i], 4);
  //   double sigma = 2*hubble[i] + 2*c3*pow(hubble[i], 3)*pow(intvar[i+1]*x[i], 3) - 15*c4*pow(hubble[i], 5)*pow(intvar[i+1]*x[i], 4) + 21*c5*pow(hubble[i], 7)*pow(intvar[i+1]*x[i], 5) + 6*cG*pow(hubble[i], 3)*pow(intvar[i+1]*x[i], 2);
  //   double lambda = 3*pow(hubble[i], 2) + orad/pow(intvar[i+1], 4) + c2/2*pow(hubble[i], 2)*pow(intvar[i+1]*x[i], 2) - 2*c3*pow(hubble[i], 4)*pow(intvar[i+1]*x[i], 3) + 7.5*c4*pow(hubble[i], 6)*pow(intvar[i+1]*x[i], 4) - 9*c5*pow(hubble[i], 8)*pow(intvar[i+1]*x[i], 5) - cG*pow(hubble[i], 4)*pow(intvar[i+1]*x[i], 2);
  //   double omega = 2*c3*pow(hubble[i], 4)*pow(intvar[i+1]*x[i], 2) - 12*c4*pow(hubble[i], 6)*pow(intvar[i+1]*x[i], 3) + 15*c5*pow(hubble[i], 8)*pow(intvar[i+1]*x[i], 4) + 4*cG*pow(hubble[i], 4)*intvar[i+1]*x[i];
	
  //   double x_prime = -intvar[i+1]*x[i]+(alpha*lambda-sigma*gamma)/(sigma*beta-alpha*omega);
  //   double h_prime = (omega*gamma-lambda*beta)/(sigma*beta-alpha*omega);
	
  //   double rho = c2/2*pow(hubble[i], 2)*pow(intvar[i+1]*x[i], 2) - 6*c3*pow(hubble[i], 4)*pow(intvar[i+1]*x[i], 3) + 22.5*c4*pow(hubble[i], 6)*pow(intvar[i+1]*x[i], 4) - 21*c5*pow(hubble[i], 8)*pow(intvar[i+1]*x[i], 5) - 9*cG*pow(hubble[i], 4)*pow(intvar[i+1]*x[i], 2);
  //   double p = c2/2*pow(hubble[i], 2)*pow(intvar[i+1]*x[i], 2) + 2*c3*pow(hubble[i], 3)*pow(intvar[i+1]*x[i], 2)*(h_prime*intvar[i+1]*x[i]+x_prime*hubble[i]) - c4*(4.5*pow(hubble[i], 6)*pow(intvar[i+1]*x[i], 4) + 12*pow(hubble[i], 6)*pow(intvar[i+1]*x[i], 3)*x_prime + 15*pow(hubble[i], 5)*pow(intvar[i+1]*x[i], 4)*h_prime) + 3*c5*pow(hubble[i], 7)*pow(intvar[i+1]*x[i], 4)*(5*hubble[i]*x_prime+7*h_prime*intvar[i+1]*x[i]+2*hubble[i]*intvar[i+1]*x[i]) + cG*(6*pow(hubble[i], 3)*pow(intvar[i+1]*x[i], 2)*h_prime + 4*pow(hubble[i], 4)*intvar[i+1]*x[i]*x_prime + 3*pow(hubble[i], 4)*pow(intvar[i+1]*x[i], 2));
	
  //   double weff = p/rho;
  //   double hubble_LCDM = sqrt(om/pow(intvar[i+1], 3)+orad/pow(intvar[i+1], 4)+(1-om-orad));
  //   double ratio = hubble[i]/hubble_LCDM;
  //   fprintf(f, "%.16f ; %.16f ; %.16f ; %.16f ; %.16f ; %.16f\n", intvar[i+1], hubble[i], x[i], hubble_LCDM, ratio, weff);
  // }

  double hubble_interp[hubble.size()];
  double x_interp[x.size()];
  double intvar_interp[hubble.size()];

  std::copy(hubble.rbegin(), hubble.rend(), hubble_interp);
  std::copy(x.rbegin(), x.rend(), x_interp);
  std::copy(intvar.rbegin(), intvar.rend(), intvar_interp);
  
  spline_h = gsl_spline_alloc(gsl_interp_cspline, hubble.size());
  spline_x = gsl_spline_alloc(gsl_interp_cspline, x.size());  
  gsl_spline_init(spline_h, intvar_interp, hubble_interp, hubble.size());
  gsl_spline_init(spline_x, intvar_interp, x_interp, x.size());

  return status;
  
}

// Interpolate the values of x and h between the stored array points
extern "C" double* handxofa_(double* point){
  if(intvar.size() == 0 || hubble.size() == 0 || x.size() == 0){
    printf("One of the global arrays is empty\n");
    exit(EXIT_FAILURE);
  }

  if((*point) < 0.0 || (*point) > 1.1){
    printf("Forbidden value of a : %.12f\n", (*point));
    exit(EXIT_FAILURE);
  }

  static double hx[2];

  // // Fill the vector of a with a geometric sequence
  double alpha = 0;
  int i = 0;

  // Spline interpolation
  hx[0] = gsl_spline_eval(spline_h, *point, acc);
  hx[1] = gsl_spline_eval(spline_x, *point, acc);

  // double OmegaP = (0.5*c2*pow(hx[0], 2)*pow((*point)*hx[1], 2) - 6*c3*pow(hx[0], 4)*pow((*point)*hx[1], 3) + 22.5*c4*pow(hx[0], 6)*pow((*point)*hx[1], 4) - 21*c5*pow(hx[0], 8)*pow((*point)*hx[1], 5) - 9*cG*pow(hx[0], 4)*pow((*point)*hx[1], 2))/(3.0*pow(hx[0], 2));
  // double OmegaM = 1 - OmegaP - orad/(pow((*point),4)*pow(hx[0], 2));
  // double OmTest = om/(pow((*point),3)*pow(hx[0], 2));
  // if ( fabs(OmegaM - OmTest)>1e-4  ) {
  //   printf("Integration error : %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", om, orad, c2, c3, c4, c5, cG, c0, OmTest, OmegaM, OmegaP, *point, pow(hx[0], 2), fabs(OmegaM - OmTest));
  // }

  double a3 = (*point)*(*point)*(*point);

  // if(*point < 3e-6){
  //   double hubble_LCDM = sqrt(om/a3+orad/(a3*(*point))+(1-om-orad));
  //   if(fabs((hx[0]-hubble_LCDM)/hx[0])>1e-3) fprintf(stderr, "Warning : no continuity between LCDM and galileon background at very early time ( a = %f, h_LCDM = %f and h_gal = %f)\n", (*point), hubble_LCDM, hx[0]);
  // }

  return hx;

}


extern "C" double grhogal_(double* point){

  double a2 = (*point)*(*point);
  double* hx = (*point >= 9.99999e-7) ? handxofa_(point) : 0;
  double h = (*point >= 9.99999e-7) ? (*hx) : sqrt(om/(a2*(*point))+orad/(a2*a2));
  double xgal = (*point >= 9.99999e-7) ? (*point)*(*(hx+1)) : 0; // here take xgal as a function of ln(a)

  // Define variables to save memory
  double xgal2 = xgal*xgal;
  double xgal3 = xgal2*xgal;
  double xgal4 = xgal3*xgal;
  double xgal5 = xgal4*xgal;
  double h2 = h*h;
  double h3 = h2*h;
  double h4 = h2*h2;
  double h5 = h4*h;
  double h6 = h4*h2;
  double h7 = h6*h;
  double h8 = h6*h2;

  if((*point) != acurr){
    acurr = (*point);
    hcurr = h*acurr; // 1/a*da/dtau
    xcurr = xgal; // dphi/dlna
    
    // Time evolution of the background (to get the derivatives of h and x)
    double alpha = c2/6*h*xgal-3*c3*h3*xgal2 + 15*c4*h5*xgal3 - 17.5*c5*h7*xgal4 - 3*cG*h3*xgal;
    double gamma = c2/3*h2*xgal-c3*h4*xgal2 + 2.5*c5*h8*xgal4 - 2*cG*h4*xgal;
    double beta = c2/6*h2 -2*c3*h4*xgal + 9*c4*h6*xgal2 - 10*c5*h8*xgal3 - cG*h4;
    double sigma = 2*h + 2*c3*h3*xgal3 - 15*c4*h5*xgal4 + 21*c5*h7*xgal5 + 6*cG*h3*xgal2;
    double lambda = 3*h2 + orad/(a2*a2) + c2/2*h2*xgal2 - 2*c3*h4*xgal3 + 7.5*c4*h6*xgal4 - 9*c5*h8*xgal5 - cG*h4*xgal2;
    double omega = 2*c3*h4*xgal2 - 12*c4*h6*xgal3 + 15*c5*h8*xgal4 + 4*cG*h4*xgal;

    dx = -xgal+(alpha*lambda-sigma*gamma)/(sigma*beta-alpha*omega); // derivative wrt ln(a)
    dh = (omega*gamma-lambda*beta)/(sigma*beta-alpha*omega);    
  }

  double grhogal = 0;

  if(*point >= 9.99999e-7) grhogal = h0*h0*a2*(c2/2*h2*xgal2 - 6*c3*h4*xgal3 + 22.5*c4*h4*h2*xgal4 - 21*c5*h4*h4*xgal5 - 9*cG*h4*xgal2);
  // if(*point >= 1e-4) grhogal = h0*h0*a2*(c2/2*h2*xgal2 - 6*c3*h4*xgal3 + 22.5*c4*h4*h2*xgal4 - 21*c5*h4*h4*xgal5 - 9*cG*h4*xgal2);

  //printf("grhogal : %.12f \t %.12f \t %.12f \t %.12f\n", *point, h, xgal, grhogal);

  return grhogal;

}


extern "C" double gpresgal_(){

  double a2 = acurr*acurr;
  double h = hcurr/acurr; // here h is 1/a*da/dt
  double xgal = xcurr;
  double xprime = dx;
  double hprime = dh;

  //  Define variables to save memory
  double xgal2 = xgal*xgal;
  double xgal4 = xgal2*xgal2;
  double h2 = h*h;
  double h3 = h2*h;
  double h4 = h3*h;
  double h6 = h4*h2;

  double gpresgal = 0;

  if(acurr >= 9.99999e-7) gpresgal = h0*h0*a2*(c2/2*h2*xgal2 + 2*c3*h3*xgal2*(hprime*xgal+xprime*h) - c4*(4.5*h6*xgal4 + 12*h6*xgal2*xgal*xprime + 15*h4*h*xgal4*hprime) + 3*c5*h6*h*xgal4*(5*h*xprime+7*hprime*xgal+2*h*xgal) + cG*(6*h3*xgal2*hprime + 4*h4*xgal*xprime + 3*h4*xgal2));
  // if(acurr >= 1e-4) gpresgal = h0*h0*a2*(c2/2*h2*xgal2 + 2*c3*h3*xgal2*(hprime*xgal+xprime*h) - c4*(4.5*h6*xgal4 + 12*h6*xgal2*xgal*xprime + 15*h4*h*xgal4*hprime) + 3*c5*h6*h*xgal4*(5*h*xprime+7*hprime*xgal+2*h*xgal) + cG*(6*h3*xgal2*hprime + 4*h4*xgal*xprime + 3*h4*xgal2));

  return gpresgal;

}


// Function that calculates the density perturbation of galileon
extern "C" double Chigal_(double* dgrho, double* eta, double* dphi, double* dphiprime, double* k){

  double a2 = acurr*acurr;
  double a4 = a2*a2;
  double a6 = a4*a2;
  double h = hcurr;
  double xgal = xcurr;

  //  Define variables to save memory
  double xgal2 = xgal*xgal;
  double xgal3 = xgal2*xgal;
  double xgal4 = xgal3*xgal;
  double xgal5 = xgal4*xgal;
  double h2 = h*h;
  double h3 = h2*h;
  double h4 = h3*h;
  double h6 = h4*h2;

  double alpha_Z = -2*c3/a2*xgal3*h2 
    + 15*c4/a4*xgal4*h4 
    - 21*c5/a6*xgal5*h6 
    - 6*cG/a2*xgal2*h2;
  double alpha_eta = 1.5*c4/a4*xgal4*h4
    - 3*c5/a6*xgal5*h6
    - cG/a2*xgal2*h2;
  double beta = 1./(1-0.5*alpha_Z);
  double ChitildeG = c2*h0*xgal*h*(*dphiprime) 
    - c3/a2*(18*h0*xgal2*h3*(*dphiprime) + 2*pow(*k, 2)*xgal2*h2*(*dphi)) 
    + c4/a4*(90*h0*xgal3*h4*h*(*dphiprime) + 12*pow(*k, 2)*xgal3*h4*(*dphi)) 
    - c5/a6*(105*h0*xgal4*h6*h*(*dphiprime) + 15*pow(*k, 2)*xgal4*h6*(*dphi)) 
    - cG/a2*(18*h0*xgal*h3*(*dphiprime) + 4*pow(*k, 2)*xgal*h2*(*dphi));

  double ChiG = 0;

  if (acurr >= 9.99999e-7) ChiG = beta*(ChitildeG + 0.5*alpha_Z*(*dgrho) + (alpha_Z - 2*alpha_eta)*(*k)*(*eta));
  // if (acurr >= 1e-4) ChiG = beta*(ChitildeG + 0.5*alpha_Z*(*dgrho) + (alpha_Z - 2*alpha_eta)*(*k)*(*eta));

  // FILE* g = fopen("chigal/chigal_q0001.dat", "a");
  // fprintf(g, "%.16f ; %.16f ; %.16f ; %.16f ; %.16f ; %.16f\n", acurr, ChiG, (*dgrho), beta*ChitildeG, 0.5*beta*alpha_Z*(*dgrho), beta*(alpha_Z - 2*alpha_eta)*(*k)*(*eta));
  // fclose(g);

  return ChiG;

}


// Perturbation of heat flux from galileon
extern "C" double qgal_(double* dgq, double* eta, double* dphi, double* dphiprime, double* k){

  double a2 = acurr*acurr;
  double a4 = a2*a2;
  double a6 = a4*a2;
  double h = hcurr;
  double xgal = xcurr;

  //  Define variables to save memory
  double xgal2 = xgal*xgal;
  double xgal3 = xgal2*xgal;
  double xgal4 = xgal3*xgal;
  double h2 = h*h;
  double h3 = h2*h;
  double h4 = h3*h;
  double h6 = h4*h2;

  double alpha = c4/a4*xgal4*h4 
    - 2*c5/a6*xgal4*xgal*h6 
    - cG/(1.5*a2)*xgal2*h2;
  double beta = 1./(1-1.5*alpha);
  double qtildeG = c2*(*k)*h0*xgal*h*(*dphi) 
    - c3/a2*(*k)*(6*h0*xgal2*h3*(*dphi) - 2*xgal2*h2*(*dphiprime)) 
    + c4/a4*(*k)*(-12*xgal3*h4*(*dphiprime) + 18*h0*xgal3*h4*h*(*dphi)) 
    - c5/a6*(*k)*(-15*xgal4*h6*(*dphiprime) + 15*h0*xgal4*h6*h*(*dphi)) 
    - cG/a2*(*k)*(-4*xgal*h2*(*dphiprime) + 6*h0*xgal*h3*(*dphi));

  double qG = 0;

  // if(-1e-5 < 1.5*alpha - 1 && 1.5*alpha - 1 < 1e-5) printf("WARNING : 1/beta_q is zero");
  if(acurr >= 9.99999e-7) qG = beta*(qtildeG + 1.5*alpha*(*dgq));
  // if(acurr >= 1e-4) qG = beta*(qtildeG + 1.5*alpha*(*dgq));

  // FILE* g = fopen("qgal/qgal_q0001.dat", "a");
  // fprintf(g, "%.16f ; %.16f ; %.16f ; %.16f ; %.16f\n", acurr, qG, (*dgq), beta*qtildeG, 1.5*beta*alpha*(*dgq));
  // fclose(g);

  return qG;

}


// Perturbation of anisotropic stress from galileon
extern "C" double Pigal_(double* dgrho, double* dgq, double* dgpi, double* eta, double* dphi, double* k){

  double a2 = acurr*acurr;
  double a4 = a2*a2;
  double a6 = a4*a2;
  double h = hcurr;
  double xgal = xcurr;

  //  Define variables to save memory
  double xgal2 = xgal*xgal;
  double xgal3 = xgal2*xgal;
  double xgal4 = xgal3*xgal;
  double xgal5 = xgal4*xgal;
  double h2 = h*h;
  double h3 = h2*h;
  double h4 = h3*h;
  double h5 = h4*h;
  double h6 = h5*h;

  // Derivatives of h and x
  double xprime = dx;
  double hprime = acurr*dh + h;
  double xhprime = xprime*h + xgal*hprime; // Derivative of the product (xgal*H)'

  double alpha_sigprime = c4/a4*xgal4*h4
    - 3*c5/a6*xgal4*h5*xhprime;
  double alpha_sig = c4/a4*(3*xgal4*h4 - 6*xgal3*h3*xhprime)
    - c5/a6*(-3*xgal5*h5*hprime + 12*xgal5*h6 - 15*xgal4*h5*xhprime)
    + 2*cG/a2*xgal*h*xhprime;
  double alpha_phi = -c4/a4*xgal4*h4
    - c5/a6*(6*xgal4*h5*xhprime - 6*xgal5*h6)
    + 2*cG/a2*xgal2*h2;
  double beta_pi = 1./(1+0.5*alpha_phi - alpha_sigprime);
  double PitildeG = pow(*k, 2)*(c4/a4*(4*xgal3*h4*(*dphi) - 6*xgal2*h3*xhprime*(*dphi))
    - c5/a6*(-12*xgal3*h5*xhprime*(*dphi) + 12*xgal4*h6*(*dphi) - 3*xgal4*h5*hprime*(*dphi))
    + 2*cG/a2*h*xhprime*(*dphi));

  double PiG = 0;

  // if(-1e-5 < (alpha_phi - 2*alpha_sigprime+2) && (alpha_phi - 2*alpha_sigprime+2) < 1e-5) printf("WARNING : 1/beta_pi is zero");
  if(acurr >= 9.99999e-7) PiG = beta_pi*(PitildeG + (alpha_sigprime - 0.5*alpha_phi)*(*dgpi) + 0.5*(2*alpha_sigprime + alpha_sig - alpha_phi)*((*dgrho) + 3*h0*h/(*k)*(*dgq)) + (alpha_sig + alpha_sigprime)*(*k)*(*eta));
  // if(acurr >= 1e-4) PiG = beta_pi*(PitildeG + (alpha_sigprime - 0.5*alpha_phi)*(*dgpi) + 0.5*(2*alpha_sigprime + alpha_sig - alpha_phi)*((*dgrho) + 3*h0*h/(*k)*(*dgq)) + (alpha_sig + alpha_sigprime)*(*k)*(*eta));

  // FILE* g = fopen("pigal/pigal_q0001.dat", "a");
  // fprintf(g, "%.16f ; %.16f ; %.16f ; %.16f ; %.16f ; %.16f ; %.16f ; %.16f\n", acurr, PiG, (*dgpi), beta_pi*PitildeG, beta_pi*(0.5*alpha_phi - alpha_sigprime)*(*dgpi), beta_pi*(0.5*alpha_sig - 0.5*alpha_phi + alpha_sigprime)*(*dgrho), 3*beta_pi*(0.5*alpha_sig - 0.5*alpha_phi + alpha_sigprime)*h0*h/(*k)*(*dgq), beta_pi*(alpha_sig + alpha_sigprime)*(*k)*(*eta));
  // fclose(g);

  return PiG;

}

// Perturbation of the galileon field
extern "C" double dphisecond_(double* dgrho, double* dgq, double* eta, double* dphi, double* dphiprime, double* k, double* deltafprime){

  double a2 = acurr*acurr;
  double a4 = a2*a2;
  double a6 = a4*a2;
  double h = hcurr;
  double xgal = xcurr;

  //  Define variables to save memory
  double xgal2 = xgal*xgal;
  double xgal3 = xgal2*xgal;
  double xgal4 = xgal3*xgal;
  double xgal5 = xgal4*xgal;
  double h2 = h*h;
  double h3 = h2*h;
  double h4 = h3*h;
  double h5 = h4*h;
  double h6 = h5*h;
  double h7 = h6*h;
  double h8 = h7*h;

  // Derivatives of h and x
  double xprime = dx;
  double hprime = acurr*dh + h;
  double xhprime = xprime*h + xgal*hprime; // Derivative of the product (xgal*H)'

  // Expression of dphisecond
  double alpha_gammasecond = c2
    - 12*c3/a2*xgal*h2
    + 54*c4/a4*xgal2*h4
    - 60*c5/a6*xgal3*h6
    - 6*cG/a2*h2;
  double alpha_gammaprime = 2*c2
    - c3/a2*(12*h*xhprime + 12*xgal*h*hprime)
    + c4/a4*(-108*xgal2*h4 + 108*xgal*h3*xhprime + 108*xgal2*h3*hprime)
    - c5/a6*(-240*xgal3*h6 + 180*xgal2*h5*xhprime + 180*xgal3*h5*hprime)
    - 12*cG/a2*h*hprime;
  double alpha_gamma = c2
    - c3/a2*(4*xgal*h2 + 4*h*xhprime)
    + c4/a4*(-10*xgal2*h4 + 24*xgal*h3*xhprime + 12*xgal2*h3*hprime)
    - c5/a6*(-36*xgal3*h6 + 36*xgal2*h5*xhprime + 24*xgal3*h5*hprime)
    - cG/a2*(2*h2 + 4*h*hprime);
  double alpha_Z = c2*xgal
    - c3/a2*(6*xgal2*h2 + 4*xgal*h*xhprime)
    + c4/a4*(-6*xgal3*h4 + 36*xgal2*h3*xhprime + 12*xgal3*h3*hprime)
    - c5/a6*(-45*xgal4*h6 + 60*xgal3*h5*xhprime + 30*xgal4*h5*hprime)
    - cG/a2*(6*xgal*h2 + 4*h*xhprime + 4*xgal*h*hprime);
  double alpha_Zprime = -2*c3/a2*xgal2*h2
    + 12*c4/a4*xgal3*h4
    - 15*c5/a6*xgal4*h6
    - 4*cG/a2*xgal*h2;
  double alpha_eta = c4/a4*(-4*xgal3*h4 + 6*xgal2*h3*xhprime)
    - c5/a6*(-12*xgal4*h6 + 3*xgal4*h5*hprime + 12*xgal3*h5*xhprime)
    - 2*cG/a2*h*xhprime;
  
  // Derivatives of density perturbations
  double dotdeltaf = (*deltafprime);

  // Expression of Z'
  double chiprimehat = c2*(-2*pow(h0, 2)*xgal*h2*(*dphiprime) + pow(h0, 2)*xgal*h*hprime*(*dphiprime) + pow(h0, 2)*h2*xprime*(*dphiprime))
    - c3/a2*(-72*pow(h0, 2)*xgal2*h4*(*dphiprime) + 54*pow(h0, 2)*xgal2*h3*hprime*(*dphiprime) + 36*pow(h0, 2)*xgal*h4*xprime*(*dphiprime)
  			 + pow(*k, 2)*(2*xgal2*h2*(*dphiprime) - 8*h0*xgal2*h3*(*dphi) + 4*h0*xgal2*h2*hprime*(*dphi) + 4*h0*xgal*h3*xprime*(*dphi)))
    + c4/a4*(-540*pow(h0, 2)*xgal3*h6*(*dphiprime) + 450*pow(h0, 2)*xgal3*h5*hprime*(*dphiprime) + 270*pow(h0, 2)*xgal2*h6*xprime*(*dphiprime)
  			 + pow(*k, 2)*(12*xgal3*h4*(*dphiprime) - 72*h0*xgal3*h5*(*dphi) + 48*h0*xgal3*h4*hprime*(*dphi) + 36*h0*xgal2*h5*xprime*(*dphi)))
    - c5/a6*(-840*pow(h0, 2)*xgal4*h8*(*dphiprime) + 735*pow(h0, 2)*xgal4*h7*hprime*(*dphiprime) + 420*pow(h0, 2)*xgal3*h8*xprime*(*dphiprime)
  			 + pow(*k, 2)*(15*xgal4*h6*(*dphiprime) - 120*h0*xgal4*h7*(*dphi) + 90*h0*xgal4*h6*hprime*(*dphi) + 60*h0*xgal3*h7*xprime*(*dphi)))
    - cG/a2*(-72*pow(h0, 2)*xgal*h4*(*dphiprime) + 54*pow(h0, 2)*xgal*h3*hprime*(*dphiprime) + 18*pow(h0, 2)*h4*xprime*(*dphiprime)
  			 + pow(*k, 2)*(4*xgal*h2*(*dphiprime) - 16*h0*xgal*h3*(*dphi) + 8*h0*xgal*h2*hprime*(*dphi) + 4*h0*h3*xprime*(*dphi)))
    ;
  double beta_gammasecond = c2*h0*xgal*h - 18*c3*h0/a2*xgal2*h3 + 90*c4*h0/a4*xgal3*h5 - 105*c5*h0/a6*xgal4*h7 - 18*cG*h0/a2*xgal*h3;
  double beta_Z = -2*c3/a2*xgal3*h2
    + 15*c4/a4*xgal4*h4
    - 21*c5/a6*xgal5*h6
    - 6*cG/a2*xgal2*h2;
  double beta_eta = 1.5*c4/a4*xgal4*h4
    - 3*c5/a6*xgal5*h6
    - cG/a2*xgal2*h2;
  double beta_Z_prime = -c3*h0/a2*(-4*xgal3*h3 + 4*xgal3*h2*hprime + 6*xgal2*h3*xprime)
    + c4*h0/a4*(-60*xgal4*h5 + 60*xgal4*h4*hprime + 60*xgal3*h5*xprime)
    - c5*h0/a6*(-126*xgal5*h7 + 126*xgal5*h6*hprime + 105*xgal4*h7*xprime)
    - cG*h0/a2*(-12*xgal2*h3 + 12*xgal2*h2*hprime + 12*xgal*h3*xprime);
  double beta_eta_prime = c4*h0/a4*(-6*xgal4*h5 + 6*xgal4*h4*hprime + 6*xgal3*h5*xprime)
    - c5*h0/a6*(-18*xgal5*h7 + 18*xgal5*h6*hprime + 15*xgal4*h7*xprime)
    - cG*h0/a2*(-2*xgal2*h3 + 2*xgal2*h2*hprime + 2*xgal*h3*xprime);


  double ksi = alpha_Zprime/(h0*h*(2-beta_Z));
  double dphisecond = 0;
  if(acurr >= 9.99999e-7) dphisecond = -(alpha_gammaprime*h0*h*(*dphiprime) + alpha_gamma*pow(*k, 2)*(*dphi) + (0.5*alpha_Z + 2*ksi*(h0*h - 0.5*(h0*hprime + beta_Z*h0*h) + 0.25*(beta_Z_prime + beta_Z*h0*hprime)))*(*dgrho) + ksi*(1-beta_eta)*(*k)*(*dgq) + ksi*dotdeltaf + ksi*chiprimehat + (alpha_Z - 2*alpha_eta + 2*ksi*(2*beta_eta*h0*h - h0*hprime - beta_Z*h0*h - beta_eta_prime + 0.5*(beta_Z_prime + beta_Z*h0*hprime)))*(*k)*(*eta))/(alpha_gammasecond + ksi*beta_gammasecond);
  // if(acurr >= 1e-4) dphisecond = -(alpha_gammaprime*h0*h*(*dphiprime) + alpha_gamma*pow(*k, 2)*(*dphi) + (0.5*alpha_Z + 2*ksi*(h0*h - 0.5*(h0*hprime + beta_Z*h0*h) + 0.25*(beta_Z_prime + beta_Z*h0*hprime)))*(*dgrho) + ksi*(1-beta_eta)*(*k)*(*dgq) + ksi*dotdeltaf + ksi*chiprimehat + (alpha_Z - 2*alpha_eta + 2*ksi*(2*beta_eta*h0*h - h0*hprime - beta_Z*h0*h - beta_eta_prime + 0.5*(beta_Z_prime + beta_Z*h0*hprime)))*(*k)*(*eta))/(alpha_gammasecond + ksi*beta_gammasecond);

  // FILE* g = fopen("dphisecond/dphisecond_q0001.dat", "a");
  // fprintf(g, "%.16f ; %.16f ; %.16f ; %.16f ; %.16f ; %.16f ; %.16f ; %.16f\n", acurr, dphisecond, alpha_gammaprime/alpha_gammasecond*h0*h*(*dphiprime), alpha_gamma/alpha_gammasecond*pow(*k, 2)*(*dphi), 0.5/alpha_gammasecond*(alpha_Z - 2*alpha_Zprime)*(*dgrho), (alpha_Z - alpha_Zprime - 2*alpha_eta)/alpha_gammasecond*(*k)*(*eta), (*dphiprime), (*dphi));
  // fclose(g);

  return dphisecond;

}

// Calculate the conformal time derivative of pigal
extern "C" double pigalprime_(double* dgrho, double* dgq, double* dgpi, double* pidot, double* eta, double* dphi, double* dphiprime, double* k, double* grho, double* gpres){

  double a2 = acurr*acurr;
  double a4 = a2*a2;
  double a6 = a4*a2;
  double h = hcurr;
  double hoft = hcurr/acurr;
  double xgal = xcurr;

  //  Define variables to save memory
  double hoft2 = hoft*hoft;
  double hoft3 = hoft2*hoft;
  double hoft4 = hoft3*hoft;
  double hoft5 = hoft4*hoft;
  double hoft6 = hoft5*hoft;
  double hoft7 = hoft6*hoft;
  double hoft8 = hoft7*hoft;
  double xgal2 = xgal*xgal;
  double xgal3 = xgal2*xgal;
  double xgal4 = xgal3*xgal;
  double xgal5 = xgal4*xgal;
  double h2 = h*h;
  double h3 = h2*h;
  double h4 = h3*h;
  double h5 = h4*h;
  double h6 = h5*h;
  double h7 = h6*h;

  // Derivatives of h and x wrt lna
  double alpha = c2/6*hoft*xgal-3*c3*hoft3*xgal2 + 15*c4*hoft5*xgal3 - 17.5*c5*hoft7*xgal4 - 3*cG*hoft3*xgal;
  double gamma = c2/3*hoft2*xgal-c3*hoft4*xgal2 + 2.5*c5*hoft8*xgal4 - 2*cG*hoft4*xgal;
  double beta = c2/6*hoft2 -2*c3*hoft4*xgal + 9*c4*hoft6*xgal2 - 10*c5*hoft8*xgal3 - cG*hoft4;
  double delta = 2*hoft + 2*c3*hoft3*xgal3 - 15*c4*hoft5*xgal4 + 21*c5*hoft7*xgal5 + 6*cG*hoft3*xgal2;
  double lambda = 3*hoft2 + orad/(a2*a2) + c2/2*hoft2*xgal2 - 2*c3*hoft4*xgal3 + 7.5*c4*hoft6*xgal4 - 9*c5*hoft8*xgal5 - cG*hoft4*xgal2;
  double omega = 2*c3*hoft4*xgal2 - 12*c4*hoft6*xgal3 + 15*c5*hoft8*xgal4 + 4*cG*hoft4*xgal;

  double xprime = dx;
  double hprime = acurr*dh + h;
  double xhprime = xprime*h + xgal*hprime; // Derivative of the product (xgal*H)'

  // Second derivatives of h and x wrt lna
  double alpha_prime = (c2/6*xgal - 9*c3*hoft2*xgal2 + 75*c4*hoft4*xgal3 - 122.5*c5*hoft6*xgal4 - 9*cG*hoft2*xgal)*(hprime-h)/acurr + (c2/6*hoft - 6*c3*hoft3*xgal + 45*c4*hoft5*xgal2 - 70*c5*hoft7*xgal3 - 3*cG*hoft3)*xprime;
  double gamma_prime = (2.*c2/3*hoft*xgal - 4*c3*hoft3*xgal2 + 20*c5*hoft7*xgal4 - 8*cG*hoft3*xgal)*(hprime-h)/acurr + (c2/3*hoft2 - 2*c3*hoft4*xgal + 10*c5*hoft8*xgal3 - 2*cG*hoft4)*xprime;
  double beta_prime = (c2/3*hoft - 8*c3*hoft3*xgal + 54*c4*hoft5*xgal2 - 80*c5*hoft7*xgal3 - 4*cG*hoft3)*(hprime-h)/acurr + (-2*c3*hoft4 + 18*c4*hoft6*xgal - 30*c5*hoft8*xgal2)*xprime;
  double delta_prime = (2 + 6*c3*hoft2*xgal3 - 75*c4*hoft4*xgal4 + 147*c5*hoft6*xgal5 + 18*cG*hoft2*xgal2)*(hprime-h)/acurr + (6*c3*hoft3*xgal2 - 60*c4*hoft5*xgal3 + 105*c5*hoft7*xgal4 + 12*cG*hoft3*xgal)*xprime;
  double lambda_prime = -4*orad/a4 + (6*hoft + c2*hoft*xgal2 - 8*c3*hoft3*xgal3 + 45*c4*hoft5*xgal4 - 72*c5*hoft7*xgal5 - 4*cG*hoft3*xgal2)*(hprime-h)/acurr + (c2*hoft2*xgal - 6*c3*hoft4*xgal2 + 30*c4*hoft6*xgal3 - 45*c5*hoft8*xgal4 - 2*cG*hoft4*xgal)*xprime;
  double omega_prime = (8*c3*hoft3*xgal2 - 72*c4*hoft5*xgal3 + 120*c5*hoft7*xgal4 + 16*cG*hoft3*xgal)*(hprime-h)/acurr + (4*c3*hoft4*xgal - 36*c4*hoft6*xgal2 + 60*c5*hoft8*xgal3 + 4*cG*hoft4)*xprime;

  double xprimedot = h0*h*(-xprime + ((alpha_prime*lambda + alpha*lambda_prime - delta_prime*gamma - delta*gamma_prime)*(delta*beta-alpha*omega) - (alpha*lambda-delta*gamma)*(delta_prime*beta + delta*beta_prime - alpha_prime*omega - alpha*omega_prime))/pow(delta*beta-alpha*omega, 2));
  double hprimedot = acurr*h0*h*((omega*gamma-lambda*beta)/(delta*beta-alpha*omega) + ((omega_prime*gamma + omega*gamma_prime - lambda_prime*beta - lambda*beta_prime)*(delta*beta-alpha*omega) - (omega*gamma-lambda*beta)*(delta_prime*beta + delta*beta_prime - alpha_prime*omega - alpha*omega_prime))/pow(delta*beta-alpha*omega, 2)) + h0*h*hprime;

  double piprimetilde = pow(*k, 2)*c4/a4*(4*xgal3*h4*(*dphiprime) - 6*xgal3*h3*hprime*(*dphiprime) - 6*xgal2*h4*xprime*(*dphiprime) - 24*h0*xgal3*h5*(*dphi) + 52*h0*xgal3*h4*hprime*(*dphi) - 18*h0*xgal3*h3*pow(hprime, 2)*(*dphi) + 48*h0*xgal2*h5*xprime*(*dphi) - 42*h0*xgal2*h4*hprime*xprime*(*dphi) - 12*h0*xgal*h5*pow(xprime, 2)*(*dphi) - 6*xgal3*h3*hprimedot*(*dphi) - 6*xgal2*h4*xprimedot*(*dphi))
    - pow(*k, 2)*c5/a6*(12*xgal4*h6*(*dphiprime) - 15*xgal4*h5*hprime*(*dphiprime) - 12*xgal3*h6*xprime*(*dphiprime) - 96*h0*xgal4*h7*(*dphi) + 192*h0*xgal4*h6*hprime*(*dphi) - 75*h0*xgal4*h5*pow(hprime, 2)*(*dphi) + 144*h0*xgal3*h7*xprime*(*dphi) - 132*h0*xgal3*h6*hprime*xprime*(*dphi) - 36*h0*xgal2*h7*pow(xprime, 2)*(*dphi) - 15*xgal4*h5*hprimedot*(*dphi) - 12*xgal3*h6*xprimedot*(*dphi))
    - pow(*k, 2)*cG/a2*(-2*xgal*h*hprime*(*dphiprime) - 2*h2*xprime*(*dphiprime) + 8*h0*xgal*h2*hprime*(*dphi) - 2*h0*xgal*h*pow(hprime, 2)*(*dphi) + 8*h0*h3*xprime*(*dphi) - 6*h0*h2*hprime*xprime*(*dphi) - 2*xgal*h*hprimedot*(*dphi) - 2*h2*xprimedot*(*dphi));

  double alpha_sig = c4/a4*(3*xgal4*h4 - 6*xgal3*h3*xhprime)
    - c5/a6*(-3*xgal5*h5*hprime + 12*xgal5*h6 - 15*xgal4*h5*xhprime)
    + 2*cG/a2*xgal*h*xhprime;
  double alpha_sigprime = c4/a4*xgal4*h4
    - 3*c5/a6*xgal4*h5*xhprime;
  double alpha_phi = -c4/a4*xgal4*h4
    - c5/a6*(6*xgal4*h5*xhprime - 6*xgal5*h6)
    + 2*cG/a2*xgal2*h2;
  
  double alpha_sig_prime = c4/a4*(-12*h0*xgal4*h5 + 36*h0*xgal4*h4*hprime - 18*h0*xgal4*h3*pow(hprime, 2) + 36*h0*xgal3*h5*xprime - 48*h0*xgal3*h4*hprime*xprime - 18*h0*xgal2*h5*pow(xprime, 2) - 6*xgal4*h3*hprimedot - 6*xgal3*h4*xprimedot)
    - c5/a6*(-72*h0*xgal5*h7 + 180*h0*xgal5*h6*hprime - 90*h0*xgal5*h5*pow(hprime, 2) + 150*h0*xgal4*h7*xprime - 180*h0*xgal4*h6*hprime*xprime - 60*h0*xgal3*h7*pow(xprime, 2) - 18*xgal5*h5*hprimedot - 15*xgal4*h6*xprimedot)
    - cG/a2*(4*h0*xgal2*h2*hprime - 2*h0*xgal2*h*pow(hprime, 2) + 4*h0*xgal*h3*xprime - 8*h0*xgal*h2*hprime*xprime - 2*h0*h3*pow(xprime, 2) - 2*xgal2*h*hprimedot - 2*xgal*h2*xprimedot);
  double alpha_sigprime_prime = c4/a4*(-4*h0*xgal4*h5 + 4*h0*xgal4*h4*hprime + 4*h0*xgal3*h5*xprime)
    - c5/a6*(-18*h0*xgal5*h6*hprime + 15*h0*xgal5*h5*pow(hprime, 2) - 18*h0*xgal4*h7*xprime + 33*h0*xgal4*h6*hprime*xprime + 12*h0*xgal3*h7*pow(xprime, 2) + 3*xgal5*h5*hprimedot + 3*xgal4*h6*xprimedot);
  double alpha_phi_prime = c4/a4*(4*h0*xgal4*h5 - 4*h0*xgal4*h4*hprime - 4*h0*xgal3*h5*xprime)
    - c5/a6*(36*h0*xgal5*h7 - 72*h0*xgal5*h6*hprime + 30*h0*xgal5*h5*pow(hprime, 2) - 66*h0*xgal4*h7*xprime + 66*h0*xgal4*h6*hprime*xprime + 24*h0*xgal3*h7*pow(xprime, 2) + 6*xgal5*h5*hprimedot + 6*xgal4*h6*xprimedot)
    - cG/a2*(4*h0*xgal2*h3 - 4*h0*xgal2*h2*hprime - 4*h0*xgal*h3*xprime);

  double piGdot = 0;

  if(acurr >= 9.99999e-7) piGdot = (piprimetilde + (alpha_sigprime - 0.5*alpha_phi)*(*pidot) + (-3.5*h0*h*alpha_sigprime + 0.5*h0*hprime*alpha_sigprime + 0.25*(alpha_phi - alpha_sigprime)*((*grho)+(*gpres))/(h0*h) + 0.5*h0*hprime*alpha_sig + 0.5*alpha_sig_prime + alpha_sigprime_prime - 2*h0*h*alpha_sig + 1.5*h0*h*alpha_phi - 0.5*alpha_phi_prime)*(*dgrho) + (alpha_sig_prime + alpha_sigprime_prime + alpha_sigprime*h0*hprime + alpha_sig*h0*hprime - 3*h0*h*alpha_sigprime - 3*h0*h*alpha_sig + 0.5*(alpha_phi - alpha_sigprime)*((*grho)+(*gpres))/(h0*h))*(*k)*(*eta) + 0.5/(*k)*(3*pow(h0, 2)*h*hprime*alpha_sigprime + 3*pow(h0, 2)*h*hprime*alpha_sig + 1.5*(alpha_phi - alpha_sigprime)*((*grho)+(*gpres)) + 3*h0*h*alpha_sig_prime - 12*pow(h0, 2)*h2*alpha_sig - 21*pow(h0, 2)*h2*alpha_sigprime - 3*h0*h*alpha_phi_prime + 9*pow(h0, 2)*h2*alpha_phi + 6*h0*h*alpha_sigprime_prime + pow(*k, 2)*(alpha_phi - alpha_sigprime))*(*dgq) + (-2*h0*h*alpha_sigprime + h0*h*alpha_phi - 0.5*alpha_phi_prime + alpha_sigprime_prime - h0*h*alpha_sig)*(*dgpi))/(1. - alpha_sigprime + 0.5*alpha_phi);
  // if(acurr >= 1e-4) piGdot = (piprimetilde + (alpha_sigprime - 0.5*alpha_phi)*(*pidot) + (-3.5*h0*h*alpha_sigprime + 0.5*h0*hprime*alpha_sigprime + 0.25*(alpha_phi - alpha_sigprime)*((*grho)+(*gpres))/(h0*h) + 0.5*h0*hprime*alpha_sig + 0.5*alpha_sig_prime + alpha_sigprime_prime - 2*h0*h*alpha_sig + 1.5*h0*h*alpha_phi - 0.5*alpha_phi_prime)*(*dgrho) + (alpha_sig_prime + alpha_sigprime_prime + alpha_sigprime*h0*hprime + alpha_sig*h0*hprime - 3*h0*h*alpha_sigprime - 3*h0*h*alpha_sig + 0.5*(alpha_phi - alpha_sigprime)*((*grho)+(*gpres))/(h0*h))*(*k)*(*eta) + 0.5/(*k)*(3*pow(h0, 2)*h*hprime*alpha_sigprime + 3*pow(h0, 2)*h*hprime*alpha_sig + 1.5*(alpha_phi - alpha_sigprime)*((*grho)+(*gpres)) + 3*h0*h*alpha_sig_prime - 12*pow(h0, 2)*h2*alpha_sig - 21*pow(h0, 2)*h2*alpha_sigprime - 3*h0*h*alpha_phi_prime + 9*pow(h0, 2)*h2*alpha_phi + 6*h0*h*alpha_sigprime_prime + pow(*k, 2)*(alpha_phi - alpha_sigprime))*(*dgq) + (-2*h0*h*alpha_sigprime + h0*h*alpha_phi - 0.5*alpha_phi_prime + alpha_sigprime_prime - h0*h*alpha_sig)*(*dgpi))/(1. - alpha_sigprime + 0.5*alpha_phi);

  return piGdot;

}

// Cross checks that independent equations are satisfied by perturbations
extern "C" double*  crosschecks_(double* dgrho, double* dgq, double* dgpi, double* eta, double* dphi, double* dphiprime, double* dphiprimeprime, double* k, double* grho, double* gpres, double* deltafprime){

  // static double eq[2];
  static double eq[7];

  double a2 = acurr*acurr;
  double a4 = a2*a2;
  double a6 = a4*a2;
  double h = hcurr;
  double xgal = xcurr;

  //  Define variables to save memory
  double xgal2 = xgal*xgal;
  double xgal3 = xgal2*xgal;
  double xgal4 = xgal3*xgal;
  double xgal5 = xgal4*xgal;
  double h2 = h*h;
  double h3 = h2*h;
  double h4 = h3*h;
  double h5 = h4*h;
  double h6 = h5*h;
  double h7 = h6*h;
  double h8 = h7*h;

  // Derivatives of h and x
  double xprime = dx;
  double hprime = acurr*dh + h;
  double xhprime = xprime*h + xgal*hprime; // Derivative of the product (xgal*H)'

  // // Derivatives of density perturbations of other fluids
  double dotdeltaf = (*deltafprime);

  // Expression of Z' and ChiP to be used in the conservation equations
  double chiprimehat = c2*(-2*pow(h0, 2)*xgal*h2*(*dphiprime) + pow(h0, 2)*xgal*h*hprime*(*dphiprime) + pow(h0, 2)*h2*xprime*(*dphiprime))
    - c3/a2*(-72*pow(h0, 2)*xgal2*h4*(*dphiprime) + 54*pow(h0, 2)*xgal2*h3*hprime*(*dphiprime) + 36*pow(h0, 2)*xgal*h4*xprime*(*dphiprime)
  			 + pow(*k, 2)*(2*xgal2*h2*(*dphiprime) - 8*h0*xgal2*h3*(*dphi) + 4*h0*xgal2*h2*hprime*(*dphi) + 4*h0*xgal*h3*xprime*(*dphi)))
    + c4/a4*(-540*pow(h0, 2)*xgal3*h6*(*dphiprime) + 450*pow(h0, 2)*xgal3*h5*hprime*(*dphiprime) + 270*pow(h0, 2)*xgal2*h6*xprime*(*dphiprime)
  			 + pow(*k, 2)*(12*xgal3*h4*(*dphiprime) - 72*h0*xgal3*h5*(*dphi) + 48*h0*xgal3*h4*hprime*(*dphi) + 36*h0*xgal2*h5*xprime*(*dphi)))
    - c5/a6*(-840*pow(h0, 2)*xgal4*h8*(*dphiprime) + 735*pow(h0, 2)*xgal4*h7*hprime*(*dphiprime) + 420*pow(h0, 2)*xgal3*h8*xprime*(*dphiprime)
  			 + pow(*k, 2)*(15*xgal4*h6*(*dphiprime) - 120*h0*xgal4*h7*(*dphi) + 90*h0*xgal4*h6*hprime*(*dphi) + 60*h0*xgal3*h7*xprime*(*dphi)))
    - cG/a2*(-72*pow(h0, 2)*xgal*h4*(*dphiprime) + 54*pow(h0, 2)*xgal*h3*hprime*(*dphiprime) + 18*pow(h0, 2)*h4*xprime*(*dphiprime)
  			 + pow(*k, 2)*(4*xgal*h2*(*dphiprime) - 16*h0*xgal*h3*(*dphi) + 8*h0*xgal*h2*hprime*(*dphi) + 4*h0*h3*xprime*(*dphi)));
  double beta_gammasecond = c2*h0*xgal*h - 18*c3*h0/a2*xgal2*h3 + 90*c4*h0/a4*xgal3*h5 - 105*c5*h0/a6*xgal4*h7 - 18*cG*h0/a2*xgal*h3;

  double beta_Z = -2*c3/a2*xgal3*h2
    + 15*c4/a4*xgal4*h4
    - 21*c5/a6*xgal5*h6
    - 6*cG/a2*xgal2*h2;
  double beta_eta = 1.5*c4/a4*xgal4*h4
    - 3*c5/a6*xgal5*h6
    - cG/a2*xgal2*h2;
  double beta_Z_prime = -c3*h0/a2*(-4*xgal3*h3 + 4*xgal3*h2*hprime + 6*xgal2*h3*xprime)
    + c4*h0/a4*(-60*xgal4*h5 + 60*xgal4*h4*hprime + 60*xgal3*h5*xprime)
    - c5*h0/a6*(-126*xgal5*h7 + 126*xgal5*h6*hprime + 105*xgal4*h7*xprime)
    - cG*h0/a2*(-12*xgal2*h3 + 12*xgal2*h2*hprime + 12*xgal*h3*xprime);
  double beta_eta_prime = c4*h0/a4*(-6*xgal4*h5 + 6*xgal4*h4*hprime + 6*xgal3*h5*xprime)
    - c5*h0/a6*(-18*xgal5*h7 + 18*xgal5*h6*hprime + 15*xgal4*h7*xprime)
    - cG*h0/a2*(-2*xgal2*h3 + 2*xgal2*h2*hprime + 2*xgal*h3*xprime);

  // kHZ'
  double khzprime = ((h0*h - 0.5*(h0*hprime + beta_Z*h0*h) + 0.25*(beta_Z_prime + beta_Z*h0*hprime))*(*dgrho) + 0.5*(1-beta_eta)*(*k)*(*dgq) - (h0*hprime + beta_Z*h0*h + beta_eta_prime - 2*beta_eta*h0*h - 0.5*(beta_Z_prime + beta_Z*h0*hprime))*(*k)*(*eta) + 0.5*dotdeltaf + 0.5*chiprimehat + 0.5*beta_gammasecond*(*dphiprimeprime))/(1-0.5*beta_Z); //Barreira

  // Derivative of ChiG
  // double dotdeltagal = chiprimehat + beta_gammasecond*(*dphiprimeprime) + 0.5*(beta_Z_prime + beta_Z*h0*hprime - 2*beta_Z*h0*h)*(*dgrho) + beta_Z*khzprime + (beta_Z_prime + beta_Z*h0*hprime - 2*beta_Z*h0*h - 2*beta_eta_prime + 4*beta_eta*h0*h)*(*k)*(*eta) - beta_eta*(*k)*(*dgq);

  // Derivative of q and Chi
  double dotdeltaq = -2./3.*((*k)*khzprime/(h0*h) + pow(*k, 2)*(*eta) + (*k)*((*dgpi) + (*dgrho)) + 6*h0*h*(*dgq));
  // double dotdelta = dotdeltaf + dotdeltagal;
  double dotdelta = 2*khzprime + (h0*hprime - 2*h0*h)*(*dgrho) + 2*h0*hprime*(*k)*(*eta) - (*k)*(*dgq);

  // Equation (49) of arXiv:1208.0600
  std::vector<double> eq49;
  eq49.push_back(dotdelta);
  eq49.push_back(((*grho)+(*gpres))/(h0*h)*(0.5*(*dgrho) + (*k)*(*eta)));
  eq49.push_back(h0*h*((*dgrho) - 2*(*k)*(*eta)) - 2*khzprime);
  eq49.push_back((*k)*(*dgq));
  // double max0 = fmax(dotdelta, fmax(((*grho)+(*gpres))/(h0*h)*(0.5*(*dgrho) + (*k)*(*eta)), fmax(h0*h*((*dgrho) - 2*(*k)*(*eta)) - 2*khzprime, (*k)*(*dgq))));
  double max0 = maxVec(eq49);
  // eq[0] = (dotdelta + ((*grho)+(*gpres))/(h0*h)*(0.5*(*dgrho) + (*k)*(*eta)) + (h0*h*((*dgrho) - 2*(*k)*(*eta)) - 2*khzprime) + (*k)*(*dgq))/max0;
  eq[0] = (dotdelta + ((*grho)+(*gpres))/(h0*h)*(0.5*(*dgrho) + (*k)*(*eta)) + (h0*h*((*dgrho) - 2*(*k)*(*eta)) - 2*khzprime) + (*k)*(*dgq))/dotdelta;
  eq[1] = (dotdelta + ((*grho)+(*gpres))/(h0*h)*(0.5*(*dgrho) + (*k)*(*eta)) + (h0*h*((*dgrho) - 2*(*k)*(*eta)) - 2*khzprime) + (*k)*(*dgq))/((*grho)+(*gpres))/(h0*h)*(0.5*(*dgrho) + (*k)*(*eta));
  eq[2] = (dotdelta + ((*grho)+(*gpres))/(h0*h)*(0.5*(*dgrho) + (*k)*(*eta)) + (h0*h*((*dgrho) - 2*(*k)*(*eta)) - 2*khzprime) + (*k)*(*dgq))/((h0*h*((*dgrho) - 2*(*k)*(*eta)) - 2*khzprime) + (*k)*(*dgq));

  // Equation (50) of arXiv:1208.0600  
  std::vector<double> eq50;
  eq50.push_back(dotdeltaq);
  eq50.push_back(4*h0*h*(*dgq));
  eq50.push_back(2./3.*(*k)*(khzprime/(h0*h) + (*k)*(*eta) + (*dgrho)));
  eq50.push_back(2./3.*(*k)*(*dgpi));
  // double max1 = fmax(dotdeltaq, fmax(4*h0*h*(*dgq), fmax(2./3.*(*k)*(khzprime/(h0*h) + (*k)*(*eta) + (*dgrho)), 2./3.*(*k)*(*dgpi))));
  double max1 = maxVec(eq50);
  // eq[1] = (dotdeltaq + 4*h0*h*(*dgq) + 2./3.*(*k)*(khzprime/(h0*h) + (*k)*(*eta) + (*dgrho)) + 2./3.*(*k)*(*dgpi))/max1;
  eq[3] = (dotdeltaq + 4*h0*h*(*dgq) + 2./3.*(*k)*(khzprime/(h0*h) + (*k)*(*eta) + (*dgrho)) + 2./3.*(*k)*(*dgpi))/dotdeltaq;
  eq[4] = (dotdeltaq + 4*h0*h*(*dgq) + 2./3.*(*k)*(khzprime/(h0*h) + (*k)*(*eta) + (*dgrho)) + 2./3.*(*k)*(*dgpi))/(4*h0*h*(*dgq));
  eq[5] = (dotdeltaq + 4*h0*h*(*dgq) + 2./3.*(*k)*(khzprime/(h0*h) + (*k)*(*eta) + (*dgrho)) + 2./3.*(*k)*(*dgpi))/(2./3.*(*k)*(khzprime/(h0*h) + (*k)*(*eta) + (*dgrho)));
  eq[6] = (dotdeltaq + 4*h0*h*(*dgq) + 2./3.*(*k)*(khzprime/(h0*h) + (*k)*(*eta) + (*dgrho)) + 2./3.*(*k)*(*dgpi))/(2./3.*(*k)*(*dgpi));


  // FILE* g = fopen("crosschecks/crosschecks_galileon1_eq49ratios_q1.dat", "a");
  // fprintf(g, "%.16f \t %.16f \t %.16f\n", acurr, eq[0], eq[1]);
  // fprintf(g, "%.16f \t %.16f \t %.16f \t %.16f \t %.16f\n", acurr, dotdelta, ((*grho)+(*gpres))/(h0*h)*(0.5*(*dgrho) + (*k)*(*eta)), (h0*h*((*dgrho) - 2*(*k)*(*eta)) - 2*khzprime) + (*k)*(*dgq), eq[0]);
  // fclose(g);

  return eq;

}

extern "C" void freegal_(){

  fflush(stdout);

  gsl_spline_free(spline_h);
  gsl_spline_free(spline_x);
  gsl_interp_accel_free(acc);

  // printf("before -> size a : %i, size hubble : %i, size x : %i\n", intvar.size(), hubble.size(), x.size());

  intvar = std::vector<double>();
  hubble = std::vector<double>();
  x = std::vector<double>();

  // printf("after -> size a : %i, size hubble : %i, size x : %i\n", intvar.size(), hubble.size(), x.size());

}

int test(){

  fflush(stdout);

  // orad = 8.063127541638e-5;
  // orad = 8e-5;
  // orad = 2.469e-5/(0.736*0.736)*(1+0.2271*3.046);
  orad = 0;

  // // Scenario 1
  // om = 0.279;
  // h0 = 60.*1000/(2.99792458e8);
  // c2 = -4.278;
  // c3 = -1.580;
  // c4 = -0.772;
  // cG = 0;

  // Scenario 2
  om = 0.275;
  h0 = 73.6*1000/(2.99792458e8);
  c2 = -4.1;
  c3 = -1.5;
  c4 = -0.78;
  cG = 0.;

  // // Scenario 3
  // om = 0.280;
  // h0 = 72.7*1000/(2.99792458e8);
  // c2 = -3.4;
  // c3 = -1.1;
  // c4 = -0.61;
  // cG = 0.15;

  // // Scenario 4
  // om = 0.275;
  // h0 = 72.7*1000/(2.99792458e8);
  // c2 = -4.1;
  // c3 = -3.375;
  // c4 = -0.775;
  // cG = 0.;

  double grhormass[] = {0,0,0,0,0};
  double nu_masses[] = {0,0,0,0,0};
  int nu_mass_eigenstates = 0;
  double a1 = 0.1;
  double a2 = 0.01;
  double a3 = 0.001;
  double a4 = 0.0001;
  double a5 = 0.00001;

  arrays_(&orad, &om, &h0, &c2, &c3, &c4, &cG, grhormass, nu_masses, &nu_mass_eigenstates);

  double* hx1 = handxofa_(&a1);
  const double y1[] = {(*hx1), *(hx1+1), 0};

  double* hx2 = handxofa_(&a2);
  const double y2[] = {(*hx2), *(hx2+1), 0};

  double* hx3 = handxofa_(&a3);
  const double y3[] = {(*hx3), *(hx3+1), 0};

  double* hx4 = handxofa_(&a4);
  const double y4[] = {(*hx4), *(hx4+1), 0};

  double* hx5 = handxofa_(&a5);
  const double y5[] = {(*hx5), *(hx5+1), 0};

  calcPertOmC2C3C4C5CGC0(a1, y1);
  calcPertOmC2C3C4C5CGC0(a2, y2);
  calcPertOmC2C3C4C5CGC0(a3, y3);
  calcPertOmC2C3C4C5CGC0(a4, y4);
  calcPertOmC2C3C4C5CGC0(a5, y5);

  freegal_();

  return 0;

}
