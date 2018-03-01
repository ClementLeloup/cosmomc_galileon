#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>

#include <vector>
#include <iostream>


// External functions from fortran
extern"C" void massivenu_mp_nu_rho_(double* am, double* rhonu);
extern"C" void massivenu_mp_nu_background_(double* am, double* rhonu, double* pnu);

int calcPertOmC2C3C4C5CG(double a, const double y[3], void* params) {

  fflush(stdout);

 const double *in_params = static_cast<const double*>(params);
  // bool store_derivs = (int)in_params[0];
  double grhormass[5];
  for(int i = 0; i<5; i++){grhormass[i] = in_params[i];}
  double nu_masses[5];
  for(int i = 0; i<5; i++){nu_masses[i] = in_params[5+i];}
  int nu_mass_eigenstates = in_params[10];
  double c2 = in_params[11];
  double c3 = in_params[12];
  double c4 = in_params[13];
  double c5 = in_params[14];
  double cG = in_params[15];
  double orad = in_params[16];
  double h0 = in_params[17];

  double h0_Mpc = h0*1000/(2.99792458e8);

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
  double alpha = c2/6.0*prod - 3*c3*h*prod2 + 15*c4*h2*prod3 - 17.5*c5*h3*prod4 - 3.0*cG*h2*prod;
  double gamma = c2/3.0*h*prod - c3*h2*prod2 + 2.5*c5*h4*prod4 - 2.0*cG*h3*prod;
  double beta = -2*c3*h3*prod + c2/6.0*h2 + 9*c4*h4*prod2 - 10*c5*h4*h*prod3 - cG*h4;
  double sigma = 2.0*h + 2.0*c3*prod3 - 15.0*c4*h*prod4 + 21.0*c5*h2*prod5 + 6.0*cG*h*prod2;
  double lambda = 3.0*h2 + orad/(a*a*a*a) - 2.0*c3*h*prod3 + c2/2.0*prod2 + 7.5*c4*h2*prod4 - 9.0*c5*h3*prod5 - cG*h2*prod2;
  double omega = 2*c3*h2*prod2 - 12*c4*h3*prod3 + 15*c5*h4*prod4 + 4.0*cG*h3*prod;
  double denom = sigma*beta - alpha*omega;

  // if(a == 1.) printf("theta_ a : %.16f \t alpha : %.16f \t gamma : %.16f \t beta : %.16f \t sigma : %.16f \t lambda : %.16f \t omega : %.16f \t h : %.16f \t x : %.16f\n", a, alpha, gamma, beta, sigma, lambda, omega, h, xgal);

  // Contribution from massive neutrinos
  if(nu_mass_eigenstates>0){
    for(int i = 0; i<nu_mass_eigenstates; i++){
      double rhonu = 0;
      double pnu = 0;
      double am = a*nu_masses[i];
      massivenu_mp_nu_background_(&am, &rhonu, &pnu);
      lambda += grhormass[i]*pnu/(a*a*a*a*h0_Mpc*h0_Mpc);
    }
  }

  // printf("theta_ a : %.16f \t alpha : %.16f \t gamma : %.16f \t beta : %.16f \t sigma : %.16f \t lambda : %.16f \t omega : %.16f \t h : %.16f \t x : %.16f\n", a, alpha, gamma, beta, sigma, lambda, omega, h, xgal);

  double dhdlna = (omega*gamma - lambda*beta)/denom;
  double dxdlna = (alpha*lambda - sigma*gamma)/denom - xgal;

  double k1 = -6*c4*h*prod2*(dhdlna*xgal + h*dxdlna + prod/3.0) + c5*h2*prod3*(12*h*dxdlna + 15*dhdlna*xgal + 3*prod) + 2.0*cG*(dhdlna*prod+h2*dxdlna+h*prod);
  double k2 = -0.5*c2 + 6*c3*h*prod - 27*c4*h2*prod2 + 30*c5*h3*prod3 + 3.0*cG*h2;
  double k3 = -1.0 - 0.5*c4*prod4 - 3*c5*h*prod4*(h*dxdlna + dhdlna*xgal) + cG*prod2;
  double k4 = -2.0 + 3.0*c4*prod4 - 6*c5*h*prod5 - 2.0*cG*prod2;
  double k5 = 2*c3*prod2 - 12*c4*h*prod3 + 15.0*c5*h2*prod4 + 4.0*cG*h*prod;
  double k6 = 0.5*c2 - 2.0*c3*( h2*dxdlna + dhdlna*prod + 2*h*prod ) + c4*( 12*h3*prod*dxdlna + 18*h*prod2*dhdlna + 13*h2*prod2 ) - c5*( 18*h2*h2*prod2*dxdlna + 30*h2*prod3*dhdlna + 12*h3*prod3 ) - cG*( 2.0*h*dhdlna + 3.0*h2 );

  double noghost = k2 + 1.5*k5*k5/k4;
  double cs2 = (4*k1*k4*k5 - 2*k3*k5*k5 - 2*k4*k4*k6)/(k4*(2*k4*k2 + 3*k5*k5));
  // Here we can't have a condition noghost>0 because with some set of
  // parameters we have noghost=1e-16, condition passes anyway
  // and then rk4 fails and gives NaN values.

  // Check instabilities for tensorial modes (from Felice & Tsujikawa)
  double noghostt = 0.5 - 0.75*c4*prod4 + 1.5*c5*prod5*h + 0.5*cG*prod2;
  double ct2 = (0.5 + 0.25*c4*prod4 + 1.5*c5*prod4*h*(h*dxdlna+dhdlna*xgal) - 0.5*cG*prod2) / (0.5 - 0.75*c4*prod4 + 1.5*c5*prod5*h + 0.5*cG*prod2);;

  // // To test theoretical constraints
  // printf("a = %.12f \t Qs = %.12f \t cs2 = %.12f \t Qt = %.12f \t cT2 = %f \t OmegaPi = %.12f\n om = %.12f \t orad = %.12f \t c2 = %.12f \t c3 = %.12f \t c4 = %.12f \t cG = %.12f\n", a, noghost, cs2, noghostt, ct2, OmegaP, om, orad, c2, c3, c4, cG);

  // printf("theta_ a : %.16f \t Qs : %.16f \t cs2 : %.16f \t QT : %.16f \t cT2 : %.16f \t h : %.16f \t x : %.16f\n", a, noghost, cs2, noghostt, ct2, h, xgal);
  // printf("arrays_ a : %.16f \t dhdlna : %.16f \t dxdlna : %.16f \t k1 : %.16f \t k2 : %.16f \t k3 : %.16f \t k4 : %.16f \t k5 : %.16f \t k6 : %.16f\n", a, dhdlna, dxdlna, k1, k2, k3, k4, k5, k6);

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
    fprintf(stderr, "Error : Laplace stability condition not verified : a = %.12f \t cT2 = %.12f\n", a, ct2);
    return 4;
  }

  return 0;

}

int calcValOmC2C3C4C5CG(double a, const double y[3], double f[3], void* params){

  const double *in_params = static_cast<const double*>(params);
  // bool store_derivs = (int)in_params[0];
  double grhormass[5];
  for(int i = 0; i<5; i++){grhormass[i] = in_params[i];}
  double nu_masses[5];
  for(int i = 0; i<5; i++){nu_masses[i] = in_params[5+i];}
  int nu_mass_eigenstates = in_params[10];
  double c2 = in_params[11];
  double c3 = in_params[12];
  double c4 = in_params[13];
  double c5 = in_params[14];
  double cG = in_params[15];
  double orad = in_params[16];
  double h0 = in_params[17];

  double h0_Mpc = h0*1000/(2.99792458e8);

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
  double alpha = c2/6.0*prod - 3*c3*h*prod2 + 15*c4*h2*prod3 - 17.5*c5*h3*prod4 - 3.0*cG*h2*prod;
  double gamma = c2/3.0*h*prod - c3*h2*prod2 + 2.5*c5*h4*prod4 - 2.0*cG*h3*prod;
  double beta = c2/6.0*h2 - 2*c3*h3*prod + 9*c4*h4*prod2 - 10*c5*h4*h*prod3 - cG*h4;
  double sigma = 2.0*h + 2.0*c3*prod3 - 15.0*c4*h*prod4 + 21.0*c5*h2*prod5 + 6.0*cG*h*prod2;
  double lambda = 3.0*h2 + orad/(a2*a2) + c2/2.0*prod2 - 2.0*c3*h*prod3 + 7.5*c4*h2*prod4 - 9.0*c5*h3*prod5 - cG*h2*prod2;

  // Contribution from massive neutrinos
  if(nu_mass_eigenstates>0){
    for(int i=0; i<nu_mass_eigenstates; i++){
      double rhonu = 0;
      double pnu = 0;
      double am = a*nu_masses[i];
      massivenu_mp_nu_background_(&am, &rhonu, &pnu);
      lambda += grhormass[i]*pnu/(a2*a2*h0_Mpc*h0_Mpc);
      // printf("i : %d \t a : %.16f \t am : %.16f \t rhonu : %.16f \t lambda_numass : %.16f\n", i, a, am, rhonu, grhormass[i]*pnu/(a2*a2*h0*h0));
    }
  }

  double omega = 2*c3*h2*prod2 - 12*c4*h3*prod3 + 15*c5*h4*prod4 + 4.0*cG*h3*prod;

  double denom = sigma*beta - alpha*omega;
  f[0] = (omega*gamma - lambda*beta) / (a*denom);
  f[1] = (alpha*lambda - sigma*gamma) / (a2*denom) - 2.0*y[1]/a;
  f[2] = y[1];

  return GSL_SUCCESS;
}

int calcHubble(double* intvar, double* hbar, double* grhormass, double* nu_masses, int* nu_mass_eigenstates, double om, double orad, double c2, double c3, double c4, double c5, double cG, double h0, int nb){

  double params[18]; // parameters for massive neutrinos
  for(int i = 0; i<5; i++){ params[i] = grhormass[i]; }
  for(int i = 0; i<5; i++){ params[5+i] = nu_masses[i]; }
  params[10] = (*nu_mass_eigenstates);
  params[11] = c2;
  params[12] = c3;
  params[13] = c4;
  params[14] = c5;
  params[15] = cG;
  params[16] = orad;
  params[17] = h0;

  double h0_Mpc = h0*1000/(2.99792458e8);

  const gsl_odeiv_step_type * T = gsl_odeiv_step_rkf45;
  gsl_odeiv_step * s = gsl_odeiv_step_alloc(T, 3);
  gsl_odeiv_control * c = gsl_odeiv_control_y_new(1e-18, 1e-16);
  gsl_odeiv_evolve * e = gsl_odeiv_evolve_alloc(3);
  gsl_odeiv_system sys;
  sys.function = calcValOmC2C3C4C5CG;
  sys.dimension = 3;
  sys.params = &params;
 
  double y[3] = { 1.0, 1.0, 0.0 };  //Inital value of integral
  double a = 1;
  // int nstep = min(hubble.size(),xgalileon.size());
  int nstep = nb+1;
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
      // fprintf(stderr, "\nFailure with om = %f, c2 = %f, c3 = %f, c4 = %f, c5 = %f, cG = %f and c0 = %f at a = %f, h = %f, x = %f and y = %f\n", om, c2, c3, c4, c5, cG, c0, acurrtarg, y[0], y[1], y[2]);
      return 8;
    }

    // testPert = calcPertOmC2C3C4C5CGC0(intvar[i], y);
    int testPert = calcPertOmC2C3C4C5CG(a, y, params);
    if(testPert != 0){
      gsl_odeiv_evolve_free(e);
      gsl_odeiv_control_free(c);
      gsl_odeiv_step_free(s);
      return 7;
    }
    
    // Variables to optimize computation
    double h2 = y[0]*y[0];
    double x1 = a*y[1];
    double x2 = x1*x1;
    double x3 = x2*x1;
    double a3 = a*a*a;

    // Friedmann equation for background
    double OmegaP = (0.5*c2*x2 - 6.0*c3*h2*x3 + 22.5*c4*h2*h2*x2*x2 - 21.0*c5*h2*h2*h2*x3*x2 - 9.0*cG*h2*x2)/3.0;
    // OmegaM = 1 - OmegaP - orad/(a3*intvar[i]*h2);
    double OmegaM = 1 - OmegaP - orad/(a3*a*h2);

    // Contribution from massive neutrinos
    if((*nu_mass_eigenstates)>0){
      for(int j=0; j<(*nu_mass_eigenstates); j++){
    	double rhonu = 0;
	double pnu = 0;
    	double am = intvar[i]*nu_masses[j];
    	massivenu_mp_nu_background_(&am, &rhonu, &pnu);
    	OmegaM -= grhormass[j]*rhonu/(a3*a*3.*h2*h0_Mpc*h0_Mpc);
	// printf("i : nu_masses = %.16f\tgrhormass = %.16f\n", nu_masses[j], grhormass[j]);
      }
    }

    double OmTest = om/(a3*h2);

    //For test
    // printf("Test on OmegaM : %f %f %f %f %f %f %f %f %f %f %f %f %f %.12f\n", om, orad, c2, c3, c4, c5, cG, OmTest, OmegaM, orad/(a3*a*h2), OmegaP, a3*h2, intvar[i], fabs((OmegaM - OmTest)/OmTest));
    // printf("lambda : %.12f\tusual : %.12f\tfrom massive neutrinos : %.12f\n", lambda, 3.0*h2 + orad/(a3*intvar[i]) - 2.0*c0*h2*x1 - 2.0*c3*h2*h2*x3 + c2/2.0*h2*x2 + 7.5*c4*h2*h2*h2*x2*x2 - 9.0*c5*h2*h2*h2*h2*x3*x2 - cG*h2*h2*x2, grhormass[0]*nuPres(intvar[i]*nu_masses[0])/(a3*intvar[i]*h0*h0));

    if(OmegaP<0){
      // fprintf(stderr, "Negative galileon energy density : a = %.12f \t %.12f\n", a, OmegaP);
      int truc = calcPertOmC2C3C4C5CG(a, y, params);
      gsl_odeiv_evolve_free(e);
      gsl_odeiv_control_free(c);
      gsl_odeiv_step_free(s);
      return 5;
    }

    if ( fabs((OmegaM - OmTest)/OmTest)>1e-4 ) {
      fprintf(stderr, "Integration error : %f %f %f %f %f %f %f %f %f %f %f %.12f %f\n", om, orad, c2, c3, c4, c5, cG, OmTest, OmegaM, OmegaP, OmegaM - 1 + OmegaP + orad/(a3*intvar[i]*h2), a, fabs((OmegaM - OmTest)/OmTest));
      
      gsl_odeiv_evolve_free(e);
      gsl_odeiv_control_free(c);
      gsl_odeiv_step_free(s);
      return 4;
    }

    hbar[i] = y[0];
  }

  gsl_odeiv_evolve_free(e);
  gsl_odeiv_control_free(c);
  gsl_odeiv_step_free(s);

  return 0;

}



// Function to replace Cosmomctheta in galileon case to allow multithreading inversion of theta(H0)
extern "C" double theta_(double* omegam, double* orad, double* ombh2, double* h0, double* astar, double* c2in, double* c3in, double* c4in, double* cGin, double* grhormass, double* nu_masses, int* nu_mass_eigenstates){

  int nb = 29999;

  double om = (*omegam);
  double omegar = (*orad);

  double intvar[nb+1];
  double hbar[nb+1];

  double c2 = (*c2in);
  double c3 = (*c3in);
  double c4 = (*c4in);
  double c5 = 0;
  double cG = (*cGin);

  if((*c3in) == 0 && (*c4in) == 0){
    c3 = 0.5*(-1. + om + omegar + c2/6. - 3*cG);

    // Add massive neutrinos
    double rhonu = 0;
    double grhom = 3*pow((*h0)*1000/2.99792458e8, 2); // critical density
    if((*nu_mass_eigenstates)>0){
      for(int i = 0; i<(*nu_mass_eigenstates); i++){
	massivenu_mp_nu_rho_(&(nu_masses[i]), &rhonu);
	// printf("i : nu_masses = %.18f\tgrhormass = %.18e\trhonu = %.18f\tgrhom = %.18f\ttest1 = %.18f\ttest2 = %.18f\ttest3 = %.18f\tc5 = %.18f\n", nu_masses[i], grhormass[i], rhonu, grhom, om + orad, c2/6., - 2*c3 + 7.5*c4 - 3*cG, c5);
	c3 += 0.5*rhonu*grhormass[i]/(grhom);
	// printf("Omeganu of mass eigenstate %d = %.16f\n", i, rhonu*grhormass[i]/grhom);
      }
    }

    c4 = 0;
    c5 = 0;

  } else if((*c4in) == 0){
    c3 = (*c3in);
    c4 = (1. - om - omegar - c2/6. + 2*c3 - 7.5*c4 + 3*cG)/7.5;

    // Add massive neutrinos
    double rhonu = 0;
    double grhom = 3*pow((*h0)*1000/2.99792458e8, 2); // critical density
    if((*nu_mass_eigenstates)>0){
      for(int i = 0; i<(*nu_mass_eigenstates); i++){
	massivenu_mp_nu_rho_(&(nu_masses[i]), &rhonu);
	// printf("i : nu_masses = %.18f\tgrhormass = %.18e\trhonu = %.18f\tgrhom = %.18f\ttest1 = %.18f\ttest2 = %.18f\ttest3 = %.18f\tc5 = %.18f\n", nu_masses[i], grhormass[i], rhonu, grhom, om + orad, c2/6., - 2*c3 + 7.5*c4 - 3*cG, c5);
	c4 += -rhonu*grhormass[i]/(7.5*grhom);
	// printf("Omeganu of mass eigenstate %d = %.16f\n", i, rhonu*grhormass[i]/grhom);
      }
    }

    c5 = 0;

  } else{
    c3 = (*c3in);
    c4 = (*c4in);
    c5 = (-1. + om + omegar + c2/6. - 2*c3 + 7.5*c4 - 3*cG)/7.;

    // Add massive neutrinos
    double rhonu = 0;
    double grhom = 3*pow((*h0)*1000/2.99792458e8, 2); // critical density
    if((*nu_mass_eigenstates)>0){
      for(int i = 0; i<(*nu_mass_eigenstates); i++){
	massivenu_mp_nu_rho_(&(nu_masses[i]), &rhonu);
	// printf("i : nu_masses = %.18f\tgrhormass = %.18e\trhonu = %.18f\tgrhom = %.18f\ttest1 = %.18f\ttest2 = %.18f\ttest3 = %.18f\tc5 = %.18f\n", nu_masses[i], grhormass[i], rhonu, grhom, om + orad, c2/6., - 2*c3 + 7.5*c4 - 3*cG, c5);
	c5 += rhonu*grhormass[i]/(7.*grhom);
	// printf("Omeganu of mass eigenstate %d = %.16f\n", i, rhonu*grhormass[i]/grhom);
      }
    }
  }

    int status = 0;

  // Where to put initial condition
  // double amin = 9.99999e-7;
  double amin = 1e-10;
  double amax = 1.;

  // Fill the vector of a with a geometric sequence
  // double nb = 30000; // number of a points using q0
  // double nb = 300000; // number of a points using q0
  double q = pow(amin/amax, 1./nb);
  for(int i = 0; i<=nb; i++){
    // intvar.push_back(amax*pow(q, i));
    intvar[i] = amax*pow(q, i);
  }

  // printf("theta_ OmegaM0 = %.18f\nOmegaR0 = %.18f\nc2 = %.18f\nc3 = %.18f\nc4 = %.18f\nc5 = %.18f\ncG = %.18f\nh0 = %.18f km/s/Mpc\n", om, (*orad), c2, c3, c4, c5, cG, (*h0));
  // fflush(stdout);

  status = calcHubble(intvar, hbar, grhormass, nu_masses, nu_mass_eigenstates, om, (*orad), c2, c3, c4, c5, cG, (*h0), nb);

  printf("status = %i\n", status);
  if(status!=0) return -1;

  // Interpolation tools
  gsl_interp_accel* acc;
  gsl_spline* spline_rs;
  gsl_spline* spline_DA;
  gsl_spline* spline_hbar;
  spline_rs = gsl_spline_alloc(gsl_interp_cspline, nb+1);
  spline_DA = gsl_spline_alloc(gsl_interp_cspline, nb+1);
  spline_hbar = gsl_spline_alloc(gsl_interp_cspline, nb+1);

  double intvar_interp[nb+1];
  double rs_interp[nb+1];
  double DA_interp[nb+1];
  double hbar_interp[nb+1];
  for(int i = 0; i<=nb; i++){
    intvar_interp[i] = intvar[nb-i];
    hbar_interp[i] = hbar[nb-i];
    rs_interp[i] = 2.99792458e8/(intvar[nb-i]*intvar[nb-i]*hbar[nb-i]*(*h0)*1000*sqrt(3*(1+3e4*intvar[nb-i]*(*ombh2))));
    DA_interp[i] = 2.99792458e8/(intvar[nb-i]*intvar[nb-i]*hbar[nb-i]*(*h0)*1000);
  }
  gsl_spline_init(spline_rs, intvar_interp, rs_interp, nb+1);
  gsl_spline_init(spline_DA, intvar_interp, DA_interp, nb+1);
  gsl_spline_init(spline_hbar, intvar_interp, hbar_interp, nb+1);

  double rs_int = gsl_spline_eval_integ(spline_rs, 1e-8, (*astar), acc);
  double DA_int = gsl_spline_eval_integ(spline_DA, (*astar), 1.0, acc);

  double theta = rs_int/DA_int;
  
  // printf("theta 2 = %.12f\n", theta);

  gsl_spline_free(spline_rs);
  gsl_spline_free(spline_DA);
  gsl_interp_accel_free(acc);

  return theta;

}
