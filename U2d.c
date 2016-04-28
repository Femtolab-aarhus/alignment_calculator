/*
   Copyright 2016 Anders Aspegren Søndergaard / Femtolab, Aarhus University

   This file is part of Alignment calculator.

   Alignment calculator is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   Alignment calculator is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with Alignment calculator. If not, see <http://www.gnu.org/licenses/>.
 */


#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <assert.h>

#include <math.h>

#define M_PI 3.14159265358979323846

#define max(a,b) (a > b) ? a : b
#define min(a,b) (a < b) ? a : b

// This library can be compiled without any GSL dependences.
// If you do this, define NO_GSL. The expand_U2d function
// for calculating expansion coefficients will not be exported
// in this case.

#ifndef NO_GSL
#include <gsl/gsl_sf.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_integration.h>


static double integrand(double phi, double u) {
     const double usq = u*u;
     double sin2phi;

     if (u == 0) return 0;

     sin2phi = sin(phi);
     sin2phi = sin2phi*sin2phi;
     return usq/(usq+(1-usq)*sin2phi);
}

struct integrand_params {
     double u;
};


static double integrand_wrap(double psi, void *params) {
     struct integrand_params *p = params;
     return integrand(psi,p->u);
}

struct inner_params {
     gsl_function integrand;
     gsl_integration_workspace *workspace;
     int l;
     size_t limit;
};

static double inner(double u, void *params) {
  
     struct inner_params *p = params;
     struct integrand_params *ip = p->integrand.params;
     const double abstol = 1e-13, reltol = 1e-13; 
     double result,abserr;

     ip->u = u;
     gsl_integration_qag(&(p->integrand),0.0,2*M_PI,abstol,reltol,p->limit,GSL_INTEG_GAUSS61,p->workspace,&result,&abserr);

     return result*gsl_sf_legendre_sphPlm(p->l, 0, u);
}

static double outer(int l, gsl_integration_workspace *ws_in, gsl_integration_workspace *ws_out, size_t limit) {

     double res, abserr;
     const double abstol = 1e-13, reltol = 1e-13; 
     gsl_function inner_F, F;
     
     inner_F.function = &integrand_wrap;
     struct integrand_params inner_params = {.u=0};
     inner_F.params = &inner_params;

     F.function = &inner;
     struct inner_params params = {.integrand=inner_F,.workspace=ws_in,.l=l,.limit=limit};
     F.params = &params;

     gsl_integration_qag(&F, -1.0, 1.0, abstol, reltol, limit, GSL_INTEG_GAUSS61, ws_out, &res, &abserr);

     return res;

}

// Calculate relevant expansion coefficients by numerical integration
void expand_U2d(int lmax, double coeff[lmax+1]) {

     const size_t limit = 409600;
     int l;
     
     memset(coeff,0,((size_t)lmax+1)*sizeof(double));
     gsl_integration_workspace *inner_workspace,*workspace;

     printf("Expanding cos^2 theta 2d in Legendre polynomials...\n");

     #pragma omp parallel default(shared) private(l,inner_workspace,workspace)
     {
     inner_workspace = gsl_integration_workspace_alloc(limit);
     workspace = gsl_integration_workspace_alloc(limit);
     #pragma omp for schedule(guided,4)
     for (l=0; l<=lmax; l+=2) {
          coeff[l] = outer(l,inner_workspace,workspace,limit);
          //printf("%i: %e\n",l,coeff[l]);
          printf("l = %i, ",l);
          fflush(stdout);
     }
     gsl_integration_workspace_free(workspace);
     gsl_integration_workspace_free(inner_workspace);
     }
     printf("\n");

}

// Export Wigner 3j symbol (using GSL)
double three(int j1, int j2, int j3, int m1, int m2, int m3) {
     return gsl_sf_coupling_3j(2*j1,2*j2,2*j3,2*m1,2*m2,2*m3);
}
// Same but with m=0
double three0(int j1, int j2, int j3) {
     // Apparently we loose some speed here. There is a faster
     // way to calculate (j1,j2,j3,0,0,0) than gsl does.
     return three(j1,j2,j3,0,0,0);
}

#endif 


// Prototype for the drc3jj netlib function.
int drc3jj(int l2, int l3, int m2, int m3, int *l1min, int *l1max, double *thrcof, int ndim);




// Calculate the matrix elements. Coeff are the expansion coefficients.
void populate_U2d(const int lmax, const int k, const int m, const int lppmax, double coeff[lppmax+1], double U2d[lmax+1][lmax+1]) {

     int l,lp,lpp;
     double w3jk[2*lmax+1];
     double w3jm[2*lmax+1];
     double fac;
     int lppmin, lppmax_;

     memset(U2d,0,sizeof(double)*((size_t)lmax+1)*((size_t)lmax+1));

     #pragma omp parallel default(shared) private(l,lp,lpp,fac,w3jk,w3jm,lppmin,lppmax_)
     {
     #pragma omp for schedule(guided,4)
     for (l = max(abs(m),abs(k)); l <= lmax; l++) {
     for (lp = l; lp <= lmax; lp += ((k==0 || m==0) ? 2 : 1)) {
     
     // Evaluate all 3j symbols in one go using recursion formula
     assert(drc3jj(l,lp,k,-k,&lppmin,&lppmax_,w3jk,2*lmax+1)==0);
     assert(lppmin==abs(l-lp));
     assert(lppmax_==l+lp);
     assert(drc3jj(l,lp,m,-m,&lppmin,&lppmax_,w3jm,2*lmax+1)==0);
     assert(lppmin==abs(l-lp));
     assert(lppmax_==l+lp);

     lppmax_ = min(lppmax_,lppmax);
          
     U2d[l][lp] = 0;
     for (lpp = abs(l-lp); lpp <= lppmax_; lpp++) {
          fac = sqrt(2.0*lpp+1);
          U2d[l][lp] += coeff[lpp]*fac*w3jk[lpp-abs(l-lp)]*w3jm[lpp-abs(l-lp)];
     }

     U2d[l][lp] *= sqrt((2.0*((double) l)+1)*(2.0*((double) lp)+1)/(4*M_PI));
     if ((m-k)&1) U2d[l][lp] = -U2d[l][lp]; // odd m-k: add phase factor.

     U2d[lp][l] = U2d[l][lp];
     }
     }
     }

}

