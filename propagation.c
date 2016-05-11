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

#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <stdbool.h>
#include <string.h>
#include <assert.h>

#ifndef NO_GSL
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#endif

#define M_PI 3.14159265358979323846

#define min(x,y) ((x < y) ? x : y)
#define max(x,y) ((x > y) ? x : y)


// Fast matrix-vector multiplication when the matrix is symmetric 5-diagonal
// Assumes the dimension is at least 4x4 !
static void fast2band(const size_t N, const double diag[N], const double band1[N-1], const double band2[N-2], double complex psi_out[N], const double complex psi[N]) {

     size_t i;

     assert(N>3);

     psi_out[0] = diag[0]*psi[0] + band1[0]*psi[1] + band2[0]*psi[2];

     psi_out[1] = diag[1]*psi[1] + band1[0]*psi[0] + band1[1]*psi[2] + band2[0]*psi[3];

     for (i = 2; i < N-2; i++) {
          psi_out[i] = diag[i]*psi[i];
          psi_out[i] += band1[i]*psi[i+1];
          psi_out[i] += band1[i-1]*psi[i-1];
          psi_out[i] += band2[i]*psi[i+2];
          psi_out[i] += band2[i-2]*psi[i-2];
     }

     psi_out[N-2] = diag[N-2]*psi[N-2] + band2[N-4]*psi[N-4] + \
                band1[N-3]*psi[N-3] + band1[N-2]*psi[N-1];
     psi_out[N-1] = diag[N-1]*psi[N-1] + band2[N-3]*psi[N-3];
}

static void matvec(size_t N, const double mat[N][N], double complex vec[N]) {

     double complex out[N];
     size_t i,j;

     for (i = 0; i < N; i++) {
          out[i] = 0;
          for (j = 0; j < N; j++) {
               out[i] += mat[i][j]*vec[j];
          }
     }
     memcpy(vec,out,N*sizeof(double complex));
}

static void matvec_nc(size_t N, const double mat[N][N], double complex out[N], const double complex vec[N]) {

     size_t i,j;

     memset(out,0,N*sizeof(double complex));
     for (i = 0; i < N; i++) {
          for (j = 0; j < N; j++) {
               out[i] += mat[i][j]*vec[j];
          }
     }
}

static void matTvec(size_t N, const double mat[N][N], double complex vec[N]) {

     double complex out[N];
     size_t i,j;
     complex double tmp;

     memset(out,0,N*sizeof(double complex));
     for (j = 0; j < N; j++) {
          tmp = vec[j];
          for (i = 0; i < N; i++) {
               out[i] += mat[j][i]*tmp;
          }
     }

     memcpy(vec,out,N*sizeof(double complex));
}

static void vecvec(size_t N, double complex v1[N], const double complex v2[N]) {
     size_t i;
     for (i = 0; i < N; i++) {
          v1[i] = v1[i]*v2[i];
     }
}
/*
static double complex dot(size_t N, const double complex a[N], const double complex b[N]) {
     complex double res;
     size_t i;

     res = 0;
     for (i = 0; i < N; i++) {
          res += conj(a[i])*b[i];
     }
     return res;
}
static double norm(size_t N, const double complex a[N]) {
     return sqrt(creal(dot(N,a,a)));
}
*/


void fieldfree_propagation(size_t Jmax, const double complex psi_0[Jmax+1], double t0, size_t num_times, const double times[num_times], const double E_rot[Jmax+1], const double Udiag[Jmax+1], const double Uband1[Jmax], const double Uband2[Jmax-1], double complex psi_result[num_times][Jmax+1], double cos2[num_times], _Bool do_cos2d, const double U2d[Jmax+1][Jmax+1], double cos2d[num_times]) {


     size_t i;
     double t, Dt;
     double complex phase, U_psi[Jmax+1], U2d_psi[Jmax+1];
     double complex phase_incr, phase_incr_base;
     double Bconst;
     size_t J;

     Bconst = (E_rot[2]-E_rot[1])/4; // also works for symmetric tops!

     for (i = 0; i < num_times; i++) {
     
          t = times[i];
          Dt = t-t0;

          /*
          for (J = 0; J <= Jmax; J++) { // Propagate wave function
               phase = cexp(-E_rot[J]*Dt*I);
               psi_result[i][J] = psi_0[J]*phase;
          }*/

          // The same as the commented out code,
          // but exploits the Bj(j+1) structure of the energy levels.
          // This avoids the very costly evaluation of cexp()
          // Note: It won't work if you want to implement centrifugal
          // distorsion.
          
          phase_incr_base = cexp(remainder(-2*Bconst*Dt,2*M_PI)*I);
          phase_incr = 1;
          phase = cexp(-E_rot[0]*Dt*I);
          psi_result[i][0] = psi_0[0]*phase;
          for (J = 1; J <= Jmax; J++) {
               phase_incr *= phase_incr_base;
               phase *= phase_incr;
               psi_result[i][J] = psi_0[J]*phase;
          }
 
          //U_psi = U dot Psi
          fast2band(Jmax+1, Udiag,Uband1,Uband2, U_psi, psi_result[i]);

          cos2[i] = 0;
          for (J = 0; J <= Jmax; J++) { // calculate cos^2 = conj(psi)*U*psi
               cos2[i] += conj(psi_result[i][J])*U_psi[J];
          }
          
          // U2d_psi = U2d dot psi, if U2d is specified
          if (do_cos2d) {
               // U2d_psi = U2d dot psi
               matvec_nc(Jmax+1, U2d, U2d_psi, psi_result[i]);
               cos2d[i] = 0;
               for (J = 0; J <= Jmax; J++)
                    cos2d[i] += conj(psi_result[i][J])*U2d_psi[J];
          }

     }
}

// Like the python version, except without all the implicit memory allocations
// and deallocations
int propagate_field(size_t Nsteps, size_t Nsteps_inner, size_t dim, double t, double dt, double E0_sq_max, double sigma, double complex psi_t[Nsteps][dim], const double eig[dim], const double eigvec[dim][dim], const double complex expRot[dim], const double complex expRot2[dim]) {
     

     size_t i, j, k;
     double fac;
     double sigma_sq = sigma*sigma;
     
     for (i = 1; i < Nsteps; i++) {
          for (j = 0; j < dim; j++)
               psi_t[i][j] = expRot2[j]*psi_t[i-1][j];
          for (k = 0; k < Nsteps_inner; k++) {
               if (k > 0)
                    vecvec(dim,psi_t[i],expRot); // psi = expRot*psi;
               matTvec(dim,eigvec,psi_t[i]); // psi = eigvec^T dot psi
               fac = -E0_sq_max*exp(-(t*t)/(2*sigma_sq))*dt;
               t = t + dt;
               for (j = 0; j < dim; j++) // psi = psi*exp(the diagonal)
                    psi_t[i][j] *= cexp(fac*eig[j]*I);
               
               matvec(dim,eigvec,psi_t[i]); // psi = eigvec dot psi

               if (fmax(pow(cabs(psi_t[i][dim-2]),2),pow(cabs(psi_t[i][dim-1]),2)) > 1e-5)
                    return 1; // Basis size too small
          }
          vecvec(dim,psi_t[i],expRot2); // psi = expRot2*expRot2
     }

     return 0;
}


#ifndef NO_GSL
struct deriv_params {
     size_t jmax;
     double peak_field_amplitude_squared, sigma;
     const double *E_rot;
     const double *V0, *V1, *V2;
     int ncalls;
};


int deriv (double t, const double psi[], double dPsidt[], void * params) {
     // Note: psi and dPsidt are complex, but GSL ode interface requires double.
     // Index 2j is the j'th real component, and 2j+1 is the j'th imaginary
     // component. This is a somewhat dangerous assumption that
     // can depend on the compiler implementation. It holds for GCC.
 
     struct deriv_params *p = params;
     size_t j;
     const size_t jmax = p->jmax;
     // Gaussian pulse:
     double E_0_squared = p->peak_field_amplitude_squared * \
                          exp(-t*t/(2*(p->sigma)*(p->sigma)));
     double tmp;

     
     // Multiply the interaction term
     fast2band(jmax+1,p->V0,p->V1,p->V2,(double complex *) dPsidt, (const double complex *) psi);
     // Then add the diagonal and multiply every component with -i.
     // I.e. -i*(a+ib) = b-ia.
     for (j = 0; j < jmax; j++) {
          dPsidt[2*j] *= E_0_squared;
          dPsidt[2*j+1] *= E_0_squared;
          dPsidt[2*j] += p->E_rot[j]*psi[2*j];
          dPsidt[2*j+1] += p->E_rot[j]*psi[2*j+1];

          // Multiplication by -i:
          tmp = dPsidt[2*j+1]; // Imaginary part
          dPsidt[2*j+1] = -dPsidt[2*j];
          dPsidt[2*j]  = tmp; // becomes real part.
     }

     //p->ncalls++;

     return GSL_SUCCESS;
}


int propagate_field_ODE(size_t Nsteps, const size_t dim, double t, double dt, double E0_sq_max, double sigma, double complex psi_t[Nsteps][dim], const double V0[], const double V1[], const double V2[], const double E_rot[], double abstol, double reltol) {

     size_t i;
     struct deriv_params p = {dim-1,E0_sq_max,sigma,E_rot,V0,V1,V2,0};
     gsl_odeiv2_system sys = {.function=deriv,.jacobian=NULL,\
          .dimension=2*dim,.params=&p};

     double initial_step_size = 2.4*sigma/150; // 150 steps per pulse
     gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new (&sys, \
               gsl_odeiv2_step_rk8pd, initial_step_size, abstol,reltol);
               // rk8pd seems fastest. Did not try the methods
               // that rely on the Jacobian. d(dpsi/dt)/dpsi is diagonal,
               // but d(dpsi/dt)/dt is less trivial.
               //gsl_odeiv2_step_rkck, initial_step_size, abstol,reltol);
               //gsl_odeiv2_step_msadams, initial_step_size, abstol,reltol);
               //gsl_odeiv2_step_rk2, initial_step_size, abstol,reltol);
               //gsl_odeiv2_step_rkf45, initial_step_size, abstol,reltol);

     double t_run = t;

     for (i = 1; i < Nsteps; i++) {
     
          memcpy(psi_t[i],psi_t[i-1],sizeof(double complex)*dim);
          gsl_odeiv2_driver_apply(d,&t_run,t+(double)i*dt,(double *) psi_t[i]);
     
          // In the truncated basis, the ODE equation interaction term
          // never transfers population to the lase expansion coefficient.
          // For K=0 (eg linear molecules) also it does not allow population
          // of the last two levels. Therefore, we check the levels
          // before that.
          if (fmax(pow(cabs(psi_t[i][dim-4]),2),pow(cabs(psi_t[i][dim-3]),2)) > 1e-5)
               return 1; // Basis size too small

     }

     gsl_odeiv2_driver_free(d);

     //printf("ncalls: %i.\n",p.ncalls);

     return 0;
}

#endif


/* Experimentation with krylov subspace method. 
// Note: Clutters psi
size_t lanczos_banded(size_t N, size_t m, double complex psi[N], complex double V[m+1][N],  double diag[N], double offdiag[N], double Adiag[N], double band1[N-1], double band2[N-2]) {

     size_t i,j;

     // First vector of V is psi:
     memcpy(V[0],psi,N*sizeof(complex double));

     for (i = 0; i < m; i++) {
          
          // Multiply vector with the matrix
          fast2band(N, Adiag,band1,band2,psi);

          if (i > 0) {
               // Project out previous vector
               for (j = 0; j < N ; j++) 
                    psi[j] -= offdiag[i-1]*V[i-1][j];
          }
          
          // Store overlap with "next previous" vector
          diag[i] = creal(dot(N,V[i],psi));
          // And then project it out
          for (j = 0; j < N ; j++)
               psi[j] -= diag[i]*V[i][j];

          offdiag[i] = norm(N,psi); // Store norm for later reconstruction

          if (offdiag[i] < 1e-6) { // Happy breakdown
               m = i+1;
               memset(V[i+1],0,N*sizeof(complex double));
               break;
          }

          // Renormalize 
          for (j = 0; j < N; j++)
               psi[j] /= offdiag[i];
          
          memcpy(V[i+1],psi,N*sizeof(complex double));

     }

     return m;
}

*/
