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

#define ODD(n) ((n)&1) // n (un)signed integer

// Fast matrix-vector multiplication when the matrix is symmetric 5-diagonal
static void fast2band(const size_t N, const double diag[N], const double band1[N-1], const double band2[N-2], double complex out[N], const double complex in[N]) {

     size_t i;

     if (N > 3) {
     out[0] = diag[0]*in[0] + band1[0]*in[1] + band2[0]*in[2];

     out[1] = diag[1]*in[1] + band1[0]*in[0] + band1[1]*in[2] + band2[1]*in[3];
     
     for (i = 2; i < N-2; i++) {
          out[i] = diag[i]*in[i]+band1[i]*in[i+1]+band1[i-1]*in[i-1]+band2[i]*in[i+2]+band2[i-2]*in[i-2];
     }

     out[N-2] = diag[N-2]*in[N-2] + band2[N-4]*in[N-4] + \
                band1[N-3]*in[N-3] + band1[N-2]*in[N-1];
     out[N-1] = diag[N-1]*in[N-1] + band1[N-2]*in[N-2] + band2[N-3]*in[N-3];

     } else { // N <= 3 
          for (i = 0; i < N; i++) out[i] = diag[i]*in[i];
          if (N == 2) {
               out[0] += in[1]*band1[0];
               out[1] += in[0]*band1[0];
          } else if (N == 3) {
               out[0] += in[1]*band1[0]+in[2]*band2[0];
               out[1] += in[0]*band1[0]+in[2]*band1[1];
               out[2] += in[0]*band2[0]+in[1]*band1[1];
          }
     }
}

// Fast matrix-vector multiplication when the matrix is symmetric 3-diagonal
void fast1band(const size_t N, const double diag[N], const double band[N-1], double complex out[N], const double complex in[N]) {

     size_t i;

     if (N > 2) {
     out[0] = diag[0]*in[0] + band[0]*in[1];
     for (i = 1; i < N-1; i++) {
          out[i] = diag[i]*in[i] + band[i]*in[i+1] + band[i-1]*in[i-1];
     }
     out[N-1] = diag[N-1]*in[N-1] + band[N-2]*in[N-2];

     } else { // N <= 2
          for (i = 0; i < N; i++) out[i] = diag[i]*in[i];
          if (N == 2) {
               out[0] += in[1]*band[0];
               out[1] += in[0]*band[0];
          }
     }
}



// nc: no clobber, i.e. no overwriting of the input vector. 
static void matvec_nc(const size_t N, const double mat[N][N], double complex out[N], const double complex vec[N]) {

     size_t i,j,k;
     // Do the matrix-vector operation num_rows matrix rows at a time.
     // This utilizes the CPU cache and the register bank better, as we only
     // need to reload the input vector N/num_rows times instead of N times.
     // At some point, we start to have to reload rows into the cache/registers
     // so there is an optimal num_rows size.
     // On my laptop, this is 4. Some other machine might have a
     // diferent optimal num_rows size, eg. machines with more registers,
     // like ones with the AVX-512 instruction set.
     const size_t num_rows=4; //important that this is known at compile time!
     double complex tmp[num_rows];

     const size_t remaining_rows = N%num_rows;
     
     for (i = 0; i < N-remaining_rows; i+=num_rows) {
          for (k = 0; k < num_rows; k++) tmp[k]=0;
          for (j = 0; j < N; j++) {
               for (k = 0; k < num_rows; k++) {
                    tmp[k] += mat[i+k][j]*vec[j];
               }
          }
          for (k = 0; k < num_rows; k++) {
               out[i+k] = tmp[k];
          }
     }
     // Do the remaining rows  
     for (k = 0; k < remaining_rows; k++) tmp[k]=0;
     for (j = 0; j < N; j++) {
          for (k = 0; k < remaining_rows; k++) {
               tmp[k] += mat[i+k][j]*vec[j];
          }
     }
     for (k = 0; k < remaining_rows; k++) {
          out[N-remaining_rows+k] = tmp[k];
     }
}

static void matvec(size_t N, const double mat[N][N], double complex vec[N]) {

     double complex out[N];
     matvec_nc(N,mat,out,vec);
     memcpy(vec,out,N*sizeof(double complex));
}
/*
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
}*/

static void vecvec(size_t N, double complex v1[N], const double complex v2[N]) {
     size_t i;
     for (i = 0; i < N; i++) {
          v1[i] = v1[i]*v2[i];
     }
}

double complex dot(size_t N, const double complex a[N], const double complex b[N]) {
     complex double res;
     size_t i;

     res = 0;
     for (i = 0; i < N; i++) {
          res += conj(a[i])*b[i];
     }
     return res;
}
double norm(size_t N, const double complex a[N]) {
     return sqrt(creal(dot(N,a,a)));
}


int detect_parity(size_t dim, const double complex psi[dim]) {
     size_t i;

     _Bool possible_parity[2] = {true,true};

     for (i = 0; i < dim; i++) {
          if (psi[i] != 0) possible_parity[i&1]=false;
     }

     if (possible_parity[0]) return -1;
     if (possible_parity[1]) return 1;
     return 0;
}

static int reduce_cos2dmat(int K, int M, size_t Jmax, size_t *dimm, const double complex psi_0[Jmax+1], void *buffer,const double U2d[Jmax+1][Jmax+1]) {

     size_t i,j;
     int parity;
     unsigned int offset;
     size_t dim;


     if (K == 0 || M == 0) {
          parity = detect_parity(Jmax+1,psi_0);
          if (parity != 0) {
               dim = (parity==1) ? (Jmax+1)/2 + ODD(Jmax+1) : Jmax/2 + ODD(Jmax);
               *dimm = dim;
               offset = (parity == 1) ? 0 : 1;
               double (*U2d_reduced)[dim][dim] = buffer;
               for (i = 0; i < dim; i++)
                    for (j = 0; j < dim; j++) {
                         (*U2d_reduced)[i][j] = U2d[2*i+offset][2*j+offset];
                    }

               return parity;
          }
     }
     return 0;
}

static inline void downscale(size_t dim, size_t jmax, unsigned int offset, double complex out[dim], double complex in[jmax+1]) {
     size_t i;
     for (i = 0; i < dim; i++) out[i] = in[2*i+offset];
}

void fieldfree_propagation(int K, int M, const size_t Jmax, const double complex psi_0[Jmax+1], double t0, size_t num_times, const double times[num_times], const double E_rot[Jmax+1], const double Udiag[Jmax+1], const double Uband1[Jmax], const double Uband2[Jmax-1], double complex psi_result[num_times][Jmax+1], double cos2[num_times], const _Bool do_cos2d, const double U2d[Jmax+1][Jmax+1], double cos2d[num_times]) {

     size_t i;
     double t, Dt;
     double complex phase, U_psi[Jmax+1];
     double complex phase_incr, phase_incr_base;
     double Bconst = (E_rot[1]-E_rot[0])/2; // also works for symmetric tops!
     size_t J;
     int parity = 0;
     unsigned int offset = 0;
     size_t dim = 0;
     double U2d_reduced[(Jmax/2+2)*(Jmax/2+2)]; //Allocate enough, if a bit much

     if (do_cos2d) {
          // Check if the wave function has parity.
          // If it does, we don't need to multiply every other row in the
          // U2d matrix.
          parity = reduce_cos2dmat(K,M,Jmax,&dim,psi_0,U2d_reduced,U2d);
          offset = (parity == 1) ? 0 : 1;
     }

     double complex U2d_psi[dim], psi_downscaled[dim];

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
          // exp(-i*B*j*(j+1)) = exp(-i*B*j*(j-1)*exp(-2Bj*i)
          // = exp(-i*B*j*(j-1)*exp(-2B*i)^j,
          // i.e. phase(j) = phase(j-1)*exp(-2B*i)^j
          // This avoids the very costly evaluation of cexp() multiple times.
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
               cos2d[i] = 0;
               // U2d_psi = U2d dot psi
               if (parity == 0) {
                    matvec_nc(Jmax+1, U2d, U_psi, psi_result[i]);
                    for (J = 0; J <= Jmax; J++)
                         cos2d[i] += conj(psi_result[i][J])*U_psi[J];
               } else {
                    // Use the downscaled U2d matrix because every other
                    // entry in psi is 0. This gives a fairly significant
                    // speedup
                    downscale(dim,Jmax,offset,psi_downscaled,psi_result[i]);
                    matvec_nc(dim,(const double (*)[dim]) &U2d_reduced[0],U2d_psi,psi_downscaled);
                    for (J=0; J<dim; J++) {
                         cos2d[i] += conj(psi_downscaled[J])*U2d_psi[J];
                    }

               }
          }
     }
}


// Like the python version, except without all the implicit memory allocations
// and deallocations
int propagate_field(size_t Nsteps, size_t Nsteps_inner, size_t dim, double t, double dt, double E0_sq_max, double sigma, double complex psi_t[Nsteps][dim], const double eig[dim], const double eigvec[dim][dim], const double eigvecT[dim][dim], const double complex expRot[dim], const double complex expRot2[dim]) {
     

     size_t i, j, k;
     double fac;
     double sigma_sq = sigma*sigma;
     
     for (i = 1; i < Nsteps; i++) {
          for (j = 0; j < dim; j++)
               psi_t[i][j] = expRot2[j]*psi_t[i-1][j];
          for (k = 0; k < Nsteps_inner; k++) {
               if (k > 0)
                    vecvec(dim,psi_t[i],expRot); // psi = expRot*psi;
               //matTvec(dim,eigvec,psi_t[i]); // psi = eigvec^T dot psi
               matvec(dim,eigvecT,psi_t[i]); // psi = eigvec^T dot psi
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
     size_t dim;
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
     const size_t dim = p->dim;
     // Gaussian pulse:
     double E_0_squared = p->peak_field_amplitude_squared * \
                          exp(-t*t/(2*(p->sigma)*(p->sigma)));
     double tmp;

     
     // Multiply the interaction term
     fast2band(dim,p->V0,p->V1,p->V2,(double complex *) dPsidt, (const double complex *) psi);

     // Then add the diagonal and multiply every component with -i.
     // I.e. -i*(a+ib) = b-ia.
     for (j = 0; j < dim; j++) {
          dPsidt[2*j] = E_0_squared*dPsidt[2*j] + p->E_rot[j]*psi[2*j];
          dPsidt[2*j+1] = E_0_squared*dPsidt[2*j+1] + p->E_rot[j]*psi[2*j+1];

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
     struct deriv_params p = {dim,E0_sq_max,sigma,E_rot,V0,V1,V2,0};
     gsl_odeiv2_system sys = {.function=deriv,.jacobian=NULL,\
          .dimension=2*dim,.params=&p};

     double initial_step_size = min(2.4*sigma/150,1/E_rot[dim-1]/5.7); // 150 steps per pulse
     gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new (&sys, \
               gsl_odeiv2_step_rk8pd, initial_step_size, abstol,reltol);
               // rk8pd seems fastest. Did not try the methods
               // that rely on the Jacobian. d(dpsi/dt)/dpsi is diagonal,
               // but d(dpsi/dt)/dt is less trivial.
               //gsl_odeiv2_step_rkck, initial_step_size, abstol,reltol);
               //gsl_odeiv2_step_msadams, initial_step_size, abstol,reltol);
               //gsl_odeiv2_step_rk2, initial_step_size, abstol,reltol);
               //gsl_odeiv2_step_rkf45, initial_step_size, abstol,reltol);

     //gsl_odeiv2_driver_set_hmax(d,max step);

     double t_run = t;

     for (i = 1; i < Nsteps; i++) {
     
          memcpy(psi_t[i],psi_t[i-1],sizeof(double complex)*dim);
          gsl_odeiv2_driver_apply(d,&t_run,t+(double)i*dt,(double *) psi_t[i]);
     
          if (fmax(pow(cabs(psi_t[i][dim-2]),2),pow(cabs(psi_t[i][dim-1]),2)) > 1e-5) {
               gsl_odeiv2_driver_free(d);
               return 1; // Basis size too small
          }

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
