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

#include <gsl/gsl_sf.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_monte.h>
// protip: you can switch between vegas and miser just by search and replace
// through the source file.
#include <gsl/gsl_monte_vegas.h>


//#define WINDOWS

#ifndef WINDOWS
#include <unistd.h>
#include <sys/mman.h>
#include <sys/file.h>
#endif

#define max(a,b) (a > b) ? a : b

double integrand(double u, double phi, int l, int lp, int m) {
     double Pl1, Pl2;
     const double usq = u*u;
     double sin2phi;

     if (u == 0) return 0;

     Pl1 = gsl_sf_legendre_sphPlm(l, m, u);
     if (l != lp) {
          Pl2 = gsl_sf_legendre_sphPlm(lp,m,u);
     } else {
          Pl2 = Pl1;
     }

     sin2phi = sin(phi);
     sin2phi = sin2phi*sin2phi;
     return Pl1*Pl2*usq/(usq+(1-usq)*sin2phi);
}

struct integrand_monte_params {
     int l, lp, m;
};

double integrand_monte_wrap(double *x, size_t dim, void *params) {
     (void)(dim); // Avoid unused warning
     struct integrand_monte_params *p = params;
     return integrand(x[0],x[1],p->l,p->lp,p->m);

}

double vegas_integral(int l, int lp, int m, gsl_monte_vegas_state *vegas_state, gsl_rng *rng) {

     double res, abserr;
     //const double abstol = 1e-14, reltol = 1e-14; 
     gsl_monte_function F;
     double xl[2] = {-1, 0};
     double xu[2] = { 1, 2*M_PI};
     size_t calls;
     
     gsl_monte_vegas_init(vegas_state);
     F.f = &integrand_monte_wrap;
     struct integrand_monte_params params = {.l=l,.lp=lp,.m=m};
     F.params = &params;
     F.dim=2;

     calls = 1000000;
     gsl_monte_vegas_integrate(&F,xl,xu,2,10000,rng,vegas_state,&res,&abserr);
     gsl_monte_vegas_integrate(&F,xl,xu,2,calls,rng,vegas_state,&res,&abserr);
     printf("%e\n",abserr);

     return res;

}


struct integrand_params {
     double phi;
     int l, lp, m;
};


double integrand_wrap(double u, void *params) {
     struct integrand_params *p = params;
     return integrand(u,p->phi,p->l,p->lp,p->m);
}

struct inner_params {
     gsl_function integrand;
     gsl_integration_workspace *workspace;
     size_t limit;
};

double inner(double phi, void *params) {
  
     struct inner_params *p = params;
     struct integrand_params *ip = p->integrand.params;
     const double abstol = 1e-14, reltol = 1e-14; 
     double result,abserr;

     ip->phi = phi;
     gsl_integration_qag(&(p->integrand),-1.0,1.0,abstol,reltol,p->limit,GSL_INTEG_GAUSS61,p->workspace,&result,&abserr);
     return result;
}

double outer(int l, int lp, int m, gsl_integration_workspace *ws_in, gsl_integration_workspace *ws_out, size_t limit) {

     double res, abserr;
     const double abstol = 1e-14, reltol = 1e-14; 
     gsl_function inner_F, F;
     
     inner_F.function = &integrand_wrap;
     struct integrand_params inner_params = {.phi=0,.l=l,.lp=lp,.m=m};
     inner_F.params = &inner_params;

     F.function = &inner;
     struct inner_params params = {.integrand=inner_F,.workspace=ws_in,.limit=limit};
     F.params = &params;

     gsl_integration_qag(&F, 0, 2*M_PI, abstol, reltol, limit, GSL_INTEG_GAUSS61, ws_out, &res, &abserr);

     return res;

}

void populate_U2d(int lmax, int lmin, int m, double U2d[lmax+1][lmax+1]) {

     int l,lp;
     const size_t limit = 409600;
     gsl_integration_workspace *inner_workspace,*workspace;

     #pragma omp parallel default(shared) private(l,lp,inner_workspace,workspace)
     {
     inner_workspace = gsl_integration_workspace_alloc(limit);
     workspace = gsl_integration_workspace_alloc(limit);
     #pragma omp for schedule(guided,4)
     for (l = m; l <= lmax; l++) {
     for (lp = max(l,lmin); lp <= lmax; lp+=2) {

     U2d[l][lp] = outer(l,lp,m,inner_workspace,workspace,limit);
     U2d[lp][l] = U2d[l][lp];

     }
     }

     gsl_integration_workspace_free(workspace);
     gsl_integration_workspace_free(inner_workspace);
     }
}


void populate_U2d_vegas(int lmax, int lmin, int m, double U2d[lmax+1][lmax+1]) {

     int l,lp;
     gsl_monte_vegas_state *vegas_state;
     gsl_rng *rng;
     unsigned long int seed;
     FILE *urnd = fopen("/dev/urandom","r");

     #pragma omp parallel default(shared) private(l,lp,vegas_state,rng,seed)
     {
     vegas_state = gsl_monte_vegas_alloc(2);
     //rng = gsl_rng_alloc(gsl_rng_mt19937);
     rng = gsl_rng_alloc(gsl_rng_taus);
     fread(&seed,sizeof(seed),1,urnd);
     gsl_rng_set(rng,seed);
     #pragma omp for schedule(guided,4)
     for (l = m; l <= lmax; l++) {
     for (lp = max(l,lmin); lp <= lmax; lp+=2) {

     U2d[l][lp] = vegas_integral(l,lp,m,vegas_state,rng);
     U2d[lp][l] = U2d[l][lp];

     }
     }

     gsl_monte_vegas_free(vegas_state);
     gsl_rng_free(rng);
     }

     fclose(urnd);
}



struct shm_state {
     void *p;
     size_t size;
#ifdef WINDOWS
     FILE *f;
#endif
};

#ifndef WINDOWS
void *shm_allocate(struct shm_state *state, const char *filename, size_t size) {

     const int prot = PROT_READ|PROT_WRITE;
     const int flags = MAP_SHARED;
     int fd;

     state->size = size;
     fd = open(filename,O_CREAT|O_RDWR,00660); 
     ftruncate(fd,(off_t) size);
     state->p = mmap(NULL,size,prot,flags,fd,0);
     close(fd);
     return state->p;
}

void shm_deallocate(struct shm_state *state) {
     munmap(state->p,state->size);
}
#else
void *shm_allocate(struct shm_state *state, const char *filename, size_t size) {
     state->size = size;
     state->f = fopen(filename,"wb"); 
     state->p = malloc(size);
     return state->p;
}

void shm_deallocate(struct shm_state *state) {
     fwrite(state->p,1,state->size,state->f);
     fclose(state->f);
     free(state->p);
}
#endif

int main(int argc, char **argv) {

     int lmax, lmin, k, m, kmsign;
     size_t U2dsize;
     char *filename;
     struct shm_state mmap_state;

     if (argc != 6 && argc != 7) {
          fprintf(stderr,"\n");
          fprintf(stderr,"Usage: %s filename lmax k m kmsign [lmin]\n", argv[0]);
          fprintf(stderr,"\n");
          fprintf(stderr,"Calculate <l'km| cos^2theta 2d |lkm> matrix elements.\n");
          fprintf(stderr,"Both k and m must be nonnegative.\n");
          fprintf(stderr,"If both k and m are nonzero, kmsign is the sign of k*m. Otherwise it is 1.\n");
          fprintf(stderr,"Result stored in the file by name filename.\n");
          fprintf(stderr,"\n");
          exit(1);
     }
     
     filename = argv[1];
     lmax = atoi(argv[2]);
     k = atoi(argv[3]);
     m = atoi(argv[4]);
     kmsign = atoi(argv[5]);
     if (argc == 7) {
          lmin = atoi(argv[6]);
     } else {
          lmin = 0;
     }


     if (k != 0) {
          // This requires the evaluation of a tripple integral
          // over all 3 Euler angles
          fprintf(stderr,"k != 0 not implemented.\n");
          exit(1);
     }
     if ((kmsign != 1 && kmsign != -1) || (k*m == 0 && kmsign != 1) ) {
          fprintf(stderr, "Invalid KMsign. Must be -1 or 1.\n");
          exit(1);
     }

     U2dsize = (size_t) ((lmax+1)*(lmax+1));
     U2dsize *= sizeof(double);

     double (*U2d)[lmax+1][lmax+1];

     U2d = shm_allocate(&mmap_state,filename,U2dsize);

     memset(U2d,0,U2dsize);
 
     //populate_U2d(lmax, lmin, m, *U2d);
     populate_U2d_vegas(lmax, lmin, m, *U2d);

     shm_deallocate(&mmap_state);


     return 0;
}
