#include <stdio.h>

#include "propagation_nostatic.c"


void test_fieldfree() {

     int K=0, M=0;
     size_t Jmax = 140;
     double complex psi_0[Jmax+1];
     double t0=0;
     size_t num_times=440000;
     double times[num_times];
     double E_rot[Jmax+1];
     double Udiag[Jmax+1];
     double Uband1[Jmax];
     double Uband2[Jmax-1];
     double complex psi_result[num_times][Jmax+1];
     double cos2[num_times];
     _Bool do_cos2d = false;
     const double U2d[Jmax+1][Jmax+1];
     double cos2d[num_times];

     fieldfree_propagation(K, M, Jmax, psi_0, t0, num_times, times, E_rot, Udiag, Uband1, Uband2, psi_result, cos2, do_cos2d, U2d, cos2d);

}

int main() {
     
     test_fieldfree();

}
