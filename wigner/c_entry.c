// Prototype netlib/slatec drc3jj fortran function:
void drc3jj_(double *L2, double *L3, double *M2, double *M3, double *L1min, double *L1max, double *thrcof, int *ndim, int *ier);

// Small C wrapper:
int drc3jj(int l2, int l3, int m2, int m3, int *l1min, int *l1max, double *thrcof, int ndim) {

     int ierr;

     double l2d, l3d, m2d, m3d, l1mind, l1maxd;
     l2d = (double) l2;
     l3d = (double) l3;
     m2d = (double) m2;
     m3d = (double) m3;


     drc3jj_(&l2d,&l3d,&m2d,&m3d,&l1mind,&l1maxd,thrcof,&ndim,&ierr);
     *l1min = (int) l1mind;
     *l1max = (int) l1maxd;
     return ierr;

}
