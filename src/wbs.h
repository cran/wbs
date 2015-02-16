#ifndef WBS_H
#define WBS_H

#include <R.h>
#include <Rmath.h>

#define IDX(i,j,ld) ((((j)-1) * (ld))+((i)-1))

/* 
  
    List of functions:
    
    wbs_rec - recursive implementation of the Wild Binary Segmentation algorithm,
    wbs_rec_wrapper - wrapper for wbs_rec; can be called from R,
    wbs_int_rec - recursive implementation of the Wild Binary Segmentation algorithm in its augmented version,
    wbs_int_rec_wrapper - wrapper for wbs_int_rec; can be called from R,
    wbs_ipi - finds the values of the CUSUM statistics on given interval [s,e]; quick implementation in O(n),
    ipi_arg_max - finds the arg_max of the CUSUM's computed over interval [s,e],
    bs_rec - ecursive implementation of the Binary Segmentation algorithm,
    bs_rec_wrapper - wrapper for bs_rec; can be called from R.
    
*/

void wbs_rec(double *x, int n, int s, int e, double *res, double *wbsres, int *index, int indexn, int M, int scale);
void wbs_rec_wrapper(double *x, int *n, double *res,  int *intervals, int *M);
void wbs_int_rec(double *x, int n, int s, int e, double *res, double *iplus, double *iminus, double *ipres, double *wbsres, int *index, int indexn, int M, double minth, int scale);
void wbs_int_rec_wrapper(double *x, int *n, double *res,  int *intervals, int *M);
void wbs_ipi(double *x, int n, double *res, double *iplus, double *iminus, int *ipargmax, double *ipmax);
void ipi_arg_max(double *res, int n, int *ipargmax, double *ipmax);
void bs_rec(double *x, int n, int s, int e,  double *res,double *iplus, double *iminus, double *ipres, double minth, int scale);
void bs_rec_wrapper(double *x, int *n, double *res);




#endif
