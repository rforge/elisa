/* Babak Naimi, July 2016
   naimi.b@gmail.com
   July 2016
   v 1.0
*/


#include <R.h>
#include <Rinternals.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <Rmath.h>
#include "Rdefines.h"
#include "R_ext/Random.h"
#include "R_ext/Rdynload.h"



///--- semivariance for raster dataset at given lag distance
SEXP semivar(SEXP v, SEXP nc, SEXP nr, SEXP rr, SEXP cc) {
  int nProtected=0;
  int c, row, col, ngb, nnr, nnc, nrow, ncol, cellnr, n;
  double xi;
  
  R_len_t i, j;
  
  SEXP ans;
  double  *xv, *xans;
  int *xrr, *xcc ;
  
  nrow=INTEGER(nr)[0];
  ncol=INTEGER(nc)[0];
  
  n=length(v);
  
  PROTECT(v = coerceVector(v, REALSXP));
  ++nProtected;
  
  
  PROTECT(ans = allocVector(REALSXP, n));
  ++nProtected;
  
  
  PROTECT(rr = coerceVector(rr, INTSXP));
  ++nProtected;
  
  PROTECT(cc = coerceVector(cc, INTSXP));
  ++nProtected;
  
  ngb=length(rr);
  
  xans=REAL(ans);
  xv=REAL(v);
  xrr=INTEGER(rr);
  xcc=INTEGER(cc);
  
  int nv=0;
  double Eij,nn;
  int *xcells=malloc(n*sizeof(int));
  for (i=0;i < n;i++) {
    if (!R_IsNA(xv[i])) {
      xcells[nv]=i;
      nv++;
    } else {
      xans[i]=R_NaReal;
    }
  }
  if (nv < n) xcells=realloc(xcells,nv*sizeof(int));
  nn=nv;
  
  for (c=0;c < nv;c++) {
    j=xcells[c];
    xi=xv[j];
    row = (j / ncol) + 1;
    col = (j + 1) - ((row - 1) * ncol);
    
    int q=-1;
    Eij=0;
    for (i=0; i < ngb; i++) {
      nnr= row + xrr[i];
      nnc = col + xcc[i];
      
      if ((nnr > 0) & (nnr <= nrow) & (nnc > 0) & (nnc <= ncol)) {
        cellnr = ((nnr - 1) * ncol) + nnc;
        if (!R_IsNA(xv[(cellnr-1)])) {
          Eij=Eij+pow(xi - xv[(cellnr-1)],2);
          q+=1;
        }
      }
    }
    //----
    xans[j] = (Eij/(double)q)/2;
  }
  free(xcells);
  UNPROTECT(nProtected);
  return(ans);
}
/////
/////
SEXP semivar_vector(SEXP v, SEXP nb) {
  int nProtected=0;
  int c, ngb, n, q;
  double xi, a;
  
  R_len_t i, j;
  
  SEXP ans;
  double  *xv, *xans;
  
  n=length(v);
  
  PROTECT(v = coerceVector(v, REALSXP));
  ++nProtected;
  
  PROTECT(ans = allocVector(REALSXP, n));
  ++nProtected;
  
  xans=REAL(ans);
  xv=REAL(v);
  
  int nv=0;
  double Eij,nn;
  int *xcells=malloc(n*sizeof(int));
  for (i=0;i < n;i++) {
    if (!R_IsNA(xv[i])) {
      xcells[nv]=i;
      nv++;
    } else {
      xans[i]=R_NaReal;
    }
  }
  if (nv < n) xcells=realloc(xcells,nv*sizeof(int));
  nn=nv;
  
  for (c=0;c < nv;c++) {
    j=xcells[c];
    xi=xv[j];
    ngb = length(VECTOR_ELT(nb,c));
    q=0;
    Eij=0;
    for (i=0;i < ngb;i++) {
      a=xv[INTEGER_POINTER(VECTOR_ELT(nb,c))[i] - 1];
      if (!R_IsNA(a)) {
        Eij=Eij+pow(xi - a,2);
        q+=1;
      }
    }
    
    //----
    xans[j] = (Eij/(double)q)/2;
  }
  free(xcells);
  UNPROTECT(nProtected);
  return(ans);
}
/////