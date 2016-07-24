/* Babak Naimi, August 2014 
   naimi.b@gmail.com
*/
#include <R.h>
#include <Rinternals.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <Rmath.h>
#include "Rdefines.h"
#include "R_ext/Rdynload.h"

SEXP elsac(SEXP v, SEXP nc, SEXP nr, SEXP nclass, SEXP rr, SEXP cc, SEXP classes, SEXP dif) {
  int nProtected=0;
  int c, row, col, ngb, q, nnr, nnc, nrow, ncol, cellnr, ncl, n, a, nw;
  double e, w, s, xi, qq, count, maxW;
  
  R_len_t i, j;
  
  SEXP ans;
  
  double *xans, *xv, *xdif;
  int *xrr, *xcc, *xcls;
  
  nrow=INTEGER(nr)[0];
  ncol=INTEGER(nc)[0];
  ncl=INTEGER(nclass)[0]; //nclass for categorical variables referes to the number of unique classes
  
  n=length(v);
  
  
  PROTECT(v = coerceVector(v, REALSXP));
  ++nProtected;
  
  PROTECT(ans = allocVector(REALSXP, n));
  ++nProtected;
  
  PROTECT(rr = coerceVector(rr, INTSXP));
  ++nProtected;

  PROTECT(cc = coerceVector(cc, INTSXP));
  ++nProtected;
  
  PROTECT(classes = coerceVector(classes, INTSXP));
  ++nProtected;
  
  PROTECT(dif = coerceVector(dif, REALSXP));
  ++nProtected;
  
  ngb=length(rr);
  
  xans=REAL(ans);
  xv=REAL(v);
  xrr=INTEGER(rr);
  xcc=INTEGER(cc);
  xcls=INTEGER(classes);
  xdif=REAL(dif);
  
  maxW=0;
  for (i=0; i < length(dif);i++) {
    if (xdif[i] > maxW) maxW=xdif[i];
  }
  
  for (c=0;c < n;c++)  {
    xi=xv[c];
    if (!R_IsNA(xi)) {
      row = (c / ncol) + 1;
      col = (c + 1) - ((row - 1) * ncol);
      
      double xn[ngb], xw[ngb];
      //------
      for (i=0; i < ncl;i++) {
        if (xcls[i] == xi) {
          nw=i;
          break;
        }
      }
      //-------
      
      q=-1;
      for (i=0; i < ngb; i++) {
        nnr= row + xrr[i];
        nnc = col + xcc[i];
        
        if ((nnr > 0) & (nnr <= nrow) & (nnc > 0) & (nnc <= ncol)) {
          cellnr = ((nnr - 1) * ncol) + nnc;
          if (!R_IsNA(xv[(cellnr-1)])) {
            q+=1;
            xn[q]=xv[(cellnr-1)];
            for (j=0;j < ncl;j++) {
              if (xcls[j] == xn[q]) {
                xw[q]=xdif[(nw*ncl)+j];
                break;
              }
            }
          }
        }
      }
       
      // sort
      for (i=0;i <= (q-1);i++) {
        for (j=i+1;j <= q;j++) {
          if (xn[i] > xn[j]) {
            a=xn[i];
            xn[i]=xn[j];
            xn[j]=a;
          }
        }
      }
      //------
      
      
      
      
      //
      a=xn[0];
      count=1;
      e=0;
      qq=q+1;
      
      for (i=1;i <= q;i++) {
        if (xn[i] != a) {
          e = e + ((count / qq) * log2(count / qq));
          a=xn[i];
          count=1;
        } else {
          count+=1;
        }
      }
      e = e + ((count / qq) * log2(count / qq));
      w=0;
      for (i=0; i <= q;i++) {
        w = w + xw[i];
      }
      w = w / ((qq - 1) * maxW);
      
      if (qq > ncl) {
        s = log2(ncl);
      } else {
        s = log2(qq);
      }
      
      xans[c] = (-e * w) / s;
      
    } else {
      xans[c]=R_NaReal;
    }
  }
  UNPROTECT(nProtected);
  return(ans);
  
}
//------

SEXP elsac_cell(SEXP v, SEXP nc, SEXP nr, SEXP nclass, SEXP rr, SEXP cc, SEXP classes, SEXP dif, SEXP cells) {
  int nProtected=0;
  int c, row, col, ngb, q, nnr, nnc, nrow, ncol, cellnr, ncl, n, a, nw, cn;
  double e, w, s, xi, qq, count, maxW;
  
  R_len_t i, j;
  
  SEXP ans;
  
  double *xans, *xv, *xdif;
  int *xrr, *xcc, *xcls, *xcells;
  
  nrow=INTEGER(nr)[0];
  ncol=INTEGER(nc)[0];
  ncl=INTEGER(nclass)[0]; //nclass for categorical variables referes to the number of unique classes
  
  n=length(cells);
  
  
  PROTECT(v = coerceVector(v, REALSXP));
  ++nProtected;
  
  PROTECT(ans = allocVector(REALSXP, n));
  ++nProtected;
  
  PROTECT(rr = coerceVector(rr, INTSXP));
  ++nProtected;

  PROTECT(cc = coerceVector(cc, INTSXP));
  ++nProtected;
  
  PROTECT(classes = coerceVector(classes, INTSXP));
  ++nProtected;
  
  PROTECT(dif = coerceVector(dif, REALSXP));
  ++nProtected;
  
  PROTECT(cells = coerceVector(cells, INTSXP));
  ++nProtected;
  
  ngb=length(rr);
  
  xans=REAL(ans);
  xv=REAL(v);
  xrr=INTEGER(rr);
  xcc=INTEGER(cc);
  xcls=INTEGER(classes);
  xdif=REAL(dif);
  xcells=INTEGER(cells);
  
  maxW=0;
  for (i=0; i < length(dif);i++) {
    if (xdif[i] > maxW) maxW=xdif[i];
  }
  
  for (c=0;c < n;c++)  {
    cn=xcells[c]-1;
    xi=xv[cn];
    if (!R_IsNA(xi)) {
      row = (cn / ncol) + 1;
      col = (cn + 1) - ((row - 1) * ncol);
      
      double xn[ngb], xw[ngb];
      //------
      for (i=0; i < ncl;i++) {
        if (xcls[i] == xi) {
          nw=i;
          break;
        }
      }
      //-------
      
      q=-1;
      for (i=0; i < ngb; i++) {
        nnr= row + xrr[i];
        nnc = col + xcc[i];
        
        if ((nnr > 0) & (nnr <= nrow) & (nnc > 0) & (nnc <= ncol)) {
          cellnr = ((nnr - 1) * ncol) + nnc;
          if (!R_IsNA(xv[(cellnr-1)])) {
            q+=1;
            xn[q]=xv[(cellnr-1)];
            for (j=0;j < ncl;j++) {
              if (xcls[j] == xn[q]) {
                xw[q]=xdif[(nw*ncl)+j];
                break;
              }
            }
          }
        }
      }
       
      // sort
      for (i=0;i <= (q-1);i++) {
        for (j=i+1;j <= q;j++) {
          if (xn[i] > xn[j]) {
            a=xn[i];
            xn[i]=xn[j];
            xn[j]=a;
          }
        }
      }
      //
      a=xn[0];
      count=1;
      e=0;
      qq=q+1;
      
      for (i=1;i <= q;i++) {
        if (xn[i] != a) {
          e = e + ((count / qq) * log2(count / qq));
          a=xn[i];
          count=1;
        } else {
          count+=1;
        }
      }
      e = e + ((count / qq) * log2(count / qq));
      w=0;
      for (i=0; i <= q;i++) {
        w = w + xw[i];
      }
      w = w / ((qq - 1) * maxW);
      
      if (qq > ncl) {
        s = log2(ncl);
      } else {
        s = log2(qq);
      }
      
      xans[c] = (-e * w) / s;
      
    } else {
      xans[c]=R_NaReal;
    }
  }
  UNPROTECT(nProtected);
  return(ans);
  
}
//////

SEXP elsac_vector(SEXP v, SEXP nb, SEXP nclass) {
  int nProtected=0;
  int  ncl, n, a, q,xi;
  double e, w, s,  qq, count;
  R_len_t i, j, c;
  
  SEXP ans;
  double *xans;
  int *xv;
  ncl=INTEGER(nclass)[0];
  n=length(v);
  
  PROTECT(v = coerceVector(v, INTSXP));
  ++nProtected;
  
  PROTECT(ans = allocVector(REALSXP, n));
  ++nProtected;
  
  xans=REAL(ans);
  xv=INTEGER(v);
  
  for (c=0;c < n;c++)  {
    R_CheckUserInterrupt();
    xi=xv[c];
    if (!R_IsNA(xi)) {
      
      q = length(VECTOR_ELT(nb,c));
      
      int xn[q+1];
      
      for (i=0;i < q;i++) {
        xn[i]=xv[INTEGER_POINTER(VECTOR_ELT(nb,c))[i] - 1];
      }
      
      xn[q]=xi;
      
      // sort
      for (i=0;i <= (q-1);i++) {
        for (j=i+1;j <= q;j++) {
          if (xn[i] > xn[j]) {
            a=xn[i];
            xn[i]=xn[j];
            xn[j]=a;
          }
        }
      }
      //------
      
      a=xn[0];
      count=1;
      e=0;
      qq=q+1;
      
      for (i=1;i <= q;i++) {
        if (xn[i] != a) {
          e = e + ((count / qq) * log2(count / qq));
          a=xn[i];
          count=1;
        } else {
          count+=1;
        }
      }
      e = e + ((count / qq) * log2(count / qq));
      w=0;
      for (i=0; i <= q;i++) {
        w = w + abs(xn[i] - xi);
      }
      w = w / ((qq - 1) * (ncl - 1));
      
      if (qq > ncl) {
        s = log2(ncl);
      } else {
        s = log2(qq);
      }
      
      xans[c] = (-e * w) / s;
      
    } else {
      xans[c]=R_NaReal;
    }
  }
  UNPROTECT(nProtected);
  return(ans);
}