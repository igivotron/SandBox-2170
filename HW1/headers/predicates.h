#ifndef _ROBUST_PREDICATES_H_
#define _ROBUST_PREDICATES_H_

// namespace necessary to avoid conflicts with predicates used by Tetgen
double exactinit();
double incircle(double *pa, double *pb, double *pc, double *pd);
double insphere(double *pa, double *pb, double *pc, double *pd, double *pe);
double orient2d(double *pa, double *pb, double *pc);
double orient3d(double *pa, double *pb, double *pc, double *pd);


#endif // _ROBUST_PREDICATES_H_