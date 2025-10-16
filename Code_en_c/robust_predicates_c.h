#ifndef ROBUST_PREDICATES_C_H
#define ROBUST_PREDICATES_C_H

#ifdef __cplusplus
extern "C" {
#endif

/* incircle: appelle la version robuste (Shewchuk).
 * Chaque point est un tableau double[2] (x,y) ou double[3] si vous utilisez 3D.
 * Retour:
 *   >0 : pd est à l'intérieur du cercle circonscrit à (pa,pb,pc) si (pa,pb,pc) sont CCW
 *   <0 : à l'extérieur
 *   =0 : sur le cercle
 */
double incircle_c(const double *pa, const double *pb, const double *pc, const double *pd);

/* Appeler une fois avant d'utiliser les prédicats. (Prépare les constantes internes). */
void exactinit_c(void);

#ifdef __cplusplus
}
#endif

#endif /* ROBUST_PREDICATES_C_H */