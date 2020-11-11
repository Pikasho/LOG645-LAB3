/* Solver
 * 
 * Fonctions permettant de solve le problème de transfert de chaleur sur une plaque 2D
 * 
* Philippe Bolduc
* Rémi Thibault
 */

#ifndef SOLVER_HPP
#define SOLVER_HPP

void solveSeq(int rows, int cols, int iterations, double td, double h, int sleep, double * matrix);
void solvePar(int rows, int cols, int iterations, double td, double h, int sleep, double * matrix);

#endif
