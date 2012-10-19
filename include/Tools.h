#ifndef TOOLS_H
#define TOOLS_H

#include <iostream>
#include <vector>

using std::ostream;

#include "Matrix.h"
#include "input.h"
#include "R.h"

/**
 * @author Brecht Verstichel
 * @date 15-04-2012\n\n
 * This class Tools is a class written to precalculate and store some of the properties needed
 * by the integral calculaters, to speed up the program.
 */
class Tools {

   public:
      
      static double gnorm(int,int);

      static double ggamma(int,int,int,int,int,int);

      static const R &gP(int,int,int,int,int,int);

      static const CI_SPPM &gS();

      static const CI_SPPM &gT();

      static const CI_SPPM_m &gU(int);

      static const CI_TPPM &gV();

      static double gprefac_overlap(int,int,int,int,int,int);

      static double gdfac(int);

      static double gfac(int);

      static double gcomb(int,int);

      static double norm_const(int,int,int,double);

      static void init();

      static void clear();

   private:

      //!easy to store this number, the number of different gaussian orbitals
      static int gaussdim;

      //!and of course a list to transform this index to core and j indices
      static vector< vector<int> > gauss2ij;

      //!and a reverse list
      static int **ij2gauss;

      //!a list that transforms the index of pairs of different cores to the different cores
      static vector< vector<int> > tc2i;

      //!the reverse list
      static int **i2tc;

      //!double array containing the norms of the different primitives
      static double **norm;

      //!matrix of "Gaussian index" dimension containing matrices of the sum of the exponents of the primitives
      static double ****gamma;

      //!matrix that stores the "P" objects of the gaussians on different cores
      static R ******P;

      //!matrix that stores the overlap prefactor for the objects on different cores
      static double *****prefac_overlap;

      //!array containg double faculties
      static double *dfac;

      //!array containing regular faculties
      static double *fac;

      //!array containing the combinations of k out of n
      static double **comb;

      //!will be used for the recursive overlapmatrix calculations
      static CI_SPPM *S;

      //!will be used for the recursive kinetic energy matrix calculations
      static CI_SPPM *T;

      //!will be used for the recursive nuclear attraction matrix element calculations
      static CI_SPPM_m **U;

      //!will be used for the evaluation of the electron repulsion integrals
      static CI_TPPM *V;

};

#endif

/* vim: set ts=3 sw=3 expandtab :*/
