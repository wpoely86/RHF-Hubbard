#include <assert.h>
#include <string>

#include "include.h"

/**
 * copy constructor
 * @param spm_copy content of this matrix will be copied into the constructed matrix
 */
DIIS::DIIS(const DIIS &diis_copy)
{
   FD = diis_copy.FD;
   F = diis_copy.F;
   b = diis_copy.b;

   B.reset(new Matrix(*diis_copy.B));
}


std::ostream &operator<<(std::ostream &output,DIIS &diis_p)
{
   return output;
}

/**
 * add a new commutator to the DIIS object
 * @param FD_i new commutator to add
 */
void DIIS::push_comm(const RSPM &FD_i)
{
   FD.push_back(FD_i);
}

/**
 * add a new Fock matrix to the DIIS object
 * @param F_i new Fock matrix to add
 */
void DIIS::push_F(const RSPM &F_i)
{
   F.push_back(F_i);
}

/** 
 * @return the vector containing the commutator matrices
 */
const std::vector<RSPM> &DIIS::gFD() const
{
   return FD;
}

/**
 * access to the individual commutators
 * @param i index of the iteration
 * @return commutator with index i
 */
const RSPM &DIIS::gFD(int i) const
{
   assert(i < FD.size());
   return FD[i];
}

/**
 * access to the individual Fock Matrices
 * @param i index of the iteration
 * @return Fock matrix with index i
 */
const RSPM &DIIS::gF(int i) const
{
   assert(i < F.size());
   return F[i];
}

/**
 * construct the DIIS matrix
 */ 
void DIIS::construct()
{
   //construct the B matrix
   B.reset(new Matrix(FD.size() + 1));

#pragma omp parallel for
   for(unsigned int i = 0;i < FD.size();++i)
   {
      for(unsigned int j = i;j < FD.size();++j)
         (*B)(i,j) = (*B)(j,i) = FD[i].ddot(FD[j]);

      (*B)(i,FD.size()) = (*B)(FD.size(),i) = 1.0;
   }

   (*B)(FD.size(),FD.size()) = 0.0;

   //and the vector on the right hand side
   b.resize(FD.size() + 1);

   for(unsigned int i = 0;i < FD.size();++i)
      b[i] = 0;

   b[FD.size()] = 1;

   //now solve!
   B->solve(b);
}

/**
 * access to the relaxation coefficients
 * @param i index of the coefficient
 */
double DIIS::gb(int i) const
{
   assert(i < b.size());

   return b[i];
}

/**
 * @return the size of the system
 */
int DIIS::size() const
{
   return FD.size();
}

/* vim: set ts=3 sw=3 expandtab :*/
