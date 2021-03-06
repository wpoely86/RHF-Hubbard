#include <iostream>

#include "include.h"

int RSPM::M;
int RSPM::N;
double RSPM::shift;

/**
 * static function that initializes the static variables 
 * @param M_i dimension of single particle space and dimension of the Matrix
 * @param N_i Nr of particles
 */
void RSPM::init(int M_i,int N_i, double shift_i)
{
   M = M_i;
   N = N_i;
   shift = shift_i;
}

RSPM::RSPM() : Matrix(M) { }

/**
 * copy constructor
 * @param spm_copy content of this matrix will be copied into the constructed matrix
 */
RSPM::RSPM(const RSPM &spm_copy) : Matrix(spm_copy) { }

/**
 * @return nr of particles
 */
int RSPM::gN() const
{
   return N;
}

/**
 * @return dimension of sp space and of matrix
 */
int RSPM::gM() const
{
   return M;
}

std::ostream &operator<<(std::ostream &output,RSPM &spm_p){

   for(int i = 0;i < spm_p.gM();++i)
      for(int j = 0;j < spm_p.gM();++j)
         output << i << "\t" << j << "\t" << spm_p(i,j) << std::endl;

   return output;
}

/**
 * fill the RSPM with the unit matrix on the first N/2 orbitals
 */
void RSPM::unit(){

   *this = 0.0;

   for(int i = 0;i < N/2;++i)
      (*this)(i,i) = 1.0;
}

/**
 * construct the Fock matrix, given the rspm from the previous iteration
 */
void RSPM::construct_Fock(const RSPM &spm)
{
   for(int a = 0;a < M;++a)
      for(int c = a;c < M;++c){

         //first the SP part
         (*this)(a,c) = Hubbard::gT(a,c);

         //the Hartree term
         for(int b = 0;b < M;++b)
            for(int d = 0;d < M;++d)
               (*this)(a,c) += 2.0 *Hubbard::gV(a,b,c,d) * spm(b,d);

         //Fock
         for(int b = 0;b < M;++b)
            for(int d = 0;d < M;++d)
               (*this)(a,c) -= Hubbard::gV(a,b,d,c) * spm(b,d);

      }

   daxpy(shift, spm);

   this->symmetrize();
}

void RSPM::setshift(double shift_i)
{
   shift = shift_i;
}

/**
 * construct the new rspm guess using the diagonalized Fock matrix
 */
void RSPM::update(const RSPM &F)
{
   for(int a = 0;a < M;++a)
      for(int b = a;b < M;++b){

         (*this)(a,b) = 0.0;

         for(int n = 0;n < N/2;++n)
            (*this)(a,b) += F(a,n) * F(b,n);

         (*this)(b,a) = (*this)(a,b);
      }
}

/**
 * construct the Fock matrix, given the rspm from the previous iteration
 */
void RSPM::construct_sp_ham()
{
   for(int a = 0;a < M;++a)
      for(int c = a;c < M;++c)
      {
         //first the SP part
         (*this)(a,c) = (*this)(c,a) = Hubbard::gT(a,c);
      }
}

/**
 * take a linear combination of the matrices Fock matrices in the DIIS object using the relaxation coefficients
 */
void RSPM::relax(const DIIS &diis)
{
   *this = 0.0;

   for(int i = 0;i < diis.size();++i)
      this->daxpy(diis.gb(i),diis.gF(i));
}

/**
 * @return the expectation value of the energy
 */
double RSPM::energy() const
{
   double energy = 0.0;

   for(int a = 0;a < M;++a)
      for(int c = 0;c < M;++c){

         energy += 2.0*(*this)(a,c) * Hubbard::gT(a,c);

         for(int b = 0;b < M;++b)
            for(int d = 0;d < M;++d)
               energy += (*this)(a,c) * ( 2.0 * Hubbard::gV(a,b,c,d) - Hubbard::gV(a,b,d,c) ) * (*this)(b,d);

      }

   return energy;
}

/* vim: set ts=3 sw=3 expandtab :*/
