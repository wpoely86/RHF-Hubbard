/**
 * @mainpage 
 * This is a implementation of ThING using Libint for the electron repulsion integrals, and using faster
 * methods for the calculation of the S, T and U matrices.
 * @author Brecht Verstichel, Sebastian Wouters
 * @date 08-06-2012
 */

#include <iostream>
#include <fstream>
#include <cmath>
#include <getopt.h>
#include <libint2.h>

using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;

//includes all important headers and defines which conditions are
//going to be used:
#include "include.h"

/**
 * Here all the libraries come together to do the actual calculations.
 */
int main(int argc, char **argv){

   //initialize the random nr generator
   srand(time(NULL));

   cout.precision(10);

   input::init("start.stp");

   LibInt::init();

   CI_SPM::init();
   CI_TPM::init();
   CI_SPPM::init();
   CI_SPPM_m::init();
   CI_TPPM::init();

   CartInt::init();

   SI_SPM::init();
   SI_TPM::init();

   SphInt::init();

   Tools::init();

   //here the cartesian integrals are calculated
   CartInt::calc_integrals();

   //have to be normalised before transforming to spherical
   CartInt::orthogonalize();

   const int M = CI_SPM::gdim();//dim of spin restricted sp hilbert space
   const int N = input::NumberOfElectrons();//nr of particles

   RSPM::init(M,N); 

   //initialize spm
   RSPM spm;
   spm.unit();

   //some variables needed in scf loop
   RSPM copy,F;
   Vector v(M);

   DIIS diis;

   RSPM commutator;

   double convergence = 1.0;

   int iter = 0;

   //basic scf loop, easy peasy
   while(convergence > 1.0e-14){

      ++iter;

      //construct the fock matrix
      F.construct_Fock(spm);

      //add it to the diis object
      diis.push_F(F);

      //calculate the commutator between F and the spm
      commutator.commute(F,spm);

      //add it to diis object
      diis.push_comm(commutator);

      //construct the B matrix and solve the system: output = b coefficients
      diis.construct();

      //now construct the "relaxed" Fock matrix
      F.relax(diis);

      v.diagonalize(F);

      copy = spm;

      spm.update(F);

      copy -= spm;

      convergence = std::sqrt(copy.ddot(copy));

      cout << iter << "\t" << convergence << endl;
   
   }

   cout << endl;
   cout << "single-particle spectrum is:" << endl;
   cout << endl;

   F.construct_Fock(spm);
   v.diagonalize(F);

   cout << v << endl;

   F.construct_sp_ham();

   double energy = 2.0*F.ddot(spm);

   for(int n = 0;n < N/2;++n)
      energy += 2.0*v[n];

   cout << "Ground-state energy:" << endl;
   cout << endl;

   cout << "E_0 = \t" << energy/2.0 + input::gNucRepEn() << endl;

   //here all the static classes are removed
   Tools::clear();

   SphInt::clear();

   SI_TPM::clear();
   SI_SPM::clear();

   CartInt::clear();

   CI_TPPM::clear();
   CI_SPPM_m::clear();
   CI_SPPM::clear();
   CI_TPM::clear();
   CI_SPM::clear();

   input::clear();
   LibInt::clear();

   return 0;

}
