#ifndef RSPM_H
#define RSPM_H

#include <iostream>

#include "Matrix.h"

class DIIS;

/**
 * @author Brecht Verstichel
 * @date 14-10-2012\n\n
 * This class RSPM was written for restricted single particle matrices. It inherits from the class Matrix and expands it with
 * specific memberfunction and a knowledge of the nr of sp orbitals and particles. The Fock Matrix and the 1DM are part of this class.
 */
class RSPM : public Matrix {

   /**
    * Output stream operator overloaded, the usage is simple, if you want to print to a file, make an
    * ifstream object and type:\n\n
    * object << spm_p << endl;\n\n
    * For output onto the screen type: \n\n
    * cout << spm_p << endl;\n\n
    * @param output The stream to which you are writing (e.g. cout)
    * @param spm_p de RSPM you want to print
    */
   friend std::ostream &operator<<(std::ostream &output,RSPM &spm_p);

   public:
      
      //constructor
      RSPM();

      //copy constructor
      RSPM(const RSPM &);

      //destructor
      virtual ~RSPM() = default;

      using Matrix::operator=;

      int gN() const;

      int gM() const;

      void unit();

      void construct_Fock(const RSPM &);

      void construct_sp_ham();

      void update(const RSPM &);

      void relax(const DIIS &);

      static void init(int,int,double);

      double energy() const;

      static void setshift(double shift_i);

   private:

      //!dimension of single particle space
      static int M;

      //!nr of particles
      static int N;

      static double shift;

};

#endif

/* vim: set ts=3 sw=3 expandtab :*/
