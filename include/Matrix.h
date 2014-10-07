#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <cstdlib>
#include <memory>
#include <vector>

/**
 * @author Brecht Verstichel
 * @date 18-02-2010\n\n
 * This is a class written for symmetric matrices. It is a wrapper around a double pointer and
 * redefines much used lapack and blas routines as memberfunctions
 */

class Vector;

class Matrix
{
   /**
    * Output stream operator overloaded, the usage is simple, if you want to print to a file, make an
    * ifstream object and type:\n\n
    * object << matrix << endl;\n\n
    * For output onto the screen type: \n\n
    * cout << matrix << endl;\n\n
    * @param output The stream to which you are writing (e.g. cout)
    * @param matrix_p de Matrix you want to print
    */
   friend std::ostream &operator<<(std::ostream &output,Matrix &matrix_p);

   public:

      //constructor
      Matrix(int n);

      //copy constructor
      Matrix(const Matrix &);

      Matrix(Matrix &&);

      //destructor
      virtual ~Matrix();

      //overload equality operator
      Matrix &operator=(const Matrix &);

      Matrix &operator=(Matrix &&);

      Matrix &operator=(double );

      //overload += operator
      Matrix &operator+=(const Matrix &);

      //overload -= operator
      Matrix &operator-=(const Matrix &);

      Matrix &daxpy(double alpha,const Matrix &);

      Matrix &operator*=(double);

      Matrix &operator/=(double);

      Matrix &mprod(const Matrix &,const Matrix &);

      //easy to change the numbers
      double &operator()(int i,int j);

      //easy to access the numbers
      double operator()(int i,int j) const;

      //get the pointer to the matrix
      double *gMatrix();

      const double *gMatrix() const;

      int gn() const;

      double trace() const;

      Vector diagonalize();

      Vector diagonalize_2x2();

      double ddot(const Matrix &) const;

      void invert();

      void invert_2x2();

      void dscal(double alpha);

      void fill_Random();

      void fill_Random(int seed);

      //positieve of negatieve vierkantswortel uit de matrix
      void sqrt(int option);

      void sqrt_2x2(int option);

      void mdiag(const Vector &);

      void L_map(const Matrix &,const Matrix &);

      void L_map_2x2(const Matrix &,const Matrix &);

      void symmetrize();

      void SaveToFile(const std::string) const;

      void ReadFromFile(const std::string);

      void sep_pm(Matrix &,Matrix &);

      void sep_pm_2x2(Matrix &,Matrix &);

      Matrix & commute(const Matrix &,const Matrix &);

      void solve(std::vector<double> &);


   private:

      //!pointer of doubles, contains the numbers, the matrix
      std::unique_ptr<double []> matrix;

      //!dimension of the matrix
      int n;
};

#endif

/* vim: set ts=3 sw=3 expandtab :*/
