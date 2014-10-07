#include <iostream>
#include <algorithm>
#include <time.h>
#include <cmath>
#include <assert.h>

#include "Vector.h"
#include "lapack.h"

/**
 * constructor that takes dimension as input
 * @param n dimension of the vector
 */
Vector::Vector(int n)
{
   vector.resize(n);
}

/**
 * Construct and initialize the Vector object by diagonalizing a Matrix object:
 */
Vector::Vector(Matrix &matrix)
{
   vector.resize(matrix.gn());

   diagonalize(matrix);
}

/**
 * copy constructor 
 * @param vec_copy The vector you want to be copied into the object you are constructing
 */
Vector::Vector(const Vector &vec_copy)
{
   vector = vec_copy.vector;
}

Vector::Vector(Vector &&vec_copy)
{
   vector = std::move(vec_copy.vector);
}


/**
 * overload the equality operator
 * @param vector_copy The vector you want to be copied into this
 */
Vector &Vector::operator=(const Vector &vec_copy)
{
   assert(vec_copy.vector.size() == vector.size());
   vector = vec_copy.vector;

   return *this;
}

Vector &Vector::operator=(Vector &&vec_copy)
{
   assert(vec_copy.vector.size() == vector.size());

   vector = std::move(vec_copy.vector);

   return *this;
}

/**
 * Make all the number in your vector equal to the number a, e.g. usefull for initialization (Vector M = 0)
 * @param a the number
 */
Vector &Vector::operator=(double a)
{
   for(unsigned int i = 0;i < vector.size();++i)
      vector[i] = a;

   return *this;
}

/**
 * overload the += operator for matrices
 * @param vector_pl The vector you want to add to this
 */
Vector &Vector::operator+=(const Vector &vector_pl)
{
   daxpy(1,vector_pl);

   return *this;
}

/**
 * overload the -= operator for matrices
 * @param vector_pl The vector you want to deduct from this
 */
Vector &Vector::operator-=(const Vector &vector_pl)
{
   daxpy(-1,vector_pl);

   return *this;
}

/**
 * add the vector vector_pl times the constant alpha to this
 * @param alpha the constant to multiply the vector_pl with
 * @param vector_pl the Vector to be multiplied by alpha and added to this
 */
Vector &Vector::daxpy(double alpha,const Vector &vector_pl)
{
   int inc = 1;
   int n = vector.size();

   daxpy_(&n,&alpha,vector_pl.vector.data(),&inc,vector.data(),&inc);

   return *this;
}

/**
 * /= operator overloaded: divide by a constant
 * @param c the number to divide your vector through
 */
Vector &Vector::operator/=(double c)
{
   dscal(1.0/c);

   return *this;
}

/**
 * *= operator overloaded: divide by a constant
 * @param c the number to divide your vector through
 */
Vector &Vector::operator*=(double c)
{
   dscal(c);

   return *this;
}

/**
 * write access to your vector, change the number on index i
 * @param i row number
 * @return the entry on place i
 */
double &Vector::operator[](int i)
{
   assert(i<vector.size());

   return vector[i];
}

/**
 * read access to your vector, change the number on index i: const version
 * @param i row number
 * @return the entry on place i
 */
double Vector::operator[](int i) const
{
   assert(i<vector.size());

   return vector[i];
}

/**
 * Diagonalize the Matrix matrix when you have allready allocated the memory of the vector
 * on the correct dimension.
 */
void Vector::diagonalize(Matrix &matrix)
{
   *this = matrix.diagonalize();
}

/**
 * @return the underlying pointer to vector, useful for mkl and lapack applications
 */
double *Vector::gVector()
{
   return vector.data();
}

/**
 * @return the underlying pointer to vector, useful for mkl and lapack applications: const version
 */
const double *Vector::gVector() const
{
   return vector.data();
}

/**
 * @return the dimension of the vector
 */
unsigned int Vector::gn() const
{
   return vector.size();
}

/**
 * @return the sum of all the elements in the vector
 */
double Vector::sum() const
{
   double ward = 0;

   for(auto &elem : vector)
      ward += elem;

   return ward;
} 

double Vector::trace() const
{
   return sum();
}

/**
 * @return the logarithm of the product of all the elements in the vector (so the sum of all the logarithms)
 */
double Vector::log_product() const
{
   double ward = 0;

   for(auto &elem : vector)
      ward += log(elem);

   return ward;
}

/**
 * @return inproduct of (*this) vector with vector_i
 * @param vector_i input vector
 */
double Vector::ddot(const Vector &vector_i) const
{
   int inc = 1;
   int n = vector.size();

   return ddot_(&n,vector.data(),&inc,vector_i.vector.data(),&inc);
}

/**
 * Scale the vector (*this) with parameter alpha
 * @param alpha scalefactor
 */
void Vector::dscal(double alpha)
{
   int inc = 1;
   int n = vector.size();

   dscal_(&n,&alpha,vector.data(),&inc);
}

/**
 * Fill the vector with random numbers.
 */
void Vector::fill_Random()
{
   srand(time(NULL));

   for(auto &elem : vector)
      elem = (double) rand()/RAND_MAX;
}

std::ostream &operator<<(std::ostream &output,const Vector &vector_p)
{
   for(int i = 0;i < vector_p.gn();++i)
      output << i << "\t" << vector_p[i] << std::endl;

   return output;
}

/**
 * @return the minimal element present in this Vector object.
 * watch out, only works when Vector is filled with the eigenvalues of a diagonalized Matrix object
 */
double Vector::min() const
{
   double min = vector[0];

   for(unsigned int i=1;i<gn();i++)
      if(vector[i]<min)
         min = vector[i];

   return min;
}

/**
 * @return the maximal element present in this Vector object.
 * watch out, only works when Vector is filled with the eigenvalues of a diagonalized Matrix object
 */
double Vector::max() const
{
   double max = vector[gn()-1];

   for(int i=gn()-2;i>=0;i--)
      if(vector[i]>max)
         max = vector[i];

   return max;
}

void Vector::invert()
{
   for(int i = 0;i < gn();++i)
      vector[i] = 1.0/vector[i];
}

void Vector::sqrt(int option)
{
   if(option == 1)
      for(int i = 0;i < gn();++i)
         vector[i] = std::sqrt(vector[i]);
   else 
      for(int i = 0;i < gn();++i)
         vector[i] = 1.0/std::sqrt(vector[i]);
}

void Vector::symmetrize()
{
}

Vector &Vector::mprod(const Vector &x,const Vector &y)
{
   assert(x.gn() == gn() && y.gn() == gn());

   for(int i = 0;i < gn();++i)
      vector[i] = x[i] * y[i];

   return *this;
}

/**
 * Sort the vector, small to large
 */
void Vector::sort()
{
   std::sort(vector.begin(), vector.end());
}

/* vim: set ts=3 sw=3 expandtab :*/
