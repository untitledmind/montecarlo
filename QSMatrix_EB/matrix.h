#ifndef __QS_MATRIX_H
#define __QS_MATRIX_H

/* need C++11 or higher to compile */

#include <vector>
#include <iostream>
#include <cmath>			//YS
#include <complex>			//YS
#include <initializer_list> //EB
using namespace std;

template <typename T> 
class QSMatrix {
 private:
 vector<vector<T> > mat;
  unsigned rows;
  unsigned cols;
  T d;	//determinant

 public:
    // QSMatrix() {}
  QSMatrix(unsigned _rows = 0, unsigned _cols = 0, const T _initial = 0); // EB added default arguments
//  QSMatrix(unsigned _rows, unsigned _cols, const T* _initial);
  QSMatrix(const QSMatrix<T>& rhs);
  QSMatrix(initializer_list<initializer_list<T> >, unsigned _rows = 0, unsigned _cols = 0, const T _initial = 0);   // EB
  virtual ~QSMatrix();

  // Operator overloading, for "standard" mathematical matrix operations
  QSMatrix<T>& operator=(const QSMatrix<T>& rhs);

  // Matrix mathematical operations
  QSMatrix<T> operator+(const QSMatrix<T>& rhs);
  QSMatrix<T>& operator+=(const QSMatrix<T>& rhs);
  QSMatrix<T> operator-(const QSMatrix<T>& rhs);
  QSMatrix<T>& operator-=(const QSMatrix<T>& rhs);
  QSMatrix<T> operator*(const QSMatrix<T>& rhs);
  QSMatrix<T>& operator*=(const QSMatrix<T>& rhs);
  bool operator== (const QSMatrix<T>& rhs);     // EB
  bool circa(const QSMatrix<T>& rhs, const T& tolerance);           // EB

//Auxiliary Matrix functions
  QSMatrix<T> transpose();
  void swapRows(int i, int j);
  void Permute(int col);						//permute rows of matrix in descending order
  void Pivot(int row, int col);						//create next privot element
  QSMatrix<T> PermutationMatrix(int Size, int row1, int row2);		//creates a permutation matrix
  int maxCol(int col);
  QSMatrix<T> Inverse();
  QSMatrix<T> block_multiply(const QSMatrix<T>& rhs);       // new EB function
  T Determinant();			//compute the determinant of this matrix
  T determinant(QSMatrix<T>& rhs, int row, int col);
  vector<T> gauss_jordan(const char& pivoting);
  vector<T> jacobi(const T& tolerance, const int& max_iter = 100);
  vector<T> gauss_seidel(const T& tolerance, const int& max_iter = 100);

  // Matrix/scalar operations
  QSMatrix<T> operator+(const T& rhs);
  QSMatrix<T> operator-(const T& rhs);
  QSMatrix<T> operator*(const T& rhs);
  QSMatrix<T> operator/(const T& rhs);

  // Matrix/vector operations
  std::vector<T> operator*(const std::vector<T>& rhs);
  std::vector<T> diag_vec();

  // Access the individual elements
  T& operator()(const unsigned& row, const unsigned& col);
  const T& operator()(const unsigned& row, const unsigned& col) const;
  vector<T>& operator[](const unsigned& el);        // EB
  const vector<T>& operator[](const unsigned& el) const;    // EB

  //set the individual elements
  void set_value(int i, int j, const T& value){mat[i][j] = value;}

  // Access the row and column sizes
  unsigned get_rows() const;
  unsigned get_cols() const;

  template<typename U>
  friend ostream& operator<<(ostream &os, QSMatrix<U> rhs);
};

#include "matrix.cpp"

#endif
