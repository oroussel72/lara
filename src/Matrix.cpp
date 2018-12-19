#include <iostream>
#include <stdlib.h>
#include <math.h>
using namespace std;
#include "Global.h"
#include "Vector.h"
#include "Matrix.h"

/*
_______________________________________________________________________________

Empty constructor
_______________________________________________________________________________

*/

Matrix::Matrix()
{
	_U = new real;
	_lines = _columns = 1;
	*_U = 0.;
}

/*
_______________________________________________________________________________

Constructor with size
_______________________________________________________________________________

*/

Matrix::Matrix(const int n, const int m)
{
	int i;

	_U = new real[n*m];
	
	for (i=0; i<(n*m); i++)
		_U[i] = 0.;

	_lines = n;
	_columns = m;
	
}

/*
_______________________________________________________________________________

Copy constructor
_______________________________________________________________________________

*/

Matrix::Matrix(const Matrix& M)
{
    int i,j;
  
    _U = new real[M.lines()*M.columns()];
    _lines = M._lines;
    _columns = M._columns;
    
    for (i = 1; i <= _lines; i++)
    for (j = 1; j <= _columns; j++)
      set(i,j,M.at(i,j));
}

/*
_______________________________________________________________________________

Distructor
_______________________________________________________________________________

*/

Matrix::~Matrix()
{
  delete[] _U;
}

/*
_______________________________________________________________________________

Reset
_______________________________________________________________________________

*/

void Matrix::reset()
{
	int i,j;

	for (i=1;i<=lines();i++)
	for (j=1;j<=columns();j++)
		set(i, j, 0.);
}

/*
_______________________________________________________________________________

Resize
_______________________________________________________________________________

*/

void Matrix::resize(const int n, const int m)
{
  int i,j;
  int minLines, minColumns;
  
  // Copy Matrix into M
  
  Matrix M(lines(), columns());
  
  for (i = 1; i <= lines(); i++)
  for (j = 1; j <= columns(); j++)
    M.set(i,j,at(i,j));
  
  // Resize matrix and reset values
    
  raise(n,m);
    
  // Copy values of M into local matrix
  
  minLines = min(lines(),M.lines());
  minColumns = min(columns(), M.columns());

  for (i = 1; i <= minLines; i++)
  for (j = 1; j <= minColumns; j++)
    set(i,j,M.at(i,j));
}


/*
_______________________________________________________________________________

Raise matrix and reset values
_______________________________________________________________________________

*/

void Matrix::raise(const int n, const int m)
{
      
  delete [] _U;
  _U = new real[n*m];
  _lines = n;
  _columns = m;
  reset();
  
}

/*
_______________________________________________________________________________

Operator = (assessment)
_______________________________________________________________________________

*/

void Matrix::operator= (const Matrix& M)
{
  int i,j;
  
  if (lines() != M.lines() || columns() != M.columns())
    raise(M.lines(),M.columns());
  
  for (i=1; i<=lines(); i++)
  for (j=1; j<=columns(); j++)
    set(i,j,M.at(i,j));
}

/*
_______________________________________________________________________________

Get a vector from a matrix using a given column
_______________________________________________________________________________

*/

Vector Matrix::vectorColumn (const int column) const
{
	Vector V(lines());
	int i;

	for (i=1; i<=lines();i++)
		V.set(i, at(i,column));
	
	return V;
}


/*
_______________________________________________________________________________

Add a matrix into the same matrix
_______________________________________________________________________________

*/

void Matrix::operator+= (const Matrix& M)
{
	int i,j;

#ifdef debug_mode
	if (lines() != M.lines() || columns() != M.columns())
	{
		cout << "in Matrix.cpp: operator += \n";
		cout << "Matrices have different sizes \n";
		exit(-1);
	}
#endif

	for (i=1; i<=lines(); i++)
	for (j=1; j<=columns(); j++)
		set(i, j, at(i,j) + M.at(i,j));	
}

/*
_______________________________________________________________________________

Add a matrix into another matrix
_______________________________________________________________________________

*/

Matrix Matrix::operator+ (const Matrix& M) const
{
	int i,j;
	Matrix P(lines(),columns());

#ifdef debug_mode
	if (lines() != M.lines() || columns() != M.columns())
	{
		cout << "in Matrix.cpp: operator + \n";
		cout << "Matrices have different sizes \n";
		exit(-1);
	}
#endif

	for (i=1; i<=lines(); i++)
	for (j=1; j<=columns(); j++)
		P.set(i, j, at(i,j) + M.at(i,j));	

	return P;
}

/*
_______________________________________________________________________________

Substract a matrix into the same matrix
_______________________________________________________________________________

*/

void Matrix::operator-= (const Matrix& M)
{
	int i,j;

#ifdef debug_mode
	if (lines() != M.lines() || columns() != M.columns())
	{
		cout << "in Matrix.cpp: operator -= \n";
		cout << "Matrices have different sizes \n";
		exit(-1);
	}
#endif

	for (i=1; i<=lines(); i++)
	for (j=1; j<=columns(); j++)
		set(i, j, at(i,j) - M.at(i,j));	
}

/*
_______________________________________________________________________________

Substract a matrix into another matrix
_______________________________________________________________________________

*/

Matrix Matrix::operator- (const Matrix& M) const
{
	int i,j;
	Matrix P(lines(),columns());

#ifdef debug_mode
	if (lines() != M.lines() || columns() != M.columns())
	{
		cout << "in Matrix.cpp: operator- \n";
		cout << "Matrices have different sizes \n";
		exit(-1);
	}
#endif

	for (i=1; i<=lines(); i++)
	for (j=1; j<=columns(); j++)
		P.set(i, j, at(i,j) - M.at(i,j));	

	return P;
}

/*
_______________________________________________________________________________

Multiply a matrix with a real
_______________________________________________________________________________

*/

void Matrix::operator*= (const real x)
{
	int i,j;

	for (i=1; i<=lines(); i++)
	for (j=1; j<=columns(); j++)
		set(i,j,at(i,j)*x);	
}

/*
_______________________________________________________________________________

Internal product
_______________________________________________________________________________

*/

Matrix Matrix::operator* (const Matrix& M) const
{
	int i,j,k;
	Matrix P(lines(),M.columns());

#ifdef debug_mode
	if (columns() != M.lines())
	{
		cout << "in Matrix.cpp: operator* (internal product) \n";
		cout << "Matrices do not match \n";
		exit(-1);
	}
#endif

	for (i=1; i<=lines(); i++)
	for (j=1; j<=M.columns(); j++)
	for (k=1; k<=columns(); k++)
		P.set(i, j, P.at(i,j) + at(i,k)*M.at(k,j));	

	return P;
}

/*
_______________________________________________________________________________

Matrix-vector product
_______________________________________________________________________________

*/

Vector Matrix::operator* (const Vector& V) const
{
	int i,k;
	Vector W(lines());

#ifdef debug_mode
	if (columns() != V.size())
	{
		cout << "in Matrix.cpp: operator* (matrix-vector product) \n";
		cout << "Matrix and vector do not match \n";
		exit(-1);
	}
#endif

	for (i=1; i<=lines(); i++)
	for (k=1; k<=columns(); k++)
		W.set(i, W.at(i) + at(i,k)*V.at(k));	

	return W;
}

/*
_______________________________________________________________________________

Swap lines1 and lines2
_______________________________________________________________________________

*/

void Matrix::swapLines(const int line1, const int line2)
{
  Vector V(columns());
  int j;
  
  // Store data of line 1 into V
  
  for (j=1; j <= columns(); j++)
    V.set(j,at(line1,j));
  
  // Copy data from line 2 to line 1
    
  for (j=1; j <= columns(); j++)
    set(line1,j,at(line2,j));
  
  // Copy data from line 1 to line 2 via stored data in V 

  for (j=1; j <= columns(); j++)
    set(line2,j,V.at(j));   
}

/*
_______________________________________________________________________________

Swap column1 and colums2
_______________________________________________________________________________

*/

void Matrix::swapColumns(const int column1, const int column2)
{
  Vector V(lines());
  int i;
  
  // Store data of column 1 into V
  
  for (i=1; i <= lines(); i++)
    V.set(i,at(i, column1));
  
  // Copy data from column 2 to column 1
    
  for (i=1; i <= lines(); i++)
    set(i,column1,at(i,column2));
  
  // Copy data from column 1 to column 2 via stored data in V 

  for (i=1; i <= lines(); i++)
    set(i,column2,V.at(i));   
}


/*
_______________________________________________________________________________

EXTERNAL: Multiply a matrix with a real (store into a new matrix)
_______________________________________________________________________________

*/

Matrix operator* (const real x, const Matrix& M)
{
	int i,j;
	Matrix P(M.lines(),M.columns());

	for (i=1; i<=M.lines(); i++)
	for (j=1; j<=M.columns(); j++)
		P.set(i, j, M.at(i,j)*x);

	return P;

}

/*
_______________________________________________________________________________

EXTERNAL: Absolute value of a matrix
_______________________________________________________________________________

*/

Matrix fabs(const Matrix& M)
{
  int i,j;
  Matrix P(M.lines(),M.columns());
  
  for (i=1; i<=M.lines(); i++)
  for (j=1; j<=M.columns(); j++)
    P.set(i, j, fabs(M.at(i,j)));

  return P;
  
}
