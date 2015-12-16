////////////////////////////////////////////////////////////////////////////////////
// Copyright © Charalambos "Charis" Poullis, charalambos@poullis.org    	  //
// This work can only be used under an exclusive license of the author.           //
////////////////////////////////////////////////////////////////////////////////////

#ifndef __MATRIX_H__
#define __MATRIX_H__

#include <stdarg.h>
#include <math.h>
#include <memory.h>
#include <iostream>
#include "Vector.h"
#include <malloc.h>
#include <opencv2/opencv.hpp>

#define MATRIX_EPSILON	1e-06

/**		Matrix
* This is a templated class for matrix manipulations.
*/

template <typename T, int Y, int X>
class Matrix	{

	public:
			///Constructor
			Matrix<T,Y,X>()	{Zero();}
			//NOTE: PAY ATTENTION!!! THIS WORKS IF ALL NUMBERS ENTERED ARE OF THE SAME TYPE
			///Constructor
// 			Matrix<T,Y,X>(T first, ...)	{
// 				va_list list;
// 				//Initialize the va_list i.e T variable 'list' by the address
// 				//of first unknown variable argument by a call of va_start() macro
// 				va_start(list,first);
// 				//In loop, retrieve each argument
// 				//second argument to va_arg is datatype
// 				//of expected argument
// 				matrix[0][0] = (T) first;
// 				for (int i=1;i<Y;i++)	{
// 					for (int j=0;j<X;j++)	{
// 						T value = va_arg(list,T);
// 						matrix[i][j] = (T) value;
// 					}
// 				}
// 				va_end(list);
// 			}
			///Constructor
			Matrix<T,Y,X>(Matrix<T,Y,X> const &mat) {memcpy(matrix,mat.matrix,sizeof(T)*X*Y);}
			///Constructor
			Matrix<T,Y,X>(T * array)		{memcpy(matrix,array,sizeof(T)*X*Y);}

			///Destructor
			~Matrix()	{}

			///Operators
			T operator() ( int elementNumberY, int elementNumberX)	const	{return(matrix[elementNumberY][elementNumberX]);}
			T& operator() ( int elementNumberY, int elementNumberX)		{return(matrix[elementNumberY][elementNumberX]);}
			Matrix<T,Y,X> &operator+= (Matrix<T,Y,X> const &mat)	{
				for (int i=0;i<Y;i++)	{
					for (int j=0;j<X;j++)	{
						matrix[i][j] += mat.matrix[i][j];
					}
				}
				return (*this);
			}
			Matrix<T,Y,X> &operator-= ( Matrix<T,Y,X> const &mat)	{
				for (int i=0;i<Y;i++)	{
					for (int j=0;j<X;j++)	{
						matrix[i][j] -= mat.matrix[i][j];
					}
				}
				return (*this);
			}
			Matrix<T,Y,X> &operator*= ( T scalar)	{
				for (int i=0;i<Y;i++)	{
					for (int j=0;j<X;j++)	{
						matrix[i][j] *= scalar;
					}
				}
				return (*this);
			}
			Matrix<T,Y,X> &operator/= ( T scalar)	{
				for (int i=0;i<Y;i++)	{
					for (int j=0;j<X;j++)	{
						matrix[i][j] /= scalar;
					}
				}
				return (*this);
			}


			///Methods
			Matrix<T,Y,X> getTranspose() const {
               Matrix<T,Y,X> ret(*this);
				ret.transpose();
				return ret;
			}

			Matrix<T,Y,X> getInverse() const {
				 Matrix<T,Y,X> ret(*this);
				ret.invert();
				return ret;
			}

			cv::Mat toCvMat() const {

				Matrix<T,Y,X> tmp(*this);
				cv::Mat m(Y,X,CV_64F);

				for(int i=0; i<Y; i++)
					for(int j=0; j<X; j++)
						m.at<double>(i,j) = tmp(i,j);
				return m;
			}

			void Zero(void)	{for (int i=0;i<Y;i++)
								for (int j=0;j<X;j++)
									matrix[i][j] = (T) 0.0;}
			void identity(void)	{for (int i=0;i<Y;i++)
									for (int j=0;j<X;j++)
										if (i==j) matrix[i][j] = (T) 1.0;}
			void normalize()	{
				T denom = (T) 0.0;
				for (int y=0;y<Y;y++)	{
					for (int x=0;x<X;x++)	{
						denom += matrix[y][x]*matrix[y][x];
					}
				}
				T k;
				if (denom > T(MATRIX_EPSILON)) k = T(1.0 /sqrt(double(denom)));
				else k = T(1.0);
				for (int y=0;y<Y;y++)	{
					for (int x=0;x<X;x++)	{
						matrix[y][x] *= k;
					}
				}
				return;
			}

			T Determinant(void)	{std::cout << "Matrix.h:GetDeterminant: not yet done\n" << std::endl;return 0.0;}
			T *Array(void)		{
								 T *mat=new T[Y*X];
								 for (int i=0;i<Y;i++)
									for (int j=0;j<X;j++)
										mat[i*X+j]= matrix[i][j];
								 return mat;
								  /*return matrix*/;}

			int transpose()	{
				 if (Y!=X)    {
					std::cout << "Not a square matrix." << std::endl;
					return -1;
				}
				for (int i = 0; i < Y; i++)	{
					for (int j = 0; j < i; j++)	{
						T temp = matrix[i][j];
						matrix[i][j] = matrix[j][i];
						matrix[j][i] = temp;
					}
				}
				return 1;
			}

			int invert()	{
				if (Y!=X)	{
					std::cout << "Not a square matrix." << std::endl;
					return -1;
				}
				double 	aext,                   /* Extreme (largest) value in sub-matrix */
					atemp,                  /* Holds temp. value of the matrix       */
					de = 1.0;               /* Initial determate value.              */
				int     ir[ 64],           /* Use to remember swaps pos. in rows    */
					ic[ 64 ],          /* Use to remember swaps pos. in columns */
					i,j,k,kk,               /* Ioop control variables                */
					itemp,idim,             /*                                       */
					iext,jext,              /* I extreme, J extreme for swapping pos */
					flag = 0;               /* Boolean test for finding correct rows */
								/*      and columns in step 5.           */
								/*---------------------------------------*/
				int size = Y;
				/* INITIALIZE SWAP ARRAYS TO CORRECT POSITIONS */
				for( j=0; j < size; j++ ) {
					ic[ j ] = j;
					ir[ j ] = j;
				}
				/* STEP 1 */

				/* FIND INDEXES  (jext and iext) OF LARGEST ELEMENT IN REMAINING MATRIX */
				for( k=0; k < size; k++ ) {
					aext = 0.0;
					for( i=k; i < size; i++ ) {
						for( j=k; j < size; j++ ) {
							if( aext < fabs( matrix[i][j] ) ) {
								iext = i;
								jext = j;
								aext = fabs( matrix[i][j] );
							}
						}
					}

					if( aext <= 0.0 ){
						printf("\nInverse_matrix; singular\n");
						//*det = 0;
						return(-100);
						/*exit( -20 );
						*/
					}

					/* STEP 2 */

					if( k != iext ) {
						/*  MOVE COLUMN WITH BIGGEST ELEMENT TO ROW K. */
						de = - de;
						for( j=0; j < size; j++ ) {
							atemp		      = matrix[k][j];
							matrix[k][j]    = matrix[iext][j];
							matrix[iext][j] = atemp;
						}
						itemp    = ic[k];
						ic[k]    = ic[iext];
						ic[iext] = itemp;
					}

					if( k != jext ) {
						/*  MOVE ROW WITH BIGGEST ELEMENT TO COLUMN K. */
						de = - de;
						for( i=0; i < size; i++ ) {
							atemp	      	      = matrix[i][k];
							matrix[i][k]     = matrix[i][jext];
							matrix[i][jext] = atemp;
						}
						itemp    = ir[k];
						ir[k]    = ir[jext];
						ir[jext] = itemp;
					}


					/* STEP 3 */
					/* PERFORM GAUSS-JORDAN ELIMINATION WITH CURRENT ROW K. */

					aext = matrix[k][k];
					de = de * aext;
					matrix[k][k] = 1.0;
					for( j=0; j < size; j++ )
						matrix[k][j] = matrix[k][j] / aext;
						for( i=0; i < size; i++ ) {
							if( k != i ) {
									aext = matrix[i][k];
									if( aext != 0.0 ) {
										matrix[i][k] = 0.0;
										for( j=0; j < size; j++ )
											matrix[i][j] = matrix[i][j] -aext * matrix[k][j];
									}
							}
						}
				}


				/* STEP 5 */
				/* INVERSE ELEMENTS HAVE BEEN CALCULATED RETURN
					MATRIX TO ORIGINAL ROW, COLUMN ORDER;         */

				idim = size - 1;
				for( k=0; k < idim; k++ ) {
					kk = k + 1;
					if( k != ic[k] ) {

				/*  LOOP UNTIL FIND CORRECT POSITON IN SWAP ARRAY.
					IF CORRECT POSITION IS NOT FOUND; ERROR OCCURED. */

						i = kk;
						flag = 0;
						while((i < size) && (!flag)) {
							if( k == ic[i] ) flag = 1;
							else i++;
						}

						if( flag ) {
							for( j=0; j < size; j++ ) {
								atemp                  = matrix[j][k];
								matrix[j][k] = matrix[j][i];
								matrix[j][i] = atemp;
							}
							itemp = ic[i];
							ic[i] = ic[k];
							ic[k] = itemp;
						}
						else {
							printf("\nInverse_matrix: 1 suspected machine error.\n");
							printf(" k = %d  ic[ %d ] = %d \n",k,i,ic[i]);
							/*exit( -15 ) ;
							*/
						}
					}
					if( k != ir[k] ) {

				/*  LOOP UNTIL FIND CORRECT POSITON IN SWAP ARRAY.
					IF CORRECT POSITION IS NOT FOUND; ERROR OCCURED. */

						j = kk;
						flag = 0;
						while((j < size) && (!flag)) {
							if( k == ir[j] ) flag = 1;
							else j++;
						}

						if( flag ) {
								for( i=0; i < size; i++ ) {
								atemp 	           = matrix[k][i];
								matrix[k][i] = matrix[j][i];
								matrix[j][i] = atemp;
								}
								itemp = ir[j];
								ir[j] = ir[k];
								ir[k] = itemp;
						}
						else {
								printf(" \nInverse_matrix: 2 suspected machine error.\n");
								printf(" k = %d  ic[ %d ] = %d \n",k,i,ic[i]);
								/*exit( -15 );
								*/
						}
					}
				}

				return( 1 );
			}

/*
			Recursive definition of determinate using expansion by minors.
*/
			double determinant(T **in_a=0x00,int in_n=-1)
			{
				int n;
				T a[Y][X];
				if (in_a==0x00 && in_n == -1)	{
					if (Y!=X)	{
						std::cerr << "Not a square matrix" << std::endl;
						return -666.0;
					}
					n = Y;
					for (unsigned int i=0;i<n;i++)	{
						for (unsigned int j=0;j<n;j++)	{
							a[i][j] = matrix[i][j];
						}
					}
				}
				else	{
					n = in_n;
					for (unsigned int i=0;i<n;i++)	{
						for (unsigned int j=0;j<n;j++)	{
							a[i][j] = in_a[i][j];
						}
					}
				}

				int i,j,j1,j2;
				double det = 0.0;
				T **m = NULL;

				if (n < 1) { /* Error */

				} else if (n == 1) { /* Shouldn't get used */
					det = (double) a[0][0];
				} else if (n == 2) {
					det =(double) a[0][0] * a[1][1] - a[1][0] * a[0][1];
				} else {
					det = 0.0;
					for (j1=0;j1<n;j1++) {
						m = (T **) malloc((n-1)*sizeof(T *));
						for (i=0;i<n-1;i++)
							m[i] = (T *) malloc((n-1)*sizeof(T));
						for (i=1;i<n;i++) {
							j2 = 0;
							for (j=0;j<n;j++) {
								if (j == j1)
									continue;
								m[i-1][j2] = a[i][j];
								j2++;
							}
						}
						det += pow(-1.0,1.0+j1+1.0) * double(a[0][j1]) * determinant(m,n-1);
						for (i=0;i<n-1;i++)
							free(m[i]);
						free(m);
					}
				}
				return(det);
			}

	private:
			//structure
			T matrix[Y][X];

			//Friend functions
			friend std::ostream &operator<< (std::ostream &os, Matrix<T,Y,X> const &matrix)	{
				for (int i=0;i<Y;i++)	{
					for (int j=0;j<X;j++)
						os<< T(matrix(i,j)) << ' ';
					os << "\n" ;
				}
				return os;
			}
			friend std::istream &operator>> (std::istream &is, Matrix<T,Y,X> const &matrix)	{
				int i,j;
				T temp;
				for (i=0;i<Y;i++)	{
					for (j=0;j<X;j++)
						is >>temp;	matrix.matrix[i][j] = temp;
				}
				return is;
			}
			friend Matrix<T,Y,X> operator- (Matrix<T,Y,X> const &mat)	{
				Matrix<T,Y,X> new_mat;
				for (int i=0;i<Y;i++)	{
					for (int j=0;j<X;j++)
						new_mat.matrix[i][j] = -mat.matrix[i][j];
				}
				return new_mat;
			}
			friend Matrix<T,Y,X> operator+ ( Matrix<T,Y,X> const &matrix1,  Matrix<T,Y,X> const &matrix2)	{
				Matrix<T,Y,X> result;
				for (int i=0;i<Y;i++)	{
					for (int j=0;j<X;j++)	{
						result.matrix[i][j] = matrix1.matrix[i][j] + matrix2.matrix[i][j];
					}
				}
				return (result);
			}
			friend Matrix<T,Y,X> operator- ( Matrix<T,Y,X> const &matrix1,  Matrix<T,Y,X> const &matrix2)	{
				Matrix<T,Y,X> result;
				for (int i=0;i<Y;i++)
					for (int j=0;j<X;j++)
						result.matrix[i][j] = matrix1.matrix[i][j] - matrix2.matrix[i][j];
				return (result);
			}
			template<int Z>
			friend Matrix<T,Y,Z> operator* ( Matrix<T,Y,X> const &matrix1,  Matrix<T,X,Z> const &matrix2)	{
				Matrix<T,Y,Z> result;
				result.Zero();
				for (int i=0;i<Y;i++)	{
					for (int j=0;j<Z;j++)	{
						for (int k=0;k<X;k++)	{
							result(i,j) += (T) matrix1(i,k) * matrix2(k,j);
						}
					}
				}
				return (result);
			}
			friend Vector<T,Y> operator* ( Matrix<T,Y,X> const &mat,  Vector<T,X> const &vec)	{
				Vector<T,Y> result;
				for (int i=0;i<Y;i++)
					for (int j=0;j<X;j++)
						result(i) += mat.matrix[i][j] * vec(j);
				return result;
			}
			friend Matrix<T,Y,X> operator* ( Matrix<T,Y,X> const &mat,  T scalar)	{
				Matrix<T,Y,X> result;
				for (int i=0;i<Y;i++)
					for (int j=0;j<X;j++)
						result.matrix[i][j] = mat.matrix[i][j] * scalar;
				return (result);
			}
			friend Matrix<T,Y,X> operator* ( T scalar,  Matrix<T,Y,X> const &mat)	{
				Matrix<T,Y,X> result;
				for (int i=0;i<Y;i++)
					for (int j=0;j<X;j++)
						result.matrix[i][j] = mat.matrix[i][j] *scalar;

				return (result);
			}
			friend Matrix<T,Y,X> operator/ ( Matrix<T,Y,X> const &mat,  T scalar)	{
				Matrix<T,Y,X> result;
				for (int i=0;i<Y;i++)	{
					for (int j=0;j<X;j++)	{
						result.matrix[i][j] = mat.matrix[i][j] / scalar;
					}
				}
				return (result);

			}
			friend bool operator== ( Matrix<T,Y,X>const &matrix1,  Matrix<T,Y,X> const &matrix2)	{
				for (int i=0;i<Y;i++)
					for (int j=0;j<X;j++)
						if (fabs(matrix1.matrix[i][j] - matrix2.matrix[i][j]) > MATRIX_EPSILON)
							return false;
				return true;
			}

};

///Instantiate the classes
typedef Matrix<int,2,2>				Matrix2i;
typedef Matrix<unsigned int,2,2>		Matrix2ui;
typedef Matrix<float,2,2>			Matrix2f;
typedef Matrix<double,2,2>			Matrix2d;

typedef Matrix<int,3,3>				Matrix3i;
typedef Matrix<unsigned int,3,3>		Matrix3ui;
typedef Matrix<float,3,3>			Matrix3f;
typedef Matrix<double,3,3>			Matrix3d;

typedef Matrix<int,4,4>				Matrix4i;
typedef Matrix<unsigned int,4,4>		Matrix4ui;
typedef Matrix<float,4,4>			Matrix4f;
typedef Matrix<double,4,4>			Matrix4d;

typedef Matrix<float,5,5>			Matrix5f;
typedef Matrix<float,6,6>			Matrix6f;
typedef Matrix<float,7,7>			Matrix7f;

typedef Matrix<double,3,4>			Matrix3x4d;

typedef Matrix<double,1,2>			Matrix1x2d;
typedef Matrix<double,1,4>			Matrix1x4d;

typedef Matrix<double,4,1>			Matrix4x1d;
typedef Matrix<double,4,2>			Matrix4x2d;


template<class T> static Matrix<T,3,3> J(Vector<T,3> const &v)	{
  Matrix<T,3,3> m;

  m(0,0) =  (T) 0.0f;   m(0,1) = (T) -v(2);  m(0,2) = (T)  v(1);
  m(1,0) =  (T) v(2);  	m(1,1) =  (T) 0.0f;  m(1,2) = (T) -v(0);
  m(2,0) =  (T) -v(1);  m(2,1) = (T) v(0);   m(2,2) =  (T) 0.0f;

  return m;
}


#endif
