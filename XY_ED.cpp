/*
Title: Eigenstates and Reduced Density Matrices of the XY Model

Author: Reyhaneh Khasseh

Description: 

The purpose of this project is to evaluate the eigenstates of the XY model Hamiltonian, which is defined as:

H = \sum_{ i = 1 }^{ N } [ ( 1 + gamma ) / 2 \sigma_i^x  \sigma_{ i + 1 }^x +  ( 1 - gamma ) / 2 \sigma_i^y  \sigma_{ i + 1 }^y ] + hz * \sum_{i=1}^{N} \sigma_i^z

The values of gamma that are commonly used in the literature are gamma = 0 and gamma = 1, which correspond to the XX and Ising models, respectively.

The obtained eigenstates can then be used to evaluate the corresponding reduced density matrix and subsystem configuration probabilities.

*/

#include<iostream>
#include<stdlib.h>
#include<cmath>
#include <stdio.h>
#include "mkl_lapacke.h"
#include<iomanip>

using namespace std;

/*defined functions*/

int* decToBinary( int n , int N ); //convert decimal number to binary

void Hamiltonian( int N , int d , double*& H , double JXX , double JYY , double hz ); //Hamiltonian 

void RDM_func( int d_r , int ns , double* V_in , double*& rho ); //Reduced density matrix of given subsystem

/*main code starts here*/

int main(){

  /*defined parameters*/	

  MKL_INT info_re;   

  int N = 2 , d = pow( 2 , N ); //N: systems size, d: Hilbert space dimension

  int N_r = 2 , d_r = pow( 2 , N_r ) , ns = pow( 2 , N - N_r ); //N_r: subsystem size, d_r, n_s: Hilbert space dimension of subsystem, rest of the system

  double gamma = 1 , JXX = ( 1. + gamma ) / 4. , JYY = ( 1. - gamma ) / 4. , hz = -0.5; //Hamiltonian tunable parameters

  int sum_h , num1 , num2 , S_t;

  /*defined arrays*/

  int* arr = new int[ N ]; //binary representation of configurations

  double* H = new double [ d * d ]; //Hamiltonian

  double* w_re = new double[ d ]; //eigenenergies

  double* V_in = new double [ d ]; //eigenstate

  double* rho = new double [ d_r * d_r ]; //reduced density matrix

  /*Hamiltonian function*/

  Hamiltonian( N , d , H , JXX , JYY , hz );

  /*ED eigensolver*/

  info_re = LAPACKE_dsyev( LAPACK_ROW_MAJOR , 'V' , 'U', d , H , d , w_re );

  for( int j = 0; j < 1; j++ ){ //j-th eigenstate loop

     cout << setprecision( 10 ) << w_re[ j ] <<endl; //eigenvalue of j-th eigenstate

     /*configuration probabilities of j-th eigenstate*/

     for( int i = 0 ; i < d ; i++ ){

       arr = decToBinary( i, N );

       cout << setprecision(5) <<pow( H[i*d+j] , 2  ) <<"  ";

       V_in[ i ] = H[ j + d * i ];

     } //end of i loop

     cout << endl;

     /*Reduced density matrix*/

     RDM_func( d_r , ns , V_in , rho );

     /*configuration probabilities of subsystem*/

    int iter=0;
    for( int jd = 0; jd < d_r; jd ++ ){

      double pr=rho[ iter++ + jd * d_r ];

      cout << pr << "  ";

    }

    cout << endl;

  } //end of j-th eigenstates loop


  return 0;

}  

int* decToBinary( int n, int N ){

  /* This function convert decimal numbers to binary
   *
   * Arguments:
   *
   * n -- Decimal number
   * N -- Length of the binary array
   * binaryNum[N] -- The array of binary
   *
   * Return:
   *
   * The binary array - Note: here we replace 0 with -1
   *
   */

  int* binaryNum = new int [N];

  for( int i = 0 ; i < N ; i++ )  binaryNum[ i ] = -1;

  //counter for binary array
  int i = 0;

  while (n > 0) {

    binaryNum[ N - 1 - i ] = 2 * ( n % 2 ) - 1;
    n = n / 2;
    i++;

  }

  return binaryNum;

}

void Hamiltonian(int N , int d , double*& H , double JXX , double JYY , double hz){

  /* This function construct the Hamiltonian of XY spin chain
   *
   * Arguments:
   *
   * N -- number of spins
   * d -- Hilbert space dimension
   * H -- d*d matrix of Hamiltonian
   * JXX, JYY -- interaction strength 
   * hz -- transverse field
   *
   *
   */	

  int* arr = new int[ N ];
  int num1 , num2;
  double sum_hz;

  for( int j = 0; j < d * d; j++ )  H[ j ] = 0.0;

  /*Intraction terms*/

  for( int j = 0; j < d; j++ ){//row loop

    arr = decToBinary( j, N );

    for( int l = 0; l < N; l++ ){//column loop

      num1 = j - arr[ l ] * pow( 2 , N - 1 - l );
      num2 = num1 - arr[ ( l + 1 ) % N ] * pow(2 ,  N - 1 - ( ( l + 1 )  % N));

      H[ num2 + j * d ] = - ( JXX - arr[ l ] * arr[ ( l + 1 ) % N ] * JYY );

    }//end of column loop

  }//end of row loop


  /*transverse field*/

  int it = 0;
  int it1 = 0;

  for( int j=0; j < d; j++ ){

    arr = decToBinary( j, N );

    sum_hz = 0.0;

    for( int l = 0; l < N; l++ ){

      sum_hz += arr[ l ];

    }

    H[ j * d + it++ ] =  0.5 * hz * sum_hz + H[ j * d + it1++ ];

  }//end of j


}	

void RDM_func(int d_r, int ns, double* V_in, double*& rho){

  /* This function evaluate reduced density matrix for corresponding pure state
   *
   * Arguments:
   *
   * d_r -- Dimension of reduced density matrix
   * ns -- Hilbert space dimension of subsystem we trace out
   * V_in -- The array of eigenstate vector
   * rho -- d_r*d_r array of reduced density matrix 
   *
   *
   */	

  int iter = 0;

  for( int k = 0 ; k < d_r ; k++ ){

    for( int l = 0 ; l < d_r; l++ ){

      double sum_vv = 0.0;

      for( int ii = 0; ii < ns; ii++ ){

        for( int j = 0; j < ns; j++ ){

          if(ii == j) sum_vv += V_in[ ii + k * ns ]*V_in[ j + l * ns ];

        }// end of j loop

      }// end of ii loop

      rho[ iter++ ] = sum_vv;

    }// end of l loop

  }// end of k loop

}
