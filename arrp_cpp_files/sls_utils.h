//sls_utils.h
//****************************************************************************

//******************************GENERAL ROUTINES******************************

//****************************************************************************

#ifndef INCLUDED_SLS_GENERAL_ROUTINES
#define INCLUDED_SLS_GENERAL_ROUTINES

#ifndef _MSC_VER //UNIX program
#else
#define _CRT_SECURE_NO_WARNINGS
#endif

#include <complex>
#include <iostream>
#include <map>
#include <vector>
#include <fstream>
#include <float.h>
#include <algorithm>
#include <errno.h>

using namespace std;

namespace Sls { 

	struct error//struct to handle exceptions
	{
		std::string st;
		error(std::string st_,long int error_code_){st=st_;error_code=error_code_;};
		long int error_code;
	};



	class utils{

		public:

		//rounding double values
		static double round(//returns nearest integer to x_
		const double &x_);

		//Gauss method
		static bool Gauss(
			std::vector<std::vector<double> > A_,//matrix n*(n+1)
			std::vector<double> &x_,//solution
			double inside_eps_,
			std::vector<std::vector<double> > *inv_A_=0);

			//Output of a matrix
			static void print_matrix(
				const std::vector<std::vector<double> > A_);

			//Multiplication of matrices
			static void multiply_matrices(
				const std::vector<std::vector<double> > &A_,
				const std::vector<std::vector<double> > &B_,
				std::vector<std::vector<double> > &res_);

			//Multiplication of a matrix and a vector
			static void multiply_matrix_and_vector(
				const std::vector<std::vector<double> > &A_,
				const std::vector<double> &y_,
				std::vector<double> &res_);

			//Matrix transposition
			static void transpose_matrix(
				const std::vector<std::vector<double> > &A_,
				std::vector<std::vector<double> > &res_);

			static bool the_value_is_double(
				std::string str_,
				double &val_);

			static bool the_value_is_double_old(
				string str_,
				double &val_);


			static bool the_value_is_long(
				std::string str_,
				long int &val_);


		static void assert_mem(void *pointer_)
		{
			if(!pointer_)
			{
				throw error("Memory allocation error\n",41);
			};
		}


		template<class T>
		static inline T Tmax(T i_, T j_)
		{
			if(i_>j_)
			{
				return i_;
			};
			return j_;
		}


		template<class T>
		static inline T Tmin(T i_, T j_)
		{
			if(i_<j_)
			{
				return i_;
			};
			return j_;
		}

	
	};
}

#endif

