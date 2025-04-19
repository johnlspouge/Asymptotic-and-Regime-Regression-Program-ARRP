
//sls_asymptotic_regression.h

//****************************************************************************

//********************ASYMPTOTIC REGRESSION ROUTINES**************************

//****************************************************************************


#ifndef INCLUDED_SLS_ARRP_REGRESSION
#define INCLUDED_SLS_ARRP_REGRESSION

#include "sls_robust_regression.h"
#include "sls_utils.h"

#include <complex>
#include <iostream>
#include <map>
#include <vector>
#include <fstream>
#include <float.h>
#include <algorithm>

using namespace std;

namespace Sls { 


	class asymptotic_regression{

	public:


		//constant asymptotic regression
		void asymptotic_regression_const_beta1_is_defined(
		long int min_length_,//minimum number of data points in the optimal interval
		long int number_of_elements_,//number of data points
		//regression for data points (x_n,y_n) with errors e_n (n numerates data points)
		//x_n are not important for constant regression
		double *values_,//y_n
		double *errors_,//e_n
		bool cut_left_tail_,//cut or not the left side (0 or 1)
		bool cut_right_tail_,//cut or not the right side (0 or 1)
		double y_,//parameter y of the regression; y=sqrt(c)
		//result is ~beta0_+beta1_*x; beta1_ is known beforehand
		double &beta0_,
		double beta1_,//known parameter
		double &beta0_error_,//error of beta0_
		double beta1_error_,//error of beta1_; known parameter
		long int &k1_opt_,//optimal interval is [k1_opt_,k2_opt_]
		long int &k2_opt_,
		bool &res_was_calculated_);//=true if the regression is successful and false otherwise

		//constant asymptotic regression; single step
		double function_for_asymptotic_regression_const_beta1_is_defined(
		double *values_,
		double *errors_,
		long int number_of_elements_,
		long int k_start_,
		double c_,
		double &beta0_,
		double beta1_,
		double &beta0_error_,
		double beta1_error_,
		bool &res_was_calculated_);

		static double median(
		long int dim_,
		double *array_);

		static double median_error(
		long int number_of_elements_,
		double *e_);

		static double standard_error(
		long int number_of_elements_,
		double *values_);

		double standard_error(
		long int number_of_elements_,
		double *args_,
		double *values_);


		//linear asymptotic regression
		void asymptotic_regression_linear(
		long int min_length_,//minimum number of data points in the optimal interval
		long int number_of_elements_,//number of data points
		//regression for data points (x_n,y_n) with errors e_n (n numerates data points)
		double *args_,//x_n
		double *values_,//y_n
		double *errors_,//e_n
		bool cut_left_tail_,//cut or not the left side (0 or 1)
		bool cut_right_tail_,//cut or not the right side (0 or 1)
		double y_,//parameter y of the regression; y=sqrt(c)
		//result is ~beta0_+beta1_*x
		double &beta0_,
		double &beta1_,
		double &beta0_error_,//error of beta0_
		double &beta1_error_,//error of beta1_
		long int &k1_opt_,//optimal interval is [k1_opt_,k2_opt_]
		long int &k2_opt_,
		bool &res_was_calculated_);//=true if the regression is successful and false otherwise

		//linear asymptotic regression; single step
		double function_for_asymptotic_regression_linear(
		double *args_,
		double *values_,
		double *errors_,
		long int number_of_elements_,
		double c_,
		double &beta0_,
		double &beta1_,
		double &beta0_error_,
		double &beta1_error_,
		bool &res_was_calculated_);

	};
}

#endif

