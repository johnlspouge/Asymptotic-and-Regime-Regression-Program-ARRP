// sls_robust_regression.h
#ifndef INCLUDED_SLS_REGRESSION
#define INCLUDED_SLS_REGRESSION

#include <complex>
#include <iostream>
#include <vector>
#include <algorithm>

#include "sls_utils.h"

using namespace std;

/*****************************************************************************
 
************************** ROBUST REGRESSION ROUTINES ************************
  
*****************************************************************************/


namespace Sls {

	typedef double w_function_type(double x_,void* func_number_);

	

	class robust_regression
	{
	public:


		//Calculaton of a median
		static double median_calculator(
			std::vector<double> x_);


		//W functions for different methods
		static double w_Least_squares(
			double z_,
			void* func_number_);//possible additional parameter

		static double w_Andrews(
			double z_,
			void* func_number_);//possible additional parameters

		static double w_Hampel17A(
			double z_,
			void* func_number_);//possible additional parameters

		static double w_Huber(
			double z_,
			void* func_number_);//possible additional parameters

		static double w_Ramsay(
			double z_,
			void* func_number_);//possible additional parameters

		//Iteration for classical robust regression
		static bool M_estimator(
			w_function_type *w_function_,
			void* func_number_,
			std::vector<std::vector<double> > &X_,
			std::vector<double> &y_,
			double eps_,
			std::vector<double> &res_beta_,
			std::vector<double> &res_beta_errors_,
			double inside_eps_,
			std::vector<double> &med_,
			double &var_,
			std::string flag_,
			//flag="Simple regression"
			//flag="Classical iteration process"
			std::vector<double> &limits_
			);

		//distance between vectors
		static double norma(
			const std::vector<double> &res_beta_,
			const std::vector<double> &previous_beta_);

		//copy of a vector
		static void copy_vector(
			std::vector<double> &res_beta_,
			const std::vector<double> &previous_beta_);

		//Single iteration for classical robust regression
		static bool M_estimator_single_step(
			w_function_type *w_function_,
			void* func_number_,
			const std::vector<std::vector<double> > &X_,
			const std::vector<double> &y_,
			const std::vector<double> &previous_beta_,
			std::vector<double> &res_beta_,
			std::vector<double> &res_beta_errors_,
			double inside_eps_,
			std::vector<double> med_,
			double &var_,
			std::vector<double> &W0_);

		//least square method
		static bool least_square_method(
			const std::vector<std::vector<double> > &X_,
			const std::vector<double> &y_,
			std::vector<double> &res_beta_,
			double inside_eps_);



	};
}
#endif //INCLUDED_SLS_REGRESSION

