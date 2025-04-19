//sls_asymptotic_regression.cpp

#include "sls_asymptotic_regression.h"

using namespace Sls;

void asymptotic_regression::asymptotic_regression_linear(
long int min_length_,
long int number_of_elements_,
double *args_,
double *values_,
double *errors_,
bool cut_left_tail_,
bool cut_right_tail_,
double y_,
double &beta0_,
double &beta1_,
double &beta0_error_,
double &beta1_error_,
long int &k1_opt_,
long int &k2_opt_,
bool &res_was_calculated_)
{

	if(number_of_elements_<2)
	{
		throw error("Unexpected error\n",4);
	};

	//minization of the function
	double c=y_*y_;


	long int k1_start,k1_end;
	long int k2_start,k2_end;

	if(cut_left_tail_&&cut_right_tail_)
	{
		k1_start=0;
		k1_end=number_of_elements_-1;

		k2_start=0;
		k2_end=number_of_elements_-1;

	}
	else
	{
		if(cut_left_tail_&&!cut_right_tail_)
		{
			k1_start=0;
			k1_end=number_of_elements_-1;

			k2_start=number_of_elements_-1;
			k2_end=number_of_elements_-1;

		}
		else
		{
			if(!cut_left_tail_&&cut_right_tail_)
			{
				k1_start=0;
				k1_end=0;

				k2_start=0;
				k2_end=number_of_elements_-1;
			}
			else
			{
				k1_start=0;
				k1_end=0;

				k2_start=number_of_elements_-1;
				k2_end=number_of_elements_-1;

			};
		};
	};
	
	long int k1_opt=0,k2_opt=0;

	double func_opt=DBL_MAX;
	double beta0_opt=0;
	double beta1_opt=0;
	double beta0_opt_error=0;
	double beta1_opt_error=0;

	long int k1,k2;

	
	res_was_calculated_=false;


	for(k1=k1_start;k1<=k1_end;k1++)
	{

		for(k2=utils::Tmax(k1+(long int)1,utils::Tmax(k1,k2_start));k2<=k2_end;k2++)
		{
			if(k2-k1+1<min_length_)
			{
				continue;
			};

			double beta0_opt_tmp,beta1_opt_tmp,beta0_opt_error_tmp,beta1_opt_error_tmp;
			bool res_was_calculated;


			double tmp=function_for_asymptotic_regression_linear(
				args_+k1,
				values_+k1,
				errors_+k1,
				k2-k1+1,
				c,
				beta0_opt_tmp,
				beta1_opt_tmp,
				beta0_opt_error_tmp,
				beta1_opt_error_tmp,
				res_was_calculated);



			if(tmp<func_opt&&res_was_calculated)
			{
				func_opt=tmp;
				beta0_opt=beta0_opt_tmp;
				beta1_opt=beta1_opt_tmp;
				beta0_opt_error=beta0_opt_error_tmp;
				beta1_opt_error=beta1_opt_error_tmp;
				k1_opt=k1;
				k2_opt=k2;
				res_was_calculated_=true;
			};

		};
	};

	if(res_was_calculated_)
	{
		beta0_=beta0_opt;
		beta1_=beta1_opt;
		beta0_error_=beta0_opt_error;
		beta1_error_=beta1_opt_error;
		k1_opt_=k1_opt;
		k2_opt_=k2_opt;
	};

}

double asymptotic_regression::standard_error(
long int number_of_elements_,
double *args_,
double *values_)
{
			long int i;
			double x_mean=0;
			double y_mean=0;

			for(i=0;i<number_of_elements_;i++)
			{
				y_mean+=values_[i];
				x_mean+=args_[i];
			};

			x_mean/=(double)number_of_elements_;
			y_mean/=(double)number_of_elements_;

			double xy=0;
			double xx=0;
			for(i=0;i<number_of_elements_;i++)
			{
				xy+=(args_[i]-x_mean)*(values_[i]-y_mean);
				xx+=(args_[i]-x_mean)*(args_[i]-x_mean);
			};

			if(xx==0)
			{
				return -1.0;
			};

			double slope=xy/xx;

			double *errors_=new double[number_of_elements_];
			utils::assert_mem(errors_);

			for(i=0;i<number_of_elements_;i++)
			{
				errors_[i]=values_[i]-slope*args_[i];
			};

			double sd=asymptotic_regression::standard_error(number_of_elements_,errors_);
			delete[]errors_;
			return sd;

}

double asymptotic_regression::function_for_asymptotic_regression_linear(
double *args_,
double *values_,
double *errors_,
long int number_of_elements_,
double c_,
double &beta0_,
double &beta1_,
double &beta0_error_,
double &beta1_error_,
bool &res_was_calculated_)
{

	long int i;
	double a11=0;
	double a12=0;
	double a21=0;
	double a22=0;
	double y1=0;
	double y2=0;

	double y1_error=0;
	double y2_error=0;

	std::vector<std::vector<double> > B(2,std::vector<double>(number_of_elements_,0));

	for(i=0;i<number_of_elements_;i++)
	{
		if(errors_[i]!=0)
		{
			double tmp=1.0/(errors_[i]*errors_[i]);

			a11+=tmp;
			a12+=args_[i]*tmp;
			a22+=args_[i]*args_[i]*tmp;
			y1+=values_[i]*tmp;
			B[0][i]=tmp;
			y1_error+=tmp*tmp*errors_[i]*errors_[i];
			y2+=args_[i]*values_[i]*tmp;
			B[1][i]=args_[i]*tmp;
			y2_error+=args_[i]*args_[i]*tmp*tmp*errors_[i]*errors_[i];
		};
	};

	a21=a12;
	y1_error=sqrt(y1_error);
	y2_error=sqrt(y2_error);

	double eps=1e-10*utils::Tmax(fabs(a11*a22),fabs(a21*a12));

	double den=a11*a22-a21*a12;
	if(fabs(den)<=eps)
	{
		res_was_calculated_=false;
		return 0;
	}
	else
	{
		res_was_calculated_=true;
	};

	beta0_=(y1*a22-a12*y2)/den;
	beta1_=(a11*y2-a21*y1)/den;

	
	beta0_error_=sqrt(y1_error*y1_error*a22*a22+a12*a12*y2_error*y2_error)/den;
	beta1_error_=sqrt(a11*a11*y2_error*y2_error+a21*a21*y1_error*y1_error)/den;

	//new error calculation
	{
		std::vector<std::vector<double> > A(2,std::vector<double>(2,0));
		A[0][0]=a22/den;
		A[0][1]=-a12/den;
		A[1][0]=-a21/den;
		A[1][1]=a11/den;


		std::vector<std::vector<double> > B_matrix;
		utils::multiply_matrices(A,B,B_matrix);
		std::vector<std::vector<double> > B_matrix_transpose;
		utils::transpose_matrix(B_matrix,B_matrix_transpose);


		std::vector<double> tmp_zero(number_of_elements_,0);
		std::vector<std::vector<double> > error_matrix(number_of_elements_,tmp_zero);

		long int i;
		for(i=0;i<number_of_elements_;i++)
		{
			error_matrix[i][i]=errors_[i]*errors_[i];
		};

		std::vector<std::vector<double> > tmp3;
		utils::multiply_matrices(B_matrix,error_matrix,tmp3);

		std::vector<std::vector<double> > var_matrix;
		utils::multiply_matrices(tmp3,B_matrix_transpose,var_matrix);


		if(var_matrix[0][0]<=0)
		{
			beta0_error_=0;
		}
		else
		{
			beta0_error_=sqrt(var_matrix[0][0]);
		};
	
		if(var_matrix[1][1]<=0)
		{
			beta1_error_=0;
		}
		else
		{
			beta1_error_=sqrt(var_matrix[1][1]);
		};

	};

	double res=0;
	for(i=0;i<number_of_elements_;i++)
	{
		if(errors_[i]!=0)
		{
			double tmp=(beta0_+beta1_*args_[i]-values_[i])/errors_[i];
			res+=tmp*tmp-c_;
		};
	};

	return res;
}

double asymptotic_regression::median_error(
long int number_of_elements_,
double *e_)
{
	double m;

	double *e_m=new double[number_of_elements_];
	utils::assert_mem(e_m);
	
	long int j;

	m=asymptotic_regression::median(number_of_elements_,e_);

	for(j=0;j<number_of_elements_;j++)
	{
		e_m[j]=fabs(e_[j]-m)/0.6745;
	};

	double res=asymptotic_regression::median(number_of_elements_,e_m);

	delete[]e_m;

	return res;
}

double asymptotic_regression::standard_error(
long int number_of_elements_,
double *values_)
{
	long int i;
	double sd=0; 
	double mean=0;
	for(i=0;i<number_of_elements_;i++)
	{
		sd+=values_[i]*values_[i];
		mean+=values_[i];
	};

	sd/=(double)number_of_elements_;
	mean/=(double)number_of_elements_;
	sd-=mean*mean;
	sd=sqrt(sd);


	return sd;
}


//arbitrary x

void asymptotic_regression::asymptotic_regression_const_beta1_is_defined(
long int min_length_,
long int number_of_elements_,
double *values_,
double *errors_,
bool cut_left_tail_,
bool cut_right_tail_,
double y_,
double &beta0_,
double beta1_,
double &beta0_error_,
double beta1_error_,
long int &k1_opt_,
long int &k2_opt_,
bool &res_was_calculated_)
{
	//minization of the function
	double c=y_*y_;


	long int k1_start,k1_end;
	long int k2_start,k2_end;

	if(cut_left_tail_&&cut_right_tail_)
	{
		k1_start=0;
		k1_end=number_of_elements_-1;

		k2_start=0;
		k2_end=number_of_elements_-1;

	}
	else
	{
		if(cut_left_tail_&&!cut_right_tail_)
		{
			k1_start=0;
			k1_end=number_of_elements_-1;

			k2_start=number_of_elements_-1;
			k2_end=number_of_elements_-1;

		}
		else
		{
			if(!cut_left_tail_&&cut_right_tail_)
			{
				k1_start=0;
				k1_end=0;

				k2_start=0;
				k2_end=number_of_elements_-1;
			}
			else
			{
				k1_start=0;
				k1_end=0;

				k2_start=number_of_elements_-1;
				k2_end=number_of_elements_-1;

			};
		};
	};
	
	long int k1_opt=0,k2_opt=0;

	double func_opt=DBL_MAX;
	double beta0_opt=0;
	double beta1_opt=beta1_;
	double beta0_opt_error=-1;
	double beta1_opt_error=beta1_error_;

	long int k1,k2;

	res_was_calculated_=false;


	for(k1=k1_start;k1<=k1_end;k1++)
	{

		for(k2=utils::Tmax(k1,k2_start);k2<=k2_end;k2++)
		{
			if(k2-k1+1<min_length_)
			{
				continue;
			};

			double beta0_opt_tmp,beta1_opt_tmp,beta0_opt_error_tmp,beta1_opt_error_tmp;
			bool res_was_calculated;

			beta1_opt_tmp=beta1_;
			beta1_opt_error_tmp=beta1_error_;

			double tmp=function_for_asymptotic_regression_const_beta1_is_defined(
				values_+k1,
				errors_+k1,
				k2-k1+1,
				k1,
				c,
				beta0_opt_tmp,
				beta1_opt_tmp,
				beta0_opt_error_tmp,
				beta1_opt_error_tmp,
				res_was_calculated);



			if(tmp<func_opt&&res_was_calculated)
			{
				func_opt=tmp;
				beta0_opt=beta0_opt_tmp;
				beta1_opt=beta1_opt_tmp;
				beta0_opt_error=beta0_opt_error_tmp;
				beta1_opt_error=beta1_opt_error_tmp;
				k1_opt=k1;
				k2_opt=k2;
				res_was_calculated_=true;
			};

		};
	};


	if(res_was_calculated_)
	{
		beta0_=beta0_opt;
		beta0_error_=beta0_opt_error;
		k1_opt_=k1_opt;
		k2_opt_=k2_opt;
	};

	
}

double asymptotic_regression::function_for_asymptotic_regression_const_beta1_is_defined(
double *values_,
double *errors_,
long int number_of_elements_,
long int k_start_,
double c_,
double &beta0_,
double beta1_,
double &beta0_error_,
double beta1_error_,
bool &res_was_calculated_)
{

	long int i;
	double a11=0;
	double y1=0;


	double y1_error=0;


	for(i=0;i<number_of_elements_;i++)
	{
		if(errors_[i]!=0)
		{
			double tmp=1.0/(errors_[i]*errors_[i]);

			a11+=tmp;
			y1+=(values_[i]-(double)(k_start_+i)*beta1_)*tmp;
			double error_tmp=errors_[i]*errors_[i]+(double)(k_start_+i)*(double)(k_start_+i)*beta1_error_*beta1_error_;
			y1_error+=tmp*tmp*error_tmp;
		};
	};

	y1_error=sqrt(y1_error);

	double eps=1e-10*fabs(a11);

	double den=a11;
	if(fabs(den)<=eps)
	{
		res_was_calculated_=false;
		return 0;
	}
	else
	{
		res_was_calculated_=true;
	};

	beta0_=y1/den;

	beta0_error_=y1_error/den;


	double res=0;
	for(i=0;i<number_of_elements_;i++)
	{
		if(errors_[i]!=0)
		{
			double tmp=(beta0_+beta1_*(i+k_start_)-values_[i])/errors_[i];
			res+=tmp*tmp-c_;
		};
	};

	return res;

}

double asymptotic_regression::median(
long int dim_,
double *array_)
{
	vector<double> array_vect(dim_);
	long int i;
	for(i=0;i<dim_;i++)
	{
		array_vect[i]=array_[i];
	};
	sort(array_vect.begin(),array_vect.end());
	if(dim_%2==0)
	{
		long int k=(long int)utils::round((double)dim_/2.0);
		return 0.5*(array_vect[k-1]+array_vect[k]);
	}
	else
	{
		long int k=(long int)utils::round((double)(dim_-1.0)/2.0);
		return array_vect[k];

	};
}

