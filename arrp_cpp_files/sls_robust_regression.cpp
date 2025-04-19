// sls_robust_regression.cpp

#include "sls_robust_regression.h"

using namespace std;
using namespace Sls;

double robust_regression::median_calculator(
std::vector<double> x_)
{
	sort(x_.begin(),x_.end());
	long int x_size=(long int)x_.size();
	if(x_size%2==1)
	{
		long int tmp=(long int)utils::round((double)(x_size+1)/2.0)-1;
		return x_[tmp];
	};
	long int tmp=(long int)utils::round((double)(x_size)/2.0)-1;
	return 0.5*(x_[tmp]+x_[tmp+1]);

}

double robust_regression::w_Andrews(
double z_,
void* func_number_)//possible additional parameters
{
	double Pi=3.14159265358979323846;
	double a=*((double*)func_number_);
	if(a==0)
	{
		cout<<"Error in parameters of robust_regression::w_Andrews\n";
		return 0;
	};
	if(z_==0)
	{
		return 1.0;
	};
	if(fabs(z_)<=a*Pi)
	{
		return sin(z_/a)/(z_/a);
	};
	return 0;

}

double robust_regression::w_Huber(
double z_,
void* func_number_)//possible additional parameters

{
	double a=*((double*)func_number_);
	if(a==0)
	{
		cout<<"Error in parameters of robust_regression::w_Huber\n";
		return 0;
	};
	if(fabs(z_)<=a)
	{
		return 1.0;
	};
	return a/fabs(z_);
}

double robust_regression::w_Ramsay(
double z_,
void* func_number_)//possible additional parameters

{
	double a=*((double*)func_number_);
	if(a==0)
	{
		cout<<"Error in parameters of robust_regression::w_Ramsay\n";
		return 0;
	};
	return exp(-a*fabs(z_));
}

bool robust_regression::M_estimator(
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
)
{
	if(!(flag_=="Simple regression"||
		flag_=="Classical iteration process"))
	{
		cout<<"Error in robust_regression::M_estimator\n";
		return false;
	};
	std::vector<double> previous_beta_;
	//res_beta_.clear();

if(flag_=="Classical iteration process")
{
	if(res_beta_.size()==0)
	{
	if(!least_square_method(
		X_,
		y_,
		res_beta_,
		inside_eps_))
	{
		return false;
	};
	};
};

if(flag_=="Simple regression")
{
	if(!least_square_method(
		X_,
		y_,
		res_beta_,
		inside_eps_))
	{
		return false;
	};
	return true;
};

		long int step=0;


		std::vector<double> W0;
		double norma_tmp;
		do
		{
			step++;
			if(step>100)
			{
				break;
			};
			copy_vector(previous_beta_,res_beta_);
			if(!M_estimator_single_step(
				w_function_,
				func_number_,
				X_,
				y_,
				previous_beta_,
				res_beta_,
				res_beta_errors_,
				inside_eps_,
				med_,
				var_,
				W0))
				{
					return false;
				};


			norma_tmp=norma(res_beta_,previous_beta_);

			{
				limits_.clear();
				long int k_left=0;
				while(W0[k_left]==0)
				{
					k_left++;
					if(k_left>=(long int)W0.size())
					{
						break;
					};
				};

				long int k_right=(long int)W0.size()-1;
				while(W0[k_right]==0)
				{
					k_right--;
					if(k_left<=-1)
					{
						break;
					};
				};

				if(k_left>k_right)
				{
					cout<<"Unexpected error in defining ranges in robust regression\n"<<endl;
					return false;
				};
				limits_.resize(2);
				limits_[0]=k_left;
				limits_[1]=k_right;

			};

			long int k;
			for(k=0;k<(long int)W0.size();k++)
			{
				limits_.push_back(W0[k]);
			};

		}
 		while(norma_tmp>=eps_);

		return true;
}

void robust_regression::copy_vector(
std::vector<double> &res_beta_,
const std::vector<double> &previous_beta_)
{
	res_beta_.clear();
	res_beta_.resize(previous_beta_.size());

	long int i;
	for(i=0;i<(long int)previous_beta_.size();i++)
	{
		res_beta_[i]=previous_beta_[i];
	};

}

double robust_regression::norma(
const std::vector<double> &res_beta_,
const std::vector<double> &previous_beta_)
{
	if(res_beta_.size()!=previous_beta_.size())
	{
		cout<<"Error in robust_regression::norma\n";
		return 0;
	};
	long int i;
	double tmp=0;
	for(i=0;i<(long int)res_beta_.size();i++)
	{
		tmp+=fabs(res_beta_[i]-previous_beta_[i]);
	};
	return tmp;
}

bool robust_regression::M_estimator_single_step(
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
std::vector<double> &W0_)
{
	long int size1=(long int)X_.size();
	if(size1==0)
	{
		cout<<"Error in parameters of robust_regression::M_estimator\n";
		return false;
	};
	W0_.clear();
	W0_.resize(size1);
	long int size2=(long int)X_[0].size();
	long int i,j;
	for(i=1;i<size1;i++)
	{
		if((long int)X_[i].size()!=size2)
		{
			cout<<"Error in parameters of robust_regression::M_estimator\n";
			return false;
		};
	};

	if((long int)y_.size()!=size1)
	{
		cout<<"Error in parameters of robust_regression::M_estimator\n";
		return false;
	};

	res_beta_.clear();
	res_beta_.resize(size2+1);
	std::vector<std::vector<double> > X_expanded(size1);
	for(i=0;i<size1;i++)
	{
		X_expanded[i].resize(size2+1);
	};
	
	for(i=0;i<size1;i++)
	{
		X_expanded[i][0]=1;
		for(j=1;j<=size2;j++)
		{
			X_expanded[i][j]=X_[i][j-1];
		};
	};

	std::vector<std::vector<double> > X_expanded_transposed;
	utils::transpose_matrix(X_expanded,X_expanded_transposed);

	std::vector<double> e;
	utils::multiply_matrix_and_vector(X_expanded,previous_beta_,e);
	for(i=0;i<size1;i++)
	{
		e[i]=y_[i]-e[i];
	};

	double median=median_calculator(e);
	std::vector<double> e_tmp(size1);
	for(i=0;i<size1;i++)
	{
		e_tmp[i]=fabs(e[i]-median)/0.6745;
	};
	double s=median_calculator(e_tmp);

	if(s==0)
	{
		if(med_.size()!=y_.size())
		{
			cout<<"Median is equal to 0 in robust_regression::M_estimator_single_step\n";
			for(i=0;i<size2;i++)
			{
				res_beta_[i]=previous_beta_[i];
			};
			return false;
		};
	};

	if(med_.size()!=y_.size())
	{
		med_.resize(y_.size());
		for(i=0;i<(long int)y_.size();i++)
		{
			med_[i]=s;
		};
	};

	e_tmp.clear();

	std::vector<std::vector<double> > W0(size1);
	for(i=0;i<size1;i++)
	{
		W0[i].resize(size1,0);
		if(fabs(e[i])<inside_eps_)
		{
			W0[i][i]=1.0;
		}
		else
		{
			W0[i][i]=(*w_function_)(e[i]/med_[i],func_number_)/(med_[i]*med_[i]);
		};

		W0_[i]=W0[i][i];
	};

	std::vector<std::vector<double> > tmp;
	std::vector<std::vector<double> > A;
	std::vector<double> b;

	utils::multiply_matrices(X_expanded_transposed,W0,tmp);
	utils::multiply_matrices(tmp,X_expanded,A);

	utils::multiply_matrix_and_vector(tmp,y_,b);

	for(i=0;i<size2+1;i++)
	{
		(A[i]).push_back(b[i]);
	};
	b.clear();

	std::vector<std::vector<double> > *inv_A = new std::vector<std::vector<double> >;


	if(!utils::Gauss(A,//matrix n*(n+1)
			res_beta_,//solution
			inside_eps_,
			inv_A))
	{
		delete inv_A;
		return false;
	};

	//calculation of variance matrix
	std::vector<std::vector<double> > beta_matrix;
	utils::multiply_matrices(*inv_A,tmp,beta_matrix);

	std::vector<double> tmp_zero(med_.size(),0);
	std::vector<std::vector<double> > median_matrix(med_.size(),tmp_zero);

	for(i=0;i<(long int)med_.size();i++)
	{
		median_matrix[i][i]=med_[i]*med_[i];
	};

	std::vector<std::vector<double> > tmp3;

	utils::multiply_matrices(beta_matrix,median_matrix,tmp3);

	std::vector<std::vector<double> > beta_matrix_transposed;
	utils::transpose_matrix(beta_matrix,beta_matrix_transposed);

	std::vector<std::vector<double> > var_matrix;
	utils::multiply_matrices(tmp3,beta_matrix_transposed,var_matrix);

	res_beta_errors_.resize(res_beta_.size());
	for(i=0;i<(long int)res_beta_.size();i++)
	{
		if(var_matrix[i][i]<=0)
		{
			res_beta_errors_[i]=0;
		}
		else
		{
			res_beta_errors_[i]=sqrt(var_matrix[i][i]);
		};
	};


	if(res_beta_.size()==1)
	{
		double num=0;
		for(i=0;i<size1;i++)
		{
			num+=W0[i][i];
		};

		double den=0;
		for(i=0;i<size1;i++)
		{
			den+=med_[i]*med_[i]*W0[i][i]*W0[i][i];
		};
		var_=sqrt(den/(num*num));

	};

	delete inv_A;
	return true;
}

bool robust_regression::least_square_method(
const std::vector<std::vector<double> > &X_,
const std::vector<double> &y_,
std::vector<double> &res_beta_,
double inside_eps_)
{
	long int size1=(long int)X_.size();
	if(size1==0)
	{
		cout<<"Error in parameters of robust_regression::least_square_method\n";
		return false;
	};
	long int size2=(long int)X_[0].size();
	long int i,j;
	for(i=1;i<size1;i++)
	{
		if((long int)X_[i].size()!=size2)
		{
			cout<<"Error in parameters of robust_regression::least_square_method\n";
			return false;
		};
	};

	if((long int)y_.size()!=size1)
	{
		cout<<"Error in parameters of robust_regression::least_square_method\n";
		return false;
	};

	res_beta_.clear();
	res_beta_.resize(size2+1);
	std::vector<std::vector<double> > X_expanded(size1);
	for(i=0;i<size1;i++)
	{
		X_expanded[i].resize(size2+1);
	};
	
	for(i=0;i<size1;i++)
	{
		X_expanded[i][0]=1;
		for(j=1;j<=size2;j++)
		{
			X_expanded[i][j]=X_[i][j-1];
		};
	};


	std::vector<std::vector<double> > X_expanded_transposed;
	utils::transpose_matrix(X_expanded,X_expanded_transposed);

	std::vector<std::vector<double> > tmp;
	std::vector<std::vector<double> > A;
	std::vector<double> b;

	utils::multiply_matrices(X_expanded_transposed,X_expanded,A);

	utils::multiply_matrix_and_vector(X_expanded_transposed,y_,b);
	tmp.clear();
	X_expanded.clear();
	X_expanded_transposed.clear();

	for(i=0;i<size2+1;i++)
	{
		(A[i]).push_back(b[i]);
	};
	b.clear();


	if(!utils::Gauss(A,//matrix n*(n+1)
			res_beta_,//solution
			inside_eps_))
	{
		return false;
	};
	
	return true;
}

