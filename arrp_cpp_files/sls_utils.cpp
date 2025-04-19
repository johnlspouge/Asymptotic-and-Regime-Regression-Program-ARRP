// sls_utils.cpp

#include "sls_utils.h"

using namespace std;
using namespace Sls;

double utils::round(//returns nearest integer to x_
const double &x_)
{
	double x_floor=floor(x_);
	double x_ceil=ceil(x_);
	if(fabs(x_-x_floor)<0.5)
	{
		return x_floor;
	};
	return x_ceil;
}

bool utils::Gauss(
std::vector<std::vector<double> > A_,//matrix n*(n+1)
std::vector<double> &x_,//solution
double inside_eps_,
std::vector<std::vector<double> > *inv_A_)
{
	

	long int i,j,jj;
	long int matr_size=(long int)A_.size();
	if(matr_size==0)
	{
		cout<<"Error in robust_regression::Gauss - sizes of matrix are wrong\n";
		return false;
	};

	std::vector<std::vector<double> > E;

	if(inv_A_)
	{
		vector<double> zero(matr_size,0);
		(*inv_A_).resize(matr_size, zero);
		
		E.resize(matr_size, zero);
		long int i;
		for(i=0;i<matr_size;i++)
		{
			E[i][i]=1.0;
		};
	};


	for(i=0;i<matr_size;i++)
	{
		if((long int)A_[i].size()!=matr_size+1)
		{
			cout<<"Error in Statistics::Gauss - sizes of matrix are wrong\n";
			return false;
		};
	};
	x_.clear();
	x_.resize(matr_size);
	//forward trace
	for(j=0;j<matr_size;j++)
	{
		long int absmax=j;
		for(i=j+1;i<matr_size;i++)
		{
			if(fabs(A_[absmax][j])<fabs(A_[i][j]))
			{
				absmax=i;
			};
		};

		if(j!=absmax)
		{
			for(jj=j;jj<matr_size+1;jj++)
			{
				double tmp=A_[absmax][jj];
				A_[absmax][jj]=A_[j][jj];
				A_[j][jj]=tmp;
			};

			if(inv_A_)
			{
				for(jj=0;jj<matr_size;jj++)
				{
					double tmp=E[absmax][jj];
					E[absmax][jj]=E[j][jj];
					E[j][jj]=tmp;
				};
			};
		};

		if(fabs(A_[j][j])<=inside_eps_)
		{
			cout<<"Error in utils::Gauss - matrix is singular\n";
			print_matrix(A_);
			x_.clear();
			return false;
		};

		for(i=j+1;i<matr_size;i++)
		{
			double tmp=A_[i][j]/A_[j][j];
			for(jj=j+1;jj<matr_size+1;jj++)
			{
				A_[i][jj]=A_[i][jj]-tmp*A_[j][jj];
			};

			if(inv_A_)
			{
				for(jj=0;jj<matr_size;jj++)
				{
					E[i][jj]=E[i][jj]-tmp*E[j][jj];
				};
			};
		};
	};

	//reverse trace
	x_[matr_size-1]=A_[matr_size-1][matr_size]/A_[matr_size-1][matr_size-1];
	for(i=matr_size-2;i>=0;i--)
	{
		x_[i]=A_[i][matr_size];
		for(j=i+1;j<matr_size;j++)
		{
			x_[i]-=A_[i][j]*x_[j];
		};
		x_[i]/=A_[i][i];
	};

	if(inv_A_)
	{
		long int k;
		for(k=0;k<matr_size;k++)
		{
			(*inv_A_)[matr_size-1][k]=E[matr_size-1][k]/A_[matr_size-1][matr_size-1];
			long int i;
			for(i=matr_size-2;i>=0;i--)
			{
				(*inv_A_)[i][k]=E[i][k];
				long int j;
				for(j=i+1;j<matr_size;j++)
				{
					(*inv_A_)[i][k]-=A_[i][j]*(*inv_A_)[j][k];
				};
				(*inv_A_)[i][k]/=A_[i][i];
			};
		};

	};

	return true;
}

void utils::print_matrix(
const std::vector<std::vector<double> > A_)
{
	long int i,j;
	for(i=0;i<(long int)A_.size();i++)
	{
		for(j=0;j<(long int)A_[i].size();j++)
		{
			if(j<(long int)A_[i].size()-1)
			{
				cout<<A_[i][j]<<"\t";
			}
			else
			{
				cout<<A_[i][j]<<"\n";
			};
		};
	};

}

void utils::multiply_matrices(
const std::vector<std::vector<double> > &A_,
const std::vector<std::vector<double> > &B_,
std::vector<std::vector<double> > &res_)
{
	long int size1=(long int)A_.size();
	if(size1==0)
	{
		cout<<"Error in robust_regression::multiply_matrices\n";
		return;
	};
	long int size2=(long int)A_[0].size();
	if(size2==0)
	{
		cout<<"Error in robust_regression::multiply_matrices\n";
		return;
	};

	long int i,j,k;
	for(i=1;i<size1;i++)
	{
		if((long int)A_[i].size()!=size2)
		{
			cout<<"Error in robust_regression::multiply_matrices\n";
			return;
		};
	};

	if(size2!=(long int)B_.size())
	{
		cout<<"Error in robust_regression::multiply_matrices\n";
		return;
	};

	long int size3=(long int)B_[0].size();
	if(size3==0)
	{
		cout<<"Error in robust_regression::multiply_matrices\n";
		return;
	};

	for(i=1;i<size2;i++)
	{
		if((long int)B_[i].size()!=size3)
		{
			cout<<"Error in robust_regression::multiply_matrices\n";
			return;
		};
	};

	res_.clear();
	res_.resize(size1);
	for(i=0;i<size1;i++)
	{
		res_[i].resize(size3,0);
	};

	for(i=0;i<size1;i++)
	{
		for(j=0;j<size3;j++)
		{
			for(k=0;k<size2;k++)
			{
				res_[i][j]+=A_[i][k]*B_[k][j];
			};
		};
	};
}

void utils::multiply_matrix_and_vector(
const std::vector<std::vector<double> > &A_,
const std::vector<double> &y_,
std::vector<double> &res_)
{
	long int size1=(long int)A_.size();
	if(size1==0)
	{
		cout<<"Error in robust_regression::multiply_matrices\n";
		return;
	};
	long int size2=(long int)A_[0].size();
	if(size2==0)
	{
		cout<<"Error in robust_regression::multiply_matrices\n";
		return;
	};

	long int i,k;
	for(i=1;i<size1;i++)
	{
		if((long int)A_[i].size()!=size2)
		{
			cout<<"Error in robust_regression::multiply_matrices\n";
			return;
		};
	};

	if(size2!=(long int)y_.size())
	{
		cout<<"Error in robust_regression::multiply_matrices\n";
		return;
	};


	res_.clear();
	res_.resize(size1,0);

	for(i=0;i<size1;i++)
	{
		for(k=0;k<size2;k++)
		{
			res_[i]+=A_[i][k]*y_[k];
		};
	};
}

void utils::transpose_matrix(
const std::vector<std::vector<double> > &A_,
std::vector<std::vector<double> > &res_)
{
	long int size1=(long int)A_.size();
	if(size1==0)
	{
		res_.clear();
		return;
	};
	long int size2=(long int)A_[0].size();

	long int i,j;
	for(i=1;i<size1;i++)
	{
		if((long int)A_[i].size()!=size2)
		{
			cout<<"Error in robust_regression::multiply_matrices\n";
			return;
		};
	};

	res_.clear();
	res_.resize(size2);
	for(i=0;i<size2;i++)
	{
		res_[i].resize(size1);
	};

	for(i=0;i<size2;i++)
	{
		for(j=0;j<size1;j++)
		{
			res_[i][j]=A_[j][i];
		};
	};
}

bool utils::the_value_is_double_old(
string str_,
double &val_)
{
	bool res=true;
	int flag=sscanf(str_.c_str(),"%lf",&val_);
	if(flag!=1)
	{
		res=false;
	}
	return res;
}

bool utils::the_value_is_double(
string str_,
double &val_)
{
	if(str_=="")
	{
		return false;
	};

	bool res=false;
	errno=0;
	char *p;
	val_=strtod(str_.c_str(),&p);
	if(errno!=0)
	{
		res=false;
	}
	else
	{
		res=(*p==0);
	};
	return res;
	
}

bool utils::the_value_is_long(
string str_,
long int &val_)
{

	if(str_.length()==0)
	{
		return false;
	};
	if(!(str_[0]=='+'||str_[0]=='-'||isdigit(str_[0])))
	{
		return false;
	};

	long int start_digit=0;

	if(str_[0]=='+'||str_[0]=='-')
	{
		start_digit=1;
	};


	long int i;
	for(i=start_digit;i<(long int)str_.size();i++)
	{
		if(!isdigit(str_[i]))
		{
			return false;
		};
	};

	if(((long int)str_.size()-start_digit)<=0)
	{
		return false;
	};

	if(((long int)str_.size()-start_digit)>1)
	{
		/*
		if(str_[start_digit]=='0')
		{
			return false;
		};
		*/

		while(str_[start_digit]=='0')
		{
			string::iterator it=str_.begin()+start_digit;


			str_.erase(it);
			if((long int)str_.size()<=start_digit+1)
			{
				break;
			};
		};
	};

	if(((long int)str_.size()-start_digit>10)||((long int)str_.size()-start_digit)<=0)
	{
		return false;
	};


	if((long int)str_.size()-start_digit==10)
	{
		if(!(str_[start_digit]=='1'||str_[start_digit]=='2'))
		{
			return false;
		};

		if(str_[start_digit]=='2')
		{

			long int val2;
			string str2=str_.substr(start_digit+1,9);
			int flag=sscanf(str2.c_str(),"%ld",&val2);
			if(flag!=1)
			{
				return false;
			};

			bool positive=true;
			if(start_digit>0)
			{
				if(str_[0]=='-')
				{
					positive=false;
				};
			};

			if(positive)
			{
				if(val2>147483647)
				{
					return false;
				};
			}
			else
			{
				if(val2>147483648)
				{
					return false;
				};
			};

		};
	};

	int flag=sscanf(str_.c_str(),"%ld",&val_);
	if(flag!=1)
	{
		return false;
	};

	return true;
}

