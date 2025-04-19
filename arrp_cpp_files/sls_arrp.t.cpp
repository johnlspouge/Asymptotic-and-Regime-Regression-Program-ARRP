// sls_arrp.t.cpp

#include "sls_utils.h"
#include "sls_asymptotic_regression.h"
#include "sls_robust_regression.h"

using namespace Sls;

void error_message()
{
	cout<<
"ARRP 1.1\n\
The program is an implementation of robust and asymptotic regressions.\n\
\n\
USAGE:\n\
\n\
arrp.exe <OPTIONS>\n\
The <OPTIONS> is a series of options separated by spaces.\n\
\n\
Possible options for ARRP (square brackets indicate that the option can be omitted):\n\
\n\
-in <Input file name>\n\
-out <Output file name>\n\
[ -append ], defines whether the output file is appended (the option is present) or rewritten\n\
[ -flat ], an option for switching to constant regression on a line of 0 slope\n\
[ -regression <type> ], type of regression (asymptotic regression is performed if the option is omitted); the possible alternatives are\n\
-regression robust: ARRP performs robust regression\n\
-regression weighted: ARRP performs weighted least-squares regression\n\
-regression asymptotic: ARRP performs asymptotic regression\n\
\n\
Additional options for asymptotic regression:\n\
[ -c <Parameter 'c'> ], the parameter 'c' of asymptotic regression\n\
[ -include <parameter> ], if the option is omitted, asymptotic regression permits removal of data points on both the left and right sides; if the option is present, then the possible alternatives are\n\
-include left: the asymptotic regression must include points on the left\n\
-include right: the asymptotic regression must include points on the right\n\
-include neither: permits removal of data points on both the left and right sides\n\
\n\
Additional options for robust regression:\n\
[ -a <Parameter 'a'> ], the parameter 'a' of robust regression\n\
\n\
The most recent version and a help page can be found here:\n\
http://www.ncbi.nlm.nih.gov/CBBresearch/Spouge/html_ncbi/html/software/program.html?uid=3\n\
";

	throw error("Error - parameters of the program are wrong\n",-1);
}

void command_line_interpreter2(//extracts parameters from the command line
	int argc, char* argv[],//arguments of the program
	string &reg_type_,
	double &c_,
	long int &cl_,
	long int &cr_,
	string &include_,
	bool &append_,
	double &a_,
	string &in_file_name_,
	string &out_file_name_)//output file name
{
	c_=4.0;
	a_=1.339;

	bool c_flag=false;
	bool a_flag=false;
	bool in_file_name_flag=false;
	bool out_file_name_flag=false;

	append_=false;
	bool append_flag=false;

	bool flat_=false;
	bool flat_flag=false;

	string regression_="asymptotic";
	bool regression_flag=false;

	include_="neither";
	bool include_flag=false;

	long int i=1;
	while(i<argc)
	{

		string tmp=argv[i];

		if(tmp=="-in")
		{
			if(in_file_name_flag)
			{
				cerr<<"Error - parameter -in has been defined already\n";
				cout<<endl;
				error_message();
			};

			i++;

			if(i>=argc)
			{
				cerr<<"Error - parameter -in is wrong\n";
				cout<<endl;
				error_message();
			};

			in_file_name_=argv[i];
			in_file_name_flag=true;
			i++;
			continue;
		};

		if(tmp=="-out")
		{
			if(out_file_name_flag)
			{
				cerr<<"Error - parameter -out has been defined already\n";
				cout<<endl;
				error_message();
			};

			i++;

			if(i>=argc)
			{
				cerr<<"Error - parameter -out is wrong\n";
				cout<<endl;
				error_message();
			};

			out_file_name_=argv[i];
			out_file_name_flag=true;
			i++;
			continue;
		};

		if(tmp=="-append")
		{
			if(append_flag)
			{
				cerr<<"Error - parameter -append has been defined already\n";
				cout<<endl;
				error_message();
			};

			append_=true;
			append_flag=true;
			i++;
			continue;
		};

		if(tmp=="-flat")
		{
			if(flat_flag)
			{
				cerr<<"Error - parameter -flat has been defined already\n";
				cout<<endl;
				error_message();
			};

			flat_=true;
			flat_flag=true;
			i++;
			continue;
		};

		if(tmp=="-regression")
		{
			if(regression_flag)
			{
				cerr<<"Error - parameter -regression has been defined already\n";
				cout<<endl;
				error_message();
			};

			i++;

			if(i>=argc)
			{
				cerr<<"Error - parameter -regression is wrong\n";
				cout<<endl;
				error_message();
			};

			regression_=argv[i];
			std::transform(regression_.begin(),regression_.end(),regression_.begin(),::tolower);
			regression_flag=true;
			i++;
			continue;
		};

		if(tmp=="-c")
		{
			if(c_flag)
			{
				cerr<<"Error - parameter -c has been defined already\n";
				cout<<endl;
				error_message();
			};

			i++;

			if(i>=argc)
			{
				cerr<<"Error - parameter -c is wrong\n";
				cout<<endl;
				error_message();
			};


			if(!utils::the_value_is_double(argv[i],c_))
			{
				cerr<<"Error - parameter -c is wrong\n";
				cout<<endl;
				error_message();
			};

			//rescaling c parameter
			c_*=2.0;

			c_flag=true;
			i++;
			continue;
		};

		if(tmp=="-include")
		{
			if(include_flag)
			{
				cerr<<"Error - parameter -include has been defined already\n";
				cout<<endl;
				error_message();
			};

			i++;

			if(i>=argc)
			{
				cerr<<"Error - parameter -include is wrong\n";
				cout<<endl;
				error_message();
			};

			include_=argv[i];
			std::transform(include_.begin(),include_.end(),include_.begin(),::tolower);
			include_flag=true;
			i++;
			continue;
		};

		if(tmp=="-a")
		{
			if(a_flag)
			{
				cerr<<"Error - parameter -a has been defined already\n";
				cout<<endl;
				error_message();
			};

			i++;

			if(i>=argc)
			{
				cerr<<"Error - parameter -a is wrong\n";
				cout<<endl;
				error_message();
			};


			if(!utils::the_value_is_double(argv[i],a_))
			{
				cerr<<"Error - parameter -a is wrong\n";
				cout<<endl;
				error_message();
			};


			a_flag=true;
			i++;
			continue;
		};


		cerr<<"Error - the command line with input parameters is incorrect\n";
		cout<<endl;
		error_message();


	};

	//checking parameters
	if(!in_file_name_flag)
	{
		cerr<<"Error - missing parameter -in\n";
		cout<<endl;
		error_message();
	};

	if(!out_file_name_flag)
	{
		cerr<<"Error - missing parameter -out\n";
		cout<<endl;
		error_message();
	};

	if(!(regression_=="robust"||regression_=="weighted"||regression_=="asymptotic"))
	{
		cerr<<"Error - the parameter -regression is incorrect\n";
		cout<<endl;
		error_message();
	};


	if(regression_=="robust")
	{
		//robust regression
		if(c_flag)
		{
			cerr<<"Error - unnecessary parameter -c for robust regression\n";
			cout<<endl;
			error_message();
		};

		if(include_flag)
		{
			cerr<<"Error - unnecessary parameter -include for robust regression\n";
			cout<<endl;
			error_message();
		};

		if(a_<=0)
		{
			cerr<<"Error - the parameter -a must be strictly positive\n";
			cout<<endl;
			error_message();
		};

		if(flat_)
		{
			reg_type_="cr0";
		}
		else
		{
			reg_type_="cr1";
		};

	};


	if(regression_=="asymptotic")
	{
		//asymptotic regression
		if(a_flag)
		{
			cerr<<"Error - unnecessary parameter -a for asymptotic regression\n";
			cout<<endl;
			error_message();
		};

		if(c_<=0)
		{
			cerr<<"Error - the parameter -c must be strictly positive\n";
			cout<<endl;
			error_message();
		};

		if(!(include_=="right"||include_=="left"||include_=="neither"))
		{
			cerr<<"Error - the parameter -include can only take the values \"left\", \"right\", \"neither\"\n";
			cout<<endl;
			error_message();
		};

		if(flat_)
		{
			reg_type_="ar0";
		}
		else
		{
			reg_type_="ar1";
		};

		if(include_=="neither")
		{
			cl_=1;
			cr_=1;
		};

		if(include_=="left")
		{
			cl_=0;
			cr_=1;
		};

		if(include_=="right")
		{
			cl_=1;
			cr_=0;
		};

	};

	if(regression_=="weighted")
	{
		//weighted least-squares regression
		cl_=0;
		cr_=0;

		if(c_flag)
		{
			cerr<<"Error - unnecessary parameter -c for weighted least-squares regression\n";
			cout<<endl;
			error_message();
		};

		if(include_flag)
		{
			cerr<<"Error - unnecessary parameter -include for weighted least-squares regression\n";
			cout<<endl;
			error_message();
		};

		if(a_flag)
		{
			cerr<<"Error - unnecessary parameter -a for weighted least-squares regression\n";
			cout<<endl;
			error_message();
		};

		if(flat_)
		{
			reg_type_="ar0";
		}
		else
		{
			reg_type_="ar1";
		};

		include_="both";

	};



}

int main_tmp(int argc, char* argv[], 
			 pair<double,double> &betas_, 
			 bool &betas_flag_, 
			 double &k_left_,
			 double &k_right_)
{


	betas_flag_=false;

	double c;
	long int cl;
	long int cr;
	string include;
	string reg_type;
	bool append;
	double a_in;
	string in_file_name;
	string out_file_name;//output file name
	asymptotic_regression tmp;


	try
	{


		command_line_interpreter2(//extracts parameters from the command line
		argc,argv,//arguments of the program
		reg_type,
		c,
		cl,
		cr,
		include,
		append,
		a_in,
		in_file_name,
		out_file_name);//output file name

 		}
		catch (error er)
		{
			return 1;
		};



		//read data begin
		long int number_of_elements;

		ifstream fin(in_file_name.data());
		if(!fin)
		{
			throw error("Error - input file "+in_file_name+" is not found\n",1);
		};


		long int current_ind=-1;

		vector<double> args_vect;
		vector<double> array_vect;
		vector<double> array_errors_vect;

		while(!fin.eof())
		{
			string st;
			getline(fin,st);

			vector<std::string> st_strings;
			istringstream st_in(st);
			while(st_in)
			{
				std::string st_tmp="";
				st_in>>st_tmp;
				if(st_tmp!="")
				{
					st_strings.push_back(st_tmp);
				};
			};

			if(st_strings.size()==0)
			{
				continue;
			};

			if(st_strings.size()!=3)
			{
				throw error("Error in the input file "+in_file_name+"\n",1);
			};

			current_ind++;

			double args_tmp,array_tmp,array_errors_tmp;

			if(!(utils::the_value_is_double(st_strings[0],args_tmp)&&
				utils::the_value_is_double(st_strings[1],array_tmp)&&
				utils::the_value_is_double(st_strings[2],array_errors_tmp)
				)
				)
			{
				throw error("Error in the input file "+in_file_name+": not valid real number\n",1);
			};

			if(array_errors_tmp<=0)
			{
				throw error("Error in the input file "+in_file_name+": values of errors must be strictly positive\n",1);
			};

			args_vect.push_back(args_tmp);
			array_vect.push_back(array_tmp);
			array_errors_vect.push_back(array_errors_tmp);

		};

		number_of_elements=current_ind+1;

		if(reg_type=="ar0"||reg_type=="cr0")
		{
			if(number_of_elements<=0)
			{
				throw error("Error- number of points cannot be less than 1 for constant regression\n",1);
			};
		};

		if(reg_type=="ar1"||reg_type=="cr1")
		{
			if(number_of_elements<=1)
			{
				throw error("Error - number of points cannot be less than 2 for linear regression\n",1);
			};
		};


		double *args=new double[number_of_elements];
		utils::assert_mem(args);

		double *array=new double[number_of_elements];
		utils::assert_mem(array);

		double *array_errors=new double[number_of_elements];
		utils::assert_mem(array_errors);

		long int j;

		for(j=0;j<number_of_elements;j++)
		{
			args[j]=args_vect[j];
			array[j]=array_vect[j];
			array_errors[j]=array_errors_vect[j];
		};

		for(j=1;j<number_of_elements;j++)
		{
			if(args[j]<args[j-1])
			{
				throw error("Error in the input file "+in_file_name+"; data points must be sorted\n",1);
			};

			if(args[j]==args[j-1])
			{
				throw error("Error in the input file "+in_file_name+"; X-coordinates must be different\n",1);
			};
		};

	

		fin.close();
		//read data end




		bool RR_res_flag=false;

		std::vector<double> res_beta;//output
		std::vector<double> res_beta_errors;//output
		std::vector<double> limits;//output

		if(reg_type=="cr0"||reg_type=="cr1")
		{
			//w_function_type *w_function=&(robust_regression::w_Andrews);
			//double a=1.339;

			//w_function_type *w_function=&(robust_regression::w_Huber);
			//double a=2.0;

		
			w_function_type *w_function=&(robust_regression::w_Ramsay);
			double a=0.3;


		void* func_number=&a;
		std::vector<std::vector<double> > X;
		std::vector<double> y;
		double eps=0.00000001;;
		double inside_eps=1e-50;
		std::vector<double> med;
		double var=0;//output
		std::string flag="Classical iteration process";
		//flag="Simple regression"
		//flag="Classical iteration process"

		//if(result[0]==3)
		//if(0)
		{
			long int j;
			for(j=0;j<number_of_elements;j++)
			{
				med.push_back(array_errors[j]);
			};
		};



		if(reg_type=="cr1")
		{
			long int i;
			for(i=0;i<number_of_elements;i++)
			{
				std::vector<double> tmp(1);
				tmp[0]=args[i];


				X.push_back(tmp);
				y.push_back(array[i]);
			};
		};

		if(reg_type=="cr0")
		{
			long int i;
			for(i=0;i<number_of_elements;i++)
			{
				X.push_back(std::vector<double>(0));
				y.push_back(array[i]);
			};

		};


	RR_res_flag=robust_regression::M_estimator(
		w_function,
		func_number,
		X,
		y,
		eps,
		res_beta,
		res_beta_errors,
		inside_eps,
		med,
		var,
		flag,
		//flag="Simple regression"
		//flag="Classical iteration process"
		limits
		);

	

		if(!RR_res_flag)
		{
			throw error("Error - regression failed\n",1);
		}
		else
		{

			
			a=a_in;
			w_function_type *w_function=&(robust_regression::w_Andrews);

			RR_res_flag=robust_regression::M_estimator(
			w_function,
			func_number,
			X,
			y,
			eps,
			res_beta,
			res_beta_errors,
			inside_eps,
			med,
			var,
			flag,
			//flag="Simple regression"
			//flag="Classical iteration process"
			limits
			);
			
			if(!RR_res_flag)
			{
				throw error("Error - regression failed\n",1);
			};
			
		};
	};

	if(reg_type=="ar0"||reg_type=="ar1")
	{
		res_beta.resize(2);//output
		res_beta_errors.resize(2);//output
		limits.resize(2+number_of_elements);//output

		double beta0=0;
		double beta0_error=0;

		double beta1=0;
		double beta1_error=0;

		bool cut_left_tail=false;
		bool cut_right_tail=false;

		if(cl)
		{
			cut_left_tail=true;
		};

		if(cr)
		{
			cut_right_tail=true;
		};


		double y=sqrt(c);

		long int k1_opt=0;
		long int k2_opt=0;


		bool res_was_calculated=false;

		long int min_length;

		if(reg_type=="ar0")
		{
			min_length=1;

			tmp.asymptotic_regression_const_beta1_is_defined(
			min_length,
			number_of_elements,
			array,
			array_errors,
			cut_left_tail,
			cut_right_tail,
			y,
			beta0,
			beta1,
			beta0_error,
			beta1_error,
			k1_opt,
			k2_opt,
			res_was_calculated);
		};

		if(reg_type=="ar1")
		{
			min_length=2;

			tmp.asymptotic_regression_linear(
			min_length,
			number_of_elements,
			args,
			array,
			array_errors,
			cut_left_tail,
			cut_right_tail,
			y,
			beta0,
			beta1,
			beta0_error,
			beta1_error,
			k1_opt,
			k2_opt,
			res_was_calculated);
		};






		RR_res_flag=res_was_calculated;

		if(res_was_calculated)
		{
			res_beta[0]=beta0;
			res_beta_errors[0]=beta0_error;
			res_beta[1]=beta1;
			res_beta_errors[1]=beta1_error;
			limits[0]=k1_opt;
			limits[1]=k2_opt;

			long int i;
			for(i=2+k1_opt;i<=2+k2_opt;i++)
			{
				limits[i]=1.0;
			};

		};

	};

	//output

	ofstream fout;
	bool appendflag=append;

	if(!append)
	{
		fout.open(out_file_name.data());
	}
	else
	{
		ifstream f_test1_in(out_file_name.data());
		if(!f_test1_in)
		{
			fout.open(out_file_name.data());
			appendflag=false;
		}
		else
		{
			f_test1_in.close();
			fout.open(out_file_name.data(),ios::app);
		};
	};

	if(!fout)
	{
		throw error("Error - file "+out_file_name+" is not found\n",1);
	};

	if(appendflag)
	{
		fout<<"---------------------------------------------------------\n";
	};

	fout<<"Input\t"<<in_file_name<<endl;
	fout<<"Flat\t";
	if(reg_type=="ar1"||reg_type=="cr1")
	{
		fout<<"no\tResult\tf(y)~beta0+beta1*y\n";
	}
	else
	{
		fout<<"yes\tResult\tf(y)~beta0\n";
	};

	fout<<"Regression\t";
	if(reg_type=="ar0"||reg_type=="ar1")
	{
		if(include=="both")
		{
			fout<<"weighted"<<endl;
		}
		else
		{
			fout<<"asymptotic\tc\t"<<c/2<<"\tInclude\t"<<include<<endl;
		};
	}
	else
	{
		fout<<"robust\ta\t"<<a_in<<endl;
	};
	

	if(limits.size()>=2)
	{
		k_left_=args[(long int)utils::round(limits[0])];
		k_right_=args[(long int)utils::round(limits[1])];

		fout<<"Range\t"<<k_left_<<"\t"<<k_right_<<"\n";
		cout<<"Range\t"<<k_left_<<"\t"<<k_right_<<"\n";
	}


	if(RR_res_flag&&res_beta.size()>0)
	{
		betas_.first=res_beta[0];

		fout<<"beta0\t"<<res_beta[0]<<"\terror\t"<<res_beta_errors[0]<<endl;
		cout<<"beta0\t"<<res_beta[0]<<"\terror\t"<<res_beta_errors[0]<<endl;
		if(reg_type=="ar1"||reg_type=="cr1")
		{
			betas_.second=res_beta[1];
			fout<<"beta1\t"<<res_beta[1]<<"\terror\t"<<res_beta_errors[1]<<endl;
			cout<<"beta1\t"<<res_beta[1]<<"\terror\t"<<res_beta_errors[1]<<endl;
		};
		cout<<endl;
		betas_flag_=true;
		fout<<endl;
		fout<<"X\tY\tError\tWeight\n";

		long int k;
		for(k=2;k<(long int)limits.size();k++)
		{
			fout<<args[k-2]<<"\t"<<array[k-2]<<"\t"<<array_errors[k-2]<<"\t"<<limits[k]<<endl;
		};

	}
	else
	{
		fout.close();
		delete[]args;
		delete[]array;
		delete[]array_errors;

		throw error("Error - regression failed",1);
	};



	delete[]args;
	delete[]array;
	delete[]array_errors;
	fout.close();

	return 0;
}

int main(int argc, char* argv[])
{

	try
	{
	try
	{

		pair<double,double> betas;
		bool betas_flag;
		double k_left;
		double k_right;

		int res=main_tmp(argc, argv, betas, betas_flag, k_left,k_right);

		if(res!=0)
		{
			return res;
		};
		return 0;
 	}
	catch (error er)
	{
		if(er.error_code>=0)
		{
			std::cerr<<er.st;
		}
		else
		{
			std::cout<<er.st;
		};

		return 0;
	};
	}
	catch (...)
	{
		std::cerr<<"Unexpected error\n";
		return 0;
	};



}


