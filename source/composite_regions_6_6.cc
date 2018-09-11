#include <iostream>
#include <iomanip>
#include <vector>
#include <math.h>
#include <DataTools.h>
 //#include "LUdcmp.h"

using namespace std;
extern "C" {
  void dgetrf_(int* M, int *N, double* A,
	       int* LDA, int* IPIV, int* INFO);
  void dgetrs_(char* TRANS, int* N, int* NRHS, double* A,
	       int* LDA, int* IPIV, double* B, int* LDB, int* INFO);
}

void initialize_matrix(vector<vector<double> > &matrix, const double value)
{
  int N=matrix.size();//matrix rows
  int M=matrix[0].size();//matrix columns
  for (int i=0; i<N; i++)
	  for (int j=0; j<M; j++)
		  matrix[i][j]=value;
}

class LapackLUdcmp
{
public:
	LapackLUdcmp(vector<vector<double> > &matrix);
	//	void solve();
	double det();
private:
	int N; //matrix rows
	int M; //matrix columns
	int *IPIV;
	int INFO;
	//	int LWORK;
	double *LU;
	//	double *WORK;
};

LapackLUdcmp::LapackLUdcmp(vector<vector<double> > &matrix)
{
	N=matrix.size();
	M=matrix[0].size();
	//	LWORK=N*N;
	//	WORK=new double[LWORK];
	LU = new double[N*M];
	IPIV=new int[N+1];
	INFO=0;
	/*
	 * express matrix in column major form for fortran
	 */
	for (int i=0; i<N; i++)
		for (int j=0; j<M; j++)
			LU[M*j+i]=matrix[i][j];

	dgetrf_(&N,&N,LU,&N,IPIV,&INFO);//lapack_ludcmp
	/*
	 * transfer result to matrix
	 */
	for (int i=0; i<N; i++)
		for (int j=0; j<M; j++)
			matrix[i][j]=LU[M*j+i];
}

double LapackLUdcmp::det()
{
	double determinant=1.;
	for (int i=0; i<N; i++)
		determinant*=LU[i*N+i];

	int sign=1;
	for (int i=0; i<N; i++)
	{
		if (IPIV[i]!=i+1)
			sign=-sign;
	}
	return sign*determinant;
}

class AnalyticCompositeRegion
{
public:
  AnalyticCompositeRegion();
  void find_roots(const unsigned int number_of_roots,
		  vector<double> &roots);
  void find_coefficients(const vector<double> &roots,
			 vector<double> &C, vector<double> &D);
  double find_value(const double position, const double time_in_seconds, const vector<double> &roots,
		    const vector<double> &C, const vector<double> &D);
  void form_matrix_simple(
			  vector<vector<double> >& matrix,
			  const double beta);
  double top_external_temperature(double time_in_seconds);

  double Psi(double A, double x);
  double Phi(double A, double x);
  double Psi_der(double A, double x);
  double Phi_der(double A, double x);

  void print_properties();
  void calculate_results();
private:
  unsigned int number_of_layers;
  //unsigned int number_of_equations;
  unsigned int number_of_roots;
  double h1;
  double h2;
  double initial_condition;
  double origin;
  double boundary_average;
  double boundary_amplitud;
  double boundary_period;
  double boundary_phase;
  string filename_properties;
  string filename_positions;
  string filename_times;
  bool active_point_source;
  double source_position;
  double source_amplitud;
  double source_period;
  double source_phase;

  vector<double> k;
  vector<double> rho;
  vector<double> cp;
  vector<double> alpha;
  vector<double> positions;
  vector<double> x;
  vector<double> layer_thickness;
  vector<double> times;
};

AnalyticCompositeRegion::AnalyticCompositeRegion():
		number_of_roots(500),
		h1(10000000.000),
		h2(0.000),
		initial_condition(0.),
		origin(0.),
		boundary_average(15.),
		boundary_amplitud(10.),
		boundary_period(86400.),
		boundary_phase(54000.),
		filename_properties("properties.txt"),
		filename_positions("positions.txt"),
		filename_times("times.txt"),
		active_point_source(false),
		source_position(1.05),
		source_amplitud(100.),
		source_period(86400.),
		source_phase(54000.)
{
	{
		DataTools data_tools;
		vector< vector<double> > data;
		data_tools.read_data(filename_properties,data);
		x.push_back(origin);
		for (unsigned int i=0; i<data.size(); i++)
		{
			layer_thickness.push_back(data[i][0]);
			k.push_back(data[i][1]);
			rho.push_back(data[i][2]);
			cp.push_back(data[i][3]);
			alpha.push_back(data[i][1]/(data[i][2]*data[i][3]));
			x.push_back(x[i]+data[i][0]);
		}
	}
	{
		DataTools data_tools;
		vector< vector<double> > data;
		data_tools.read_data(filename_times,data);
		for (unsigned int i=0; i<data.size(); i++)
		{
			times.push_back(data[i][0]);
		}
	}
	{
		DataTools data_tools;
		vector< vector<double> > data;
		data_tools.read_data(filename_positions,data);
		for (unsigned int i=0; i<data.size(); i++)
		{
			positions.push_back(data[i][0]);
		}
	}
	/*
	 * Check if the requested positions are within
	 * the layers rage
	 */
//	for (unsigned int i=0; i<positions.size(); i++)
//	{
//		if (positions[i]<x[0] ||
//			positions[i]>x[x.size()-1])
//		{
//			cout << "Error. Position "
//					<< positions[i] << " is out of range: "
//					<< x[0] << " --- " << x[x.size()-1] << "\n";
//		}
//	}

	number_of_layers=layer_thickness.size();
//	number_of_equations=2*number_of_layers;
}

void AnalyticCompositeRegion::print_properties()
{
	cout << "Solving a domain with " << number_of_layers << " layers\n";
	cout << "Dz(m)\tk(W/mK)\trho(kg/m3)\tcp(J/kgK)\n";
	for (unsigned int i=0; i<layer_thickness.size(); i++)
	{
		cout << layer_thickness[i] << "\t"
				<< k[i] << "\t"
				<< rho[i] << "\t"
				<< cp[i] << "\t"
				<< alpha[i] << "\n";
	}

	cout << "\nPositions (m)\n";
	for (unsigned int i=0; i<positions.size(); i++)
	{
		cout << positions[i] << "\n";
	}

	cout << "\nTimes (seconds)\n";
	for (unsigned int i=0; i<times.size(); i++)
	{
		cout << times[i] << "\n";
	}

	cout << "\nLayers limits (m)\n";
	for (unsigned int i=0; i<x.size(); i++)
	{
		cout << x[i] << "\n";
	}
}

void AnalyticCompositeRegion::calculate_results()
{
	vector<double> roots;
	vector<double> coeff_C;
	vector<double> coeff_D;

	find_roots(number_of_roots,roots);
	find_coefficients(roots,coeff_C,coeff_D);

	cout << "\nX";
	for (unsigned int t=0; t<times.size(); t++)
		cout << "\tt: " << times[t]/3600. << "h";
	cout << endl;

	for (unsigned int p=0; p<positions.size(); p++)
	{
		cout << positions[p];
		for (unsigned int t=0; t<times.size(); t++)
			cout << "\t" << find_value(positions[p],times[t],roots,coeff_C,coeff_D);
		cout << endl;
	}

}

double AnalyticCompositeRegion::top_external_temperature(const double time_in_seconds)
{
  return (boundary_average
		  +boundary_amplitud*cos((2.*M_PI/boundary_period)*(time_in_seconds-boundary_phase)));
}

void AnalyticCompositeRegion::find_roots(const unsigned int number_of_roots,
		vector<double> &roots)
{
	int Eq=2*number_of_layers;
	vector<vector<double> > matrix(Eq,vector<double>(Eq));

	double delta=0.00001;
	double epsilon=1.E-13;
	double a=1.E-11;
	double b=a+delta;
	double fa=0.;
	double fb=0.;
	/*
	 * Initial guesses fa, fb
	 */
	form_matrix_simple(matrix,a);
	LapackLUdcmp lapack_ludcmp(matrix);
	fa=lapack_ludcmp.det();

	form_matrix_simple(matrix,b);
	lapack_ludcmp=LapackLUdcmp(matrix);
	fb=lapack_ludcmp.det();
	/*
	 * Find the roots (betas)
	 */
	for (unsigned int i=0; i<number_of_roots; i++)
	{
		while (fa*fb>0)// find next interval where the root lies
		{
			a=b;
			b=b+delta;

			form_matrix_simple(matrix,a);
			lapack_ludcmp=LapackLUdcmp(matrix);
			fa=lapack_ludcmp.det();

			form_matrix_simple(matrix,b);
			lapack_ludcmp=LapackLUdcmp(matrix);
			fb=lapack_ludcmp.det();
		}
		unsigned int step=0;
		while (fabs(a-b)>epsilon)
		{
			double c=(a+b)/2.;
			form_matrix_simple(matrix,c);
			lapack_ludcmp=LapackLUdcmp(matrix);
			double fc=lapack_ludcmp.det();

			if (fc*fa<0.)
			{
				b=c;
				fb=fc;
			}
			else
			{
				a=c;
				fa=fc;
			}
			if (step>90 && step<=100)
				cout << "\troot interval f(" << a << "): " << fa << "\tf(" << b << "): "
				<< fb << "\tepsilon: " << fabs(a-b)  << "\tf(" << c << "): " << fc << endl;
			step++;
		}
		a=b;
		fa=fb;
		/*
		 * We found a root, calculate the coefficients
		 */
		{
			form_matrix_simple(matrix,a);
			lapack_ludcmp=LapackLUdcmp(matrix);

			vector<double>X(Eq);
			X[Eq-1]=1.;// this is the free variable and we give a value arbitrarily X3, D2
			for (unsigned int i=1; i<Eq; i++)
			{
				for (unsigned int j=Eq-i; j<=Eq-1; j++)
					X[Eq-i-1]-=matrix[Eq-i-1][j]*X[j];
				X[Eq-i-1]/=matrix[Eq-i-1][Eq-i-1];
			}
			roots.push_back(b);
		}
	}
}

void AnalyticCompositeRegion::find_coefficients(const vector<double> &roots,
		vector<double> &C, vector<double> &D)
{
	int Eq=2*number_of_layers;
	vector<vector<double> > matrix(Eq,vector<double>(Eq));
	LapackLUdcmp lapack_ludcmp(matrix);
	for (unsigned int beta=0; beta<roots.size(); beta++)
	{
		form_matrix_simple(matrix,roots[beta]);
		lapack_ludcmp=LapackLUdcmp(matrix);
		vector<double>X(Eq);
		X[Eq-1]=1.;// this is the free variable and we give a value arbitrarily X3, D2
		for (unsigned int i=1; i<Eq; i++)
		{
			for (unsigned int j=Eq-i; j<=Eq-1; j++)
				X[Eq-i-1]-=matrix[Eq-i-1][j]*X[j];
			X[Eq-i-1]/=matrix[Eq-i-1][Eq-i-1];
		}
		for (unsigned int i=0; i<Eq; i++)
		{
			if (i%2==0)
				C.push_back(X[i]);
			else
				D.push_back(X[i]);
		}
	}
}

double AnalyticCompositeRegion::find_value(const double position, const double time_in_seconds,
		const vector<double> &roots,
		const vector<double> &C, const vector<double> &D)
{
	double temperature_top=top_external_temperature(time_in_seconds);
	double temperature_top_initial=top_external_temperature(0.);
	unsigned int layer=0;
	for (unsigned int i=0; i<x.size(); i++)
		if (position>=x[i] && position<=x[i+1])
		{
			layer=i;
			break;
		}

	double solution=0.;
	double period_source=1.;
	for (unsigned int i=0; i<roots.size(); i++)
	{
		double norm=0.;
		double In=0.;
		double Vn=0.;
		double fn=0.;
		double Gn=0;
		double GGnn=0.;
		for (unsigned int l=0; l<number_of_layers; l++)
		{
			double Cl=C[number_of_layers*i+l];
			double Dl=D[number_of_layers*i+l];

			double A=roots[i]/sqrt(alpha[l]);

			fn+=(k[l]/alpha[l])*(
					(Cl/A)*(sin(A*x[l+1])-sin(A*x[l]))-
					(Dl/A)*(cos(A*x[l+1])-cos(A*x[l])));
			norm+=(k[l]/alpha[l])*(
					pow(Cl,2.)*(x[l+1]/2.+sin(2.*A*x[l+1])/(4.*A)-x[l]/2.-sin(2.*A*x[l])/(4.*A))
					+2.*Cl*Dl*(pow(sin(A*x[l+1]),2.)/(2.*A)-pow(sin(A*x[l]),2.)/(2.*A))
					+pow(Dl,2.)*(x[l+1]/2.-sin(2.*A*x[l+1])/(4.*A)-x[l]/2.+sin(2.*A*x[l])/(4.*A)));

			In+=(k[l]/alpha[l])*(
					(Cl/A)*(sin(A*x[l+1])-sin(A*x[l]))-
					(Dl/A)*(cos(A*x[l+1])-cos(A*x[l])));

			GGnn+=(
					(Cl/A)*(sin(A*x[l+1])-sin(A*x[l]))-
					(Dl/A)*(cos(A*x[l+1])-cos(A*x[l])));

			if (active_point_source==true)
				for (unsigned int g=1; g<=950; g++)
				{
					//double period_source=1.;
					double q=g*(2.*M_PI)/period_source;
					Gn+=cos(q*source_position)*
							(Cl*(sin((A-q)*x[l+1])/(2.*(A-q))+sin((A+q)*x[l+1])/(2.*(A+q))
								-sin((A-q)*x[l  ])/(2.*(A-q))-sin((A+q)*x[l  ])/(2.*(A+q)))+
								Dl*(-cos((A-q)*x[l+1])/(2.*(A-q))-cos((A+q)*x[l+1])/(2.*(A+q))
									+cos((A-q)*x[l  ])/(2.*(A-q))+cos((A+q)*x[l  ])/(2.*(A+q))))
								+
								sin(q*source_position)*
								(Cl*(-cos((q-A)*x[l+1])/(2.*(q-A))-cos((q+A)*x[l+1])/(2.*(q+A))
								+cos((q-A)*x[l  ])/(2.*(q-A))+cos((q+A)*x[l  ])/(2.*(q+A)))+
								Dl*(sin((q-A)*x[l+1])/(2.*(q-A))-sin((q+A)*x[l+1])/(2.*(q+A))
								-sin((q-A)*x[l  ])/(2.*(q-A))+sin((q+A)*x[l  ])/(2.*(q+A))));
				}
		}
		{
			double boundary_term=0.;
			{
				double P=2.*M_PI/boundary_period;
				boundary_term=
						(-P*boundary_amplitud)*(1./(pow(roots[i],4.)+pow(P,2.)))*
						(cos(P*boundary_phase)*(pow(roots[i],2.)*sin(P*time_in_seconds)-P*cos(P*time_in_seconds)+exp(-pow(roots[i],2.)*time_in_seconds)*P)
								-sin(P*boundary_phase)*(pow(roots[i],2.)*cos(P*time_in_seconds)+P*sin(P*time_in_seconds)-exp(-pow(roots[i],2.)*time_in_seconds)*pow(roots[i],2.)));
				boundary_term*=(In/norm);
			}

			double initial_term=0.;
			{
				initial_term=(initial_condition-temperature_top_initial)*fn*exp(-pow(roots[i],2.)*time_in_seconds)/norm;
			}

			double source_term=0.;
			if (active_point_source==true)
			{
				double P=2.*M_PI/source_period;
				//periodic source
				double Gt=-(1./(pow(roots[i],4.)+pow(P,2.)))*
						(cos(P*source_phase)*(pow(roots[i],2.)*sin(P*time_in_seconds)-P*cos(P*time_in_seconds)+exp(-pow(roots[i],2.)*time_in_seconds)*P)
						-sin(P*source_phase)*(pow(roots[i],2.)*cos(P*time_in_seconds)+P*sin(P*time_in_seconds)-exp(-pow(roots[i],2.)*time_in_seconds)*pow(roots[i],2.)));
				//constant source
				//double Gt=(1./pow(roots[i],2.))*(1.-exp(-pow(roots[i],2.)*time_in_seconds));

				double Gx=source_amplitud*
						(Gn/norm + 0.5*GGnn/norm)*(2./period_source);

				source_term=Gx*Gt;
			}
			double Cl=C[number_of_layers*i+layer];
			double Dl=D[number_of_layers*i+layer];
			double A=roots[i]/sqrt(alpha[layer]);

			double time_coefficient=initial_term-boundary_term+source_term;

			double x_coeff=Cl*cos(A*position)+Dl*sin(A*position);

			solution+=
					x_coeff*time_coefficient;
		}
	}
	return solution+temperature_top;
}

void AnalyticCompositeRegion::form_matrix_simple(vector<vector<double> >& matrix,
		const double beta)
{
	unsigned int m=number_of_layers;
	vector<double> A(m);
	for (unsigned int i=0; i<m; i++)
		A[i]=beta/sqrt(k[i]/(rho[i]*cp[i]));

	initialize_matrix(matrix, 0.);//fill matrix with zeroes

	matrix[0][0]=h1*Phi(A[0], x[0])-k[0]*Phi_der(A[0], x[0]); //C0
	matrix[0][1]=h1*Psi(A[0], x[0])-k[0]*Psi_der(A[0], x[0]); //D0
	// matrix[0][2]=0.; //C1
	// matrix[0][3]=0.; //D1
	// matrix[0][4]=0.; //C2
	// matrix[0][5]=0.; //D2

	for (unsigned int i=1; i<m; i++)
	{
		matrix[2*i-1][2*(i-1)+0]= Phi(A[i-1], x[i]); //C0
		matrix[2*i-1][2*(i-1)+1]= Psi(A[i-1], x[i]); //D0
		matrix[2*i-1][2*(i-1)+2]=-Phi(A[i  ], x[i]); //C1
		matrix[2*i-1][2*(i-1)+3]=-Psi(A[i  ], x[i]); //D1

		matrix[2*i  ][2*(i-1)+0]= (k[i-1]/k[i])*Phi_der(A[i-1], x[i]); //C0
		matrix[2*i  ][2*(i-1)+1]= (k[i-1]/k[i])*Psi_der(A[i-1], x[i]); //D0
		matrix[2*i  ][2*(i-1)+2]=              -Phi_der(A[i  ], x[i]); //C1
		matrix[2*i  ][2*(i-1)+3]=              -Psi_der(A[i  ], x[i]); //D1
	}

	// matrix[2*m-1][0]= 0.; //C0
	// matrix[2*m-1][1]= 0.; //D0
	// matrix[2*m-1][2]= 0.; //C1
	// matrix[2*m-1][3]= 0.; //D1
	matrix[2*m-1][2*(m-1)+0]= h2*Phi(A[m-1], x[m])+k[m-1]*Phi_der(A[m-1], x[m]); //C2
	matrix[2*m-1][2*(m-1)+1]= h2*Psi(A[m-1], x[m])+k[m-1]*Psi_der(A[m-1], x[m]); //D2
}

double AnalyticCompositeRegion::Psi(double A, double x)
{
  return (sin(A*x));
}

double AnalyticCompositeRegion::Phi(double A, double x)
{
  return (cos(A*x));
}

double AnalyticCompositeRegion::Psi_der(double A, double x)
{
  return (A*cos(A*x));
}

double AnalyticCompositeRegion::Phi_der(double A, double x)
{
  return (-A*sin(A*x));
}

int main()
{
	AnalyticCompositeRegion analytic_composite_region;
	analytic_composite_region.print_properties();
	analytic_composite_region.calculate_results();

	return 1;
}
