#include <cstdlib>
#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>

using namespace std;

const double c_pre = 1; //pi*r_e^2 / m_e / c^2
const double m_e = 511.0; //kev
const double s_resol = 0.01;
const double s_conv_resol = 0.001;
const double sigma_start = 0;
const double sigma_end = 1;
const double sigma_resol = 10;
const double E1 = 511;
const double E2 = 1275;

class Point {
public:
	double x;
	double y;
public:
	Point(double a, double b) : x(a), y(b) {}
	void print(ofstream* out) {
		*out<<x<<"\t"<<y<<endl;
	}
};

inline double real_compton( double s, double gamma ) {
	return ( s < 1.0 ? c_pre / gamma / gamma * ( 2.0 + s*s/gamma/gamma/(1.0-s)/(1.0-s) + s/(1.0-s)*(s-2.0/gamma) ) : 0 );
}

inline double gaussian ( double x, double mu, double sigma ) {
	return (sigma > 0 ? exp ( -0.5 * pow( (x-mu)/sigma , 2) ) / sigma / 2.5077 : ( ( fabs( x-mu ) < s_conv_resol/2 ) ? 1.0 : 0 ) );
//	return ( ( fabs (x-mu) < s_conv_resol/2) ? 1.0 : 0 ) ;
}

double broad_compton ( double s, double gamma, double sigma ) {
//	return real_compton(s,gamma);
	double sum = 0;
	for( double sp = -3.; sp < 3.0; sp += s_conv_resol) {
		sum += ( real_compton(s+sp,gamma) * gaussian(-sp,0,sigma) );
	}
	return sum;
}

vector<Point*>* simulate(int E, int sigma) {
	double Emax = 2 * E * E / (m_e + 2*E);
	double gamma = E / 511; // e / m_e c^2
	vector<Point*>* v = new vector<Point*>();
	for(double s = 0; s < 1.2; s+= s_resol)
		v->push_back( new Point( s*E, broad_compton( s,gamma, sigma ) ) );
	return v;
}

void print_max_half(vector<Point*>* v, ofstream* out) {
	double max = 0;
	double emax = 0;
	double ehalf = -1;
	double eA,eC,xA,xC;
	for( int i=0; i< v->size() ; i++) {
		if ( v->at(i)->y > max ) {
			max = v->at(i)->y;
			emax = v->at(i)->x;
			ehalf = -1;
		} else if ( v->at(i)->y < max/2 && ehalf < 0 ) {
			eA = v->at(i-1)->y;
			eC = v->at(i)->y;
			xA = v->at(i-1)->x;
			xC = v->at(i)->x;
			ehalf = xA - ( eA - eC ) / ( eA - max/2 ) * ( xA - xC );
		}
	}
	cout<<emax<<"\t"<<ehalf<<endl;
	*out <<  ehalf - emax <<endl;
}

void print(vector<Point*>* v, char* name) {
	ofstream out(name);
	for(int i=0; i<v->size(); i++) v->at(i)->print(&out);
}

int main() {
	ofstream out1("output_1.txt");
	ofstream out2("output_2.txt");
	for (int i=0; i < (sigma_end-sigma_start)/sigma_resol; i++) cout<<'-';
	cout<<endl;
	for (double sigma = sigma_start; sigma < sigma_end; sigma += sigma_resol) {
		auto v1 = simulate (E1, sigma);
		auto v2 = simulate (E2, sigma);
		print(v1,"comp1.txt");
		print(v2,"comp2.txt");
		out1<<sigma<<"\t";
		out2<<sigma<<"\t";
		print_max_half(v1,&out1);
		print_max_half(v2,&out2);
		cout<<'*'<<flush;
	}
	cout<<"Done";
	return 0;
}
