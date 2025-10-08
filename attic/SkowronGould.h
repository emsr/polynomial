#pragma once

#define MR 8
#define MT 10
#define MAXIT (MT*MR)
#define MAXM 30

class complex;

class complex {
public:
	double re;
	double im;
	complex(double, double);
	complex(double);
	complex(void);
};

double abs(complex);
complex conj(complex);
complex sqrt(complex);
double real(complex);
double imag(complex);
complex cbrt(complex);
complex expcmplx(complex);
complex operator+(complex, complex);
complex operator-(complex, complex);
complex operator*(complex, complex);
complex operator/(complex, complex);
complex operator+(complex, double);
complex operator-(complex, double);
complex operator*(complex, double);
complex operator/(complex, double);
complex operator+(double, complex);
complex operator-(double, complex);
complex operator*(double, complex);
complex operator/(double, complex);
complex operator+(int, complex);
complex operator-(int, complex);
complex operator*(int, complex);
complex operator/(int, complex);
complex operator+(complex, int);
complex operator-(complex, int);
complex operator*(complex, int);
complex operator/(complex, int);
complex operator-(complex);
bool operator==(complex, complex);
bool operator!=(complex, complex);


void cmplx_roots_gen(complex *, complex *, int, bool, bool);
void cmplx_laguerre(complex *, int, complex *, int &, bool &);
void cmplx_newton_spec(complex *, int, complex *, int &, bool &);
void cmplx_laguerre2newton(complex *, int, complex *, int &, bool &, int);
void solve_quadratic_eq(complex &, complex &, complex *);
void solve_cubic_eq(complex &, complex &, complex &, complex *);
void divide_poly_1(complex *, complex, complex, complex *, int);
void cmplx_roots_5(complex *, bool, complex *, bool);
void sort_5_points_by_separation(complex *);
void sort_5_points_by_separation_i(int *, complex *);
void find_2_closest_from_5(int *, int *, double *, complex *);
