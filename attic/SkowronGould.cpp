//Copyright 2012 Jan Skowron & Andrew Gould
// This translation to C++ was performed by Ava Hoag and Tyler Heintz under the supervision of Valerio Bozza at the University of Salerno during the month of July 2017.
//
//Licensed under the Apache License, Version 2.0 (the "License");
//you may not use this file except in compliance with the License.
//You may obtain a copy of the License at
//
//http://www.apache.org/licenses/LICENSE-2.0
//
//Unless required by applicable law or agreed to in writing, software
//distributed under the License is distributed on an "AS IS" BASIS,
//WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//See the License for the specific language governing permissions and
//limitations under the License.
//
//------------------------------------------------------------------ - //
//
//The authors also make this file available under the terms of
//GNU Lesser General Public License version 2 or any later version.
//(text of the LGPL licence 2 in NOTICE file)
//
//------------------------------------------------------------------ - //
//
//A custom in the scientific comunity is(regardless of the licence
	//you chose to use or distribute this software under)
	//that if this code was important in the scientific process or
	//for the results of your scientific work, we kindly ask you for the
	//appropriate citation of the Paper(Skowron & Gould 2012), and
	//we would be greatful if you pass the information about
	//the proper citation to anyone whom you redistribute this software to.
	//
	//------------------------------------------------------------------ - //
	//
	//No    Subroutine
	//
	//1   cmplx_roots_gen - general polynomial solver, works for random degree, not as fast or robust as cmplx_roots_5
	//2   cmplx_roots_5 - complex roots finding algorithm taylored for 5th order polynomial(with failsafes for polishing)
	//3   sort_5_points_by_separation - sorting of an array of 5 points, 1st most isolated, 4th and 5th - closest
	//4   sort_5_points_by_separation_i - sorting same as above, returns array of indicies rather than sorted array
	//5   find_2_closest_from_5 - finds closest pair of 5 points
	//6   cmplx_laguerre - Laguerre's method with simplified Adams' stopping criterion
	//7   cmplx_newton_spec - Newton's method with stopping criterion calculated every 10 steps
	//8   cmplx_laguerre2newton - three regime method : Laguerre's, Second-order General method and Newton's
	//9   solve_quadratic_eq - quadratic equation solver
	//10   solve_cubic_eq - cubic equation solver based on Lagrange's method
	//11   divide_poly_1 - division of the polynomial by(x - p)
	//
	//fortran 90 code
	//
	//Paper:  Skowron & Gould 2012
	//"General Complex Polynomial Root Solver and Its Further Optimization for Binary Microlenses"
	//
	//for a full text see :
    //http://www.astrouw.edu.pl/~jskowron/cmplx_roots_sg/
	//or http://arxiv.org/find/astro-ph
	//or http://www.adsabs.harvard.edu/abstract_service.html
	//see also file NOTICE and LICENSE
	//
	//ver. 2012.03.03 initial
	//ver. 2014.03.12 bug fix
	//ver. 2016.01.21 bug fix
	//ver. 2016.04.28 bug fix
	//


#define _USE_MATH_DEFINES
#define _CRT_SECURE_NO_WARNINGS

//#include "stdafx.h"
#include <math.h>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "SkowronGould.h"


void cmplx_roots_gen(complex *roots, complex *poly, int degree, bool polish_roots_after, bool use_roots_as_starting_points) {
	//roots - array which will hold all roots that had been found.
	//If the flag 'use_roots_as_starting_points' is set to
	//.true., then instead of point(0, 0) we use value from
	//this array as starting point for cmplx_laguerre

	//poly - is an array of polynomial cooefs, length = degree + 1,
	//poly[0] x ^ 0 + poly[1] x ^ 1 + poly[2] x ^ 2 + ...

	//degree - degree of the polynomial and size of 'roots' array

	//polish_roots_after - after all roots have been found by dividing
	//original polynomial by each root found,
	//you can opt in to polish all roots using full
	//polynomial

	//use_roots_as_starting_points - usually we start Laguerre's 
	//method from point(0, 0), but you can decide to use the
	//values of 'roots' array as starting point for each new
	//root that is searched for.This is useful if you have
	//very rough idea where some of the roots can be.
	//

	complex poly2[MAXM];
	static int i, j, n, iter;
	bool success;
	complex coef, prev;

	if (!use_roots_as_starting_points) {
		for (int jj = 0; jj < degree; jj++) {
			roots[jj] = complex(0, 0);
		}
	}

	for (j = 0; j <= degree; j++) poly2[j] = poly[j];

	// Don't do Laguerre's for small degree polynomials
	if (degree <= 1) {
		if (degree == 1) roots[0] = -poly[0] / poly[1];
		return;
	}

	for (n = degree; n >= 3; n--) {
		cmplx_laguerre2newton(poly2, n, &roots[n - 1], iter, success, 2);
		if (!success) {
			printf("Not a success!\n");
			roots[n - 1] = complex(0, 0);
			cmplx_laguerre(poly2, n, &roots[n - 1], iter, success);
		}

		// Divide by root
		//for (i = 0; i <= degree; i++) {
		//printf("Before Division...\npoly coef%d = (%f, %f)\n", i, poly[i].re, poly[i].im);
		//}
		coef = poly2[n];
		for (i = n - 1; i >= 0; i--) {
			prev = poly2[i];
			poly2[i] = coef;
			coef = prev + roots[n - 1] * coef;
		}
		//for (i = 0; i <= degree; i++) {
		//printf("After Division...\npoly coef%d = (%f, %f)\n", i, poly[i].re, poly[i].im);
		//}
	}

	//Find the to last 2 roots
	cmplx_laguerre2newton(poly2, 2, &roots[1], iter, success, 2);
	if (!success) {
		printf("Not a success!\n");
		solve_quadratic_eq(roots[1], roots[0], poly2);
	}
	else {
		roots[0] = -(roots[1] + poly2[1] / poly2[2]); // Viete's Formula for the last root
	}
	if (polish_roots_after) {
		for (n = 0; n < degree; n++) {
			cmplx_newton_spec(poly, degree, &roots[n], iter, success); // Polish roots with full polynomial
			if (!success) {
				printf("not a success\n");
			}
		}
	}

	return;
}

void cmplx_roots_5(complex *roots, bool first_3_roots_order_changed, complex *poly, bool polish_only) {
	//Subroutine finds or polishes roots of a complex polynomial
	//(degree = 5)
	//This routine is especially tailored for solving binary lens
	//equation in form of 5th order polynomial.
	//
	//Use of this routine, in comparission to 'cmplx_roots_gen' can yield
	//consideribly faster code, because it makes polishing of the roots
	//(that come in as a guess from previous solutions) secure by
	//implementing additional checks on the result of polishing.
	//If those checks are not satisfied then routine reverts to the
	//robust algorithm.These checks are designed to work for 5th order
	//polynomial originated from binary lens equation.
	//
	//Usage:
	//
	//polish_only == false - I do not know the roots, routine should
	//find them from scratch.At the end it
	//sorts roots from the most distant to closest.
	//Two last roots are the closest(in no particular
	//order).
	//polish_only = true - I do know the roots pretty well, for example
	//I have changed the coefficiens of the polynomial
	//only a bit, so the two closest roots are
	//most likely still the closest ones.
	//If the output flag 'first_3_roots_order_changed'
	//is returned as 'false', then first 3 returned roots
	//are in the same order as initialy given to the
	//routine.The last two roots are the closest ones,
	//but in no specific order(!).
	//If 'first_3_roots_order_changed' is 'true' then
	//it means that all roots had been resorted.
	//Two last roots are the closest ones.First is most
	//isolated one.
	//
	//
	//If you do not know the position of the roots just use flag
	//polish_only = .false.In this case routine will find the roots by
	//itself.
	//
	//poly - is an array of polynomial cooefs, length = degree + 1
	//poly(1) x ^ 0 + poly(2) x ^ 1 + poly(3) x ^ 2 + poly(4) x ^ 3 + ...
	//roots - roots of the polynomial('out' and optionally 'in')
	//
	//
	int const degree = 5;
	complex remainder, roots_robust[degree];
	double d2min;
	int iter, loops, go_to_robust, m, root4, root5, i, i2;
	complex poly2[MAXM], coef, prev;
	complex zero = complex(0, 0);
	bool succ;
	for (int j = 0; j < degree; j++) {
		roots_robust[j] = roots[j];
	}
	go_to_robust = 0;
	if (!polish_only) {
		go_to_robust = 1;
	}

	first_3_roots_order_changed = false;

	for (loops = 1; loops <= 3; loops++) {
		//ROBUST (we do not know the roots)
		if (go_to_robust > 0) {
			if (go_to_robust > 2) {
				for (int j = 0; j < degree; j++) {
					roots[j] = roots_robust[j];   // something is wrong, polishing creates errors so return unpolished roots
				}
				return;
			}
			for (int j = 0; j <= degree; j++) {
				poly2[j] = poly[j];
			}
			for (m = degree; m >= 4; m--) {
				cmplx_laguerre2newton(poly2, m, &roots[m - 1], iter, succ, 2);
				if (!succ) {
					roots[m - 1] = zero;
					cmplx_laguerre(poly2, m, &roots[m - 1], iter, succ);
				}
				coef = poly2[m];
				for (i = m - 1; i >= 0; i--) {
					prev = poly2[i];
					poly2[i] = coef;
					coef = prev + roots[m - 1] * coef;
				}
				//divide_poly_1(poly2, remainder, roots[m-1], poly2, m);
			}
			solve_cubic_eq(roots[0], roots[1], roots[2], poly2); // Find last 3 roots
																 // all roots found

																 //sort roots - first will be most isolated, last two will be closest
			sort_5_points_by_separation(roots);
			//copy roots in case something goes wrong when polishing
			for (int j = 0; j < degree; j++) {
				roots_robust[j] = roots[j];
			}

			//set flag that roots have been resorted
			first_3_roots_order_changed = true;
		}

		//POLISH (we know roots approximately and guess last two are closest
		for (int j = 0; j <= degree; j++) {
			poly2[j] = poly[j];
		}

		for (m = 1; m <= degree - 2; m++) {
			cmplx_newton_spec(poly2, degree, &roots[m - 1], iter, succ);

			if (!succ) {
				//go back to robust
				go_to_robust += 1;
				break;
			}
		}

		if (succ) {
			for (m = 1; m <= degree - 2; m++) {
				divide_poly_1(poly2, remainder, roots[m - 1], poly2, degree - m + 1);
			}
			// last two roots are found with quadratic solver
			// this is faster and more robust, alsthough little less accurate
			solve_quadratic_eq(roots[degree - 2], roots[degree - 1], poly2);
			//all roots found and polished

			//TEST ORDER
			//test closest roots if they are the same pair as given in polish
			//d2min = 0;
			find_2_closest_from_5(&root4, &root5, &d2min, roots);

			//check if the closest roots are not too close, this could happen
			//when using polishing with Newton only, when two roots erroneously
			//colapsed to the same root.This check is not needed for polishing
			//3 roots by Newton and using quadratic for the remaining two.
			//If the real roots are so close indeed(very low probability), this will just
			//take more time and the unpolished result be returned at the end
			//but algorithm will work, and will return accurate enough result
			//if (d2min<1d - 18) then             !POWN - polish only with Newton
			//go_to_robust = go_to_robust + 1    !POWN - polish only with Newton
			//else                             !POWN - polish only with Newton
			if ((root4 < degree - 1) || (root5 < degree - 1)) {
				// after polishing some of the 3 far roots become one of the closest ones
				// go back to robust
				if (go_to_robust > 0) {
					// if came from robust
					// copy two most isolated roots as starting points for new robust
					for (i = 1; i <= degree - 3; i++) {
						roots[degree - i] = roots_robust[i - 1];
					}
				}
				else {
					// came from users intial guess
					// copy some 2 roots (except closest ones)
					i2 = degree;
					for (i = 1; i <= degree; i++) {
						if ((i != root4) && (i != root5)) {
							roots[i2 - 1] = roots[i - 1];
							i2 -= 1;
						}
						if (i2 <= 3) break; // do not copy those that will be done by cubic in robust
					}
				}
				go_to_robust += 1;
			}
			else {
				// root4 and root5 come from the initial closest pair
				// most common case
				return;
			}
		}
	}
	return;
}

void  sort_5_points_by_separation(complex *points) {
	//Sort array of five points
	//Most isolated point will become the first point in the array
	//The closest points will be the last two points in the array

	//Algorithm works well for all dimensions. We put n = 5 as
	//a hardcoded value just for optimization purposes

	const int n = 5; //works for different n as well, but is fster for n as constant (optimization)
	int sorted_points[n];


	complex savepoints[n];
	int i;

	sort_5_points_by_separation_i(sorted_points, points);
	for (i = 0; i < n; i++) {
		savepoints[i] = points[i];
	}
	for (i = 1; i <= n; i++) {
		points[i - 1] = savepoints[sorted_points[i - 1] - 1];
	}
	return;
}

void sort_5_points_by_separation_i(int *sorted_points, complex *points) {
	//Return index array that sorts array of five points 
	//Index of the most isolated point will appear on the firts place
	//of the output array.
	//The indices of the closest 2 points will be at the last two
	//places in the 'sorted_points' array

	//Algorithm works well for all dimensions. We put n=5 as 
	//a hardcoded value just for optimization purposes.
	const int n = 5;

	double d1, d2, d;
	double distances2[n][n];
	int ki, kj, ind2, put;
	double neigh1st[n], neigh2nd[n];
	complex p;

	for (kj = 1; kj <= n; kj++) {
		for (ki = 1; ki <= kj - 1; ki++) {
			p = points[ki - 1] - points[kj - 1];
			d = real(conj(p) * p);
			distances2[ki - 1][kj - 1] = d;
			distances2[kj - 1][ki - 1] = d;
		}
	}
	//find neighbours

	for (kj = 1; kj <= n; kj++) {
		for (ki = 1; ki <= kj - 1; ki++) {
			d = distances2[kj - 1][ki - 1];
			if (d < neigh2nd[kj - 1]) {
				if (d < neigh1st[kj - 1]) {
					neigh2nd[kj - 1] = neigh1st[kj - 1];
					neigh1st[kj - 1] = d;
				}
				else {
					neigh2nd[kj - 1] = d;
				}
			}
		}

	}

	//initialize sorted_points

	for (ki = 1; ki <= n; ki++) {
		sorted_points[ki - 1] = ki;
	}

	//sort the rest 1..n-2
	for (kj = 2; kj <= n; kj++) {
		d1 = neigh1st[kj - 1];
		d2 = neigh2nd[kj - 1];
		put = 1;
		for (ki = kj - 1; ki >= 1; ki--) {
			ind2 = sorted_points[ki - 1];
			d = neigh1st[ind2 - 1];
			if (d >= d1) {
				if (d == d1) {
					if (neigh2nd[ind2 - 1] > d2) {
						put = ki + 1;
						break;
					}
				}
				else {
					put = ki + 1;
					break;
				}
			}
			sorted_points[ki] = sorted_points[ki - 1];
		}
		sorted_points[put - 1] = kj;
	}
	return;
}

void find_2_closest_from_5(int *i1, int *i2, double *d2min, complex *points) {
	//Returns indices of the two cloest points out of array of 5
	//d2min is the square of minimal distance
	int n = 5;
	//double distances2[n][n];

	double d2min1, d2;
	int i, j;
	complex p;
	d2min1 = 1.0;
	for (j = 1; j <= n; j++) {
		for (i = 1; i <= j - 1; i++) {
			p = points[i - 1] - points[j - 1];
			d2 = real(conj(p) * p);
			//distances2[i][j] = d2;
			//distances2[j][i] = d2;
			if (d2 <= d2min1) {
				*i1 = i;
				*i2 = j;
				d2min1 = d2;
			}
		}
	}
	*d2min = d2min1;
}

void cmplx_laguerre(complex *poly, int degree, complex *root, int &iter, bool &success) {
	//Subroutine finds one root of a complex polynomial using
	//Laguerre's method. In every loop it calculates simplified 
	//Adams' stopping criterion for the value of the polynomial.
	//
	//Uses 'root' value as a starting point(!!!!!)
	//Remember to initialize 'root' to some initial guess or to
	//point(0, 0) if you have no prior knowledge.
	//
	//poly - is an array of polynomial cooefs
	//
	//length = degree + 1, poly(1) is constant
	//	1              2				3
	//poly(1) x ^ 0 + poly(2) x ^ 1 + poly(3) x ^ 2 + ...
	//
	//degree - a degree of the polynomial
	//
	//root - input: guess for the value of a root
	//output : a root of the polynomial
	//iter - number of iterations performed(the number of polynomial
	//evaluations and stopping criterion evaluation)
	//
	//success - is false if routine reaches maximum number of iterations
	//
	//For a summary of the method go to :
	//http://en.wikipedia.org/wiki/Laguerre's_method
	//
	static int FRAC_JUMP_EVERY = 10;
	const int FRAC_JUMP_LEN = 10;
	double FRAC_JUMPS[FRAC_JUMP_LEN] = { 0.64109297,
		0.91577881, 0.25921289, 0.50487203,
		0.08177045, 0.13653241, 0.306162,
		0.37794326, 0.04618805, 0.75132137 }; // some random numbers

	double faq; //jump length
	double FRAC_ERR = 2.0e-15; //Fractional Error for double precision
	complex p, dp, d2p_half; //value of polynomial, 1st derivative, and 2nd derivative
	static int i, j, k;
	bool good_to_go;
	complex denom, denom_sqrt, dx, newroot;
	double ek, absroot, abs2p;
	complex fac_newton, fac_extra, F_half, c_one_nth;
	double one_nth, n_1_nth, two_n_div_n_1;
	complex c_one = complex(1, 0);
	complex zero = complex(0, 0);
	double stopping_crit2;

	//--------------------------------------------------------------------------------------------

	//EXTREME FAILSAFE! not usually needed but kept here just to be on the safe side. Takes care of first coefficient being 0
	if (true) {
		if (degree < 0) {
			printf("Error: cmplx_laguerre: degree<0");
			return;
		}
		if (poly[degree] == complex(0, 0)) {
			if (degree == 0) return;
			cmplx_laguerre(poly, degree - 1, root, iter, success);
		}
		if (degree <= 1) {
			if (degree == 0) {
				success = false; // we just checked if poly[0] is zero and it isnt
				printf("Warning: cmplx_laguerre: degree = 0 and poly[0] does not equal zero, no roots");
				return;
			}
			else {
				*root = -poly[0] / poly[1];
				return;
			}
		}
	} // End of EXTREME failsafe

	good_to_go = false;
	one_nth = 1.0 / degree;
	n_1_nth = (degree - 1.0)*one_nth;
	two_n_div_n_1 = 2.0 / n_1_nth;
	c_one_nth = complex(one_nth, 0.0);
	for (i = 1; i <= MAXIT; i++) {
		ek = abs(poly[degree]); // Preparing stopping criterion
		absroot = abs(*root);
		// Calculate the values of polynomial and its first and second derivatives
		p = poly[degree];
		dp = zero;
		d2p_half = zero;
		for (k = degree - 1; k >= 0; k--) {
			d2p_half = dp + d2p_half*(*root);
			dp = p + dp * *root;
			p = poly[k] + p*(*root); // b_k
									 //Adams, Duane A., 1967, "A stopping criterion for polynomial root finding",
									 //Communications of the ACM, Volume 10 Issue 10, Oct. 1967, p. 655
									 //ftp://reports.stanford.edu/pub/cstr/reports/cs/tr/67/55/CS-TR-67-55.pdf
									 //Eq 8.
			ek = absroot*ek + abs(p);
		}
		iter += 1;

		abs2p = real(conj(p)*p);
		if (abs2p == 0) return;
		stopping_crit2 = pow(FRAC_ERR*ek, 2.0);
		if (abs2p < stopping_crit2) {
			//(simplified a little Eq. 10 of Adams 1967)
			//do additional iteration if we are less than 10x from stopping criterion
			if (abs2p < 0.01*stopping_crit2) {
				return; // we are at a good place!
			}
			else {
				good_to_go = true;
			}
		}
		else {
			good_to_go = false;
		}

		faq = 1.0;
		denom = zero;
		if (dp != zero) {
			fac_newton = p / dp;
			fac_extra = d2p_half / dp;
			F_half = fac_newton*fac_extra;
			denom_sqrt = sqrt(c_one - two_n_div_n_1*F_half);

			//NEXT LINE PROBABLY CAN BE COMMENTED OUT. Check if compiler outputs positive real
			if (real(denom_sqrt) >= 0.0) {
				denom = c_one_nth + n_1_nth*denom_sqrt;
			}
			else {
				denom = c_one_nth - n_1_nth*denom_sqrt;
			}
		}

		if (denom == 0) {
			dx = (absroot + 1.0)*expcmplx(complex(0.0, FRAC_JUMPS[i % FRAC_JUMP_LEN] * 2 * M_PI));
		}
		else {
			dx = fac_newton / denom;
		}


		newroot = *root - dx;
		if (newroot == *root) return; //nothing changes so return
		if (good_to_go) {
			*root = newroot;
			return;
		}
		if (i % FRAC_JUMP_EVERY == 0) { //decide whether to do a jump of modified length (to break cycles)
			faq = FRAC_JUMPS[(i / FRAC_JUMP_EVERY - 1) % FRAC_JUMP_LEN];
			newroot = *root - faq*dx; // do jump of semi-random length
		}
		*root = newroot;
	}
	success = false; // too many iterations here
	return;
}

void cmplx_newton_spec(complex *poly, int degree, complex *root, int &iter, bool &success) {
	//Subroutine finds one root of a complex polynomial
	//Newton's method. It calculates simplified Adams' stopping 
	//criterion for the value of the polynomial once per 10 iterations (!),
	//after initial iteration. This is done to speed up calculations
	//when polishing roots that are known preety well, and stopping
	// criterion does significantly change in their neighborhood.

	//Uses 'root' value as a starting point (!!!!!)
	//Remember to initialize 'root' to some initial guess.
	//Do not initilize 'root' to point (0,0) if the polynomial 
	//coefficients are strictly real, because it will make going 
	//to imaginary roots impossible.

	// poly - is an array of polynomial cooefs
	//	length = degree+1, poly(1) is constant 
	//0					1				2
	//poly[0] x^0 + poly[1] x^1 + poly[2] x^2 + ...
	//degree - a degree of the polynomial
	// root - input: guess for the value of a root
	//		  output: a root of the polynomial
	//iter - number of iterations performed (the number of polynomial evaluations)
	//success - is false if routine reaches maximum number of iterations

	//For a summary of the method go to: 
	//http://en.wikipedia.org/wiki/Newton's_method

	int FRAC_JUMP_EVERY = 10;
	const int FRAC_JUMP_LEN = 10;
	double FRAC_JUMPS[FRAC_JUMP_LEN] = { 0.64109297, 0.91577881, 0.25921289, 0.50487203, 0.08177045, 0.13653241, 0.306162, 0.37794326, 0.04618805, 0.75132137 }; //some random numbers
	double faq; //jump length
	double FRAC_ERR = 2e-15;
	complex p; //value of polynomial
	complex dp; //value of 1st derivative
	int i, k;
	bool good_to_go;
	complex dx, newroot;
	double ek, absroot, abs2p;
	complex zero = complex(0, 0);
	double stopping_crit2;

	iter = 0;
	success = true;

	//the next if block is an EXTREME failsafe, not usually needed, and thus turned off in this version
	if (true) { //change false to true if you would like to use caustion about haveing first coefficient == 0
		if (degree < 0) {
			printf("Error: cmplx_newton_spec: degree<0");
			return;
		}
		if (poly[degree] == zero) {
			if (degree == 0) return;
			cmplx_newton_spec(poly, degree, root, iter, success);
			return;
		}
		if (degree <= 1) {
			if (degree == 0) {
				success = false;
				printf("Warning: cmplx_newton_spec: degree=0 and poly[0]!=0, no roots");
				return;
			}
			else {
				*root = -poly[0] / poly[1];
				return;
			}
		}
	}
	//end EXTREME Failsafe
	good_to_go = false;

	stopping_crit2 = 0.0; //value not important, will be initialized anyway on the first loop
	for (i = 1; i <= MAXIT; i++) {
		faq = 1.0;
		//prepare stoping criterion
		//calculate value of polynomial and its first two derivatives
		p = poly[degree];
		dp = zero;
		if (i % 10 == 1) { //calculate stopping criterion every tenth iteration
			ek = abs(poly[degree]);
			absroot = abs(*root);
			for (k = degree - 1; k >= 0; k--) {
				dp = p + dp * (*root);
				p = poly[k] + p * (*root); //b_k
										   //Adams, Duane A., 1967, "A stopping criterion for polynomial root finding",
										   //Communications of ACM, Volume 10 Issue 10, Oct. 1967, p. 655
										   //ftp://reports.stanford.edu/pub/cstr/reports/cs/tr/67/55/CS-TR-67-55.pdf
										   //Eq. 8
				ek = absroot * ek + abs(p);
			}
			stopping_crit2 = pow(FRAC_ERR * ek, 2);
		}
		else { // calculate just the value and derivative
			for (k = degree - 1; k >= 0; k--) { //Horner Scheme, see for eg. Numerical Recipes Sec. 5.3 how to evaluate polynomials and derivatives
				dp = p + dp * (*root);
				p = poly[k] + p * (*root);
			}
		}

		iter = iter + 1;

		abs2p = real(conj(p) * p);
		if (abs2p == 0.0) return;
		if (abs2p < stopping_crit2) { //simplified a little Eq. 10 of Adams 1967
			if (dp == zero) return; //if we have problem with zero, but we are close to the root, just accept
									//do additional iteration if we are less than 10x from stopping criterion
			if (abs2p < 0.01 * stopping_crit2) return; //return immediatley because we are at very good place
			else {
				good_to_go = true; //do one iteration more
			}
		}

		else {
			good_to_go = false; //reset if we are outside the zone of the root
		}
		if (dp == zero) {
			//problem with zero
			dx = (abs(*root) + 1.0) * expcmplx(complex(0.0, FRAC_JUMPS[i% FRAC_JUMP_LEN] * 2 * M_PI));
		}
		else {
			dx = p / dp; // Newton method, see http://en.wikipedia.org/wiki/Newton's_method
		}
		newroot = *root - dx;
		if (newroot == *root) return; //nothing changes -> return
		if (good_to_go) {//this was jump already after stopping criterion was met
			*root = newroot;
			return;
		}
		if (i % FRAC_JUMP_EVERY == 0) { // decide whether to do a jump of modified length (to break cycles)
			faq = FRAC_JUMPS[(i / FRAC_JUMP_EVERY - 1) % FRAC_JUMP_LEN];
			newroot = *root - faq * dx;
		}
		*root = newroot;
	}
	success = false;
	return;
	//too many iterations here
}

void cmplx_laguerre2newton(complex *poly, int degree, complex *root, int &iter, bool &success, int starting_mode) {
	//Subroutine finds one root of a complex polynomial using
	//Laguerre's method, Second-order General method and Newton's
	//method - depending on the value of function F, which is a 
	//combination of second derivative, first derivative and
	//value of polynomial [F=-(p"*p)/(p'p')].

	//Subroutine has 3 modes of operation. It starts with mode=2
	//which is the Laguerre's method, and continues until F
	//becames F<0.50, at which point, it switches to mode=1,
	//i.e., SG method (see paper). While in the first two
	//modes, routine calculates stopping criterion once per every
	//iteration. Switch to the last mode, Newton's method, (mode=0)
	//happens when becomes F<0.05. In this mode, routine calculates
	//stopping criterion only once, at the beginning, under an
	//assumption that we are already very close to the root.
	//If there are more than 10 iterations in Newton's mode,
	//it means that in fact we were far from the root, and
	//routine goes back to Laguerre's method (mode=2).

	//Uses 'root' value as a starting point (!!!!!)
	//Remember to initialize 'root' to some initial guess or to 
	//point (0,0) if you have no prior knowledge.

	//poly - is an array of polynomial cooefs
	//	0					1				2
	//	poly[0] x^0 + poly[1] x^1 + poly[2] x^2
	//degree - a degree of the polynomial
	//root - input: guess for the value of a root
	//		output: a root of the polynomial
	//iter - number of iterations performed (the number of polynomial
	//		 evaluations and stopping criterion evaluation)
	//success - is false if routine reaches maximum number of iterations
	//starting_mode - this should be by default = 2. However if you  
	//				  choose to start with SG method put 1 instead.
	//				  Zero will cause the routine to
	//				  start with Newton for first 10 iterations, and
	//				  then go back to mode 2.

	//For a summary of the method see the paper: Skowron & Gould (2012)

	int FRAC_JUMP_EVERY = 10;
	const int FRAC_JUMP_LEN = 10;
	double FRAC_JUMPS[FRAC_JUMP_LEN] = { 0.64109297, 0.91577881, 0.25921289, 0.50487203, 0.08177045, 0.13653241, 0.306162, 0.37794326, 0.04618805, 0.75132137 }; //some random numbers

	double faq; //jump length
	double FRAC_ERR = 2.0e-15;

	complex p; //value of polynomial
	complex dp; //value of 1st derivative
	complex d2p_half; //value of 2nd derivative
	int i, j, k;
	bool good_to_go;
	//complex G, H, G2;
	complex denom, denom_sqrt, dx, newroot;
	double ek, absroot, abs2p, abs2_F_half;
	complex fac_netwon, fac_extra, F_half, c_one_nth;
	double one_nth, n_1_nth, two_n_div_n_1;
	int mode;
	complex c_one = complex(1, 0);
	complex zero = complex(0, 0);
	double stopping_crit2;

	iter = 0;
	success = true;
	stopping_crit2 = 0; //value not important, will be initialized anyway on the first loop

						//next if block is an EXTREME failsafe, not usually needed, and thus turned off in this version.
	if (false) {//change false to true if you would like to use caution about having first coefficent == 0
		if (degree < 0) {
			printf("Error: cmplx_laguerre2newton: degree < 0");
			return;
		}
		if (poly[degree] == zero) {
			if (degree == 0) return;
			cmplx_laguerre2newton(poly, degree, root, iter, success, starting_mode);
			return;
		}
		if (degree <= 1) {
			if (degree == 0) {//// we know from previous check that poly[0] not equal zero
				success = false;
				printf("Warning: cmplx_laguerre2newton: degree = 0 and poly[0] = 0, no roots");
				return;
			}
			else {
				*root = -poly[0] / poly[1];
				return;
			}
		}
	}
	//end EXTREME failsafe

	j = 1;
	good_to_go = false;

	mode = starting_mode; // mode = 2 full laguerre, mode = 1 SG, mode = 0 newton

	for (;;) { //infinite loop, just to be able to come back from newton, if more than 10 iteration there

			   ////////////
			   ///mode 2///
			   ////////////

		if (mode >= 2) {//Laguerre's method
			one_nth = 1.0 / (degree); ///
			n_1_nth = (degree - 1) * one_nth; ////
			two_n_div_n_1 = 2.0 / n_1_nth;
			c_one_nth = complex(one_nth, 0.0);

			for (i = 1; i <= MAXIT; i++) {
				faq = 1.0;

				//prepare stoping criterion
				ek = abs(poly[degree]);
				absroot = abs(*root);
				//calculate value of polynomial and its first two derivative
				p = poly[degree];
				dp = zero;
				d2p_half = zero;
				for (k = degree; k >= 1; k--) {//Horner Scheme, see for eg.  Numerical Recipes Sec. 5.3 how to evaluate polynomials and derivatives
					d2p_half = dp + d2p_half * (*root);
					dp = p + dp * (*root);
					p = poly[k - 1] + p * (*root); // b_k
												   //Adams, Duane A., 1967, "A stopping criterion for polynomial root finding",
												   //Communications of the ACM, Volume 10 Issue 10, Oct. 1967, p. 655
												   //ftp://reports.stanford.edu/pub/cstr/reports/cs/tr/67/55/CS-TR-67-55.pdf
												   //Eq 8.
					ek = absroot * ek + abs(p);
				}
				abs2p = real(conj(p) * p); // abs(p)
				iter = iter + 1;
				if (abs2p == 0) return;

				stopping_crit2 = pow(FRAC_ERR * ek, 2);
				if (abs2p < stopping_crit2) {//(simplified a little Eq. 10 of Adams 1967)
											 //do additional iteration if we are less than 10x from stopping criterion
					if (abs2p < 0.01*stopping_crit2) return; // ten times better than stopping criterion
															 //return immediately, because we are at very good place
					else {
						good_to_go = true; //do one iteration more
					}
				}
				else {
					good_to_go = false; //reset if we are outside the zone of the root
				}

				denom = zero;
				if (dp != zero) {
					fac_netwon = p / dp;
					fac_extra = d2p_half / dp;
					F_half = fac_netwon * fac_extra;

					abs2_F_half = real(conj(F_half) * F_half);
					if (abs2_F_half <= 0.0625) {//F<0.50, F/2<0.25
												//go to SG method
						if (abs2_F_half <= 0.000625) {//F<0.05, F/2<0.02
							mode = 0; //go to Newton's
						}
						else {
							mode = 1; //go to SG
						}
					}

					denom_sqrt = sqrt(c_one - two_n_div_n_1*F_half);

					//NEXT LINE PROBABLY CAN BE COMMENTED OUT 
					if (real(denom_sqrt) > 0.0) {
						//real part of a square root is positive for probably all compilers. You can ù
						//test this on your compiler and if so, you can omit this check
						denom = c_one_nth + n_1_nth * denom_sqrt;
					}
					else {
						denom = c_one_nth - n_1_nth * denom_sqrt;
					}
				}
				if (denom == zero) {//test if demoninators are > 0.0 not to divide by zero
					dx = (abs(*root) + 1.0) + expcmplx(complex(0.0, FRAC_JUMPS[i% FRAC_JUMP_LEN] * 2 * M_PI)); //make some random jump
				}
				else {
					dx = fac_netwon / denom;
				}
				newroot = *root - dx;
				if (newroot == *root) return; // nothing changes -> return
				if (good_to_go) {//this was jump already after stopping criterion was met
					*root = newroot;
					return;
				}
				if (mode != 2) {
					*root = newroot;
					j = i + 1; //remember iteration index
					break; //go to Newton's or SG
				}
				if ((i% FRAC_JUMP_EVERY) == 0) {//decide whether to do a jump of modified length (to break cycles)
					faq = FRAC_JUMPS[((i / FRAC_JUMP_EVERY - 1) % FRAC_JUMP_LEN)];
					newroot = *root - faq * dx; // do jump of some semi-random length (0 < faq < 1)
				}
				*root = newroot;
			} //do mode 2

			if (i >= MAXIT) {
				success = false;
				return;
			}
		}

		////////////
		///mode 1///
		////////////

		if (mode == 1) {//SECOND-ORDER GENERAL METHOD (SG)

			for (i = j; i <= MAXIT; i++) {
				faq = 1.0;
				//calculate value of polynomial and its first two derivatives
				p = poly[degree];
				dp = zero;
				d2p_half = zero;
				if ((i - j) % 10 == 0) {
					//prepare stopping criterion
					ek = abs(poly[degree]);
					absroot = abs(*root);
					for (k = degree; k >= 1; k--) {//Horner Scheme, see for eg.  Numerical Recipes Sec. 5.3 how to evaluate polynomials and derivatives
						d2p_half = dp + d2p_half * (*root);
						dp = p + dp * (*root);
						p = poly[k - 1] + p * (*root); //b_k
													   //Adams, Duane A., 1967, "A stopping criterion for polynomial root finding",
													   //Communications of the ACM, Volume 10 Issue 10, Oct. 1967, p. 655
													   //ftp://reports.stanford.edu/pub/cstr/reports/cs/tr/67/55/CS-TR-67-55.pdf
													   //Eq 8.
						ek = absroot * ek + abs(p);
					}
					stopping_crit2 = pow(FRAC_ERR*ek, 2);
				}
				else {
					for (k = degree; k >= 1; k--) {//Horner Scheme, see for eg.  Numerical Recipes Sec. 5.3 how to evaluate polynomials and derivatives
						d2p_half = dp + d2p_half * (*root);
						dp = p + dp * (*root);
						p = poly[k - 1] + p * (*root); //b_k
					}
				}
				abs2p = real(conj(p) * p); //abs(p)**2
				iter = iter + 1;
				if (abs2p == 0.0) return;

				if (abs2p < stopping_crit2) {//(simplified a little Eq. 10 of Adams 1967)
					if (dp == zero) return;
					//do additional iteration if we are less than 10x from stopping criterion
					if (abs2p < 0.01*stopping_crit2) return; //ten times better than stopping criterion
															 //ten times better than stopping criterion
					else {
						good_to_go = true; //do one iteration more
					}
				}
				else {
					good_to_go = false; //reset if we are outside the zone of the root
				}
				if (dp == zero) {//test if denominators are > 0.0 not to divide by zero
					dx = (abs(*root) + 1.0) * expcmplx(complex(0.0, FRAC_JUMPS[i% FRAC_JUMP_LEN] * 2 * M_PI)); //make some random jump
				}
				else {
					fac_netwon = p / dp;
					fac_extra = d2p_half / dp;
					F_half = fac_netwon * fac_extra;

					abs2_F_half = real(conj(F_half) * F_half);
					if (abs2_F_half <= 0.000625) {//F<0.05, F/2<0.025
						mode = 0; //set Newton's, go there after jump
					}
					dx = fac_netwon * (c_one + F_half); //SG
				}
				newroot = *root - dx;
				if (newroot == *root) return; //nothing changes -> return
				if (good_to_go) {
					*root = newroot; //this was jump already after stopping criterion was met
					return;
				}
				if (mode != 1) {
					*root = newroot;
					j = i + 1; //remember iteration number
					break; //go to Newton's
				}
				if ((i% FRAC_JUMP_EVERY) == 0) {// decide whether to do a jump of modified length (to break cycles)
					faq = FRAC_JUMPS[(i / FRAC_JUMP_EVERY - 1) % FRAC_JUMP_LEN];
					newroot = *root - faq * dx; //do jump of some semi random lenth (0 < faq < 1)		
				}
				*root = newroot;
			}
			if (i >= MAXIT) {
				success = false;
				return;
			}

		}


		////////////
		///mode 0///
		////////////


		if (mode == 0) { // Newton's Method

			for (i = j; i <= j + 10; i++) { // Do only 10 iterations the most then go back to Laguerre
				faq = 1.0;

				//calc polynomial and first two derivatives
				p = poly[degree];
				dp = zero;
				if (i == j) { // Calculating stopping criterion only at the beginning
					ek = abs(poly[degree]);
					absroot = abs(*root);
					for (k = degree; k >= 1; k--) {
						dp = p + dp*(*root);
						p = poly[k - 1] + p*(*root);
						ek = absroot*ek + abs(p);
					}
					stopping_crit2 = pow(FRAC_ERR*ek, 2.0);
				}
				else {
					for (k = degree; k >= 1; k--) {
						dp = p + dp*(*root);
						p = poly[k - 1] + p*(*root);
					}
				}
				abs2p = real(conj(p)*p);
				iter = iter + 1;
				if (abs2p == 0.0) return;

				if (abs2p < stopping_crit2) {
					if (dp == zero) return;
					// do additional iteration if we are less than 10x from stopping criterion
					if (abs2p < 0.01*stopping_crit2) {
						return; // return immediately since we are at a good place
					}
					else {
						good_to_go = true; // do one more iteration
					}
				}
				else {
					good_to_go = false;
				}

				if (dp == zero) {
					dx = (abs(*root) + 1.0)*expcmplx(complex(0.0, 2 * M_PI*FRAC_JUMPS[i % FRAC_JUMP_LEN])); // make a random jump
				}
				else {
					dx = p / dp;
				}

				newroot = *root - dx;
				if (newroot == *root) return;
				if (good_to_go) {
					*root = newroot;
					return;
				}
				*root = newroot;
			}
			if (iter >= MAXIT) {
				//too many iterations
				success = false;
				return;
			}
			mode = 2; //go back to Laguerre's. Happens when could not converge with 10 steps of Newton
		}

	}/// end of infinite loop
}

void solve_quadratic_eq(complex &x0, complex &x1, complex *poly) {
	complex a, b, c, b2, delta;
	a = poly[2];
	b = poly[1];
	c = poly[0];
	b2 = b*b;
	delta = sqrt(b2 - 4 * a*c);
	if (real(conj(b)*delta) >= 0) {
		x0 = -0.5*(b + delta);
	}
	else {
		x0 = -0.5*(b - delta);
	}
	if (x0 == complex(0., 0.)) {
		x1 = complex(0., 0.);
	}
	else { //Viete's formula
		x1 = c / x0;
		x0 = x0 / a;
	}
	return;

}

void solve_cubic_eq(complex &x0, complex &x1, complex &x2, complex *poly) {
	//Cubic equation solver for comples polynomial (degree=3)
	//http://en.wikipedia.org/wiki/Cubic_function   Lagrange's method
	// poly is an array of polynomial cooefs, length = degree+1, poly[0] is constant
	//	0				1				2			3
	//poly[0] x^0 + poly[1] x^1 + poly[2] x^2 + poly[3] x^3
	complex zeta = complex(-0.5, 0.8660254037844386);
	complex zeta2 = complex(-0.5, -0.8660254037844386);
	double third = 0.3333333333333333;
	complex s0, s1, s2;
	complex E1; //x0+x1+x2
	complex E2; //x0*x1+x1*x2+x2*x0
	complex E3; //x0*x1*x2
	complex A, B, a_1, E12, delta, A2;

	complex val, x;
	a_1 = 1 / poly[3];
	E1 = -poly[2] * a_1;
	E2 = poly[1] * a_1;
	E3 = -poly[0] * a_1;

	s0 = E1;
	E12 = E1*E1;
	A = 2.0 * E1 * E12 - 9.0 * E1 * E2 + 27.0 * E3;
	B = E12 - 3.0 * E2;
	//quadratic equation z^2 - A * z + B^3 where roots are equal to s1^3 and s2^3
	A2 = A * A;
	delta = sqrt(A2 - 4.0 * (B * B * B));
	if (real(conj(A) * delta) >= 0.0) { // scalar product to decide the sign yielding bigger magnitude
		s1 = cbrt(0.5 * (A + delta));
	}
	else
	{
		s1 = cbrt(0.5 * (A - delta));
	}
	if (s1.re == 0.0 && s1.im == 0.0) {
		s2 = complex(0, 0);
	}
	else {
		s2 = B / s1;
	}

	x0 = third * (s0 + s1 + s2);
	x1 = third * (s0 + s1 * zeta2 + s2 * zeta);
	x2 = third * (s0 + s1 * zeta + s2 * zeta2);

	return;

}

void divide_poly_1(complex *polyout, complex remainder, complex p, complex *polyin, int degree) {
	//Subroutine will divide polynomial 'polyin' by(x - p)
	//results will be returned in polynomial 'polyout' of degree - 1
	//The remainder of the division will be returned in 'remainder'
	//
	//You can provide same array as 'polyin' and 'polyout' - this
	//routine will work fine, though it will not set to zero the
	//unused, highest coefficient in the output array.You just have
	//remember the proper degree of a polynomial.
	//poly - is an array of polynomial cooefs, length = degree + 1, poly(1) is constant
	//1              2             3
	//poly(1) x ^ 0 + poly(2) x ^ 1 + poly(3) x ^ 2 + ...

	int i;
	complex coef, prev;
	coef = polyin[degree];
	for (int k = 1; k <= degree; k++) {
		polyout[k - 1] = polyin[k - 1];
	}
	for (i = degree; i >= 1; i--) {
		prev = polyout[i - 1];
		polyout[i - 1] = coef;
		coef = prev + p*coef;
	}
	remainder = coef;
	return;
}

complex::complex(double a, double b) {
	re = a;
	im = b;
}

complex::complex(double a) {
	re = a;
	im = 0;
}

complex::complex(void) {
	re = 0;
	im = 0;
}

double abs(complex z) {
	return sqrt(z.re*z.re + z.im*z.im);
}

complex conj(complex z) {
	return complex(z.re, -z.im);
}

complex sqrt(complex z) {
	double md = sqrt(z.re*z.re + z.im*z.im);
	return (md>0) ? sqrt(md)*complex((sqrt((1 + z.re / md) / 2)*((z.im>0) ? 1 : -1)), sqrt((1 - z.re / md) / 2)) : 0.0;
}

complex cbrt(complex z) {
	complex zout;
	double r, r_cube, theta, theta_cube;
	r = abs(z);
	r_cube = pow(r, 1. / 3.);
	theta = atan2(z.im, z.re);
	theta_cube = theta / 3.;
	zout = complex(r_cube*cos(theta_cube), r_cube*sin(theta_cube));
	return zout;
}

double real(complex z) {
	return z.re;
}

double imag(complex z) {
	return z.im;
}

bool operator==(complex p1, complex p2) {
	if (p1.re == p2.re && p1.im == p2.im) return true;
	return false;
}

bool operator!=(complex p1, complex p2) {
	if (p1.re == p2.re && p1.im == p2.im) return false;
	return true;
}

complex expcmplx(complex p1) {
	double r = exp(p1.re);
	double theta = atan2(p1.im, p1.re);
	return complex(r*cos(theta), r*sin(theta));
}

complex operator+(complex p1, complex p2) {
	return complex(p1.re + p2.re, p1.im + p2.im);
}

complex operator-(complex p1, complex p2) {
	return complex(p1.re - p2.re, p1.im - p2.im);
}

complex operator*(complex p1, complex p2) {
	return complex(p1.re*p2.re - p1.im*p2.im, p1.re*p2.im + p1.im*p2.re);
}

complex operator/(complex p1, complex p2) {
	double md = p2.re*p2.re + p2.im*p2.im;
	return complex((p1.re*p2.re + p1.im*p2.im) / md, (p1.im*p2.re - p1.re*p2.im) / md);
}

complex operator+(complex z, double a) {
	return complex(z.re + a, z.im);
}

complex operator-(complex z, double a) {
	return complex(z.re - a, z.im);
}

complex operator*(complex z, double a) {
	return complex(z.re*a, z.im*a);
}

complex operator/(complex z, double a) {
	return complex(z.re / a, z.im / a);
}

complex operator+(double a, complex z) {
	return complex(z.re + a, z.im);
}

complex operator-(double a, complex z) {
	return complex(a - z.re, -z.im);
}

complex operator*(double a, complex z) {
	return complex(a*z.re, a*z.im);
}

complex operator/(double a, complex z) {
	double md = z.re*z.re + z.im*z.im;
	return complex(a*z.re / md, -a*z.im / md);
}

complex operator+(complex z, int a) {
	return complex(z.re + a, z.im);
}

complex operator-(complex z, int a) {
	return complex(z.re - a, z.im);
}

complex operator*(complex z, int a) {
	return complex(z.re*a, z.im*a);
}

complex operator/(complex z, int a) {
	return complex(z.re / a, z.im / a);
}

complex operator+(int a, complex z) {
	return complex(z.re + a, z.im);
}

complex operator-(int a, complex z) {
	return complex(a - z.re, -z.im);
}

complex operator*(int a, complex z) {
	return complex(a*z.re, a*z.im);
}

complex operator/(int a, complex z) {
	double md = z.re*z.re + z.im*z.im;
	return complex(a*z.re / md, -a*z.im / md);
}

complex operator-(complex z) {
	return complex(-z.re, -z.im);
}