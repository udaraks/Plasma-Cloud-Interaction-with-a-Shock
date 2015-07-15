// PH789 Final Project - Udara Senanayake

#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstring>
#include "silo.h"

using namespace std;
#define max1(x,y) (x > y ? x : y)

template <class T> T **Create2D(int n, int m)
{
   long int size = (long int)n * (long int)m + 1;
   T **array = new T *[n + 1];
   array[0] = new T[size]; array[1] = array[0];
   for(int i = 2; i <= n; i++) array[i] = array[i - 1] + m;
   return array;
}
template <class T> void Delete2D(T **array) {delete[] array[0]; delete[] array;}

void WriteSilo(double *data1, double *data2, double *data3, double *data4, double *data5, int cycle);

const double c = 0.2; //cfl factor
const double tmax = 0.045;
const int Nx = 256;
const int Ny = Nx;
double x[Nx], y[Ny];
int i, j, k, l, m, b ,z;
double dx = 1.0 / double(Nx - 1);
double dy = dx;
double dt;
const double rho1 = 3.9;
const double rho2 = 1.0847336750612;
const double p1 = 167.0;
const double p2 = 5.05;
const double u1 = 0.0;
const double u2 = -10.381410256410;
const double v1 = 0.0;
const double v2 = 0.0;
const double w1 = 0.0;
const double w2 = 0.0;
const double gama = 5.0 / 3.0;
const double epsilon = 1E-06;
double **af = Create2D<double>(Nx, Ny);
double **umax0 = Create2D<double>(Nx, Ny);
double **Ax = Create2D<double>(5, 5);
double **Ay = Create2D<double>(5, 5);
double **U[6],**Ux[6],**Uy[6], **ULx[6], **URx[6],**ULy[6], **URy[6], **Fx[6], **Fy[6], **Fy1[6],**FLx[6], **FRx[6], **FLy[6], **FRy[6], **deltax[6], **deltaAx[6], **deltay[6], **deltaAy[6], **var[6], **var0[6], **varLx[6], **varRx[6], **varLy[6], **varRy[6];
double e1, e2, umax, umax1, umax2, af1, af2, t, uh, vh, wh, hh, hL, hR, eL, eR, rho[Nx*Ny], u[Nx*Ny], v[Nx*Ny], w[Nx*Ny], p[Nx*Ny];

void init()
{
	for(i = 1; i <= 5; i++)
	{
		U[i] =  Create2D<double>(Nx, Ny);
		Ux[i] =  Create2D<double>(Nx, Ny);
		Uy[i] =  Create2D<double>(Nx, Ny);
		ULx[i] =  Create2D<double>(Nx, Ny);
		URx[i] =  Create2D<double>(Nx, Ny);
		ULy[i] =  Create2D<double>(Nx, Ny);
		URy[i] =  Create2D<double>(Nx, Ny);
		Fx[i] =  Create2D<double>(Nx, Ny);
		Fy[i] =  Create2D<double>(Nx, Ny);
		Fy1[i] =  Create2D<double>(Nx, Ny);
		FLx[i] =  Create2D<double>(Nx, Ny);
		FRx[i] =  Create2D<double>(Nx, Ny);
		FLy[i] =  Create2D<double>(Nx, Ny);
		FRy[i] =  Create2D<double>(Nx, Ny);
		deltax[i] =  Create2D<double>(Nx, Ny);
		deltaAx[i] =  Create2D<double>(Nx, Ny);
		deltay[i] =  Create2D<double>(Nx, Ny);
		deltaAy[i] =  Create2D<double>(Nx, Ny);
		var[i] =  Create2D<double>(Nx, Ny);
		var0[i] =  Create2D<double>(Nx, Ny);
		varLx[i] =  Create2D<double>(Nx, Ny);	
		varRx[i] =  Create2D<double>(Nx, Ny);
		varLy[i] =  Create2D<double>(Nx, Ny);	
		varRy[i] =  Create2D<double>(Nx, Ny);	
	}
}

void del()
{
	Delete2D(Ax); Delete2D(Ay); Delete2D(af); Delete2D(umax0);
	for(i = 1; i <= 5; i++)
	{
		Delete2D(U[i]); Delete2D(Ux[i]); Delete2D(Uy[i]); Delete2D(ULx[i]); Delete2D(URx[i]); Delete2D(ULy[i]); Delete2D(URy[i]);
		Delete2D(Fx[i]); Delete2D(Fy[i]); Delete2D(Fy1[i]); Delete2D(FLx[i]); Delete2D(FRx[i]); Delete2D(FLy[i]); Delete2D(FRy[i]);
		Delete2D(deltax[i]); Delete2D(deltaAx[i]); Delete2D(deltay[i]); Delete2D(deltaAy[i]); Delete2D(var[i]); Delete2D(var0[i]);
		Delete2D(varLx[i]); Delete2D(varRx[i]); Delete2D(varLy[i]); Delete2D(varRy[i]);	
	}
}

// Calculate |A| matrix at midpoints 
void calA(double u, double v, double w, double h, double **A)
{
	double V2, a, gama1, lamda[5][5], A1[5][5];
	int i,j,k;

	V2 = u * u + v * v + w * w;
	gama1 = gama - 1.0;
	a = sqrt(gama1 * (h - 0.5 * V2));

	memset(lamda, 0.0, (5 * 5) * sizeof(lamda[0][0]));
	//diagonal elements
	lamda[0][0] = fabs(u - a); lamda[1][1] = fabs(u); lamda[2][2] = fabs(u); 
	lamda[3][3] = fabs(u); lamda[4][4] = fabs(u + a);

	//right eigenvector
	double K[5][5] = {
					{1.0, 1.0, 0.0, 0.0, 1.0},
					{u - a, u, 0.0, 0.0, u + a},
					{v, v, 1.0, 0.0, v},
					{w, w, 0.0, 1.0, w},
					{h - u * a, 0.5 * V2, v, w, h + u * a},
				 };
	//left eigenvector
	double Kinv[5][5] = {
					{h + (a / gama1) * (u - a), -(u + (a / gama1)), -v, -w, 1.0},
					{-2.0 * h + (4.0 / gama1) * a * a, 2.0 * u, 2.0 * v, 2.0 * w, -2.0},
					{-2.0 * v * a * a / gama1, 0.0, 2.0 * a * a / gama1, 0.0, 0.0},
					{-2.0 * w * a * a / gama1, 0.0, 0.0, 2.0 * a * a / gama1, 0.0},
					{h - (a / gama1) * (u + a), -u + (a / gama1), -v, -w, 1.0},
				 };

	//lamda * const * K^(-1)
	for(i = 0; i <= 4 ;i++)
	{
		for(j = 0; j <= 4; j++)
		{
			A1[i][j] = 0.0;
			for(k = 0; k <= 4; k++)
			{
				A1[i][j] += lamda[i][k] * (gama1 / (2.0 * a * a)) * Kinv[k][j];
			}
		}
	}

	// A(at mid point) = K * A1
	for(i = 0; i <= 4 ;i++)
	{
		for(j = 0; j <= 4; j++)
		{
			A[i + 1][j + 1] = 0.0;
			for(k = 0; k <= 4; k++)
			{
				A[i + 1][j + 1] += K[i][k] * A1[k][j]; // make A index 1-5
			}
		}
	}
}

double calA_var(double varl, double varr, double rhol, double rhor)
{
	return (sqrt(rhol) * varl + sqrt(rhor) * varr) / (sqrt(rhol) + sqrt(rhor));
}

void calc_U_F(double ***U,  double ***var, double ***F, int n)
{
	for(i = 1; i <= Nx; i++) {
		for(j = 1; j <= Ny; j++) {
			// calculate conservative variables & F
			U[1][i][j] = var[1][i][j]; 
			if(n == 1)
			{		
				// in y direction, (replace u -> v, v -> w, w -> u)
				U[4][i][j] = var[1][i][j] * var[2][i][j];
				U[2][i][j] = var[1][i][j] * var[3][i][j];
				U[3][i][j] = var[1][i][j] * var[4][i][j];
			}
			else
			{
				U[2][i][j] = var[1][i][j] * var[2][i][j];
				U[3][i][j] = var[1][i][j] * var[3][i][j];
				U[4][i][j] = var[1][i][j] * var[4][i][j];
			}
 			U[5][i][j] = var[5][i][j] / (gama - 1.0) + 0.5 * var[1][i][j] *
					 (var[2][i][j] * var[2][i][j] + var[3][i][j] * var[3][i][j] + var[4][i][j] * var[4][i][j]);

			F[1][i][j] = U[2][i][j];
			F[2][i][j] = (U[2][i][j] * U[2][i][j]) / U[1][i][j] + var[5][i][j];
			F[3][i][j] = (U[2][i][j] * U[3][i][j]) / U[1][i][j];
			F[4][i][j] = (U[2][i][j] * U[4][i][j]) / U[1][i][j];
			F[5][i][j] = (var[5][i][j]  + U[5][i][j] ) * (U[2][i][j] / U[1][i][j]);	
		}
	}
}

int main(void)
{
	init();

// initial conditions 
   for(i = 1; i <= Nx; i++) {
      x[i] =  double(i - 1) * dx;

		for(j = 1; j <= Ny; j++) {
			y[j] =  double(j - 1) * dy;

		   if(x[i] < 0.6) 
			{	
				e1 = p1 / (gama - 1.0) + 0.5 * rho1 * (u1*u1 + v1*v1 + w1*w1);

				U[1][i][j] = rho1;
				U[2][i][j] = rho1 * u1;
				U[3][i][j] = rho1 * v1;
				U[4][i][j] = rho1 * w1;
				U[5][i][j] = e1;

				var[1][i][j] = var0[1][i][j] = rho1;
				var[2][i][j] = var0[2][i][j] = u1;
				var[3][i][j] = var0[3][i][j] = v1;
				var[4][i][j] = var0[4][i][j] = w1;
				var[5][i][j] = var0[5][i][j] = p1;
	 
				// find af 
				af1 = sqrt(gama * p1 / rho1);
				umax1 = fabs(u1) + af1;
			}

		   else 
			{
				// plasma cloud
				if( ((x[i] - 0.8) * (x[i] - 0.8) + (y[j] - 0.5) * (y[j] - 0.5)) <= 0.15 * 0.15 ) var[1][i][j] = 5.0;
				else var[1][i][j] = rho2;
				var0[1][i][j] = var[1][i][j];

				e2 = p2 / (gama - 1.0) + 0.5 * var[1][i][j] * (u2*u2 + v2*v2 + w2*w2);

				var[2][i][j] = var0[2][i][j] = u2;
				var[3][i][j] = var0[3][i][j] = v2;
				var[4][i][j] = var0[4][i][j] = w2;
				var[5][i][j] = var0[5][i][j] = p2;


				U[1][i][j] = var[1][i][j];
				U[2][i][j] = var[1][i][j] * u2;
				U[3][i][j] = var[1][i][j] * v2;
				U[4][i][j] = var[1][i][j] * w2;
				U[5][i][j] = e2;

				// find af 
				af2 = sqrt(gama * p2 / var[1][i][j]);
				umax2 = fabs(u2) + af2;
			}

		}
			// max (\u\+af, \v\+af) = \u\+af (since \v\ = 0)
			if(umax1 > umax2) umax = umax1;
			else umax = umax2;
			//initial timestep
			dt = c * dx / umax; 
	}

	t = 0.0;
	z = 0;

   while(t < tmax) {
		umax = 0.0;
		for(i = 1; i <= Nx - 1; i++) {
			for(j = 1; j <= Ny; j++) {
				for(k = 1; k <= 5; k++) {
					deltax[k][i][j] = var[k][i + 1][j] - var[k][i][j];
				}
			}
		}

		// boundaries 
		for(j = 1; j <= Ny; j++) {
				for(k = 1; k <= 5; k++) {
					deltax[k][0][j] = deltax[k][1][j];
					deltax[k][Nx][j] = deltax[k][Nx - 1][j];
				}
		}
		
		for(i = 1; i <= Nx; i++) {
			for(j = 1; j <= Ny - 1; j++) {
				for(k = 1; k <= 5; k++) {
					deltay[k][i][j] = var[k][i][j + 1] - var[k][i][j];
				}
			}
		}

		// boundaries 
		for(i = 1; i <= Nx; i++) {
				for(k = 1; k <= 5; k++) {
					deltay[k][i][0] = deltay[k][i][1];
					deltay[k][i][Ny] = deltay[k][i][Ny - 1];
				}
		}

		for(i = 1; i <= Nx; i++) {
			for(j = 1; j <= Ny; j++) {
				for(k = 1; k <= 5; k++) {
					deltaAx[k][i][j] = ( deltax[k][i - 1][j] * (deltax[k][i][j] * deltax[k][i][j] + epsilon) 
											+ deltax[k][i][j] * (deltax[k][i - 1][j] * deltax[k][i - 1][j] + epsilon) ) 
											/ (deltax[k][i][j] * deltax[k][i][j] + deltax[k][i - 1][j] * deltax[k][i - 1][j] + 2.0 * epsilon);

					deltaAy[k][i][j] = ( deltay[k][i][j - 1] * (deltay[k][i][j] * deltay[k][i][j] + epsilon) 
											+ deltay[k][i][j] * (deltay[k][i][j - 1] * deltay[k][i][j - 1] + epsilon) ) 
											/ (deltay[k][i][j] * deltay[k][i][j] + deltay[k][i][j - 1] * deltay[k][i][j - 1] + 2.0 * epsilon);
				}
			}
		}

		for(i = 1; i <= Nx; i++) {
			for(j = 1; j <= Ny; j++) {
				for(k = 1; k <= 5; k++) {
					varLx[k][i][j] = var[k][i][j] + 0.5 * deltaAx[k][i][j];
					varLy[k][i][j] = var[k][i][j] + 0.5 * deltaAy[k][i][j];
				}
			}
		}

		for(i = 1; i <= Nx - 1; i++) {
			for(j = 1; j <= Ny; j++) {
				for(k = 1; k <= 5; k++) {
					varRx[k][i][j] = var[k][i + 1][j] - 0.5 * deltaAx[k][i + 1][j];
				}
			}
		}

		// boundaries 
		for(j = 1; j <= Ny; j++) {
				for(k = 1; k <= 5; k++) {
					varRx[k][Nx][j] = varRx[k][Nx - 1][j];
				}
		}

		for(i = 1; i <= Nx; i++) {
			for(j = 1; j <= Ny - 1; j++) {
				for(k = 1; k <= 5; k++) {
					varRy[k][i][j] = var[k][i][j + 1] - 0.5 * deltaAy[k][i][j + 1];
				}
			}
		}

		// boundaries 
		for(i = 1; i <= Nx; i++) {
				for(k = 1; k <= 5; k++) {
					varRy[k][i][Ny] = varRy[k][i][Ny - 1];
				}
		}

		// calculate L,R conservative variables and fluxes
		calc_U_F(ULx, varLx, FLx, 0); calc_U_F(URx, varRx, FRx, 0);
		calc_U_F(ULy, varLy, FLy, 1); calc_U_F(URy, varRy, FRy, 1);

		// calculate matrix in i + 1/2,j
		for(i = 1; i <= Nx - 1; i++) {
			for(j = 1; j <= Ny - 1; j++) {
				uh = calA_var(varLx[2][i][j], varRx[2][i][j], varLx[1][i][j], varRx[1][i][j]);
				vh = calA_var(varLx[3][i][j], varRx[3][i][j], varLx[1][i][j], varRx[1][i][j]);
				wh = calA_var(varLx[4][i][j], varRx[4][i][j], varLx[1][i][j], varRx[1][i][j]);

				// calculate L,R  h 
				hL = (ULx[5][i][j] + varLx[5][i][j]) / varLx[1][i][j];
				hR = (URx[5][i][j] + varRx[5][i][j]) / varRx[1][i][j];

				hh = calA_var(hL, hR, varLx[1][i][j], varRx[1][i][j]);

				calA(uh, vh, wh, hh, Ax); // calculate |A| matrix

				for(k = 1; k <= 5; k++) Ux[k][i][j] = 0.0;

				// calculate \A\ * (UR - UL)
				for(k = 1; k <= 5; k++) {
					for(l = 1; l <= 5; l++) {
						Ux[k][i][j] += Ax[k][l] *  (URx[l][i][j] - ULx[l][i][j]);
					}
					Fx[k][i][j] = 0.5 * (FLx[k][i][j] + FRx[k][i][j] - Ux[k][i][j]);
				}
			}
		}

		// calculate matrix in i,j + 1/2
		for(i = 1; i <= Nx - 1; i++) {
			for(j = 1; j <= Ny - 1; j++) {
				uh = calA_var(varLy[2][i][j], varRy[2][i][j], varLy[1][i][j], varRy[1][i][j]);
				vh = calA_var(varLy[3][i][j], varRy[3][i][j], varLy[1][i][j], varRy[1][i][j]);
				wh = calA_var(varLy[4][i][j], varRy[4][i][j], varLy[1][i][j], varRy[1][i][j]);

				// calculate L,R h 
				hL = (ULy[5][i][j] + varLy[5][i][j]) / varLy[1][i][j];
				hR = (URy[5][i][j] + varRy[5][i][j]) / varRy[1][i][j];

				hh = calA_var(hL, hR, varLy[1][i][j], varRy[1][i][j]);

				calA(vh, wh, uh, hh, Ay); // calculate |A| matrix

				for(k = 1; k <= 5; k++) Uy[k][i][j] = 0.0;

				// calculate \A\ * (UR - UL)
				for(k = 1; k <= 5; k++) {
					for(l = 1; l <= 5; l++) {
						Uy[k][i][j] += Ay[k][l] *  (URy[l][i][j] - ULy[l][i][j]);
					}
					Fy[k][i][j] = 0.5 * (FLy[k][i][j] + FRy[k][i][j] - Uy[k][i][j]);
				}				
			}
		}

		for(i = 1; i <= Nx - 1; i++) {
			for(j = 1; j <= Ny - 1; j++) {
				for(k = 1; k <= 5; k++) {
					// change flux back in y direction
					if(k == 2) Fy1[k][i][j] = Fy[4][i][j];
					else if(k == 3) Fy1[k][i][j] = Fy[2][i][j];
					else if(k == 4) Fy1[k][i][j] = Fy[3][i][j];
					else Fy1[k][i][j] = Fy[k][i][j];

					U[k][i][j] = U[k][i][j] - (dt / dx) * (Fx[k][i][j] - Fx[k][i - 1][j]) 
													- (dt / dy) * (Fy1[k][i][j] - Fy1[k][i][j - 1]);
				}

				// calculate all the variables
				var[1][i][j] = U[1][i][j];
				var[2][i][j] = U[2][i][j] / U[1][i][j];
				var[3][i][j] = U[3][i][j] / U[1][i][j];
				var[4][i][j] = U[4][i][j] / U[1][i][j];
				var[5][i][j] = (gama - 1.0) * (U[5][i][j] - 0.5 * var[1][i][j] * (var[2][i][j] * var[2][i][j] + var[3][i][j] * var[3][i][j] + var[4][i][j] * var[4][i][j])); 	
			}
		}

		// y boundaries
		for(i = 1; i <= Nx; i++) {
				for(k = 1; k <= 5; k++) {
					var[k][i][1] = var[k][i][2];
					var[k][i][Ny] = var[k][i][Ny - 1];
				}
		}

		// x boundaries
		for(j = 1; j <= Ny; j++) {
				for(k = 1; k <= 5; k++) {
					var[k][1][j] = var0[k][2][j];
					var[k][Nx][j] = var0[k][Nx - 1][j];
				}
		}

		for(i = 1; i <= Nx ; i++) {
			for(j = 1; j <= Ny; j++) {
			   af[i][j] = sqrt(gama * var[5][i][j] / var[1][i][j]);

				umax0[i][j] = max1(fabs(var[2][i][j]) + af[i][j] , fabs(var[3][i][j]) + af[i][j]);
				if(umax0[i][j] > umax) umax = umax0[i][j];		
			}
		}
cout<<dt<<endl;
		//new timestep
		dt = c * dx / umax; // dt = dx / fabs(umax)
		t = t + dt;

		b = 0;
/*
		// calculate variables at each step for movie
		for(i = 1; i <= Nx; i++) {
			for(j = 1; j <= Ny; j++) {
				rho[b] = var[1][j][i]; u[b] = var[2][j][i]; v[b] = var[3][j][i]; w[b] = var[4][j][i]; p[b] = var[5][j][i]; b++;
			}
		}
*/
	//	WriteSilo(rho, u ,v, w, p, z);
		z = z++;
	}
/*
	// plot variables at tmax
	for(i = 1; i <= Nx; i++) {
		for(j = 1; j <= Ny; j++) {
			rho[b] = var[1][j][i]; u[b] = var[2][j][i]; v[b] = var[3][j][i]; w[b] = var[4][j][i]; p[b] = var[5][j][i]; b++;
		}
	}
	WriteSilo(rho, u ,v, w, p, z);
*/

/*
// 1D plots at x[i] = 0.55 (i ~ 141 ==> x[i] = 0.549) 
   for(j = 1; j <= Ny; j++) {
		 cout << setw(12) << setprecision(6) << y[j]<< setw(20) << setprecision(10) << var[1][141][j]
				<< setw(20) << setprecision(10) << var[2][141][j]<< setw(20) << setprecision(10) << var[3][141][j]
				<< setw(20) << setprecision(10) << var[4][141][j]<< setw(20) << setprecision(10) << var[5][141][j]<<endl;
	}
*/
	del();
}

// generate silo file for movie
void WriteSilo(double *data1, double *data2, double *data3, double *data4, double *data5, int cycle)
{
	DBfile *dbfile = NULL;
	int err,i,j,k;
	double x1[Nx + 1], y1[Ny + 1], var1[Nx * Ny], var2[Nx * Ny], var3[Nx * Ny], var4[Nx * Ny], var5[Nx * Ny];
	char silo[100];
	sprintf(silo, "output%04d.silo", cycle);
	dbfile = DBCreate(silo, DB_CLOBBER, DB_LOCAL, NULL, DB_PDB);

	for(i = 1; i <= Nx; i++) x1[i] = x[i];
	for(j = 1; j <= Ny; j++) y1[j] = y[j];

	int dims[] = {Nx, Ny};
	int ndims = 2;
	double *coords[] = {x1, y1};
	err = DBPutQuadmesh(dbfile, "quadmesh", NULL, coords, dims, ndims,DB_DOUBLE, DB_COLLINEAR, NULL);

	for(k = 0; k <= Nx * Ny; k++) 
	{
		var1[k] = data1[k]; var2[k] = data2[k], var3[k] = data3[k]; var4[k] = data4[k]; var5[k] = data5[k];
	}

	DBPutQuadvar1(dbfile, "rho", "quadmesh",var1, dims,ndims, NULL, 0, DB_DOUBLE, DB_NODECENT, NULL);
	DBPutQuadvar1(dbfile, "u", "quadmesh",var2, dims,ndims, NULL, 0, DB_DOUBLE, DB_NODECENT, NULL);
	DBPutQuadvar1(dbfile, "v", "quadmesh",var3, dims,ndims, NULL, 0, DB_DOUBLE, DB_NODECENT, NULL);
	DBPutQuadvar1(dbfile, "w", "quadmesh",var4, dims,ndims, NULL, 0, DB_DOUBLE, DB_NODECENT, NULL);
	DBPutQuadvar1(dbfile, "p", "quadmesh",var5, dims,ndims, NULL, 0, DB_DOUBLE, DB_NODECENT, NULL);
   DBClose(dbfile);
}
