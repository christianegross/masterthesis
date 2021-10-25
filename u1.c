/**
 * Christiane, Oktober 2021
 * Path Integral MC of U(1) lattice gauge theory in 3+1 dimensions
 * like in paper from Loan et al.
 * **/
 
//convention lattice numbering: temporal, spatial, spatial, spatial, direction
// t, x, y, z, mu 


#include <stdio.h>
#include "math.h"
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_block.h>

/**
 * calculates plaquette as given in paper, uses a* b* =(ab)* for complex numbers
 * **/
double plaquette(double * lattice, int * neighbour, int position, int mu, int nu){
    //~ gsl_complex p=gsl_complex_mul(lattice[position+mu], lattice[position+nu+neighbour[mu]]); //U_mu(r)*U_nu(r+mu)
    //~ gsl_complex q=gsl_complex_mul(lattice[position+neighbour[nu]+mu], lattice[position+nu]); //U_mu(r+vu)*U_nu(r)
    //~ return 1.0-GSL_REAL(gsl_complex_mul(p, gsl_complex_conjugate(q)));
    return 1.0-cos(lattice[position+mu]+lattice[position+nu+neighbour[mu]]-lattice[position+neighbour[nu]+mu]-lattice[position+nu]);
}

double action(double *lattice, double beta, double deltatau, int *neighbour, int Nt, int Ns){
    int position;
    int spaceaction=0;
    int timeaction=0;
    for (int t=0;t<Nt;t+=1){
         neighbour[0]=(t==Nt-1)?-(Nt-1)*Ns*Ns*Ns*4:Ns*Ns*Ns*4;
         for (int x=0;x<Ns;x+=1){
            neighbour[1]=(x==Ns-1)?-(Ns-1)*Ns*Ns*4:Ns*Ns*4;
            for (int y=0;y<Ns;y+=1){
                neighbour[2]=(y==Ns-1)?-(Ns-1)*Ns*4:Ns*4;
                for (int z=0;z<Ns;z+=1){
                    neighbour[3]=(z==Ns-1)?-(Ns-1)*4:4;
                    position=Ns*Ns*Ns*4*t+Ns*Ns*4*x+Ns*4*y+4*z;
                    for (int mu=1;mu<4;mu+=1){
                        for (int nu=1;nu<mu;nu+=1){
                            spaceaction+=plaquette(lattice, neighbour, position, mu, nu);
                        }
                        timeaction+=plaquette(lattice, neighbour, position, mu, 0);
                    }
                }
            }
        }
    }
    return beta*(deltatau*spaceaction+timeaction/deltatau);
}


int main(int argc, char **argv){
    double beta=2;
    double deltatau=0.8;   
    int Ns=2;       //number of sites for the spatial directions
    int Nt=2;       //number of sites for the temporal directions
    int latticesites=Nt*Ns*Ns*Ns*4;
    double lattice[latticesites];
    //~ gsl_complex lattice[latticesites];
    int neighbour[4];
    double averageplaquette=0;
    
    /** initialize lattice **/
    for (int i=0; i<latticesites; i+=1){
        //~ lattice[i]=gsl_complex_polar(1,i);
        lattice[i]=i;
    }
    
    /**neighbours: need two neighbours for every direction, one forward and one backward
	 * 0: t-forward
	 * 1: x-forward
	 * 2: y-forward
	 * 3: z-forward
	 * relative to position
	 * mu			0	1	2	3
	 * direction	t	x	y	z
     * implement periodic boundary conditions
	 * **/
    /**
     * go through the entire lattice, calculate average plaquette
     * **/ 
    int position;
    for (int t=0;t<Nt;t+=1){
         neighbour[0]=(t==Nt-1)?-(Nt-1)*Ns*Ns*Ns*4:Ns*Ns*Ns*4;
         for (int x=0;x<Ns;x+=1){
            neighbour[1]=(x==Ns-1)?-(Ns-1)*Ns*Ns*4:Ns*Ns*4;
            for (int y=0;y<Ns;y+=1){
                neighbour[2]=(y==Ns-1)?-(Ns-1)*Ns*4:Ns*4;
                for (int z=0;z<Ns;z+=1){
                    neighbour[3]=(z==Ns-1)?-(Ns-1)*4:4;
                    position=Ns*Ns*Ns*4*t+Ns*Ns*4*x+Ns*4*y+4*z;
                    for (int mu=0;mu<4;mu+=1){
                        for (int nu=0;nu<4;nu+=1){
                            averageplaquette+=plaquette(lattice, neighbour, position, mu, nu);
                        }
                    }
                }
            }
        }
    }
    averageplaquette/=latticesites*4;
    printf("%f\n", action(lattice, beta, deltatau, neighbour, Nt, Ns)); 
	return 0;
}

