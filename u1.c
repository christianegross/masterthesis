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
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

/**
 * calculates plaquette as given in paper, uses a* b* =(ab)* for complex numbers
 * **/
double plaquette(gsl_complex * lattice, int * neighbour, int position, int mu, int nu){
    gsl_complex p=gsl_complex_mul(lattice[position+mu], lattice[position+nu+neighbour[mu]]); //U_mu(r)*U_nu(r+mu)
    gsl_complex q=gsl_complex_mul(lattice[position+neighbour[nu]+mu], lattice[position+nu]); //U_mu(r+vu)*U_nu(r)
    return 1.0-GSL_REAL(gsl_complex_mul(p, gsl_complex_conjugate(q)));
    //~ return 1.0-cos(lattice[position+mu]+lattice[position+nu+neighbour[mu]]-lattice[position+neighbour[nu]+mu]-lattice[position+nu]);
}

/**
 * goes through the entire lattice to calculate the action, keeps time and space sums separate to account for anisotropy 
 * neighbour implements periodic boundary conditions
 * position keeps the code more legible
 * **/
double action(gsl_complex *lattice, double beta, double deltatau, int *neighbour, int Nt, int Ns){
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

/** 
 * calculates the change in action if the link at position, mu is changed
 * one link is involved in six plaquettes, if the link is changed, the rest does not change. 
 * uses that only the real part is needed later, so uses Re(ab*)=Re(a*b)
 * takes into account anisotropy by giving a different weight to temporal contribution
 * **/
double deltas(gsl_complex * lattice, int * neighbour, int position, double beta, double deltatau, int mu, gsl_complex change){
    gsl_complex sum=gsl_complex_rect(0,0);
    gsl_complex p=gsl_complex_rect(0,0), q=gsl_complex_rect(0,0);
    for (int nu=0;nu<4;nu+=1){if(nu!=mu){
        //contribution at position
        p=gsl_complex_mul(lattice[position+neighbour[nu]+mu], lattice[position+nu]); //U_mu(r+vu)*U_nu(r)
        p=gsl_complex_mul(lattice[position+nu+neighbour[mu]], gsl_complex_conjugate(p)); //U_nu(r+mu)
        //contribution at position-nu, use conjugate of plaquette to be able to use unconjugated missing link
        q=gsl_complex_mul(lattice[position+neighbour[4+nu]+mu], lattice[position+neighbour[4+nu]+neighbour[mu]+nu]); //U_mu(r-vu)*U_nu(r-nu+mu)
        q=gsl_complex_mul(lattice[position+nu+neighbour[4+nu]], gsl_complex_conjugate(q)); //U_nu(r-nu)
        p=gsl_complex_add(p,q);
        if (mu==0||nu==0){
            p=gsl_complex_mul_real(p, 1.0/deltatau/deltatau);
        }
        p=gsl_complex_mul_real(p, beta*deltatau);
        sum=gsl_complex_add(sum, p);
    }}
    //~ double constantsummand=4*deltatau+2/deltatau;
    //~ if(mu==0){constantsummand=6/deltatau;}
    //~ beta*constantsummand
    return -1.0*GSL_REAL(gsl_complex_mul(sum, change));
}


int main(int argc, char **argv){
    double beta=2;
    double deltatau=0.8;   
    int Ns=4;       //number of sites for the spatial directions
    int Nt=4;       //number of sites for the temporal directions
    int latticesites=Nt*Ns*Ns*Ns*4;
    //~ double lattice[latticesites];
    gsl_complex lattice[latticesites];
    int neighbour[8];
    
    gsl_rng *generator;
	generator=gsl_rng_alloc(gsl_rng_mt19937);	
    gsl_rng_set(generator, 0);
    double delta=M_PI*0.68;
    
    
    double averageplaquette=0;
    int accept=0;
    double calc_action=0;
    gsl_complex change;
    
    /** initialize lattice **/
    for (int i=0; i<latticesites; i+=1){
        lattice[i]=gsl_complex_polar(1,i);
        //~ lattice[i]=i;
    }
    
    /**neighbours: need two neighbours for every direction, one forward and one backward
	 * 0: t-forward     4: t-backward
	 * 1: x-forward     5: x-backward
	 * 2: y-forward     6: y-backward
	 * 3: z-forward     7: z-backward
     * forward mu, backward 4+mu
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
         neighbour[0]=(t==Nt-1) ?-(Nt-1)*Ns*Ns*Ns*4:Ns*Ns*Ns*4;
         neighbour[4]=(t==0)    ?(Nt-1)*Ns*Ns*Ns*4:-Ns*Ns*Ns*4;
         for (int x=0;x<Ns;x+=1){
            neighbour[1]=(x==Ns-1)  ?-(Ns-1)*Ns*Ns*4:Ns*Ns*4;
            neighbour[5]=(x==0)     ?(Ns-1)*Ns*Ns*4:-Ns*Ns*4;
            for (int y=0;y<Ns;y+=1){
                neighbour[2]=(y==Ns-1)  ?-(Ns-1)*Ns*4:Ns*4;
                neighbour[6]=(y==0)     ?(Ns-1)*Ns*4:-Ns*4;
                for (int z=0;z<Ns;z+=1){
                    neighbour[3]=(z==Ns-1)  ?-(Ns-1)*4:4;
                    neighbour[7]=(z==0)     ?(Ns-1)*4:-4;
                    position=Ns*Ns*Ns*4*t+Ns*Ns*4*x+Ns*4*y+4*z;
                    for (int mu=0;mu<4;mu+=1){
                        change=gsl_complex_polar(1, gsl_ran_flat(generator, -delta, delta));
                        //~ printf("%.1f %.1f\n", deltas(lattice, neighbour, position, beta, deltatau, mu, change), GSL_REAL(change));
                        //~ printf("%.1f %.1f\n", deltas(lattice, neighbour, position, beta, deltatau, mu, change), GSL_REAL(change));
                        if(exp(-deltas(lattice, neighbour, position, beta, deltatau, mu, change))>gsl_rng_uniform(generator)){
                        lattice[position+mu]=gsl_complex_add(lattice[position+mu], change);
                        accept+=1;
                        }
                        //~ for (int nu=0;nu<4;nu+=1){
                            //~ averageplaquette+=plaquette(lattice, neighbour, position, mu, nu);
                        //~ }
                    }
                }
            }
        }
    }
    averageplaquette/=latticesites*4;
    printf("%f\n", accept/(1.0*latticesites)); 
    gsl_rng_free(generator);
	return 0;
}

