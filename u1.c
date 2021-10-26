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
#include <gsl/gsl_vector.h>

/** for a given position in the lattice, determines the single coordinates and determines neighbours with periodic boundary conditions
 * position can also include direction, even if not taken by other functions
 * **/
void get_neighbours(int *neighbour, int position, int Nt, int Ns){
    int pos=position;
    int direction=pos%4;
    int t, x, y, z;
    pos-=direction;
    z=(pos%(4*Ns))/4;
    pos-=4*z;
    y=(pos%(4*Ns*Ns))/(4*Ns);
    pos-=4*Ns*y;
    x=(pos%(4*Ns*Ns*Ns))/(4*Ns*Ns);
    pos-=4*Ns*Ns*x;
    t=(pos%(4*Ns*Ns*Ns*Nt))/(4*Ns*Ns*Ns);
    pos-=4*Ns*Ns*Ns*t;
    if (pos!=0){printf("check needed!\n");}
    neighbour[0]=(t==Nt-1)  ?-(Nt-1)*Ns*Ns*Ns*4 :Ns*Ns*Ns*4 ;
    neighbour[4]=(t==0)     ?(Nt-1)*Ns*Ns*Ns*4  :-Ns*Ns*Ns*4;
    neighbour[1]=(x==Ns-1)  ?-(Ns-1)*Ns*Ns*4    :Ns*Ns*4    ;
    neighbour[5]=(x==0)     ?(Ns-1)*Ns*Ns*4     :-Ns*Ns*4   ;
    neighbour[2]=(y==Ns-1)  ?-(Ns-1)*Ns*4       :Ns*4       ;
    neighbour[6]=(y==0)     ?(Ns-1)*Ns*4        :-Ns*4      ;
    neighbour[3]=(z==Ns-1)  ?-(Ns-1)*4          :4          ;
    neighbour[7]=(z==0)     ?(Ns-1)*4           :-4         ;
}

//simple mean of a vector
double mean(gsl_vector *vector, int length){
    double mean=0;
    for (int i=0; i<length; i+=1){
        mean+=gsl_vector_get(vector, i);
    }
    return mean/length;
}

//simple standard deviation of a vector with given mean
double deviation(gsl_vector *vector, int length, double mean){
    double var=0;
    for (int i=0; i<length; i+=1){
        var+=pow(gsl_vector_get(vector, i)-mean, 2);
    }
    return sqrt(var/length);
}

/**
 * calculates plaquette as given in paper, uses a* b* =(ab)* for complex numbers
 * **/
double plaquette(double * lattice, int * neighbour, int position, int mu, int nu){
    //~ gsl_complex p=gsl_complex_mul(lattice[position+mu], lattice[position+nu+neighbour[mu]]); //U_mu(r)*U_nu(r+mu)
    //~ gsl_complex q=gsl_complex_mul(lattice[position+neighbour[nu]+mu], lattice[position+nu]); //U_mu(r+vu)*U_nu(r)
    //~ return 1.0-GSL_REAL(gsl_complex_mul(p, gsl_complex_conjugate(q)));
    return 1.0-cos(lattice[position+mu]+lattice[position+nu+neighbour[mu]]-lattice[position+neighbour[nu]+mu]-lattice[position+nu]);
}

/**
 * goes through the entire lattice to calculate the action, keeps time and space sums separate to account for anisotropy 
 * neighbour implements periodic boundary conditions
 * position keeps the code more legible
 * **/
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

//~ /** 
 //~ * uses complex numbers in the lattice
 //~ * calculates the change in action if the link at position, mu is changed
 //~ * one link is involved in six plaquettes, if the link is changed, the rest does not change. 
 //~ * uses that only the real part is needed later, so uses Re(ab*)=Re(a*b)
 //~ * takes into account anisotropy by giving a different weight to temporal contribution
 //~ * **/
//~ double deltas(gsl_complex * lattice, int * neighbour, int position, double beta, double deltatau, int mu, gsl_complex change){
    //~ gsl_complex sum=gsl_complex_rect(0,0);
    //~ gsl_complex p=gsl_complex_rect(0,0), q=gsl_complex_rect(0,0);
    //~ for (int nu=0;nu<4;nu+=1){if(nu!=mu){
        //~ //contribution at position
        //~ p=gsl_complex_mul(lattice[position+neighbour[nu]+mu], lattice[position+nu]); //U_mu(r+vu)*U_nu(r)
        //~ p=gsl_complex_mul(lattice[position+nu+neighbour[mu]], gsl_complex_conjugate(p)); //U_nu(r+mu)
        //~ //contribution at position-nu, use conjugate of plaquette to be able to use unconjugated missing link
        //~ q=gsl_complex_mul(lattice[position+neighbour[4+nu]+mu], lattice[position+neighbour[4+nu]+neighbour[mu]+nu]); //U_mu(r-vu)*U_nu(r-nu+mu)
        //~ q=gsl_complex_mul(lattice[position+nu+neighbour[4+nu]], gsl_complex_conjugate(q)); //U_nu(r-nu)
        //~ p=gsl_complex_add(p,q);
        //~ if (mu==0||nu==0){
            //~ p=gsl_complex_mul_real(p, 1.0/deltatau/deltatau);
        //~ }
        //~ p=gsl_complex_mul_real(p, beta*deltatau);
        //~ sum=gsl_complex_add(sum, p);
    //~ }}
    //~ //double constantsummand=4*deltatau+2/deltatau;
    //~ //if(mu==0){constantsummand=6/deltatau;}
    //~ //beta*constantsummand
    //~ return -1.0*GSL_REAL(gsl_complex_mul(sum, gsl_complex_sub(lattice[position+mu], gsl_complex_mul(lattice[position+mu], change))));
//~ }

/** 
 * calculates change in the action if link variable at position, mu is changed by change
 * uses only the angle in the lattice
 * one link contributes to six plaquettes, these plaquettes are all calculated, and the difference between the current and the proposed new state is calculated explicitly
 * **/
double deltas(double *lattice, int *neighbour, int position, double beta, double deltatau, int mu, double change){
    double sum=0, oneloop;
    double value=lattice[position+mu];
    double propose=lattice[position+mu]+change;
    //calculate S(value)-S(propose)
    for (int nu=0;nu<4;nu+=1){if(nu!=mu){
        //at position
        oneloop=cos(value+lattice[position+nu+neighbour[mu]]-lattice[position+neighbour[nu]+mu]-lattice[position+nu]);
        oneloop-=cos(propose+lattice[position+nu+neighbour[mu]]-lattice[position+neighbour[nu]+mu]-lattice[position+nu]);
        //at position-nu
        oneloop+=cos(lattice[position+neighbour[4+nu]+mu]+lattice[position+neighbour[4+nu]+neighbour[mu]+nu]-value-lattice[position+nu+neighbour[4+nu]]);
        oneloop-=cos(lattice[position+neighbour[4+nu]+mu]+lattice[position+neighbour[4+nu]+neighbour[mu]+nu]-propose-lattice[position+nu+neighbour[4+nu]]);
        if (mu==0||nu==0){  //account for anisotropy
            oneloop/=deltatau*deltatau;
        }
        oneloop*=beta*deltatau;
        sum+=oneloop;
    }}
    return -sum; //- sign also in definition of plaquette
}

/** goes through the entire lattice and attempts to change links
 * value of the average plaquette is returned, average plaquette and acceptane rate are printed out for plotting
 * can switch between ordered sweep through the lattice and random choosing of latticesites points per sweep
 * **/
double sweep(double *lattice, double beta, double deltatau, gsl_rng *generator, FILE * stream, int Nt, int Ns, double delta){
    double change;  
    int neighbour[8], accept=0, averageplaquette=0, latticesites=Ns*Ns*Ns*Nt*4;  
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
    int position, mu;
    //~ for (int t=0;t<Nt;t+=1){
         //~ neighbour[0]=(t==Nt-1) ?-(Nt-1)*Ns*Ns*Ns*4:Ns*Ns*Ns*4;
         //~ neighbour[4]=(t==0)    ?(Nt-1)*Ns*Ns*Ns*4:-Ns*Ns*Ns*4;
         //~ for (int x=0;x<Ns;x+=1){
            //~ neighbour[1]=(x==Ns-1)  ?-(Ns-1)*Ns*Ns*4:Ns*Ns*4;
            //~ neighbour[5]=(x==0)     ?(Ns-1)*Ns*Ns*4:-Ns*Ns*4;
            //~ for (int y=0;y<Ns;y+=1){
                //~ neighbour[2]=(y==Ns-1)  ?-(Ns-1)*Ns*4:Ns*4;
                //~ neighbour[6]=(y==0)     ?(Ns-1)*Ns*4:-Ns*4;
                //~ for (int z=0;z<Ns;z+=1){
                    //~ neighbour[3]=(z==Ns-1)  ?-(Ns-1)*4:4;
                    //~ neighbour[7]=(z==0)     ?(Ns-1)*4:-4;
                    //~ position=Ns*Ns*Ns*4*t+Ns*Ns*4*x+Ns*4*y+4*z;
                    //~ get_neighbours(neighbour, position, Nt, Ns);
                    //~ for (int mu=0;mu<4;mu+=1){
                        //~ change=gsl_ran_flat(generator, -delta, delta);
                        //~ if(exp(-deltas(lattice, neighbour, position, beta, deltatau, mu, change))>gsl_rng_uniform(generator)){
                        //~ //lattice[position+mu]=gsl_complex_mul(lattice[position+mu], gsl_complex_polar(1,change));
                        //~ lattice[position+mu]+=change;
                        //~ accept+=1;
                        //~ }
                        //~ for (int nu=0;nu<4;nu+=1){
                            //~ averageplaquette+=plaquette(lattice, neighbour, position, mu, nu);
                        //~ }
                    //~ }
                //~ }
            //~ }
        //~ }
    //~ }
    for (int i=0; i<latticesites; i+=1){
        position=gsl_rng_uniform_int(generator, latticesites);
        get_neighbours(neighbour, position, Nt, Ns);
        mu=position%4;
        //have to account: mu is not in position as expected by plaquette, deltas, but is also randomly generated, so faster to put into same variable as position
        position-=mu;
        change=gsl_ran_flat(generator, -delta, delta);
        if(exp(-deltas(lattice, neighbour, position, beta, deltatau, mu, change))>gsl_rng_uniform(generator)){
            //lattice[position+mu]=gsl_complex_mul(lattice[position+mu], gsl_complex_polar(1,change));
            lattice[position+mu]+=change;
            accept+=1;
        }
        for (int nu=0;nu<4;nu+=1){
            averageplaquette+=plaquette(lattice, neighbour, position, mu, nu);
        }
    }
    fprintf(stream, "%f\t%f\n", averageplaquette/(4.0*latticesites), accept/(1.0*latticesites));
    return averageplaquette/(4.0*latticesites);
}



int main(int argc, char **argv){
    double beta=1.4;
    double deltatau=1;   
    int Ns=8;       //number of sites for the spatial directions
    int Nt=8;       //number of sites for the temporal directions
    int latticesites=Nt*Ns*Ns*Ns*4;
    double measurements=15000;
    double lattice[latticesites];
    //~ gsl_complex lattice[latticesites];
    int neighbour[8];
    
    gsl_rng *generator;
	generator=gsl_rng_alloc(gsl_rng_mt19937);	
    gsl_rng_set(generator, 0);
    double delta=M_PI*0.25;     //change delta to change acceptance rate
    
    //save acceptance rate and average plaquette for plotting
    char filename[100];
    sprintf(filename, "measbeta%.1fNt%dNs%drandom.txt", beta, Nt, Ns);
    FILE * stream=fopen(filename, "w+");
    
    double averageplaquette=0;
    int accept=0;
    double change;  //stores proposed changes to link
    gsl_vector * plaq=gsl_vector_alloc(measurements);   //stores measured plaquette values to take average
    
    /** initialize lattice **/
    for (int i=0; i<latticesites; i+=1){
        //~ lattice[i]=gsl_complex_polar(1,i);
        lattice[i]=i;
    }
    
    for(int i=0; i<1000;i+=1){
        sweep(lattice, beta, deltatau, generator, stream, Nt, Ns,  delta);
    }
    printf("thermalization done\n");
    for(int i=0; i<measurements;i+=1){
        gsl_vector_set(plaq, i, sweep(lattice, beta, deltatau, generator, stream, Nt, Ns,  delta));
    }
    
    printf("%f\t%f\t%d\t%d\t%f\t%f\n", beta, deltatau,Ns, Nt, mean(plaq, measurements), deviation(plaq, measurements, mean(plaq, measurements)));
    
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
                        change=gsl_ran_flat(generator, -delta, delta);
                        if(exp(-deltas(lattice, neighbour, position, beta, deltatau, mu, change))>gsl_rng_uniform(generator)){
                        //~ lattice[position+mu]=gsl_complex_mul(lattice[position+mu], gsl_complex_polar(1,change));
                        lattice[position+mu]+=change;
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
    fclose(stream);
    gsl_vector_free(plaq);
	return 0;
}

