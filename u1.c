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
#include <gsl/gsl_sf_bessel.h>

//code for binning, bootstrap, autocorrelation from compphys project https://github.com/christianegross/CompPhys_2021/blob/main/Project/src/main/auxiliary.c
/**
 * @brief calculates the autocorrelation of measurements
 */
void autocorrelation(gsl_vector *measurements, gsl_vector *results, double mean){
	/**
	 * make block, calculate gamma with tau=0, write (0,1) to results
	 * for range of tau: calculate gamma, divide by gamma(0), write (tau, gamma/gamma) to reaults
	 */
	double czero=0;
	double ctime=0;
	for (int i=0; i<int(measurements->size); i+=1){
		 czero+=(gsl_vector_get(measurements,i)-mean)*(gsl_vector_get(measurements,i)-mean);
	}
	czero/=measurements->size;
	gsl_vector_set (results, 0,1);
	for (int time=1; time<int(results->size); time+=1){
		ctime=0;
		for (int i=0; i<int(measurements->size-time); i+=1){
		 ctime+=(gsl_vector_get(measurements,i)-mean)*(gsl_vector_get(measurements,i+time)-mean);
		}
		ctime/=(measurements->size-time);
		gsl_vector_set (results, time, ctime/czero);
	}
}

//calculates integrated autocorrelation time for given list of normalized times and window size
//neglect factor 1/2, negative times because C(tau)=C(-tau)
double tauint(gsl_vector * autocorrtimes, int W){
    double result=0.5; //account for Gamma(0)=1
    for (int i=1;i<=W;i++){
        result+=gsl_vector_get(autocorrtimes, i);
    }
    return result;        
}



/**
 * @fn 		inline void binning(gsl_block *measurements, gsl_block *binneddata, int lengthofbin)
 * @brief 	takes data points in measurements and writes bins over lengthofbin in binneddata
 */
extern inline void binning(gsl_vector *measurements, gsl_vector *binneddata, int lengthofbin){
	double bincontent=0;
	//assertion if lengths of bins, measurements fit together
	int numberofbins=measurements->size/lengthofbin;
	//printf("numberofbins=%d\tnumberofelemtns=%lu\tlengthofbin=%d\n", numberofbins, binneddata->size, lengthofbin);
	if (numberofbins!=int(binneddata->size)){fprintf(stderr, "Gave wrong lengths! %d number of bins to be calculated, but storage allocated for %lu!", numberofbins, binneddata->size); return;}
	//calculate content of single bin as arithmetic mean oover lengthofbin datapoints
	for (int bin=0; bin<numberofbins; bin+=1){
		bincontent=0;
		for (int datapoint=0; datapoint<lengthofbin; datapoint+=1){
			bincontent+=gsl_vector_get(measurements,bin*lengthofbin+datapoint);
		}
		gsl_vector_set(binneddata, bin, bincontent/lengthofbin);
	}
}

/**
 * @fn makebootstrapreplica(gsl_block * measurements, gsl_rng * generator)
 * @brief makes one bootstrapreplica out of the data in measurements
 * 
 * @param measurements contains all the measurements which can be used to calculate the replica
 * @param generator random number generator which determines which elements of measurements are used to calaculate the replica
 * 
 * @return replica arithmetic mean of chosen measurements, one bootstrapreplica
 */
double makebootstrapreplica(gsl_vector * measurements, gsl_rng * generator){
	double replica=0;
	int randomnumber;
	for (int datapoint=0; datapoint<int(measurements->size); datapoint+=1){
		randomnumber=gsl_rng_uniform_int(generator, measurements->size);
		replica+=gsl_vector_get(measurements, randomnumber);
	}
	return replica/measurements->size;
}
/**
 * @fn bootstrap(gsl_block *measurements, gsl_rng *generator, int R, double *mean, double *variance)
 * @brief uses the bootstrapmethod to calculate mean and variance of the data in measurements
 * 
 * @param measurements contains all the datappoints which are used
 * @param generator random number generator which determines which elements of measurements are used to calaculate the replica
 * @param R number of replicas which are used to calculate mean and variance
 * @param calculated mean over all replicas
 * @param variance calculated variance over all replica
 */ 
void bootstrap(gsl_vector *measurements, gsl_rng *generator, int R, double *mean, double *variance){
	*mean=0;
	*variance=0;
	gsl_block *replicalist=gsl_block_alloc(R);
	for (int i=0; i<R; i+=1){
		replicalist->data[i]=makebootstrapreplica(measurements, generator);
		*mean+=replicalist->data[i];
	}
	//calculate arithmetic mean over replicas
	*mean/=R;
	//calculate standard deviation of replicas
	for (int i=0; i<R; i+=1){
		*variance+=(*mean-replicalist->data[i])*(*mean-replicalist->data[i]);
	}
	*variance/=R-1;
	gsl_block_free(replicalist);
	}

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
    if (pos!=0){fprintf(stderr, "check needed!\n");}
    neighbour[0]=(t==Nt-1)  ?-(Nt-1)*Ns*Ns*Ns*4 :Ns*Ns*Ns*4 ;
    neighbour[4]=(t==0)     ?(Nt-1)*Ns*Ns*Ns*4  :-Ns*Ns*Ns*4;
    neighbour[1]=(x==Ns-1)  ?-(Ns-1)*Ns*Ns*4    :Ns*Ns*4    ;
    neighbour[5]=(x==0)     ?(Ns-1)*Ns*Ns*4     :-Ns*Ns*4   ;
    neighbour[2]=(y==Ns-1)  ?-(Ns-1)*Ns*4       :Ns*4       ;
    neighbour[6]=(y==0)     ?(Ns-1)*Ns*4        :-Ns*4      ;
    neighbour[3]=(z==Ns-1)  ?-(Ns-1)*4          :4          ;
    neighbour[7]=(z==0)     ?(Ns-1)*4           :-4         ;
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
        var+=(gsl_vector_get(vector, i)-mean)*(gsl_vector_get(vector, i)-mean);
    }
    return var/(length-1);
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
 //~ * 
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
 * calculates change in the action if link variable at (position, mu) is changed by change
 * uses only the angle in the lattice
 * one link contributes to six plaquettes, these plaquettes are all calculated, and the difference between the current and the proposed new state is calculated explicitly
 * **/
double deltas(double *lattice, int *neighbour, int position, double beta, double deltatau, int mu, double change){
    double sum=0, oneloop, current;
    double value=lattice[position+mu];
    //~ double propose=lattice[position+mu]+change;
    //calculate S(value)-S(propose)
    for (int nu=0;nu<4;nu+=1){if(nu!=mu){
        //at position
        current=value+lattice[position+nu+neighbour[mu]]-lattice[position+neighbour[nu]+mu]-lattice[position+nu];
        oneloop=cos(current);
        oneloop-=cos(current+change);
        //at (position-nu): current link enters backwards, so have to take minus change
        current=lattice[position+neighbour[4+nu]+mu]+lattice[position+neighbour[4+nu]+neighbour[mu]+nu]-value-lattice[position+nu+neighbour[4+nu]];
        oneloop+=cos(current);
        oneloop-=cos(current-change);
        if (mu==0||nu==0){  //account for anisotropy
            oneloop/=deltatau*deltatau;
        }
        sum+=oneloop;
    }}
    return -sum*beta*deltatau; //- sign also in definition of plaquette
}

double averageplaquette(double * lattice, int Nt, int Ns){
    int neighbour[8];
    double averageplaquette=0;
    for (int i=0;i<4*Nt*Ns*Ns*Ns;i+=4){
        get_neighbours(neighbour, i, Nt, Ns);
        for(int mu=0;mu<4;mu+=1){
            for(int nu=mu+1;nu<4;nu+=1){
                averageplaquette+=plaquette(lattice, neighbour, i, mu, nu);
            }
        }
    }
    return averageplaquette/(6*Nt*Ns*Ns*Ns);
}

/** goes through the entire lattice and attempts to change links
 * value of the average plaquette is returned, average plaquette and acceptane rate are printed out for plotting
 * can switch between ordered sweep through the lattice and random choosing of latticesites points per sweep
 * **/
double sweep(double *lattice, double beta, double deltatau, gsl_rng *generator, FILE * stream, int Nt, int Ns, double delta){
    double change;  
    int neighbour[8], accept=0, latticesites=Ns*Ns*Ns*Nt*4;  
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
    int position;//, mu;
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
                    //~ get_neighbours(neighbour, position, Nt, Ns);
                    for (int mu=0;mu<4;mu+=1){
                        change=gsl_ran_flat(generator, -delta, delta);
                        if(exp(deltas(lattice, neighbour, position, beta, deltatau, mu, change))>gsl_rng_uniform(generator)){
                        lattice[position+mu]+=change;
                        accept+=1;
                        }
                    }
                }
            }
        }
    }
    //~ /**choose site at random **/
    //~ for (int i=0; i<latticesites; i+=1){
        //~ position=gsl_rng_uniform_int(generator, latticesites);
        //~ get_neighbours(neighbour, position, Nt, Ns);
        //~ mu=position%4;
        //~ //have to account: mu is not in position as expected by plaquette, deltas, but is also randomly generated, so faster to put into same variable as position
        //~ position-=mu;
        //~ change=gsl_ran_flat(generator, -delta, delta);
        //~ if(exp(deltas(lattice, neighbour, position, beta, deltatau, mu, change))>gsl_rng_uniform(generator)){
            //~ //lattice[position+mu]=gsl_complex_mul(lattice[position+mu], gsl_complex_polar(1,change));
            //~ lattice[position+mu]+=change;
            //~ accept+=1;
        //~ }
    //~ }
    //~ double averageplaquette=0;
    //~ for (int i=0;i<latticesites;i+=4){
        //~ get_neighbours(neighbour, i, Nt, Ns);
        //~ for(int mu=0;mu<4;mu+=1){
            //~ for(int nu=mu+1;nu<4;nu+=1){
                //~ averageplaquette+=plaquette(lattice, neighbour, i, mu, nu);
            //~ }
        //~ }
    //~ }
    //fprintf(stream, "%f\t%f\n", averageplaquette/(1.5*latticesites), accept/(1.0*latticesites));
    return accept/double(latticesites);
}

void smear(double * lattice, int * neighbour, const int Nt,  const int Ns, double alpha){
    const int latticesites=4*Nt*Ns*Ns*Ns;
    double helplattice[latticesites];
    gsl_complex oneloop; 
    for (int pos=0;pos<latticesites;pos+=4){
        get_neighbours(neighbour, pos, Nt, Ns);
        //if all links should be smeared: start mu at 0, only spatial: start at 1
        //also save current state of temporal links
        helplattice[pos+0]=lattice[pos+0];
        for(int mu=1;mu<4;mu+=1){
            oneloop=gsl_complex_rect(0,0); 
            //calculate "staples"
            //question: include timelike staples(nu=0 permitted) or not?
            for(int nu=1;nu<4;nu+=1){if(mu!=nu){
                oneloop=gsl_complex_add(oneloop, gsl_complex_polar(1,lattice[pos+nu+neighbour[mu]]-lattice[pos+neighbour[nu]+mu]-lattice[pos+nu]));
                //orientation of staple important? now: chosen such that addition of link itself would provide closed plaquette
                oneloop=gsl_complex_add(oneloop, gsl_complex_polar(1,-lattice[pos+neighbour[4+nu]+mu]-lattice[pos+neighbour[4+nu]+neighbour[mu]+nu]+lattice[pos+nu+neighbour[4+nu]]));
            }}
            oneloop=gsl_complex_mul_real(oneloop, (1-alpha)/2);
            //oneloop=(1-alpha)/2 sum(staples), +=alpha U
            oneloop=gsl_complex_add(oneloop, gsl_complex_mul_real(gsl_complex_polar(1, lattice[pos+mu]), alpha));
            //project onto U(1) again: oneloop=oneloop/|oneloop|, so new abs value=1
            oneloop=gsl_complex_mul_real(oneloop, 1.0/gsl_complex_abs(oneloop));
            //~ printf("%f\t", gsl_complex_arg(oneloop));
            //store in intermediate lattice so smearing of further links is not affected by smearing already done
            helplattice[pos+mu]=gsl_complex_arg(oneloop);
        }
    }
    //copy results into lattice
    for (int i=0;i<latticesites;i+=1){
        lattice[i]=helplattice[i];
    }
}

/**calculates the loop forward direction[0]*t+direction[1]*x+direction[2]*y+direction[3]*z
 * backward direction[0]*t+direction[1]*x+direction[2]*y+direction[3]*z
 * first extracts coordinates, then goes forward, then backward, returns Re(Product(U))
 * uses explicit periodic boundary conditions
 * **/
double wilson(double * lattice, int * neighbour, int Nt, int Ns, int pos, int *direction){
    double value=0;
    int position[4]; //store t, x, y, z;
    int helppos=pos;
    
    //extract positions
    position[3]=(helppos%(4*Ns))/4;
    helppos-=4*position[3];
    position[2]=(helppos%(4*Ns*Ns))/(4*Ns);
    helppos-=4*Ns*position[2];
    position[1]=(helppos%(4*Ns*Ns*Ns))/(4*Ns*Ns);
    helppos-=4*Ns*Ns*position[1];
    position[0]=(helppos%(4*Ns*Ns*Ns*Nt))/(4*Ns*Ns*Ns);
    
    //forward
    for(int i=0;i<direction[0];i+=1){
        value+=lattice[((position[0]+i)%Nt)*4*Ns*Ns*Ns+position[1]*4*Ns*Ns+position[2]*4*Ns+position[3]*4];
    }
    for(int i=0;i<direction[1];i+=1){
        value+=lattice[((position[0]+direction[0])%Nt)*4*Ns*Ns*Ns+((position[1]+i)%Ns)*4*Ns*Ns+position[2]*4*Ns+position[3]*4+1];
    }
    for(int i=0;i<direction[2];i+=1){
        value+=lattice[((position[0]+direction[0])%Nt)*4*Ns*Ns*Ns+((position[1]+direction[1])%Ns)*4*Ns*Ns+((position[2]+i)%Ns)*4*Ns+position[3]*4+2];
    }
    for(int i=0;i<direction[3];i+=1){
        value+=lattice[((position[0]+direction[0])%Nt)*4*Ns*Ns*Ns+((position[1]+direction[1])%Ns)*4*Ns*Ns+((position[2]+direction[2])%Ns)*4*Ns+((position[3]+i)%Ns)*4+3];
    }
    
    //back
    for(int i=0;i<direction[0];i+=1){
        value-=lattice[((position[0]+direction[0]-i-1)%Nt)*4*Ns*Ns*Ns+((position[1]+direction[1])%Ns)*4*Ns*Ns+((position[2]+direction[2])%Ns)*4*Ns+((position[3]+direction[3])%Ns)*4];
    }
    for(int i=0;i<direction[1];i+=1){
        value-=lattice[position[0]*4*Ns*Ns*Ns+((position[1]+direction[1]-1-i)%Ns)*4*Ns*Ns+((position[2]+direction[2])%Ns)*4*Ns+((position[3]+direction[3])%Ns)*4+1];
    }
    for(int i=0;i<direction[2];i+=1){
        value-=lattice[position[0]*4*Ns*Ns*Ns+position[1]*4*Ns*Ns+((position[2]+direction[2]-1-i)%Ns)*4*Ns+((position[3]+direction[3])%Ns)*4+2];
    }
    for(int i=0;i<direction[3];i+=1){
        value-=lattice[position[0]*4*Ns*Ns*Ns+position[1]*4*Ns*Ns+position[2]*4*Ns+((position[3]+direction[3]-i-1)%Ns)*4+3];
    }
    return cos(value);
}

void writeconfig(double * lattice, FILE * file, int Nt, int Ns){
    int latticesites=Nt*Ns*Ns*Ns*4;
    for(int i=0;i<latticesites;i+=1){
        fprintf(file, "%f\t", lattice[i]);
    }
    fprintf(file, "\n");
}

void readconfig(double *lattice, FILE * file, int Nt, int Ns){
    int error;
    int latticesites=Nt*Ns*Ns*Ns*4;
    for(int i=0;i<latticesites;i+=1){
        error=fscanf(file, "%lf\t", &lattice[i]);
        if(error<0){fprintf(stderr, "Mistake when reading configuration!\n");}
    }
    error=fscanf(file, "\n");
    if(error<0){fprintf(stderr, "Mistake when reading configuration!\n");}
}

void gaugeinvariance(double *lattice, gsl_rng *generator, int Nt, int Ns){
    /** adds a random potential to each lattice point
     * lattice links transform like U_mu(x)->V^dagger(x)U_mu(x)V(x+mu)
     * => U_mu(x-mu)->V^dagger(x-mu)U_mu(x-mu)V(x)
     * **/
    int neighbour[8], position;
    double potential;
    for (int t=0;t<Nt;t+=1){
         for (int x=0;x<Ns;x+=1){
            for (int y=0;y<Ns;y+=1){
                for (int z=0;z<Ns;z+=1){
                    position=Ns*Ns*Ns*4*t+Ns*Ns*4*x+Ns*4*y+4*z;
                    get_neighbours(neighbour, position, Nt, Ns);
                    //V(x)
                    potential=2*M_PI*gsl_rng_uniform(generator)-M_PI;
                    for (int mu=0;mu<4;mu+=1){
                        //U_mu(x)->V^dagger(x)U_mu(x)
                        lattice[position+mu]-=potential;
                        //U_mu(x-mu)->U_mu(x-mu)V(x)
                        lattice[position+neighbour[4+mu]+mu]+=potential;
                    }
                }
            }
        }
    }
}

int main(int argc, char **argv){
    int measure=0;
    double beta;
    double deltatau=1;   
    const int Ns=10;       //number of sites for the spatial directions
    const int Nt=Ns;       //number of sites for the temporal directions
    const int latticesites=Nt*Ns*Ns*Ns*4;
    //~ printf("%d\n", latticesites);
    int thermalizations=0;
    int measurements=5000;
    double lattice[latticesites];
    //~ gsl_complex lattice[latticesites];
    //~ double potential[latticesites/4];//test for gauge invariance
    //~ int neighbour[8];
    //~ int direction[4];
    //~ printf("%d\n", latticesites);
    
    gsl_rng *generator;
	generator=gsl_rng_alloc(gsl_rng_mt19937);	
    gsl_rng_set(generator, 0);
    double delta=0.3;//M_PI*0.15;     //change delta to change acceptance rate
    double accept;
    
    //save average plaquette for plotting
    char measname[100];
    //~ sprintf(measname, "meas/measbeta%.1fNt%dNs%dmeas%d.txt", beta, Nt, Ns, measurements);
    FILE * measfile;
    char thermname[100];
    //~ sprintf(thermname, "meas/thermbeta%.1fNt%dNs%dtherm%d.txt", beta, Nt, Ns, thermalizations);
    FILE * thermfile;
    char configsavename[100];
    //~ sprintf(configsavename, "config/config%dbeta%fdeltatau%fNt%dNs%d",measurements, beta, deltatau, Nt, Ns);
    FILE * configsave; 
    
    //~ double averageplaquette=0;//, wilsonplaquette=0;
    //~ double averagelattice=0;
    //~ int plaquettes=0;
    gsl_vector * plaq=gsl_vector_alloc(measurements);   //stores measured plaquette values to take average
    gsl_vector * plaqtherm=gsl_vector_alloc(thermalizations);   //stores measured plaquette values to take average
        
    double mean_plaquette, var_plaquette, autocorrint;
	//~ gsl_vector * plaquette=gsl_vector_alloc(numberofmeasurements);
	gsl_vector * binned_plaquette_mem=gsl_vector_alloc(measurements);
	gsl_vector_view binned_plaquette;
	gsl_vector * plaquette_correlation_binned=gsl_vector_alloc(measurements/32); //longer times ain autocorrelation are uninteresting
    
    /** initialize lattice **/
    for (int i=0; i<latticesites; i+=1){
        //~ lattice[i]=gsl_complex_polar(1,i);
        lattice[i]=0;//2*M_PI*gsl_rng_uniform(generator)-M_PI;
    }

  
    double betaarray[22]={0.05,0.1,0.2,0.3,0.4,0.6,0.8,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5, 1.2, 1.4, 1.6, 1.8, 2.1};
    //~ double deltatauarray[10]={1, 0.8, 0.6, 0.5, 0.4, 0.3, 0.25, 0.2, 0.15, 0.1};
    if(measure==1){
    for(int j=21;j<22;j+=1){ 
        accept=0; 
        beta=betaarray[j];
        //~ deltatau=deltatauarray[j];
        //~ if (beta<0.8){delta=M_PI*0.8;} 
        //~ if (beta>=0.8&&beta<2.0){delta=M_PI*0.3;} 
        //~ if (beta>=2.0){delta=M_PI*0.15;} 
        sprintf(measname, "meas/measbeta%fNt%dNs%dmeas%d.txt", beta, Nt, Ns, measurements);
        measfile=fopen(measname, "w+");
        sprintf(thermname, "meas/thermbeta%fNt%dNs%dtherm%d.txt", beta, Nt, Ns, thermalizations);
        thermfile=fopen(thermname, "w+");
        sprintf(configsavename, "config/config%dbeta%fdeltatau%fNt%dNs%d",measurements, beta, deltatau, Nt, Ns);
        configsave=fopen(configsavename, "w+");
        for(int i=0; i<thermalizations;i+=1){
            sweep(lattice, beta, deltatau, generator, measfile, Nt, Ns,  delta);
            gsl_vector_set(plaqtherm, i, averageplaquette(lattice, Nt, Ns));
        }
        
        for(int i=0; i<measurements;i+=1){
            accept+= sweep(lattice, beta, deltatau, generator, measfile, Nt, Ns,  delta); 
            gsl_vector_set(plaq, i, averageplaquette(lattice, Nt, Ns));
        }
        
        fprintf(stdout, "%f\t%f\n", beta, accept/double(measurements));
        gsl_vector_fprintf(measfile, plaq, "%e");
        gsl_vector_fprintf(thermfile, plaqtherm, "%e");
        fclose(configsave);
        fclose(measfile);
        fclose(thermfile);
    }}
    
    //~ printf("%f\t%f\t%d\t%d\t%f\t%f\n", beta, deltatau,Ns, Nt, mean(plaq, measurements), deviation(plaq, measurements, mean(plaq, measurements)));
    int binarray[16]={1,2,4,8,16,32,52,72,92,112,132, 152, 172, 192, 212,232};
    int binsize;
    for(int j=0;j<21;j+=1){ 
        beta=betaarray[j];
        //~ printf("%f,", beta);
        //~ deltatau=deltatauarray[j];
        //~ sprintf(measname, "meas/measbeta%fNt%dNs%dmeas%d.txt", beta, Nt, Ns, measurements);
        sprintf(measname, "meas/measbeta%fNt%dNs%d.txt", beta, Nt, Ns);
        measfile=fopen(measname, "r");
        gsl_vector_fscanf(measfile, plaq);
        
        mean_plaquette=mean(plaq, measurements);
        var_plaquette=deviation(plaq, measurements, mean_plaquette);
        autocorrelation(plaq, plaquette_correlation_binned, mean_plaquette);
        autocorrint=tauint(plaquette_correlation_binned, 10);
        fprintf(stdout, "%f\t%f\t%.2d\t%.2d\t%.4d\t%e\t%e\t%e\n", beta, deltatau, Nt, Ns, 0, mean_plaquette, sqrt(var_plaquette), autocorrint);
        
        for (int bin=0;bin<16;bin++){
            binsize=binarray[bin];
            binned_plaquette=gsl_vector_subvector(binned_plaquette_mem, 0, plaq->size/binsize);
            binning(plaq, &binned_plaquette.vector, binsize);
            bootstrap(&binned_plaquette.vector, generator, 1000, &mean_plaquette, &var_plaquette);
            autocorrelation(&binned_plaquette.vector, plaquette_correlation_binned, mean_plaquette);
            autocorrint=tauint(plaquette_correlation_binned, 10);
            
            //~ fprintf(plaquette_data, "\nbinsize %d\n", binsize);
            //~ fprintf(plaquette_autocorrelation, "\nbinsize %d\n", binsize);
            //~ gsl_vector_fprintf(plaquette_data, &binned_plaquette.vector, "%f");
            //~ gsl_vector_fprintf(plaquette_autocorrelation, plaquette_correlation_binned, "%f");
            
            
            fprintf(stdout, "%f\t%f\t%.2d\t%.2d\t%.4d\t%e\t%e\t%e\n", beta, deltatau, Nt, Ns, binsize, mean_plaquette, sqrt(var_plaquette), autocorrint);
        }
        fclose(measfile);
        //~ printf("%f\n", beta);
    }
        //~ printf("%f\n", beta);
    //check if functions work
    //~ printf("%f\n", averageplaquette(lattice, Nt, Ns));
    //~ for(int i=0;i<20;i+=1){
        //~ gaugeinvariance(lattice, generator, Nt, Ns);
        //~ printf("%f\n", averageplaquette(lattice, Nt, Ns));
    //~ }
    //~ smear(lattice, neighbour, Nt, Ns, 0.7);
    //~ printf("%f\n", averageplaquette(lattice, Nt, Ns));
    
    

		
    //optional: write to configuration, read and do measurements after all sweeps have been done
    //takes lots of disk space
    //~ sprintf(configsavename, "config/config%dbeta%fdeltatau%fNt%dNs%d",measurements, beta, deltatau, Nt, Ns);
    //~ configsave=fopen(configsavename, "w+");
    //~ for(int i=0; i<thermalizations;i+=1){
        //~ sweep(lattice, beta, deltatau, generator, measfile, Nt, Ns,  delta);
    //~ }
    //~ for(int i=0; i<measurements;i+=1){     
        //~ gsl_vector_set(plaq, i, sweep(lattice, beta, deltatau, generator, measfile, Nt, Ns,  delta));
        //~ writeconfig(lattice, configsave, Nt, Ns);
    //~ }
    //~ fclose(configsave);
    
    //~ printf("%f\t%f\t%d\t%d\t%f\t%f\n", beta, deltatau,Ns, Nt, mean(plaq, measurements), deviation(plaq, measurements, mean(plaq, measurements)));
    
    //~ configsave=fopen(configsavename, "r+");
    //~ rewind(configsave);
    //~ for(int i=0; i<measurements;i+=1){ 
        //~ averageplaquette=0;
        //~ readconfig(lattice, configsave, Nt, Ns); 
        //~ for (int i=0;i<latticesites;i+=4){
        //~ get_neighbours(neighbour, i, Nt, Ns);
        //~ for(int mu=0;mu<4;mu+=1){
            //~ for(int nu=mu+1;nu<4;nu+=1){
                //~ averageplaquette+=plaquette(lattice, neighbour, i, mu, nu);
            //~ }
        //~ }
    //~ } 
    //~ averageplaquette/=latticesites*1.5;  
    //~ gsl_vector_set(plaq, i, averageplaquette);
    //~ }
    //~ fclose(configsave);
    //~ printf("%f\t%f\t%d\t%d\t%f\t%f\n", beta, deltatau,Ns, Nt, mean(plaq, measurements), deviation(plaq, measurements, mean(plaq, measurements)));
        
        
    gsl_rng_free(generator);
    //~ fclose(measfile);
    //~ fclose(configsave);
    gsl_vector_free(plaq);
    gsl_vector_free(plaqtherm);
    gsl_vector_free(binned_plaquette_mem);
    gsl_vector_free(plaquette_correlation_binned);
	return 0;
}

