/* 

This program computes average effective population sizes from a temporal series of allele frequencies under the assumptions of:

    1) No sampling variance on the part of the investigator;

    2) A fixed fraction of annual retention of zygotes in the egg bank from year to year. 

    3) No influence from selection or mutation.


The grand average Ne is determined from the ratio of the observed variance in allele-frequency change relative to the expectations under binomial sampling. 

Numerous independent sampling bouts of starting temporal series are invoked, and the grand average Ne is estimated by averaging over many such series.  

Computations are done in parallel for an array of different starting allele frequencies.

*/



/* ********************************************************************************************************************** */

#define ne		    1000000 			/* effective population size, Ne, assumed to be diploid */

#define phatch	    0.5					/* fraction of resting eggs hatching in egg bank per year; assumed to decline with constant probability over time */

#define nseries		40000000			/* number of independent sequential series to sample to get mean Ne estimate */

#define burnin		13  				/* number of initial burn-in increments, before starting analyses; also equal to the number of subsequent events in the temporal series actually used  */
                                        /* needs to be set high enough so that survival to the end point is arbitrarily small, say 0.0001 */

#define nreps       1000000             /* number of replicate series between printout to slurm */



#include	<stdio.h>
#include 	<math.h>
#include 	<sys/types.h>
#include	<stdlib.h>
#include	<time.h>
#include    <gsl/gsl_rng.h>
#include    <gsl/gsl_randist.h>
#include    <string.h>


/* ********************************************************************************************************** */



/* point to the output file */

FILE *stream;
char filename[100];



/* Set up the parallel runs. */

void main(int argc, char *argv[])
{
    
    

int f0, f1;
if (argc > 1)
{
    f0 = atoi(argv[1]);
    f1=f0; 
    sprintf(filename, "dataout_%d.txt", f0); 
}
else
{
    f0 = 1;
    f1 = 9;
    sprintf(filename, "dataout.txt");
}


/* Set up the binomial random number generator. */

static gsl_rng* rand_new;                                                       
gsl_rng_env_setup();                                                            
if(!rand_new)                                                                   
{                                                                               
    rand_new = gsl_rng_alloc(gsl_rng_taus2);                                    
    gsl_rng_set(rand_new, time(NULL));                                          
}     


/* Set up the normal random number generator. */

static gsl_rng* rand_new2;
gsl_rng_env_setup();
if (!rand_new2)
{
	rand_new2 = gsl_rng_alloc(gsl_rng_taus2);
	gsl_rng_set(rand_new2, time(NULL));
}



   

/* ***************************************************************************************** */



/* MAIN BODY OF PROGRAM. */

double initialfreq[40];									/* initial allele frequency */
double freq[2000];										/* temporal series of allele frequencies used to compute the final series */

double weight[100];                                     /* fraction of resting eggs surviving to hatch each year */

long igen;												/* generation counter */

int itera;												/* counter for initial allele-frequency runs */

int ig, jg, kg;											/* counters for the year-specific frequencies */

long draw;												/* new count drawn from the binomial */

double nextfreq;										/* numerator used in estimating expected allele frequencies based on sampling from egg bank */
double sumweight;										/* denominator used in estimating expected allele frequencies based on sampling from egg bank*/

double pp;											    /* probability associated with the binomial for drift */

double phi;												/* statistic used in Ne estimation */
double pmean;											/* denomintor used in phi calculation */
double sumphi, totphi;									/* summations for grand mean phi estimate */

int tint;												/* counter for the number of temporal series */

int nsamp;                                              /* counter between print to slurm */

double nestimate;										/* grand average Ne estimate */



/* Open the output file. */

remove("dataout.txt ");


/* Initial allele frequencies */

initialfreq[9] = 0.50;
initialfreq[8] = 0.40;
initialfreq[7] = 0.35;
initialfreq[6] = 0.30;
initialfreq[5] = 0.25;
initialfreq[4] = 0.20;
initialfreq[3] = 0.15;
initialfreq[2] = 0.10;
initialfreq[1] = 0.05;


weight[1] = phatch;

for (ig = 2; ig <= (burnin + 1); ++ig) {
	weight[ig] = weight[ig-1] * (1.0 - phatch); }                               /* set the surviving egg bank for each contributing year in the final series */


double start, stop, time;                                                   



for (itera = f0; itera <= f1; ++itera) {							            /* Start iterations over the set of population sizes and mutation rates. */

stream = fopen(filename, "a");		



	/* Set the initial genotype frequencies, and zero out the counters. */

	sumphi = 0.0;
	totphi = 0.0;
	
	tint = 0;
	nsamp = 0;
	

	/* ******************************************************************************************************************************************* */


	/* Iterate phi calculations over nseries of independent strings of ancestral allele frequencies. */

	while (tint < nseries)  												    /* iterate until the stopping criterion has been met. */
	{
		nsamp = nsamp + 1; 
		
		/* Sample the population for new ancestral annual genotype frequencies, and the subsequent time series allowing for egg-bank retention. */

		for (ig = 1; ig <= (4 * burnin); ++ig) {								/* zero out the initial allele frequencies */
			freq[ig] = 0.0; }

		for (ig = 1; ig <= burnin; ++ig) {										/* sets frequencies for an initial burnin generations equal to the time-zero value */
			freq[ig] = initialfreq[itera]; }

		for (ig = (burnin + 1); ig <= (4 * burnin); ++ig) {						/* computes the next series of frequencies used in time-averaged draws to be used in final computations */

			nextfreq = 0.0;
			sumweight = 0.0;

			for (jg = 1; jg <= burnin; ++jg) {
				nextfreq = nextfreq + (weight[jg] * freq[ig - jg]);
				sumweight = sumweight + weight[jg]; }

			pp = nextfreq / sumweight;                                          /* this is the expected gene frequency in hatchlings for the year */
			draw = gsl_ran_binomial_tpe(rand_new, pp, (2*ne));
			freq[ig] = 0.5 * ((double)draw) / ((double)ne);                     /* this is the new frequency after drift */
		}




		/* Get the phi statistic used in estimation of Ne. */

			for (ig = burnin + 2; ig <= ((2*burnin)+1); ++ig) {
				if ( (freq[ig] * freq[ig + 1]) > 0.0000000001 ) {                               /* don't do computations if an allele frequency hits 0.0 */
				    pmean = 0.5 * (freq[ig] + freq[ig + 1]);
				    phi = pow((freq[ig] - freq[ig + 1]), 2.0) / (pmean * (1.0 - pmean));        /* this is the scaled variance of allele-frequency change */
				    sumphi = sumphi + phi;
				    totphi = totphi + 1.0; 	}
			}

			tint = tint + 1; 

		  if (nsamp == nreps) {
            printf("%9d, %9d, %12.5f, %12.5f, %9d, %6d, %10.4f, %10.1f\n", itera, ne, phatch, initialfreq[itera], nseries, burnin, (totphi / ((double)(burnin*nseries)) ), (0.5 / (sumphi / totphi)) );
            nsamp = 0;}

	}


	/* Get the grand estimate of Ne */

	nestimate = 0.5 / (sumphi / totphi);


	/* Output the results */
	
	    /* Ne; annual hatching fraction; initial frequency; number of temporal series run; number of generations contributing to annual hatch; fraction of intervals with a non-0.0 frequency; */
	    /* number of interval attempts for estimation; mean Ne estimate; ratio of estimated to true Ne */
	
	fprintf(stream, " %9d, %12.5f, %12.5f, %9d, %6d, %10.4f, %10.1f, %10.1f,, %5.4f\n", ne, phatch, initialfreq[itera], nseries, burnin, (totphi / ((double)(burnin*nseries)) ), totphi, nestimate, (nestimate /((double) ne)) );
	
	printf("\n");

	fclose(stream);

}									                                            /* End the set of iterations over all allele frequencies. */


exit(0);

}





