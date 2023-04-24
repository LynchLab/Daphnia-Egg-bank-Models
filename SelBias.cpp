/* 

This program computes the mean and variances of estimated selection coefficients from a temporal series of allele frequencies under the assumptions of:

    1) No sampling variance on the part of the investigator, or random genetic drift;

    2) A fixed fraction of annual retention of zygotes in the egg bank from year to year. 

    3) No influence from mutation.


The grand average s and SD(s) is determined from from a series of independent sampling bouts of starting temporal series.  

Estimates of mean s are obtained in three ways: 1) logit two-point estimates; direct discrete-generation two-point estimates; and regression of logits.


Computations are done in parallel for an array of different starting allele frequencies.

*/



/* ********************************************************************************************************************** */

#define sco			0.0001	            /* mean selection coefficient alleles */

#define sds			0.01	           /* temporal standard deviation of s  */

#define phatch	    0.5					/* fraction of resting eggs hatching in egg bank per year; assumed to decline with constant probability over time */

#define ngens       9                   /* number of generations used in logit regression estimate of s */

#define nseries		50000000			/* number of independent sequential series to sample to get mean Ne estimate */

#define burnin		11 					/* number of initial burn-in increments, before starting analyses; also equal to the number of subsequent events in the temporal series actually used  */
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
	gsl_rng_set(rand_new2, time(NULL) + (20*f0));
}









/* ***************************************************************************************** */



/* MAIN BODY OF PROGRAM. */

double initialfreq[200];									/* initial allele frequency */
double freq[2000];										/* temporal series of allele frequencies used to compute the final series */

double weight[100];                                     /* fraction of resting eggs surviving to hatch each year */

long igen;												/* generation counter */

int itera;												/* counter for initial allele-frequency runs */

int ig, jg, kg;											/* counters for the year-specific frequencies */

long draw;												/* new count drawn from the binomial */

double nextfreq;										/* numerator used in estimating expected allele frequencies based on sampling from egg bank */
double sumweight;										/* denominator used in estimating expected allele frequencies based on sampling from egg bank*/

double pp;											    /* probability associated with the binomial for drift */

double selco;											/* random selection coefficient for the generation */
double selcosum, selconum;

double meanfit;											/* mean fitness */

double eta0, eta1;										/* log components used to estimate s by the logit method */

double sest;											/* estimated s for logit method */
double sumsest, sumsest2, totphi;						/* summations for grand mean estimate of s mean and SD */
double smean, ssqmean, sdsest;							/* grand means for mean s, mean-square s, and SD(s) estimates */

double sestd;											/* estimated s for direct, discrete-generation method */
double sumsestd, sumsest2d;						        /* summations for grand mean estimate of s mean and SD */
double smeand, ssqmeand, sdsestd;						/* grand means for mean s, mean-square s, and SD(s) estimates */

int tint;												/* counter for the number of temporal series */

int nsamp;                                              /* counter between print to slurm */

int check;                                              /* counters and summations for logit regression estimates of s */
double tmean, tmeansq, vart;
double cmean, cmeansq, cmeancp;
double sestr, sumsestr, numr, smeanr;

double start, stop, time;                                                   





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




for (itera = f0; itera <= f1; ++itera) {							            /* Start iterations over the set of population sizes and mutation rates. */

stream = fopen(filename, "a");		


	/* Set the initial genotype frequencies, and zero out the counters. */

	sumsest = 0.0;
	sumsest2 = 0.0;
	
	sumsestd = 0.0;
	sumsest2d = 0.0;
	sumsestr = 0.0;

    numr = 0.0;
	
	totphi = 0.0;
	
	selcosum = 0.0;
	selconum = 0.0;
	
	tint = 0;
	nsamp = 0;


    tmean = 0.0;						                                        /* Mean and variance of time for regression s */	
    tmeansq = 0.0;

    for (ig = 1; ig <= ngens; ++ig) {			
	    tmean = tmean + ((double) ig);
	    tmeansq = tmeansq + pow( ((double) ig), 2.0); }

    tmean = tmean / ((double) ngens);
    tmeansq = tmeansq / ((double) ngens);

    vart = tmeansq - (tmean * tmean);




	/* ******************************************************************************************************************************************* */


	/* Iterate phi calculations over nseries of independent strings of ancestral allele frequencies. */

	while (tint < nseries)  												    /* iterate until the stopping criterion has been met. */
	{
		nsamp = nsamp + 1; 
		
		
		
		/* Sample the population for new ancestral annual genotype frequencies, and the subsequent time series allowing for egg-bank retention. */

		for (ig = 1; ig <= (6 * burnin); ++ig) {								/* zero out the initial allele frequencies */
			freq[ig] = 0.0; }

		for (ig = 1; ig <= burnin; ++ig) {										/* sets frequencies for initial burnin generations equal to the time-zero value */
			freq[ig] = initialfreq[itera]; }

		for (ig = (burnin + 1); ig <= (6 * burnin); ++ig) {						/* computes the next series of frequencies used in time-averaged draws to be used in final computations */

			nextfreq = 0.0;
			sumweight = 0.0;

			for (jg = 1; jg <= burnin; ++jg) {
				nextfreq = nextfreq + (weight[jg] * freq[ig - jg]);
				sumweight = sumweight + weight[jg]; }

			pp = nextfreq / sumweight;                                          /* this is the expected gene frequency in hatchlings for the year */
			
			/* Draw the random selection coefficient. */
			
			selco = gsl_ran_gaussian(rand_new2, sds);
			selco = selco + sco;
			
			selcosum = selcosum + selco;
			selconum = selconum + 1.0;

			meanfit = 1.0 + (pp * selco);
			
			freq[ig] = pp * (1.0 + selco) / meanfit;				            /* this is the new frequency after selection */
		}



		/* Get the estimated selection coefficients. */

			for (ig = burnin + 2; ig <= ((2 * burnin) + 1); ++ig) {
				if ( (freq[ig] * freq[ig - 1]) > 0.0000000001 ) {                               /* don't do computations if an allele frequency hits 0.0 */

					eta0 = log((1.0 - freq[ig-1]) / freq[ig-1]);
					eta1 = log((1.0 - freq[ig]) / freq[ig]);
					
					sest = eta0 - eta1;															/* this is the estimated selection coefficient with the logit approach */
				    sumsest = sumsest + sest;
					sumsest2 = sumsest2 + (sest * sest);
					
					sestd = (freq[ig] - freq[ig-1]) / (freq[ig-1] * (1.0 - freq[ig]));			/* this is the estimated selection coefficient with the direct, discrete-generatin approach */
				    sumsestd = sumsestd + sestd;
					sumsest2d = sumsest2d + (sestd * sestd);
					
				    totphi = totphi + 1.0; 	}
			}
			
		
		    for (ig = burnin + 2; ig <= ((2 * burnin) + 1); ++ig) {                             /* calculate the logit regression estimate of s */

	            check = 1;
	            for (jg = ig; jg <= (ig + ngens); ++jg) {
	                if (freq[jg] < 0.000001) {
	                check = 0; }     }
			
	           if (check == 1) {  
		        
		            cmean = 0.0;
		            cmeansq = 0.0;
		            cmeancp = 0.0;

		            for (jg = ig; jg <= (ig + ngens - 1); ++jg) {
			            eta1 = log((1.0 - freq[jg]) / freq[jg]);
			            cmean = cmean + eta1;
			            cmeansq = cmeansq + pow(eta1, 2.0);
			            cmeancp = cmeancp + (eta1 * ( ((double) jg) - ((double) ig) + 1.0)  )      ; }
		
		            cmean = cmean / ((double) ngens);
		            cmeansq = cmeansq / ((double) ngens);
		            cmeancp = cmeancp / ((double) ngens);

		            sestr =	( cmeancp - (cmean * tmean) ) / vart;
		            sumsestr = sumsestr + sestr;
		            numr = numr + 1.0;   }
            }



			tint = tint + 1; 

		  if (nsamp == nreps) {
            printf("%9d, %16.11f, %12.5f, %12.5f, %9d, %6d, %10.1f, %16.11f, %16.11f\n", itera, sco, phatch, initialfreq[itera], nseries, burnin, totphi, sestr, (sumsestr/numr) ); 
		      
		    nsamp = 0; }
            

	}


	/* Get the grand estimate of s and its SD for logit, direct, and regression methods */

	smean = sumsest / totphi;
	ssqmean = sumsest2 / totphi;
	sdsest = ssqmean - (smean * smean);
	sdsest = pow(sdsest, 0.5);
	
	smeand = sumsestd / totphi;
	ssqmeand = sumsest2d / totphi;
	sdsestd = ssqmeand - (smeand * smeand);
	sdsestd = pow(sdsestd, 0.5);
	
	smeanr = -sumsestr / numr;
	
	

	/* Output the results */
	
	    /* Mean selection coefficient; temporal SD of s; annual hatching fraction; initial frequency; number of temporal series run; number of generations contributing to annual hatch; mean s estimate; SD(s) estimate; */
	    /* ratio of estimated mean s to true mean; ratio of estimated SD(s) to true SD */
	
	fprintf(stream, " %12.11f, %12.11f, %12.5f, %12.5f, %9d, %6d, %10.4f, %10.1f,, %12.11f, %12.11f,, %12.11f, %12.11f,, %12.11f, %12.11f,, %12.11f, %12.11f,, %12.11f,, %12.11f,, %12.11f \n", 
	    sco, sds, phatch, initialfreq[itera], nseries, burnin, (totphi / ((double)(burnin*nseries)) ), totphi, smean, sdsest, (smean/sco), (sdsest/sds), smeand, sdsestd, (smeand/sco), (sdsestd/sds), smeanr, (smeanr/sco), (selcosum/selconum)   );
	
	printf("\n");

	fclose(stream);

}									/* End the set of iterations over all allele frequencies. */


exit(0);

}





