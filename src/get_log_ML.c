#include "header.h"

SEXP C_get_log_ML(SEXP logpermsSEXP, SEXP nSEXP, SEXP SSEXP, SEXP debugSEXP){
	PROTECT(logpermsSEXP);
	PROTECT(nSEXP);
	PROTECT(SSEXP);
	PROTECT(debugSEXP);

	double * logperms = (double*)REAL(logpermsSEXP);
	int n = *(INTEGER(nSEXP));
	int S = *(INTEGER(SSEXP));
	int debug = *(INTEGER(debugSEXP));

	SEXP resultSEXP = PROTECT(allocVector(REALSXP, 1));
	double *result = REAL(resultSEXP);
	result[0] = 0;


	// as permutation numbers are always integers, their log values are always non-negative:
  	double maxval = -1;

  	for (int i = 0; i < S; ++i)
  	{
  		if(!ISNAN(logperms[i])){
  			if(logperms[i]> maxval){
  				maxval = logperms[i];
  			}
  		}
  		
  	}
  	

  	if(maxval<=-1){
  		// all logperms are NAN
  		*result = NA_REAL;
  		UNPROTECT(5);
  		//error("Error! No non-zero perms in logperms");
  		return resultSEXP;
  	}
  	*result = Clog_sum_exp(logperms, S, maxval) - log((double)S);

  	// compute log factorials

	double * log_factorials = (double * ) R_Calloc(n+1, double);

	log_factorials[0]=0.0;
	for (int i = 1; i <= n; ++i)
	{
		log_factorials[i] = log_factorials[i-1] +log((double)(i));
	}

	*result = *result -log_factorials[n];
  	
  	R_Free(log_factorials);
	UNPROTECT(5);
  	return resultSEXP;
}


SEXP C_get_log_ML_bioassay(SEXP logpermsSEXP, SEXP successesSEXP, SEXP trialsSEXP, SEXP nSEXP, 
	SEXP num_trialsSEXP, SEXP SSEXP, SEXP debugSEXP){
	PROTECT(logpermsSEXP);
	PROTECT(successesSEXP);
	PROTECT(trialsSEXP);
	PROTECT(nSEXP);
	PROTECT(num_trialsSEXP);
	PROTECT(SSEXP);
	PROTECT(debugSEXP);

	double * logperms = (double*)REAL(logpermsSEXP);
	int n = *(INTEGER(nSEXP));
	int num_trials = *(INTEGER(num_trialsSEXP));
	int S = *(INTEGER(SSEXP));
	int debug = *(INTEGER(debugSEXP));
	int * successes = INTEGER(successesSEXP);
	int * trials = INTEGER(trialsSEXP);


	SEXP resultSEXP = PROTECT(allocVector(REALSXP, 1));
	double *result = REAL(resultSEXP);
	result[0] = 0;


  	double maxval = -1;

  	for (int i = 0; i < S; ++i)
  	{
  		if(!ISNAN(logperms[i])){
  			if(logperms[i]> maxval){
  				maxval = logperms[i];
  			}
  		}
  	}
  
 	

  	if(maxval<=-1){
  		*result = NA_REAL;
  		UNPROTECT(8);
  		//error("Error! No non-zero perms in logperms");
  		return resultSEXP;
  	}
  	*result = Clog_sum_exp(logperms, S, maxval) - log((double)S);

  	// compute log factorials


	double * log_factorials = (double*) R_Calloc(n+1, double);



	log_factorials[0]=0.0;
	for (int i = 1; i <= n; ++i)
	{
		log_factorials[i] = log_factorials[i-1] +log((double)(i));
	}

	*result = *result -log_factorials[n];
  	
  	for (int j = 0; j < num_trials; ++j)
  	{
  		*result = *result + log_factorials[trials[j]] - log_factorials[successes[j]] - log_factorials[trials[j] - successes[j]];
  	}

  	R_Free(log_factorials);
  	UNPROTECT(8);
  	return resultSEXP;
  	
}