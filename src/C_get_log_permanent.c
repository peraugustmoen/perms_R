#include "header.h"

SEXP C_get_log_perms(SEXP XSEXP, SEXP tSEXP, SEXP ySEXP, SEXP nSEXP, SEXP SSEXP, 
	SEXP debugSEXP, SEXP constant_tsSEXP){

	PROTECT(XSEXP);
	PROTECT(tSEXP);
	PROTECT(ySEXP);
	PROTECT(nSEXP);
	PROTECT(SSEXP);
	PROTECT(debugSEXP);
	PROTECT(constant_tsSEXP);
	
	double * X = REAL(XSEXP);
	double * t = REAL(tSEXP);
	int * y = INTEGER(ySEXP);
	int n = *INTEGER(nSEXP);
	int S = *INTEGER(SSEXP);
	int debug = *INTEGER(debugSEXP);
	int constant_ts = *INTEGER(constant_tsSEXP);

	if(debug){
		Rprintf("Debug = %d\n",debug);
	}
	double * a = (double*) R_Calloc(n, double);
  	double * b = (double*) R_Calloc(n, double);
	
	if (constant_ts)
    {
	    for (int i = 0; i < n; ++i)
	    {
	    	if(y[i]==1){
	    		b[i] = t[i];
	    		a[i] = -DBL_MAX;
	    	}else{
	    		b[i] = DBL_MAX;
	    		a[i] = t[i];

	    	}
	    }

    	R_qsort(a, 1, n);
		R_qsort(b, 1, n);
		
    }
	
	
    if(debug){
    	Rprintf("a:\n");
    	print_float_vector(n, a);
    	Rprintf("b:\n");
    	print_float_vector(n, b);
    }
	

	SEXP logpermsSEXP = PROTECT(allocVector(REALSXP, S));
	double * logperms = REAL(logpermsSEXP);
	memset(logperms, 0, sizeof(double)*S);


	SEXP a_union_bSEXP = PROTECT(allocVector(REALSXP, 2*n));
	double * a_union_b = REAL(a_union_bSEXP);
	memset(a_union_b, 0, sizeof(double)*2*n);

	int len_a_union_b =0;

	get_union(n, a, b, &len_a_union_b, a_union_b);

	SEXP alphaSEXP = PROTECT(allocVector(INTSXP, n));
	SEXP betaSEXP = PROTECT(allocVector(INTSXP, n));
	SEXP gammaSEXP = PROTECT(allocVector(INTSXP, n));
	SEXP logfactorialsSEXP = PROTECT(allocVector(REALSXP, n+1));
	SEXP mSEXP = PROTECT(allocVector(INTSXP, 1));
	SEXP k_SEXP = PROTECT(allocVector(INTSXP, 1));

	
	int * alpha = INTEGER(alphaSEXP);
	int * beta = INTEGER(betaSEXP);
	int * gamma = INTEGER(gammaSEXP);
	double * log_factorials = REAL(logfactorialsSEXP);
	int * m = INTEGER(mSEXP);
	int * k = INTEGER(k_SEXP);


	dictionary * new_log_subperms = init_dictionary(n);
	dictionary * old_log_subperms = init_dictionary(n);

	
	memset(alpha, 0, sizeof(int)*n);
	memset(beta, 0, sizeof(int)*n);
	memset(gamma, 0, sizeof(int)*n);
	memset(log_factorials, 0, sizeof(double)*(n+1));
	memset(m, 0, sizeof(int));
	memset(k, 0, sizeof(int));

	log_factorials[0]=0.0;
	for (int i = 1; i <= n; ++i)
	{
		log_factorials[i] = log_factorials[i-1] +log((double)(i));
	}

	
	
	SEXP historySEXP = PROTECT(allocVector(INTSXP, 3*n));
	SEXP amount_historySEXP = PROTECT(allocVector(INTSXP, 6*n));

	int * history = INTEGER(historySEXP);
	int * amount_history = INTEGER(amount_historySEXP);

	memset(history, 0, sizeof(int)*3*n);
	memset(amount_history, 0, sizeof(int)*6*n);


	for (int s = 0; s < S; ++s)
	{
		
		double * x = &(X[s*n]);
		R_qsort(x, 1, n);
		
		memset(alpha, 0, sizeof(int)*n);
		memset(beta, 0, sizeof(int)*n);
		memset(gamma, 0, sizeof(int)*n);
		memset(m, 0, sizeof(int));
		memset(k, 0, sizeof(int));

		if(constant_ts){
		    if(!nonzero_perm(x, a,  b, n)){
				logperms[s] = NA_REAL;
				continue;
			}
			get_alphabetagamma(x, n, a, b, a_union_b, len_a_union_b, alpha, 
		    beta, gamma,  k, m, debug);
			if((*k)==1){
				logperms[s] = log_factorials[n];
				continue;
			}
		}else{
			for (int i = 0; i < n; ++i)
		    {
		    	if(y[i]==1){
		    		b[i] = t[i + s*n];
		    		a[i] = -DBL_MAX;
		    	}else{
		    		b[i] = DBL_MAX;
		    		a[i] = t[i + s*n];

		    	}
		    }

	    	R_qsort(a, 1, n);
			R_qsort(b, 1, n);

			if(!nonzero_perm(x, a,  b, n)){
				logperms[s] = NA_REAL;
				continue;
			}
			memset(a_union_b, 0, sizeof(double)*2*n);

			len_a_union_b =0;

			get_union(n, a, b, &len_a_union_b, a_union_b);

			get_alphabetagamma(x, n, a, b, a_union_b, len_a_union_b, alpha, 
		    beta, gamma,  k, m, debug);
		    if(!nonzero_perm(x, a,  b, n)){
				logperms[s] = NA_REAL;
				continue;
			}
			if((*k)==1){
				logperms[s] = log_factorials[n];
				continue;
			}
		}
		
		/*get_alphabetagamma(x, n, a, b, a_union_b, len_a_union_b, alpha, 
	    beta, gamma,  k, m, debug);*/


	    if(debug){
	    	Rprintf("S=%d, t=%d\n", S, s);
	    	Rprintf("len_a_union_b = %d\n", len_a_union_b);
	    	Rprintf("x:\n");
	    	print_float_vector(n,x);
	    	Rprintf("a:\n");
	    	print_float_vector(n,a);
	    	Rprintf("b:\n");
	    	print_float_vector(n,b);
	    	Rprintf("a_union_b:\n");
	    	print_float_vector(2*n,a_union_b);
	    	Rprintf("len a_union_b:%d\n", len_a_union_b);
	    	Rprintf("alpha:\n");
	    	print_int_vector(n,  alpha);
	    	Rprintf("beta:\n");
	    	print_int_vector(n,  beta);
	    	Rprintf("gamma:\n");
	    	print_int_vector(n,  gamma);
	    	Rprintf("m:%d\n", *m);
	    	Rprintf("k:%d\n", *k);
	    	
	    }

		int history_len = 0;

	
		memset(history, 0, sizeof(int)*3*n);
		memset(amount_history, 0, sizeof(int)*6*n);


		if(debug){
			Rprintf("REDUCING NOW\n");
		}
		
		//return(XSEXP);
		int result = reduction(alpha,  beta,  gamma, m, n, k, history,
				   amount_history, &history_len, debug);

		if(result != 0){

			

			memset(alpha, 0, sizeof(int)*n);
			memset(beta, 0, sizeof(int)*n);
			memset(gamma, 0, sizeof(int)*n);
			memset(m, 0, sizeof(int));
			memset(k, 0, sizeof(int));
			debug = 1;
			get_alphabetagamma(x, n, a, b, a_union_b, len_a_union_b, alpha, 
		    beta, gamma,  k, m, debug);


		    if(debug){
		    	Rprintf("len_a_union_b = %d\n", len_a_union_b);
		    	Rprintf("x:\n");
		    	print_float_vector(n,x);
		    	Rprintf("a:\n");
		    	print_float_vector(n,a);
		    	Rprintf("b:\n");
		    	print_float_vector(n,b);
		    	Rprintf("a_union_b:\n");
		    	print_float_vector(2*n,a_union_b);
		    	Rprintf("len a_union_b:%d\n", len_a_union_b);
		    	Rprintf("alpha:\n");
		    	print_int_vector(n,  alpha);
		    	Rprintf("beta:\n");
		    	print_int_vector(n,  beta);
		    	Rprintf("gamma:\n");
		    	print_int_vector(n,  gamma);
		    	Rprintf("m:%d\n", *m);
		    	Rprintf("k:%d\n", *k);
		    	
		    }

			int history_len = 0;

		
			memset(history, 0, sizeof(int)*3*n);
			memset(amount_history, 0, sizeof(int)*6*n);


			if(debug){
				Rprintf("REDUCING NOW\n");
			}
			
			//return(XSEXP);
			reduction(alpha,  beta,  gamma, m, n, k, history,
					   amount_history, &history_len, debug);


			free_dictionary(new_log_subperms);
			free_dictionary(old_log_subperms);

			error("Error in computing permanent. See output from code for debugging.");
			UNPROTECT(17);
			return NULL;
		}

		if(debug){
			Rprintf("history len = %d\n", history_len);

			Rprintf("REDUCED SUBPERMS\n");
		}
		sparse_get_reduced_log_subperms( new_log_subperms,  alpha, beta, gamma,
						log_factorials, n,  m, k);
		//Rprintf("RESULT:\n");
		//print_matrix(n+1, n+1, new_log_subperms);

		dictionary * tmp  = old_log_subperms;
		old_log_subperms = new_log_subperms;
		new_log_subperms = tmp;



		if(debug){
			Rprintf("==========\nReverse reduction:\n==========\n");
		}
		//Rprintf("old = %d\n", old_log_subperms);
		//Rprintf("new = %d\n", new_log_subperms);
		dictionary * the_log_subperms = sparse_reverse_reduction(old_log_subperms, new_log_subperms, alpha,
						   beta,  gamma, m,  n, k,  history,
				           amount_history, &history_len, log_factorials);

		


		
		double logperm =  Csparse_log_sum_exp(the_log_subperms);
		logperms[s] = logperm;
		if(debug){
			Rprintf("logperm = %f\n", logperm);

		}



	}
	free_dictionary(new_log_subperms);
	free_dictionary(old_log_subperms);

	R_Free(a);
	R_Free(b);
	UNPROTECT(17);
	return(logpermsSEXP);

}

SEXP C_get_log_perms_bioassay(SEXP XSEXP, SEXP levelsSEXP, SEXP successesSEXP, SEXP trialsSEXP,
	SEXP nSEXP, SEXP num_trialsSEXP, SEXP SSEXP, 
	SEXP debugSEXP){

	PROTECT(XSEXP);
	PROTECT(levelsSEXP);
	PROTECT(successesSEXP);
	PROTECT(trialsSEXP);
	PROTECT(nSEXP);
	PROTECT(num_trialsSEXP);
	PROTECT(SSEXP);
	PROTECT(debugSEXP);
	
	double * X = REAL(XSEXP);
	double * levels = REAL(levelsSEXP);
	int * successes = INTEGER(successesSEXP);
	int * trials = INTEGER(trialsSEXP);
	int n = *INTEGER(nSEXP);
	int num_trials = *INTEGER(num_trialsSEXP);
	int S = *INTEGER(SSEXP);
	int debug = *INTEGER(debugSEXP);

	if(debug){
		Rprintf("Debug = %d\n", debug);
	}

	double * a = (double*) R_Calloc(n, double);
  	double * b = (double*) R_Calloc(n, double);

  	int succ=0;
	int trial=0;
	int totcount = 0;
    
    for (int j = 0; j < num_trials; ++j)
    {
    	succ = successes[j];
    	trial = trials[j];
    	//printf("succ = %d, trial = %d\n", succ, trial);

    	for (int i = 0; i < succ; ++i)
    	{
    		b[totcount] = levels[j];
    		a[totcount++] = -DBL_MAX;
    		//printf("%d\n", totcount);
    	}

    	for (int i = succ; i < trial; ++i)
    	{
    		a[totcount] = levels[j];
    		b[totcount++] = DBL_MAX;
    		//printf("%d\n", totcount);
    	}


    }

    R_qsort(a, 1, n);
	R_qsort(b, 1, n);

	if(debug){
    	Rprintf("a:\n");
    	print_float_vector(n, a);
    	Rprintf("b:\n");
    	print_float_vector(n, b);
    }
	

	

	SEXP logpermsSEXP = PROTECT(allocVector(REALSXP, S));
	double * logperms = REAL(logpermsSEXP);
	memset(logperms, 0, sizeof(double)*S);


	SEXP a_union_bSEXP = PROTECT(allocVector(REALSXP, 2*n));
	double * a_union_b = REAL(a_union_bSEXP);
	memset(a_union_b, 0, sizeof(double)*2*n);

	int len_a_union_b =0;

	get_union(n, a, b, &len_a_union_b, a_union_b);

	SEXP alphaSEXP = PROTECT(allocVector(INTSXP, n));
	SEXP betaSEXP = PROTECT(allocVector(INTSXP, n));
	SEXP gammaSEXP = PROTECT(allocVector(INTSXP, n));
	SEXP logfactorialsSEXP = PROTECT(allocVector(REALSXP, n+1));
	SEXP mSEXP = PROTECT(allocVector(INTSXP, 1));
	SEXP k_SEXP = PROTECT(allocVector(INTSXP, 1));

	
	int * alpha = INTEGER(alphaSEXP);
	int * beta = INTEGER(betaSEXP);
	int * gamma = INTEGER(gammaSEXP);
	double * log_factorials = REAL(logfactorialsSEXP);
	int * m = INTEGER(mSEXP);
	int * k = INTEGER(k_SEXP);


	dictionary * new_log_subperms = init_dictionary(n);
	dictionary * old_log_subperms = init_dictionary(n);

	
	memset(alpha, 0, sizeof(int)*n);
	memset(beta, 0, sizeof(int)*n);
	memset(gamma, 0, sizeof(int)*n);
	memset(log_factorials, 0, sizeof(double)*(n+1));
	memset(m, 0, sizeof(int));
	memset(k, 0, sizeof(int));

	log_factorials[0]=0.0;
	for (int i = 1; i <= n; ++i)
	{
		log_factorials[i] = log_factorials[i-1] +log((double)(i));
	}

	
	
	SEXP historySEXP = PROTECT(allocVector(INTSXP, 3*n));
	SEXP amount_historySEXP = PROTECT(allocVector(INTSXP, 6*n));

	int * history = INTEGER(historySEXP);
	int * amount_history = INTEGER(amount_historySEXP);

	memset(history, 0, sizeof(int)*3*n);
	memset(amount_history, 0, sizeof(int)*6*n);


	for (int s = 0; s < S; ++s)
	{
		
		double * x = &(X[s*n]);
		R_qsort(x, 1, n);
		
		memset(alpha, 0, sizeof(int)*n);
		memset(beta, 0, sizeof(int)*n);
		memset(gamma, 0, sizeof(int)*n);
		memset(m, 0, sizeof(int));
		memset(k, 0, sizeof(int));


	    if(!nonzero_perm(x, a,  b, n)){
			logperms[s] = NA_REAL;
			continue;
		}
		get_alphabetagamma(x, n, a, b, a_union_b, len_a_union_b, alpha, 
	    beta, gamma,  k, m, debug);
		if((*k)==1){
			logperms[s] = log_factorials[n];
			continue;
		}
		
		
		/*get_alphabetagamma(x, n, a, b, a_union_b, len_a_union_b, alpha, 
	    beta, gamma,  k, m, debug);*/


	    if(debug){
	    	Rprintf("S=%d, t=%d\n", S, s);
	    	Rprintf("len_a_union_b = %d\n", len_a_union_b);
	    	Rprintf("x:\n");
	    	print_float_vector(n,x);
	    	Rprintf("a:\n");
	    	print_float_vector(n,a);
	    	Rprintf("b:\n");
	    	print_float_vector(n,b);
	    	Rprintf("a_union_b:\n");
	    	print_float_vector(2*n,a_union_b);
	    	Rprintf("len a_union_b:%d\n", len_a_union_b);
	    	Rprintf("alpha:\n");
	    	print_int_vector(n,  alpha);
	    	Rprintf("beta:\n");
	    	print_int_vector(n,  beta);
	    	Rprintf("gamma:\n");
	    	print_int_vector(n,  gamma);
	    	Rprintf("m:%d\n", *m);
	    	Rprintf("k:%d\n", *k);
	    	
	    }

		int history_len = 0;

	
		memset(history, 0, sizeof(int)*3*n);
		memset(amount_history, 0, sizeof(int)*6*n);


		if(debug){
			Rprintf("REDUCING NOW\n");
		}
		
		//return(XSEXP);
		int result = reduction(alpha,  beta,  gamma, m, n, k, history,
				   amount_history, &history_len, debug);

		if(result != 0){

			

			memset(alpha, 0, sizeof(int)*n);
			memset(beta, 0, sizeof(int)*n);
			memset(gamma, 0, sizeof(int)*n);
			memset(m, 0, sizeof(int));
			memset(k, 0, sizeof(int));
			debug = 1;
			get_alphabetagamma(x, n, a, b, a_union_b, len_a_union_b, alpha, 
		    beta, gamma,  k, m, debug);


		    if(debug){
		    	Rprintf("len_a_union_b = %d\n", len_a_union_b);
		    	Rprintf("x:\n");
		    	print_float_vector(n,x);
		    	Rprintf("a:\n");
		    	print_float_vector(n,a);
		    	Rprintf("b:\n");
		    	print_float_vector(n,b);
		    	Rprintf("a_union_b:\n");
		    	print_float_vector(2*n,a_union_b);
		    	Rprintf("len a_union_b:%d\n", len_a_union_b);
		    	Rprintf("alpha:\n");
		    	print_int_vector(n,  alpha);
		    	Rprintf("beta:\n");
		    	print_int_vector(n,  beta);
		    	Rprintf("gamma:\n");
		    	print_int_vector(n,  gamma);
		    	Rprintf("m:%d\n", *m);
		    	Rprintf("k:%d\n", *k);
		    	
		    }

			int history_len = 0;

		
			memset(history, 0, sizeof(int)*3*n);
			memset(amount_history, 0, sizeof(int)*6*n);


			if(debug){
				Rprintf("REDUCING NOW\n");
			}
			
			//return(XSEXP);
			reduction(alpha,  beta,  gamma, m, n, k, history,
					   amount_history, &history_len, debug);


			free_dictionary(new_log_subperms);
			free_dictionary(old_log_subperms);

			error("Error in computing permanent. See output from code for debugging.");
			UNPROTECT(18);
			return NULL;
		}

		if(debug){
			Rprintf("history len = %d\n", history_len);

			Rprintf("REDUCED SUBPERMS\n");
		}
		sparse_get_reduced_log_subperms( new_log_subperms,  alpha, beta, gamma,
						log_factorials, n,  m, k);
		//Rprintf("RESULT:\n");
		//print_matrix(n+1, n+1, new_log_subperms);

		dictionary * tmp  = old_log_subperms;
		old_log_subperms = new_log_subperms;
		new_log_subperms = tmp;



		if(debug){
			Rprintf("==========\nReverse reduction:\n==========\n");
		}
		//Rprintf("old = %d\n", old_log_subperms);
		//Rprintf("new = %d\n", new_log_subperms);
		dictionary * the_log_subperms = sparse_reverse_reduction(old_log_subperms, new_log_subperms, alpha,
						   beta,  gamma, m,  n, k,  history,
				           amount_history, &history_len, log_factorials);

		


		
		double logperm =  Csparse_log_sum_exp(the_log_subperms);
		logperms[s] = logperm;
		if(debug){
			Rprintf("logperm = %f\n", logperm);

		}



	}
	free_dictionary(new_log_subperms);
	free_dictionary(old_log_subperms);

	R_Free(a);
	R_Free(b);
	UNPROTECT(18);
	return(logpermsSEXP);

}


int nonzero_perm(double * x, double * a, double * b, int n){

	for (int i = 0; i < n; ++i)
	{
		if(x[i]< a[i] || x[i] > b[i]){
			return 0;
		}
	}
	return 1;



}