

vector<vector<int> > updatehiddenstates_unequalvar(int K, int T, int *lengths, double **Pimat, double **mu ,
						double **sigma2, double **P, double **logbeta, double *temp,
						HMM *phmm_temp, HMM *phmm_old, vector<vector<double> > data, 
						vector<vector<int> > z, int i)
{
	int t, k, j;
	for(t=0; t<T; t++){
		phmm_old->N = K;
		phmm_old->M = lengths[t];
		for(k=1; k<=K; k++){
			phmm_old->pi[k] = Pimat[t][k];
			if(phmm_old->pi[k]>0.5)
				z[t][0] = k;
			phmm_old->mu[k] = mu[t][k];
			phmm_old->sigma2[k] = sigma2[t][k];
			for(j=1; j<=K; j++)
				phmm_old->A[k][j] = P[k][j];
		}
		calculateBmatrix(phmm_old, &data[t][0]);
		samplehiddenstates(phmm_old, &data[t][0], &z[t][0], t+i, temp, logbeta, phmm_temp);
	}
	return(z);
}

vector<vector<int> > updateparas_unequalvar(int seed, int K, int T, double **P, double *mu0, double *sigma20,
				 double **mu, double **sigma2, vector<vector<double> > data, 
				 vector<vector<int> > z, int *lengths, double **count, double *sum,
				 int *num, double *tempA, double s2, double nu, gsl_vector *vv, 
				 gsl_permutation *perm, gsl_permutation *rank)
{
	gsl_rng * rng = gsl_rng_alloc (gsl_rng_taus);
    gsl_rng_set(rng, seed);

	// update transition matrix

	int i, j, t, n, k;

	for(i=1; i<=K; i++)
		for(j=1; j<=K; j++)
			count[i-1][j-1] = 1.0;
	for(t=0; t<T; t++){
		for(n=1; n<lengths[t]; n++){
			count[z[t][n-1]-1][z[t][n]-1] += 1.0;
		}
	}

	for(i=1; i<=K; i++){
		gsl_ran_dirichlet (rng, K, count[i-1] , tempA);
		for(j=1; j<=K; j++){
			P[i][j] = tempA[j-1];
		}
	}

	// update grand mean mu0

	for(k=1; k<=K; k++){
		mu0[k] = 0.0;
		for(t=0; t<T; t++){
			mu0[k] += mu[t][k];
		}
		mu0[k] = mu0[k]/T +  gsl_ran_gaussian(rng, sqrt(sigma20[k]/T));
	}

	// update corresponding variance sigma20

	for(k=1; k<=K; k++){
		sigma20[k] = 0.0;
		for(t=0; t<T; t++){
			sigma20[k] += pow((mu[t][k]-mu0[k]),2.0);
		}
		sigma20[k] = sigma20[k]/gsl_ran_chisq(rng, 1.0*(T-2));
	}

	// update means for individual traces, mu

	for(t=0; t<T; t++){
		for(k=1; k<=K; k++){
			sum[k] = 0.0;
			num[k] = 0;
		}
		for(n=0; n<lengths[t]; n++){
			sum[z[t][n]] += data[t][n];
			num[z[t][n]] += 1;
		}
		for(k=1; k<=K; k++){
		mu[t][k] = (mu0[k]*sigma2[t][k] + sigma20[k] * sum[k])/(sigma2[t][k] + sigma20[k] * num[k]) + 
			sqrt((sigma20[k] * sigma2[t][k])/(sigma20[k] * num[k] + sigma2[t][k]))* gsl_ran_gaussian(rng, 1.0);

		}

	for(k=0; k<K; k++){
		gsl_vector_set(vv,k,mu[t][k+1]);
	}
    gsl_sort_vector_index (perm, vv);
    gsl_permutation_inverse (rank, perm);
	for(k=0; k<K; k++){
		if(rank->data[k] != k){
		for(n=0; n<lengths[t]; n++){
			if(z[t][n]==k+1){
				z[t][n] = rank->data[k]+1;
			}
		}
		}
	}
	sort(mu[t]+1,mu[t]+1+K);
	}

	// update variances for individual traces, sigma2

	for(t=0; t<T; t++){
		for(k=1; k<=K; k++){
		sum[k] = s2 * nu;
		num[k] = 1;
		}
		for(n=0; n< lengths[t]; n++){
			sum[z[t][n]] += pow((data[t][n] - mu[t][z[t][n]]), 2.0);
			num[z[t][n]] += 1;
		}
		for(k=1; k<=K; k++){
			sigma2[t][k] = sum[k]/gsl_ran_chisq(rng, 1.0*(num[k] + nu));
		}
	}
	return(z);
}

void hiergibbsunequalvar(double **P, double **mu, double **sigma2, int K, int T,
			   vector<vector<double> > data, vector<vector<int> > z, int seed,
			   double **count, int *lengths, double *mu0, double *sigma20,
			   double *sum, int *num, double *tempA, double s2, double nu,
			   double **Pimat, double **logbeta, HMM *phmm_temp, HMM *phmm_old,
			   double **logalpha, double **loggamma, double ***logxi,
			   double *plogprobinit, double *plogprobfinal, int num_gibbs,
			   int burnin, int thin, double ***P_gibbs, double ***mu_gibbs,
			   double ***sigma2_gibbs, double **mu0_gibbs, double **sigma20_gibbs)
{
	int t, k, j, i;

	gsl_vector *vv;
	vv = gsl_vector_alloc(K);
	gsl_permutation * perm = gsl_permutation_alloc(K);
    gsl_permutation * rank = gsl_permutation_alloc(K);
     

	// fit individual EMs

	for(t=0; t<T; t++){
		initializeHMM_num(phmm_temp, lengths[t], K);
		calculateBmatrix(phmm_temp, &data[t][0]);
		BaumWelch(phmm_temp, &data[t][0], lengths[t], logalpha, logbeta, loggamma, 
			     plogprobinit, plogprobfinal, logxi, tempA, vv, perm, rank, nu, &s2);

		for(k=1; k<=K; k++){
			mu[t][k] = phmm_temp->mu[k];
			Pimat[t][k] = phmm_temp->pi[k];
			sigma2[t][k] = phmm_temp->sigma2[k];
		}
		z[t][0] = 1;
		for(k=1; k<=K; k++){
			if(phmm_temp->pi[k]>0.5)
				z[t][0] = k;
		}
	}

	for(k=1; k<=K; k++){
		for(j=1; j<=K; j++)
			P[k][j] = 1.0/K;
		sigma20[k] = 0.01;
	}

	int Ngibbs = burnin + (num_gibbs - 1) * thin + 1;
	i = 0;
	j = 0;

	while(i < Ngibbs){
	
		z = updatehiddenstates_unequalvar(K, T, lengths, Pimat, mu, sigma2, P, logbeta, 
		                   tempA, phmm_temp, phmm_old, data, z, i);
	
		z = updateparas_unequalvar(seed+i, K, T, P, mu0, sigma20, mu, sigma2, data, z, lengths, 
		        count, sum, num, tempA, s2, nu, vv, perm, rank);

		if(i >= burnin){
			if ((i-burnin)%thin == 0){
				for(k=1; k<=K; k++){
				mu0_gibbs[j][k] = mu0[k];
				sigma20_gibbs[j][k] = sigma20[k];
				for(t=1; t<=K; t++)
					P_gibbs[j][k][t] = P[k][t];
				for(t=0; t<T; t++){
					mu_gibbs[j][t][k] = mu[t][k];
					sigma2_gibbs[j][t][k] = sigma2[t][k];
				}
				}
				j++;
			}
		}
		i++;
	}
	gsl_vector_free(vv);
	gsl_permutation_free(perm);
	gsl_permutation_free(rank);
}

void hierEM_uneqvar(int *HQSlength, int T, vector<vector<double> > HQS, double *mu0, double *sigma20, 
			double *** logbeta, double *** logalpha, double *** loggamma, double **** logxi, 
			double * sum_mu, double ** sum_P, double * row_P, double nu, double s2, int K, 
			HMM **hmm, const char *filename_Pmu0sigam20, const char *filename_indtraces)
{
	
	int i,j,k,t;
	double sum_numerator, sum_denom;
	double epsilon = 10e-4;
	double LLnew, LLold;

	LLnew = -1000000000.0;

	for(i=1; i<=K; i++)
		sigma20[i] = 0.01;
	
	for(i=0; i < T; i++){
		initializeHMM_num(hmm[i], HQSlength[i], K);
	}

	for(k=1; k<=K; k++){
		mu0[k] = hmm[0]->mu[k];
	}

	do{
		LLold = LLnew;
		LLnew = 0.0;

/*************************** E Step ****************************/

	for(i=0; i < T; i++){
        calculateBmatrix(hmm[i], &HQS[i][0]);
		Backwardlog(hmm[i], HQSlength[i], logbeta[i], row_P);
	    Forwardlog(hmm[i], HQSlength[i], logalpha[i], row_P);
        ComputeGamma(hmm[i], HQSlength[i], logalpha[i], logbeta[i], loggamma[i]);
	    ComputeXi(hmm[i], HQSlength[i], logalpha[i], logbeta[i], logxi[i]);
	    LLnew += hmm[i]->logprobf[1];
	}

	printf("%f\n",LLnew);
/*********************************** M Step ************************************/

	for(i=1; i<= hmm[1]->N; i++)
		sum_mu[i]=0.0;

	for(i=0; i<T; i++){

		/* Update pi, mu, sigma^2 */

		for(j=1; j<= hmm[i]->N; j++){
		    hmm[i]->pi[j] = exp(loggamma[i][1][j]);
		    sum_numerator=0.0;
		    sum_denom=0.0;
			for(k=1; k<= HQSlength[i]; k++){
			   sum_denom += exp(loggamma[i][k][j]);
			   sum_numerator += exp(loggamma[i][k][j])* HQS[i][k-1];
			}
			hmm[i]->mu[j]=(sum_numerator*sigma20[j] + hmm[i]->sigma2[j]*mu0[j])/
				         (sum_denom*sigma20[j] + hmm[i]->sigma2[j]);
		}
		sort(hmm[i]->mu+1, hmm[i]->mu+1+hmm[i]->N);

		for(j=1; j<=hmm[i]->N; j++){
			sum_numerator = nu * s2;
			sum_denom = nu;
			for(k=1; k<= HQSlength[i]; k++){			  
				sum_denom += exp(loggamma[i][k][j]);
    			sum_numerator += pow((HQS[i][k-1]-hmm[i]->mu[j]),2)* exp(loggamma[i][k][j]);
			}
			hmm[i]->sigma2[j] = sum_numerator/sum_denom;
		}

	}

	for(j=1; j<= hmm[1]->N; j++)
	for(i=0; i<T; i++)
		sum_mu[j] += hmm[i]->mu[j];

        /* Update mu0, sigma0^2 */

	for(j=1; j<= hmm[1]->N; j++){
		mu0[j] = sum_mu[j]/T;
		printf("%f\t",mu0[j]);
	}

	for(j=1; j<= hmm[1]->N; j++){
		sum_numerator = 0.0;
		for(i=0; i< T; i++)
			sum_numerator += pow(hmm[i]->mu[j]-mu0[j],2.0);
		sigma20[j] = sum_numerator/T;
	}  

		/* Update Transition Matrix */

		for(i=1; i<= hmm[1]->N; i++){
			row_P[i]=0.0;
			for(j=1; j<= hmm[1]->N; j++){
				sum_P[i][j]=0;
				for(k=0; k<T; k++){
					for(t=1; t< HQSlength[k]; t++)
						sum_P[i][j]+= exp(logxi[k][t][i][j]);
				}
				row_P[i] += sum_P[i][j];
			}
			for(j=1; j<= hmm[1]->N; j++)
			hmm[0]->A[i][j]= sum_P[i][j]/row_P[i];
		}

		for(t=1; t<T; t++)
			hmm[t]->A = hmm[0]->A;  
	} while(LLnew-LLold > epsilon); 

	printf("EMfinished\n");

	// output to file

	 ofstream fout;
     fout.open(filename_Pmu0sigam20,ios::app);    
	 // open file for appending mu0, sigma20, P (global parameters)
     assert (!fout.fail( )); 
	 for(j=1; j<=K; j++)
		 fout<<mu0[j]<<"\t";
	 fout<<"\n";
 	 for(j=1; j<=K; j++)
		 fout<<sigma20[j]<<"\t";
	 fout<<"\n";
	 for(j=1; j<=K; j++){
		 for(i=1; i<=K; i++){
			 fout<<hmm[0]->A[j][i]<<"\t";
		 }
		 fout<<"\n";
	 }
	 fout.close( );       //close file
     assert(!fout.fail( ));

     fout.open(filename_indtraces,ios::app);    
	 // open file for appending pi, mu, sigma2 for individual traces
     assert (!fout.fail( )); 
	 for(i=0; i<T; i++){
		 for(j=1; j<=K; j++)
			 fout<<hmm[i]->pi[j]<<"\t";
		 for(j=1; j<=K; j++)
			 fout<<hmm[i]->mu[j]<<"\t";
		 for(j=1; j<=K; j++)
			 fout<<hmm[i]->sigma2[j]<<"\t";
		 fout<<"\n";
	 }
	 fout.close( );       //close file
     assert(!fout.fail( ));
}

