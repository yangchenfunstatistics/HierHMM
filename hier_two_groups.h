
void hier_two_groups_shortcut(vector<vector<double> > data, vector<vector<int> > z, int *state_index, int K,
					 double **mu, double **sigma2, double *mu0, double *sigma20, double **mu0_gibbs,
					 double ***P_gibbs, double **sigma20_gibbs, int T, HMM *phmm_temp, int *lengths,
					 double **logalpha, double **logbeta, double **loggamma, double *plogprobinit,
					 double *plogprobfinal, double ***logxi, double *tempA, double nu, double s2,
					 double **P, int seed, int num_gibbs, int burnin, int thin, double **Pimat,
					 HMM *phmm_old, double **count)
{	
	gsl_rng * rng = gsl_rng_alloc (gsl_rng_taus);
    gsl_rng_set(rng,seed);

    // initialize the parameters mu0,sigma20 after hierEM_uneqvar() function
	// initialize mu, sigma2 by doing Baum-Welch for each individual trace

	int t, k, j, n;
		gsl_vector *vv = gsl_vector_alloc(K);
	    gsl_permutation * perm = gsl_permutation_alloc(K);
        gsl_permutation * rank = gsl_permutation_alloc(K);
		gsl_vector *vv2 = gsl_vector_alloc(K-1);
	    gsl_permutation * perm2 = gsl_permutation_alloc(K-1);
        gsl_permutation * rank2 = gsl_permutation_alloc(K-1);

	// initialize all the paremeters and fix state indicator of first observation

	for(t=0; t<T; t++){
		initializeHMM_num(phmm_temp, lengths[t], state_index[t]);
		calculateBmatrix(phmm_temp, &data[t][0]);
//		printf("%d\n",state_index[t]);
		if(state_index[t]==K)
		BaumWelch(phmm_temp, &data[t][0], lengths[t], logalpha, logbeta, loggamma, 
			     plogprobinit, plogprobfinal, logxi, tempA, vv, perm, rank, nu, &s2);
		else if(state_index[t]==K-1)
		BaumWelch(phmm_temp, &data[t][0], lengths[t], logalpha, logbeta, loggamma, 
			     plogprobinit, plogprobfinal, logxi, tempA, vv2, perm2, rank2, nu, &s2);
		for(k=1; k<=state_index[t]; k++){
			mu[t][k] = phmm_temp->mu[k];
//			printf("%d\t%d\t%f\n",t,k,mu[t][k]);
			sigma2[t][k] = phmm_temp->sigma2[k];
//			printf("%d\t%d\t%f\n",t,k,sigma2[t][k]);
			// fix z[0] with value of pi from EM
			Pimat[t][k] = phmm_temp->pi[k];
			if(phmm_temp->pi[k]>0.5)
				z[t][0] = k;
		}
//		printf("%d\t",z[t][0]);
	}
	for(t=1; t <= K; t++){
		for(k=1; k <= K; k++)
			P[t][k] = 1.0/K;
		sigma20[t] = 0.001;
	}

	int gibbs_count = 0;
	int gibbs_save = 0;
	double sum;
	int num;
	int Ngibbs = burnin + (num_gibbs - 1) * thin + 1;

	printf("Before Gibbs\n");

	while(gibbs_count < Ngibbs){

		printf("gibbs%d\n",gibbs_count);
		// update mu0

		for(k=1; k <= K; k++){
			sum = 0.0;
			num = 0;
			if(k == 2){
				for(t=0; t < T; t++){
					if(state_index[t]==K){
//						printf("%f\t",mu[t][k]);
						sum += mu[t][k];
						num ++;
					}
				}
				mu0[k] = sum/num + gsl_ran_gaussian (rng, sqrt(sigma20[k]/num));
//				printf("mu0%d\t%f\n",k,mu0[k]);
			}
			else{
				for(t=0; t< T; t++){
					if(state_index[t]==K)
						sum += mu[t][k];
					else if(state_index[t]==K-1)
						sum += mu[t][(k+1)/2];
					else{
						printf("state index error\n");
						exit(0);
					}
				}
				mu0[k] = sum/T + gsl_ran_gaussian (rng, sqrt(sigma20[k]/T));
//				printf("mu0%d\t%f\n",k,mu0[k]);
			}
		}
		// update sigma20
		for(k=1; k <= K; k++){
			sum = 0.0;
			num = 0;
			for(t=0; t < T; t++){
				if(k==2){
					if(state_index[t]==K){
						sum += pow(mu[t][k]-mu0[k],2.0);
						num ++;
					}
				}
				else{
					if(state_index[t]==K)
						sum += pow(mu[t][k]-mu0[k],2.0);
					else if(state_index[t]==K-1)
						sum += pow(mu[t][(k+1)/2]-mu0[k],2.0);
					else{
						printf("state index error\n");
						exit(0);
					}
					num ++;
				}
			}
			sigma20[k] = sum / gsl_ran_chisq(rng, 1.0 * (num - 2));
//			printf("sigma20%d\t%f\n",k,sigma20[k]);
		}
		// update z
			
		for(t=0; t<T; t++){
		phmm_old->N = state_index[t];
		phmm_old->M = lengths[t];
		for(k=1; k<=state_index[t]; k++){
			phmm_old->pi[k] = Pimat[t][k];
			phmm_old->mu[k] = mu[t][k];
			phmm_old->sigma2[k] = sigma2[t][k];
		}
		if(state_index[t]==K){
			for(k=1; k<=state_index[t]; k++)
				for(j=1; j<=state_index[t]; j++)
					phmm_old->A[k][j] = P[k][j];
		}
		else if(state_index[t]==K-1){
			for(k=1; k<=state_index[t]; k++){
				sum = 0.0;
				for(j=1; j<=state_index[t]; j++)
					sum += P[2*k-1][2*j-1];
				for(j=1; j<=state_index[t]; j++)
					phmm_old->A[k][j] = P[2*k-1][2*j-1]/sum;
			}
		}
		else{
			printf("state index error\n");
			exit(0);
		}
		calculateBmatrix(phmm_old, &data[t][0]);
		samplehiddenstates(phmm_old, &data[t][0], &z[t][0], t+gibbs_count, tempA, logbeta, phmm_temp);
//		printf("%d%d%d%d%d\n",z[t][0],z[t][1],z[t][2],z[t][3],z[t][4]);
		}
		// update parameters for individual traces
		// update mu
		for(t=0; t<T; t++){
			if(state_index[t]==K){
				for(k=1; k<=K; k++){
					sum = 0.0;
					num = 0;
					for(j=0; j<lengths[t]; j++){
						if(z[t][j]==k){
							sum += data[t][j];
							num ++;
						}
					}
					mu[t][k] = (mu0[k]/sigma20[k] + sum/sigma2[t][k])/(1.0/sigma20[k]+1.0*num/sigma2[t][k]) + 
						gsl_ran_gaussian (rng, sqrt(1.0/(1.0/sigma20[k]+1.0*num/sigma2[t][k])));
//					printf("mu%d\t%d\t%f\n",t,k,mu[t][k]);
				}
			}
			else if(state_index[t]==K-1){
				for(k=1; k<K; k++){
					sum = 0.0;
					num = 0;
					for(j=0; j<lengths[t]; j++){
						if(z[t][j]==k){
							sum += data[t][j];
							num ++;
						}
					}
					mu[t][k] = (mu0[2*k-1]/sigma20[2*k-1] + sum/sigma2[t][k])/(1.0/sigma20[2*k-1]+1.0*num/sigma2[t][k]) + 
						gsl_ran_gaussian (rng, sqrt(1.0/(1.0/sigma20[2*k-1]+1.0*num/sigma2[t][k])));
//					printf("mu%d\t%d\t%f\n",t,k,mu[t][k]);
				}
			}
			else{
				printf("state index error\n");
				exit(0);
			}
			// sort mu[t]
			if(state_index[t]==K){	
				for(k=0; k<K; k++)
					gsl_vector_set(vv,k,mu[t][k+1]);
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
//				for(k=1; k<=K; k++)
//				printf("mu%d\t%d\t%f\n",t,k,mu[t][k]);
			}
			else if(state_index[t]==K-1){	
				for(k=0; k<K-1; k++)
					gsl_vector_set(vv2,k,mu[t][k+1]);
				gsl_sort_vector_index (perm2, vv2);
				gsl_permutation_inverse (rank2, perm2);
				for(k=0; k<K-1; k++){
					if(rank2->data[k] != k){
						for(n=0; n<lengths[t]; n++){
							if(z[t][n]==k+1){
								z[t][n] = rank2->data[k]+1;
							}
						}
					}
				}
				sort(mu[t]+1,mu[t]+K);
//				for(k=1; k<K; k++)
//				printf("mu%d\t%d\t%f\n",t,k,mu[t][k]);
			}
		}
		// update sigma2
		double num2;
		for(t=0; t<T; t++){
			for(k=1; k<=state_index[t]; k++){
			num2 = nu;
			sum = s2 * nu;
				for(j=0; j<lengths[t]; j++){
					if(z[t][j]==k){
						num2 += 1.0;
						sum += pow(data[t][j]-mu[t][k],2.0);
					}
				}
				sigma2[t][k] = sum/gsl_ran_chisq(rng,num2);
//				printf("sig%d\t%d\t%f\n",t,k,sigma2[t][k]);
			}
		}

		// update transition matrix
		
	for(k=1; k<=K; k++)
		for(j=1; j<=K; j++)
			count[k-1][j-1] = 1.0;
	for(t=0; t<T; t++){
		for(n=1; n<lengths[t]; n++){
			count[z[t][n-1]-1][z[t][n]-1] += 1.0;
		}
	}

	for(k=1; k<=K; k++){
		gsl_ran_dirichlet (rng, K, count[k-1] , tempA);
		for(j=1; j<=K; j++){
			P[k][j] = tempA[j-1];
		}
	}

		// update state_index
		for(t=0; t<T; t++){
			phmm_temp->N = state_index[t];
			phmm_old->N = 5-state_index[t];
			phmm_temp->M = lengths[t];
			phmm_old->M = lengths[t];
			for(k=1; k<=state_index[t]; k++){
				phmm_temp->pi[k] = Pimat[t][k];
				phmm_temp->mu[k] = mu[t][k];
				phmm_temp->sigma2[k] = sigma2[t][k];
			}
			if(state_index[t]==K){
				for(k=1; k<=K; k++)
					for(j=1; j<=K; j++)
						phmm_temp->A[k][j] = P[k][j];
				for(k=1; k<=K-1; k++){
					phmm_old->pi[k] = Pimat[t][2*k-1];
					phmm_old->mu[k] = mu[t][2*k-1];
					phmm_old->sigma2[k] = sigma2[t][2*k-1];
					sum = 0.0;
					for(j=1; j<=K-1; j++)
						sum += P[2*k-1][2*j-1];
					for(j=1; j<=K-1; j++)
						phmm_old->A[k][j] = P[2*k-1][2*j-1]/sum;
				}
			}
			else if(state_index[t]==K-1){
				for(k=1; k<=K-1; k++){
					sum = 0.0;
					for(j=1; j<=K-1; j++)
						sum += P[2*k-1][2*j-1];
					for(j=1; j<=K-1; j++)
						phmm_temp->A[k][j] = P[2*k-1][2*j-1]/sum;
					phmm_old->pi[2*k-1] = Pimat[t][k];
					phmm_old->mu[2*k-1] = mu[t][k];
					phmm_old->sigma2[2*k-1] = sigma2[t][k];
				}
				for(k=1; k<=K; k++)
					for(j=1; j<=K; j++)
						phmm_old->A[k][j] = P[k][j];
				phmm_old->pi[2] = 0.0;					
				phmm_old->mu[2] = mu0[2] + gsl_ran_gaussian (rng, sqrt(sigma20[2]));
				phmm_old->sigma2[2] = nu*s2/gsl_ran_chisq(rng,nu);
			}
			else{
				printf("state index error\n");
				exit(0);
			}
			// calculate marginal likelihood by doing forward step
			calculateBmatrix(phmm_temp, &data[t][0]);
			Forwardlog(phmm_temp, lengths[t], logalpha, tempA);
			calculateBmatrix(phmm_old, &data[t][0]);
			Forwardlog(phmm_old, lengths[t], logalpha, tempA);
			double diff = - phmm_temp->logprobf[1] + phmm_old->logprobf[1];
			sum = 1.0/(1.0 + exp(diff)); //probability of staying in original one
			num2 = gsl_ran_flat(rng, 0.0, 1.0);
			if(num2 > sum){
				if(state_index[t]==K){
					Pimat[t][2] = 1.0 - Pimat[t][1];
					mu[t][2] = mu[t][3];
					sigma2[t][2] = sigma2[t][3];
				}
				else if(state_index[t]==K-1){
					Pimat[t][3] = Pimat[t][2];
					Pimat[t][2] = 0.0;
					mu[t][3] = mu[t][2];
					mu[t][2] = phmm_old->mu[2];
					sigma2[t][3] = sigma2[t][2];
					sigma2[t][2] = phmm_old->sigma2[2];
				}
				state_index[t] = phmm_old->N;
				printf("%dchangeindex%d%d\n",t,phmm_temp->N,phmm_old->N);
			}
		}
		if (gibbs_count >= burnin){
			if ((gibbs_count-burnin)%thin == 0){
				for(k=1; k<=K; k++){
				mu0_gibbs[gibbs_save][k] = mu0[k];
				printf("%f\t",mu0[k]);
				sigma20_gibbs[gibbs_save][k] = sigma20[k];
				printf("%f\t",sigma20[k]);
				for(j=1; j<=K;j++){
					P_gibbs[gibbs_save][k][j] = P[k][j];
					printf("%f\t",P[k][j]);
				}
				}
				gibbs_save ++;
				printf("\n");
			}
		}
		gibbs_count ++;
	}
}



void hier_two_groups_mid_indicator(vector<vector<double> > data, vector<vector<int> > z, int *state_index, int K,
					 double **mu, double **sigma2, double *mu0, double *sigma20, double **mu0_gibbs,
					 double ***P_gibbs, double **sigma20_gibbs, int T, HMM *phmm_temp, int *lengths,
					 double **logalpha, double **logbeta, double **loggamma, double *plogprobinit,
					 double *plogprobfinal, double ***logxi, double *tempA, double nu, double s2_init,
					 double **P, int seed, int num_gibbs, int burnin, int thin, double **Pimat,
					 HMM *phmm_old, double **count)
{	
	gsl_rng * rng = gsl_rng_alloc (gsl_rng_taus);
    gsl_rng_set(rng,seed);

    // initialize the parameters mu0,sigma20 after hierEM_uneqvar() function
	// initialize mu, sigma2 by doing Baum-Welch for each individual trace

	double *s2;
	s2 = new double[K+1];
	int t, k, j, n;
	for(k=1; k<=K; k++)
		s2[k] = s2_init;
	gsl_vector *vv = gsl_vector_alloc(K);
	gsl_permutation * perm = gsl_permutation_alloc(K);
    gsl_permutation * rank = gsl_permutation_alloc(K);
	gsl_vector *vv2 = gsl_vector_alloc(K-1);
	gsl_permutation * perm2 = gsl_permutation_alloc(K-1);
    gsl_permutation * rank2 = gsl_permutation_alloc(K-1);
	double LL_old, LL_new;

	// initialize all the paremeters and fix state indicator of first observation

	for(t=0; t<T; t++){
		initializeHMM_num(phmm_temp, lengths[t], state_index[t]);
		calculateBmatrix(phmm_temp, &data[t][0]);
//		printf("%d\n",state_index[t]);
		if(state_index[t]==K)
		BaumWelch(phmm_temp, &data[t][0], lengths[t], logalpha, logbeta, loggamma, 
			     plogprobinit, plogprobfinal, logxi, tempA, vv, perm, rank, nu, &s2_init);
		else if(state_index[t]==K-1)
		BaumWelch(phmm_temp, &data[t][0], lengths[t], logalpha, logbeta, loggamma, 
			     plogprobinit, plogprobfinal, logxi, tempA, vv2, perm2, rank2, nu, &s2_init);
		for(k=1; k<=state_index[t]; k++){
			mu[t][k] = phmm_temp->mu[k];
//			printf("%d\t%d\t%f\n",t,k,mu[t][k]);
			sigma2[t][k] = phmm_temp->sigma2[k];
//			printf("%d\t%d\t%f\n",t,k,sigma2[t][k]);
			// fix z[0] with value of pi from EM
			Pimat[t][k] = phmm_temp->pi[k];
			if(phmm_temp->pi[k]>0.5)
				z[t][0] = k;
		}
//		printf("%d\t",z[t][0]);
	}
	for(t=1; t <= K; t++){
		for(k=1; k <= K; k++)
			P[t][k] = 1.0/K;
		sigma20[t] = 0.001;
	}

	int gibbs_count = 0;
	int gibbs_save = 0;
	double sum;
	double diff;
	int num;
	int Ngibbs = burnin + (num_gibbs - 1) * thin + 1;

	printf("Before Gibbs\n");

	while(gibbs_count < Ngibbs){

		printf("gibbs%d\n",gibbs_count);
		// update mu0

		for(k=1; k <= K; k++){
			sum = 0.0;
			num = 0;
			if(k == 2){
				for(t=0; t < T; t++){
					if(state_index[t]==K){
//						printf("%f\t",mu[t][k]);
						sum += mu[t][k];
						num ++;
					}
				}
				mu0[k] = sum/num + gsl_ran_gaussian (rng, sqrt(sigma20[k]/num));
//				printf("mu0%d\t%f\n",k,mu0[k]);
			}
			else{
				for(t=0; t< T; t++){
					if(state_index[t]==K)
						sum += mu[t][k];
					else if(state_index[t]==K-1)
						sum += mu[t][(k+1)/2];
					else{
						printf("state index error\n");
						exit(0);
					}
				}
				mu0[k] = sum/T + gsl_ran_gaussian (rng, sqrt(sigma20[k]/T));
//				printf("mu0%d\t%f\n",k,mu0[k]);
			}
		}
		// update sigma20
		for(k=1; k <= K; k++){
			sum = 0.0;
			num = 0;
			for(t=0; t < T; t++){
				if(k==2){
					if(state_index[t]==K){
						sum += pow(mu[t][k]-mu0[k],2.0);
						num ++;
					}
				}
				else{
					if(state_index[t]==K)
						sum += pow(mu[t][k]-mu0[k],2.0);
					else if(state_index[t]==K-1)
						sum += pow(mu[t][(k+1)/2]-mu0[k],2.0);
					else{
						printf("state index error\n");
						exit(0);
					}
					num ++;
				}
			}
			sigma20[k] = sum / gsl_ran_chisq(rng, 1.0 * (num - 2));
//			printf("sigma20%d\t%f\n",k,sigma20[k]);
		}

		// update s2
		for(k=1; k<=K; k++){
			num = 2;
			sum = 0.0;
			if(k==2){
				for(t=0; t<T; t++){
					if(state_index[t]==K){
						num ++;
						sum += 1.0/sigma2[t][k];
//						printf("%f\t",sigma2[t][k]);
					}
				}
				s2[k] = gsl_ran_chisq(rng,nu*num)/(sum*nu);
//				printf("sum%f\ts2%d\t%f\n",sum,k,s2[k]);
			}
			else{
				for(t=0; t<T; t++){
					if(state_index[t]==K){
						sum += 1.0/sigma2[t][k];
//						printf("%f\t",sigma2[t][k]);
					}
					else{
						sum += 1.0/sigma2[t][(k+1)/2];
//						printf("%f\t",sigma2[t][(k+1)/2]);
					}
				}
				s2[k] = gsl_ran_chisq(rng,1.0*(T+2))/(sum*nu);
//				printf("sum%f\ts2%d\t%f\n",sum,k,s2[k]);
			}
		}

		// update z
			
		for(t=0; t<T; t++){
		phmm_old->N = state_index[t];
		phmm_old->M = lengths[t];
		for(k=1; k<=state_index[t]; k++){
			phmm_old->pi[k] = Pimat[t][k];
			phmm_old->mu[k] = mu[t][k];
			phmm_old->sigma2[k] = sigma2[t][k];
		}
		if(state_index[t]==K){
			for(k=1; k<=state_index[t]; k++)
				for(j=1; j<=state_index[t]; j++)
					phmm_old->A[k][j] = P[k][j];
		}
		else if(state_index[t]==K-1){
			for(k=1; k<=state_index[t]; k++){
				sum = 0.0;
				for(j=1; j<=state_index[t]; j++)
					sum += P[2*k-1][2*j-1];
				for(j=1; j<=state_index[t]; j++)
					phmm_old->A[k][j] = P[2*k-1][2*j-1]/sum;
			}
		}
		else{
			printf("state index error\n");
			exit(0);
		}
		calculateBmatrix(phmm_old, &data[t][0]);
		for(k=1; k<=state_index[t]; k++){
			if(phmm_old->pi[k]>0.5)
				z[t][0] = k;
		}
//		printf("%d%d%d%d%d\n",z[t][0],z[t][1],z[t][2],z[t][3],z[t][4]);
		samplehiddenstates(phmm_old, &data[t][0], &z[t][0], t+gibbs_count, tempA, logbeta, phmm_temp);
//		printf("%d%d%d%d%d\n",z[t][0],z[t][1],z[t][2],z[t][3],z[t][4]);
		}
		// update parameters for individual traces
		// update mu
		for(t=0; t<T; t++){
			if(state_index[t]==K){
				for(k=1; k<=K; k++){
					sum = 0.0;
					num = 0;
					for(j=0; j<lengths[t]; j++){
						if(z[t][j]==k){
							sum += data[t][j];
							num ++;
						}
					}
					mu[t][k] = (mu0[k]/sigma20[k] + sum/sigma2[t][k])/(1.0/sigma20[k]+1.0*num/sigma2[t][k]) + 
						gsl_ran_gaussian (rng, sqrt(1.0/(1.0/sigma20[k]+1.0*num/sigma2[t][k])));
//					printf("mu%d\t%d\t%f\n",t,k,mu[t][k]);
				}
			}
			else if(state_index[t]==K-1){
				for(k=1; k<K; k++){
					sum = 0.0;
					num = 0;
					for(j=0; j<lengths[t]; j++){
						if(z[t][j]==k){
							sum += data[t][j];
							num ++;
						}
					}
					mu[t][k] = (mu0[2*k-1]/sigma20[2*k-1] + sum/sigma2[t][k])/(1.0/sigma20[2*k-1]+1.0*num/sigma2[t][k]) + 
						gsl_ran_gaussian (rng, sqrt(1.0/(1.0/sigma20[2*k-1]+1.0*num/sigma2[t][k])));
//					printf("mu%d\t%d\t%f\n",t,k,mu[t][k]);
				}
			}
			else{
				printf("state index error\n");
				exit(0);
			}
			// sort mu[t]
			if(state_index[t]==K){	
				for(k=0; k<K; k++)
					gsl_vector_set(vv,k,mu[t][k+1]);
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
				for(k=1; k<=K; k++)
				printf("mu%d\t%f\t",t,mu[t][k]);
				printf("\n");
			}
			else if(state_index[t]==K-1){	
				for(k=0; k<K-1; k++)
					gsl_vector_set(vv2,k,mu[t][k+1]);
				gsl_sort_vector_index (perm2, vv2);
				gsl_permutation_inverse (rank2, perm2);
				for(k=0; k<K-1; k++){
					if(rank2->data[k] != k){
						for(n=0; n<lengths[t]; n++){
							if(z[t][n]==k+1){
								z[t][n] = rank2->data[k]+1;
							}
						}
					}
				}
				sort(mu[t]+1,mu[t]+K);
				for(k=1; k<K; k++)
				printf("mu%d\t%f\t",t,mu[t][k]);
				printf("\n");
			}
		}
		// update sigma2
		double num2;
		for(t=0; t<T; t++){
			for(k=1; k<=state_index[t]; k++){
				num2 = nu;
				if(state_index[t]==K)
					sum = s2[k] * nu;
				else
					sum = s2[2*k-1] * nu;
				for(j=0; j<lengths[t]; j++){
					if(z[t][j]==k){
						num2 += 1.0;
						sum += pow(data[t][j]-mu[t][k],2.0);
					}
				}
				sigma2[t][k] = sum/gsl_ran_chisq(rng,num2);
//				printf("sig%d\t%d\t%f\n",t,k,sigma2[t][k]);
			}
		}

		// update transition matrix

	for(k=1; k<=K; k++)
		for(j=1; j<=K; j++)
			count[k-1][j-1] = 1.0;
	for(t=0; t<T; t++){
		for(n=1; n<lengths[t]; n++){
			count[z[t][n-1]-1][z[t][n]-1] += 1.0;
		}
	}

	for(k=1; k<=K; k++){
		gsl_ran_dirichlet (rng, K, count[k-1] , tempA);
		for(j=1; j<=K; j++){
			P[k][j] = tempA[j-1];
		}
	}

		// update state_index
		for(t=0; t<T; t++){
//			printf("%d\t%d\n",t,state_index[t]);
			if(state_index[t]==K){
				LL_old = log(gsl_ran_gaussian_pdf(data[t][0]-mu[t][z[t][0]],sqrt(sigma2[t][z[t][0]])));
				for(n=1; n<lengths[t]; n++){
					LL_old += log(gsl_ran_gaussian_pdf(data[t][n]-mu[t][z[t][n]],sqrt(sigma2[t][z[t][n]])));
					LL_old += log(P[z[t][n-1]][z[t][n]]);
				}
				phmm_old->N = K-1;
				phmm_old->M = lengths[t];
				phmm_old->pi[1] = Pimat[t][1];
				phmm_old->pi[2] = 1.0 - phmm_old->pi[1];
				for(j=1; j<K; j++){
					phmm_old->mu[j] = mu[t][2*j-1];
					phmm_old->sigma2[j] = sigma2[t][2*j-1];
					for(k=1; k<K; k++)
					phmm_old->A[j][k] = P[2*j-1][2*k-1]/(P[2*j-1][1]+P[2*j-1][3]);
				}
				calculateBmatrix(phmm_old, &data[t][0]);
				samplehiddenstates(phmm_old, &data[t][0], &z[t][0], t+gibbs_count, tempA, logbeta, phmm_temp);
				LL_new = log(gsl_ran_gaussian_pdf(data[t][0]-phmm_old->mu[z[t][0]],sqrt(phmm_old->sigma2[z[t][0]])));
				for(n=1; n<lengths[t]; n++){
					LL_new += log(gsl_ran_gaussian_pdf(data[t][n]-phmm_old->mu[z[t][n]],sqrt(phmm_old->sigma2[z[t][n]])));
					LL_new += log(phmm_old->A[z[t][n-1]][z[t][n]]);
				}
				for(k=1; k<K; k++){
					LL_new += log(gsl_ran_gaussian_pdf(phmm_old->mu[k]-mu0[2*k-1],sqrt(sigma20[2*k-1])));
					LL_new += log(gsl_ran_chisq_pdf(nu*s2[2*k-1]/phmm_old->sigma2[k],nu));
				}
//				printf("%f\t%f\n",LL_old,LL_new);
			}
			else if(state_index[t]=K-1){
				LL_old = log(gsl_ran_gaussian_pdf(data[t][0]-mu[t][z[t][0]],sqrt(sigma2[t][z[t][0]])));
				for(n=1; n<lengths[t]; n++){
					LL_old += log(gsl_ran_gaussian_pdf(data[t][n]-mu[t][z[t][n]],sqrt(sigma2[t][z[t][n]])));
					LL_old += log(P[2*z[t][n-1]-1][2*z[t][n]-1]/(P[2*z[t][n]-1][1]+P[2*z[t][n]-1][3]));
				}
				phmm_old->N = K;
				phmm_old->M = lengths[t];
				for(k=1; k<K; k++){
					phmm_old->pi[2*k-1] = Pimat[t][k];
					phmm_old->mu[2*k-1] = mu[t][k];
					phmm_old->sigma2[2*k-1] = sigma2[t][k];
				}
				phmm_old->pi[2] = 0.0;
				phmm_old->mu[2] = mu0[2] + gsl_ran_gaussian (rng, sqrt(sigma20[2]));
				phmm_old->sigma2[2] = nu*s2[2]/gsl_ran_chisq(rng,nu);
				for(k=1; k<=K; k++)
					for(j=1; j<=K; j++)
						phmm_old->A[k][j] = P[k][j];
				calculateBmatrix(phmm_old, &data[t][0]);
				samplehiddenstates(phmm_old, &data[t][0], &z[t][0], t+gibbs_count, tempA, logbeta, phmm_temp);
				LL_new = log(gsl_ran_gaussian_pdf(data[t][0]-phmm_old->mu[z[t][0]],sqrt(phmm_old->sigma2[z[t][0]])));
				for(n=1; n<lengths[t]; n++){
					LL_new += log(gsl_ran_gaussian_pdf(data[t][n]-phmm_old->mu[z[t][n]],sqrt(phmm_old->sigma2[z[t][n]])));
					LL_new += log(phmm_old->A[z[t][n-1]][z[t][n]]);
				}
				for(k=1; k<=K; k++){
					LL_new += log(gsl_ran_gaussian_pdf(phmm_old->mu[k]-mu0[k],sqrt(sigma20[k])));
					LL_new += log(gsl_ran_chisq_pdf(nu*s2[k]/phmm_old->sigma2[k],nu));
				}
//				printf("%f\t%f\n",LL_old,LL_new);
			}
			else{
				printf("state index error\n");
				exit(0);
			}
			diff = LL_new - LL_old;
//			printf("diff%f\n",diff);
			sum = 1.0/(1.0 + exp(diff)); //probability of staying in original one
//			printf("%f\t",sum);
			num2 = gsl_ran_flat(rng, 0.0, 1.0);
//			printf("%f\n",num2);
			if(num2 > sum){
//				printf("1\n");
				if(state_index[t]==K){
					Pimat[t][2] = 1.0 - Pimat[t][1];
					mu[t][2] = mu[t][3];
					sigma2[t][2] = sigma2[t][3];
				}
				else if(state_index[t]==K-1){
					Pimat[t][3] = Pimat[t][2];
					Pimat[t][2] = 0.0;
					mu[t][3] = mu[t][2];
					mu[t][2] = phmm_old->mu[2];
					sigma2[t][3] = sigma2[t][2];
					sigma2[t][2] = phmm_old->sigma2[2];
				}
				else{
					printf("state index error\n");
					exit(0);
				}
//				printf("%dchangeindex%d%d\n",t,state_index[t],phmm_old->N);
				state_index[t] = phmm_old->N;
			}
		}

		if (gibbs_count >= burnin){
			if ((gibbs_count-burnin)%thin == 0){
				for(k=1; k<=K; k++){
				mu0_gibbs[gibbs_save][k] = mu0[k];
				printf("%f\t",mu0[k]);
				sigma20_gibbs[gibbs_save][k] = sigma20[k];
				printf("%f\t",sigma20[k]);
				for(j=1; j<=K;j++){
					P_gibbs[gibbs_save][k][j] = P[k][j];
					printf("%f\t",P[k][j]);
				}
				}
				gibbs_save ++;
				printf("\n");
			}
		}
		gibbs_count ++;
	}
}


void hier_two_groups_old(vector<vector<double> > data, vector<vector<int> > z, int *state_index, int K,
					 double **mu, double **sigma2, double *mu0, double *sigma20, double **mu0_gibbs,
					 double ***P_gibbs, double **sigma20_gibbs, int T, HMM *phmm_temp, int *lengths,
					 double **logalpha, double **logbeta, double **loggamma, double *plogprobinit,
					 double *plogprobfinal, double ***logxi, double *tempA, double nu, double s2_init,
					 double **P, int seed, int num_gibbs, int burnin, int thin, double **Pimat,
					 HMM *phmm_old, double **count)
{	
	gsl_rng * rng = gsl_rng_alloc (gsl_rng_taus);
    gsl_rng_set(rng,seed);

    // initialize the parameters mu0,sigma20 after hierEM_uneqvar() function
	// initialize mu, sigma2 by doing Baum-Welch for each individual trace

	double *s2;
	s2 = new double[K+1];
	int t, k, j, n;
	for(k=1; k<=K; k++)
		s2[k] = s2_init;
	gsl_vector *vv = gsl_vector_alloc(K);
	gsl_permutation * perm = gsl_permutation_alloc(K);
    gsl_permutation * rank = gsl_permutation_alloc(K);
	gsl_vector *vv2 = gsl_vector_alloc(K-1);
	gsl_permutation * perm2 = gsl_permutation_alloc(K-1);
    gsl_permutation * rank2 = gsl_permutation_alloc(K-1);
	int **which_two;
	which_two = new int *[T];
	for(t=0; t<T; t++)
		which_two[t] = new int[2];
	double LL_K;
	double sum_temp;
	double *LL_temp;
	double LL_new, LL_old;
	LL_temp = new double[K];
	int kk;
	int **idxmat;
	idxmat = new int*[3];
	for(k=0; k<K; k++)
		idxmat[k] = new int[2];
	idxmat[0][0] = 1;
	idxmat[0][1] = 2;
	idxmat[1][0] = 1;
	idxmat[1][1] = 3;
	idxmat[2][0] = 2;
	idxmat[2][1] = 3;
	double **mu_temp;
	double **sig_temp;
	double **Pi_temp;
	mu_temp = new double*[K];
	for(k=0; k<K; k++)
		mu_temp[k] = new double[2];
	sig_temp = new double*[K];
	for(k=0; k<K; k++)
		sig_temp[k] = new double[2];
	Pi_temp = new double*[K];
	for(k=0; k<K; k++)
		Pi_temp[k] = new double[2];
	int num_para_iter = 100;
	int count_para;

	// initialize all the paremeters and fix state indicator of first observation

	for(t=0; t<T; t++){
		initializeHMM_num(phmm_temp, lengths[t], state_index[t]);
		calculateBmatrix(phmm_temp, &data[t][0]);
//		printf("%d\n",state_index[t]);
		if(state_index[t]==K){
		BaumWelch(phmm_temp, &data[t][0], lengths[t], logalpha, logbeta, loggamma, 
			     plogprobinit, plogprobfinal, logxi, tempA, vv, perm, rank, nu, &s2_init);
		which_two[t][0] = 0;
		which_two[t][1] = 0;
		}
		else if(state_index[t]==K-1){
		BaumWelch(phmm_temp, &data[t][0], lengths[t], logalpha, logbeta, loggamma, 
			     plogprobinit, plogprobfinal, logxi, tempA, vv2, perm2, rank2, nu, &s2_init);
		which_two[t][0] = 1;
		which_two[t][1] = 3;
		}
		for(k=1; k<=state_index[t]; k++){
			mu[t][k] = phmm_temp->mu[k];
//			printf("%d\t%d\t%f\n",t,k,mu[t][k]);
			sigma2[t][k] = phmm_temp->sigma2[k];
//			printf("%d\t%d\t%f\n",t,k,sigma2[t][k]);
			// fix z[0] with value of pi from EM
			Pimat[t][k] = phmm_temp->pi[k];
//			printf("%f\t",Pimat[t][k]);
			if(phmm_temp->pi[k]>0.5)
				z[t][0] = k;
		}
//		printf("\n");
//		printf("%d\t",z[t][0]);
	}
	for(t=1; t <= K; t++){
		for(k=1; k <= K; k++)
			P[t][k] = 1.0/K;
		sigma20[t] = 0.001;
	}

	int gibbs_count = 0;
	int gibbs_save = 0;
	double sum;
	int num;
	double num2;
	int Ngibbs = burnin + (num_gibbs - 1) * thin + 1;

	printf("Before Gibbs\n");

	while(gibbs_count < Ngibbs){

		printf("gibbs%d\n",gibbs_count);
		// update mu0

		count_para = 0;
		while(count_para < num_para_iter){
		for(k=1; k <= K; k++){
			sum = 0.0;
			num = 0;
			for(t=0; t< T; t++){					
				if(state_index[t]==K){
					sum += mu[t][k];
					num ++;
				}
				else if(state_index[t]==K-1){
					for(j=0; j<2; j++){
						if(which_two[t][j]==k){
							sum += mu[t][j+1];
							num ++;
						}
					}
				}
				else{
					printf("state index error\n");
					exit(0);
				}
			}
			mu0[k] = sum/num + gsl_ran_gaussian (rng, sqrt(sigma20[k]/num));
//			printf("%f\t%d\tmu0%d\t%f\n",sum,num,k,mu0[k]);
		}
	
		// update sigma20
		for(k=1; k <= K; k++){
			sum = 0.0;
			num = 0;
			for(t=0; t < T; t++){
				if(state_index[t]==K){	
					sum += pow(mu[t][k]-mu0[k],2.0);
					num ++;
				}
				else if(state_index[t]==K-1){
					for(j=0; j<2; j++){
						if(which_two[t][j]==k){
						sum += pow(mu[t][j+1]-mu0[k],2.0);
						num ++;
						}
					}
				}
				else{
					printf("state index error\n");
					exit(0);
				}
			}
			sigma20[k] = sum / gsl_ran_chisq(rng, 1.0 * (num - 2));
//			printf("sigma20%d\t%f\n",k,sigma20[k]);
		}

		// update s2
		for(k=1; k<=K; k++){
			num = 2;
			sum = 0.0;	
			for(t=0; t<T; t++){	
				if(state_index[t]==K){
					num ++;
					sum += 1.0/sigma2[t][k];
//					printf("%f\t",sigma2[t][k]);
				}
				else{
					for(j=0; j<2; j++){
						if(which_two[t][j]==k){
							sum += 1.0/sigma2[t][j+1];
							num ++;
						}
					}
				}
			}
			s2[k] = gsl_ran_chisq(rng,nu*num)/(sum*nu);
//			printf("sum%f\ts2%d\t%f\n",sum,k,s2[k]);
		}

		// update z
			
		for(t=0; t<T; t++){
		phmm_old->N = state_index[t];
		phmm_old->M = lengths[t];
		for(k=1; k<=state_index[t]; k++){
			phmm_old->pi[k] = Pimat[t][k];
			phmm_old->mu[k] = mu[t][k];
			phmm_old->sigma2[k] = sigma2[t][k];
		}
		if(state_index[t]==K){
			for(k=1; k<=state_index[t]; k++)
				for(j=1; j<=state_index[t]; j++)
					phmm_old->A[k][j] = P[k][j];
		}
		else if(state_index[t]==K-1){
			for(k=0; k<state_index[t]; k++){
				sum = 0.0;
				for(j=0; j<state_index[t]; j++)
					sum += P[which_two[t][k]][which_two[t][j]];
				for(j=0; j<state_index[t]; j++)
					phmm_old->A[k+1][j+1] = P[which_two[t][k]][which_two[t][j]]/sum;
			}
		}
		else{
			printf("state index error\n");
			exit(0);
		}
		calculateBmatrix(phmm_old, &data[t][0]);
		sum = 0.0;
		for(k=1; k<=state_index[t]; k++){
			if(phmm_old->pi[k]>=0.5)
				z[t][0] = k;
			if(phmm_old->pi[k]<-0.01 || phmm_old->pi[k]>1.01){
				printf("unreasonablePI\n");
				exit(0);
			}
			sum += phmm_old->pi[k];
		}
		if(sum < 0.99 || sum > 1.01){
			printf("unreasonablePIsum\n");
			exit(0);
		}
//		printf("%d%d%d%d%d\n",z[t][0],z[t][1],z[t][2],z[t][3],z[t][4]);
		samplehiddenstates(phmm_old, &data[t][0], &z[t][0], t+gibbs_count, tempA, logbeta, phmm_temp);
//		printf("%d%d%d%d%d\n",z[t][0],z[t][1],z[t][2],z[t][3],z[t][4]);
		}
		// update parameters for individual traces
		// update mu
		for(t=0; t<T; t++){
			if(state_index[t]==K){
				for(k=1; k<=K; k++){
					sum = 0.0;
					num = 0;
					for(j=0; j<lengths[t]; j++){
						if(z[t][j]==k){
							sum += data[t][j];
							num ++;
						}
					}
					mu[t][k] = (mu0[k]/sigma20[k] + sum/sigma2[t][k])/(1.0/sigma20[k]+1.0*num/sigma2[t][k]) + 
						gsl_ran_gaussian (rng, sqrt(1.0/(1.0/sigma20[k]+1.0*num/sigma2[t][k])));
//					printf("mu%d\t%d\t%f\n",t,k,mu[t][k]);
				}
			}
			else if(state_index[t]==K-1){
				for(k=1; k<K; k++){
					sum = 0.0;
					num = 0;
					for(j=0; j<lengths[t]; j++){
						if(z[t][j]==k){
							sum += data[t][j];
							num ++;
						}
					}
					mu[t][k] = (mu0[which_two[t][k-1]]/sigma20[which_two[t][k-1]] + 
						sum/sigma2[t][k])/(1.0/sigma20[which_two[t][k-1]]+1.0*num/sigma2[t][k]) + 
						gsl_ran_gaussian (rng, sqrt(1.0/(1.0/sigma20[which_two[t][k-1]]+1.0*num/sigma2[t][k])));
//					printf("mu%d\t%d\t%f\n",t,k,mu[t][k]);
				}
			}
			else{
				printf("state index error\n");
				exit(0);
			}
			// sort mu[t]
			if(state_index[t]==K){	
				for(k=0; k<K; k++)
					gsl_vector_set(vv,k,mu[t][k+1]);
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
//				for(k=1; k<=K; k++)
//				printf("mu%d\t%f\t",t,mu[t][k]);
//				printf("\n");
			}
			else if(state_index[t]==K-1){	
				for(k=0; k<K-1; k++)
					gsl_vector_set(vv2,k,mu[t][k+1]);
				gsl_sort_vector_index (perm2, vv2);
				gsl_permutation_inverse (rank2, perm2);
				for(k=0; k<K-1; k++){
					if(rank2->data[k] != k){
						for(n=0; n<lengths[t]; n++){
							if(z[t][n]==k+1){
								z[t][n] = rank2->data[k]+1;
							}
						}
					}
				}
				sort(mu[t]+1,mu[t]+K);
//				for(k=1; k<K; k++)
//				printf("mu%d\t%f\t",t,mu[t][k]);
//				printf("\n");
			}
		}
		// update sigma2
		for(t=0; t<T; t++){
			for(k=1; k<=state_index[t]; k++){
				num2 = nu;
				if(state_index[t]==K)
					sum = s2[k] * nu;
				else
					sum = s2[which_two[t][k-1]] * nu;
				for(j=0; j<lengths[t]; j++){
					if(z[t][j]==k){
						num2 += 1.0;
						sum += pow(data[t][j]-mu[t][k],2.0);
					}
				}
				sigma2[t][k] = sum/gsl_ran_chisq(rng,num2);
//				printf("sig%d\t%d\t%f\n",t,k,sigma2[t][k]);
			}
		}

		// update transition matrix

	for(k=1; k<=K; k++)
		for(j=1; j<=K; j++)
			count[k-1][j-1] = 1.0;
	for(t=0; t<T; t++){	
		if(state_index[t]==K){
		for(n=1; n<lengths[t]; n++){
			count[z[t][n-1]-1][z[t][n]-1] += 1.0;
		}
		}
		else{
//			printf("%d\t%f\t%f\t",t,Pimat[t][1],Pimat[t][2]);
		for(n=1; n<lengths[t]; n++){
//			printf("%d\twhichtwo%d\n",z[t][n-1],which_two[t][z[t][n-1]-1]-1);
			count[which_two[t][z[t][n-1]-1]-1][which_two[t][z[t][n]-1]-1] += 1.0;
		}
		}
	}

	for(k=1; k<=K; k++){
		gsl_ran_dirichlet (rng, K, count[k-1] , tempA);
		for(j=1; j<=K; j++){
			P[k][j] = tempA[j-1];
		}
	}
	count_para ++;
}
	// update state_index
	for(t=0; t<T; t++){
//		printf("%d\t%d\n",t,state_index[t]);
		if(state_index[t]==K){
			LL_K = log(gsl_ran_gaussian_pdf(data[t][0]-mu[t][z[t][0]],sqrt(sigma2[t][z[t][0]])));
			for(n=1; n<lengths[t]; n++){
				LL_K += log(gsl_ran_gaussian_pdf(data[t][n]-mu[t][z[t][n]],sqrt(sigma2[t][z[t][n]])));
				LL_K += log(P[z[t][n-1]][z[t][n]]);
			}
			for(k=1; k<=K; k++){
				LL_K += log(gsl_ran_gaussian_pdf(mu[t][k]-mu0[k],sqrt(sigma20[k])));
				LL_K += log(s2[k]*nu/2.0)*nu/2.0 - gsl_sf_lngamma(nu/2.0) - 
					s2[k]*nu/(2.0*sigma2[t][k]) -log(sigma2[t][k])*(nu/2.0+1.0);
			}
			phmm_old->N = K-1;
			phmm_old->M = lengths[t];
			for(kk=1; kk<=K; kk++){
				if(Pimat[t][idxmat[kk-1][0]]+Pimat[t][idxmat[kk-1][1]]!=0){
				    phmm_old->pi[1] = Pimat[t][idxmat[kk-1][0]]/(Pimat[t][idxmat[kk-1][0]]+Pimat[t][idxmat[kk-1][1]]);
					phmm_old->pi[2] = Pimat[t][idxmat[kk-1][1]]/(Pimat[t][idxmat[kk-1][0]]+Pimat[t][idxmat[kk-1][1]]);
				}
				else{
//					printf("warning: pi is not appropriate");
					if(Pimat[t][1]+Pimat[t][2]>0.5){
						phmm_old->pi[1] = 1.0;
						phmm_old->pi[2] = 0.0;
					}
					else{
						phmm_old->pi[1] = 0.0;
						phmm_old->pi[2] = 1.0;
					}
				}
				for(j=1; j<K; j++){
					phmm_old->mu[j] = mu[t][idxmat[kk-1][j-1]];
					phmm_old->sigma2[j] = sigma2[t][idxmat[kk-1][j-1]];
					for(k=1; k<K; k++)
					phmm_old->A[j][k] = P[idxmat[kk-1][j-1]][idxmat[kk-1][k-1]]/
					                    (P[idxmat[kk-1][j-1]][idxmat[kk-1][0]]+P[idxmat[kk-1][j-1]][idxmat[kk-1][1]]);
				}
				calculateBmatrix(phmm_old, &data[t][0]);
				samplehiddenstates(phmm_old, &data[t][0], &z[t][0], t+gibbs_count, tempA, logbeta, phmm_temp);
				LL_new = log(gsl_ran_gaussian_pdf(data[t][0]-phmm_old->mu[z[t][0]],sqrt(phmm_old->sigma2[z[t][0]])));
				for(n=1; n<lengths[t]; n++){
					LL_new += log(gsl_ran_gaussian_pdf(data[t][n]-phmm_old->mu[z[t][n]],sqrt(phmm_old->sigma2[z[t][n]])));
					LL_new += log(phmm_old->A[z[t][n-1]][z[t][n]]);
				}
				for(k=1; k<K; k++){
					LL_new += log(gsl_ran_gaussian_pdf(phmm_old->mu[k]-mu0[idxmat[kk-1][k-1]],sqrt(sigma20[idxmat[kk-1][k-1]])));
					LL_new +=  log(s2[idxmat[kk-1][k-1]]*nu/2.0)*nu/2.0 - gsl_sf_lngamma(nu/2.0) - 
						s2[idxmat[kk-1][k-1]]*nu/(2.0*phmm_old->sigma2[k]) -log(phmm_old->sigma2[k])*(nu/2.0+1);
				}
				LL_temp[kk-1] = LL_new;
				for(k=1; k<K; k++){
					Pi_temp[kk-1][k-1] = phmm_old->pi[k];
					mu_temp[kk-1][k-1] = phmm_old->mu[k];
					sig_temp[kk-1][k-1] = phmm_old->sigma2[k];
				}
				if(abs(Pi_temp[kk-1][0]+Pi_temp[kk-1][1]-1.0)>0.01){
					printf("%f\t%f\n",phmm_old->pi[1],phmm_old->pi[2]);
					printf("Pi_temp_sum_3_to_2_state_error\n");
					exit(0);
				}
			}
		}
		else if(state_index[t]==K-1){
			LL_old = log(gsl_ran_gaussian_pdf(data[t][0]-mu[t][z[t][0]],sqrt(sigma2[t][z[t][0]])));
			for(n=1; n<lengths[t]; n++){
				LL_old += log(gsl_ran_gaussian_pdf(data[t][n]-mu[t][z[t][n]],sqrt(sigma2[t][z[t][n]])));
				LL_old += log(P[which_two[t][z[t][n-1]-1]][which_two[t][z[t][n]-1]]/
					(P[which_two[t][z[t][n-1]-1]][which_two[t][0]]+P[which_two[t][z[t][n-1]-1]][which_two[t][1]]));
			}
			for(k=1; k<K; k++){
				LL_old += log(gsl_ran_gaussian_pdf(phmm_old->mu[k]-mu0[which_two[t][k-1]],sqrt(sigma20[which_two[t][k-1]])));
				LL_old +=  log(s2[which_two[t][k-1]]*nu/2.0)*nu/2.0 - gsl_sf_lngamma(nu/2.0) - 
						s2[which_two[t][k-1]]*nu/(2.0*phmm_old->sigma2[k]) -log(phmm_old->sigma2[k])*(nu/2.0+1.0);
			}
			phmm_old->N = K-1;
			phmm_old->M = lengths[t];
			kk = 0;
			while(kk < K){
				if(which_two[t][0]==idxmat[kk][0] && which_two[t][1]==idxmat[kk][1]){
					LL_temp[kk] = LL_old;
					for(k=1; k<K; k++){
						Pi_temp[kk][k-1] = Pimat[t][k];
						mu_temp[kk][k-1] = mu[t][k];
						sig_temp[kk][k-1] = sigma2[t][k];
					}
					if(abs(Pi_temp[kk][0]+Pi_temp[kk][1]-1.0)>0.01){
						printf("2_stay_Pi_temp_sum error\n");
						exit(1);
					}
					kk ++;
				}
				else{
					for(j=1; j<K; j++){
						for(k=1; k<K; k++){
							if(which_two[t][k-1]==idxmat[kk][j-1]){
								phmm_old->pi[j] = Pimat[t][k];
								phmm_old->pi[3-j] = 1.0 - phmm_old->pi[j];
								phmm_old->mu[j] = mu[t][k];
								phmm_old->sigma2[j] = sigma2[t][k];
							}
						}
						if(which_two[t][0]!=idxmat[kk][j-1] && which_two[t][1]!=idxmat[kk][j-1]){
							phmm_old->mu[j] = mu0[idxmat[kk][j-1]] + gsl_ran_gaussian (rng, sqrt(sigma20[idxmat[kk][j-1]]));
							phmm_old->sigma2[j] = nu*s2[idxmat[kk][j-1]]/gsl_ran_chisq(rng,nu);
						}
						for(k=1; k<K; k++)
							phmm_old->A[j][k] = P[idxmat[kk][j-1]][idxmat[kk][k-1]]/
					                    (P[idxmat[kk][j-1]][idxmat[kk][0]]+P[idxmat[kk][j-1]][idxmat[kk][1]]);
					}				
					calculateBmatrix(phmm_old, &data[t][0]);
					samplehiddenstates(phmm_old, &data[t][0], &z[t][0], t+gibbs_count, tempA, logbeta, phmm_temp);
				    LL_new = log(gsl_ran_gaussian_pdf(data[t][0]-phmm_old->mu[z[t][0]],sqrt(phmm_old->sigma2[z[t][0]])));				
					for(n=1; n<lengths[t]; n++){
					    LL_new += log(gsl_ran_gaussian_pdf(data[t][n]-phmm_old->mu[z[t][n]],sqrt(phmm_old->sigma2[z[t][n]])));
					    LL_new += log(phmm_old->A[z[t][n-1]][z[t][n]]);
					}
					for(k=1; k<K; k++){
					    LL_new += log(gsl_ran_gaussian_pdf(phmm_old->mu[k]-mu0[idxmat[kk][k-1]],sqrt(sigma20[idxmat[kk][k-1]])));
					    LL_new +=  log(s2[idxmat[kk][k-1]]*nu/2.0)*nu/2.0 - gsl_sf_lngamma(nu/2.0) - 
						           s2[idxmat[kk][k-1]]*nu/(2.0*phmm_old->sigma2[k]) -log(phmm_old->sigma2[k])*(nu/2.0+1);
					}
					LL_temp[kk] = LL_new;
					for(k=1; k<K; k++){
						Pi_temp[kk][k-1] = phmm_old->pi[k];
						mu_temp[kk][k-1] = phmm_old->mu[k];
						sig_temp[kk][k-1] = phmm_old->sigma2[k];
					}
					if(abs(Pi_temp[kk][0]+Pi_temp[kk][1]-1.0)>0.01){
						printf("pi_temp_sum error\n");
						exit(1);
					}
					kk ++;
				}
			}
			phmm_old->N = K;
			for(k=1; k<=K; k++){
				for(j=0; j<K-1; j++){
					if(which_two[t][j]==k){
						phmm_old->pi[k] = Pimat[t][j+1];
//						printf("%f\t",Pimat[t][j+1]);
						phmm_old->mu[k] = mu[t][j+1];
						phmm_old->sigma2[k] = sigma2[t][j+1];
					}
				}
				if(which_two[t][0]!=k && which_two[t][1]!=k){
					phmm_old->mu[k] = mu0[k] + gsl_ran_gaussian(rng,sqrt(sigma20[k]));
					phmm_old->sigma2[k] = nu*s2[k]/gsl_ran_chisq(rng,nu);
				}
				for(j=1; j<=K; j++)
					phmm_old->A[k][j] = P[k][j];
			}			
			for(k=1; k<=K; k++){
				if(which_two[t][0]!=k && which_two[t][1]!=k){
					phmm_old->pi[k] = 1.0;
					for(j=1; j<=K; j++)
						if(j!=k)
							phmm_old->pi[k] = phmm_old->pi[k] - phmm_old->pi[j];
				}
			}
			if(abs(phmm_old->pi[1]+phmm_old->pi[2]+phmm_old->pi[3]-1.0)>0.01){
				printf("Phmm_old_pi_sum_error\n");
				exit(0);
			}
			calculateBmatrix(phmm_old, &data[t][0]);
			samplehiddenstates(phmm_old, &data[t][0], &z[t][0], t+gibbs_count, tempA, logbeta, phmm_temp);
			LL_K = log(gsl_ran_gaussian_pdf(data[t][0]-phmm_old->mu[z[t][0]],sqrt(phmm_old->sigma2[z[t][0]])));
			for(n=1; n<lengths[t]; n++){
				LL_K += log(gsl_ran_gaussian_pdf(data[t][n]-phmm_old->mu[z[t][n]],sqrt(phmm_old->sigma2[z[t][n]])));
				LL_K += log(phmm_old->A[z[t][n-1]][z[t][n]]);
			}
			for(k=1; k<=K; k++){
				LL_K += log(gsl_ran_gaussian_pdf(phmm_old->mu[k]-mu0[k],sqrt(sigma20[k])));
				LL_K += log(gsl_ran_chisq_pdf(nu*s2[k]/phmm_old->sigma2[k],nu));
			}
			}
			else{
				printf("state index error\n");
				exit(0);
			}

			// determine whether to change number of states/ change which two states show
            // multinomial distribution propto exp(LL_K,LL_temp)
			num2 = gsl_ran_flat(rng, 0.0, 1.0);
			LL_new = LL_K;
			for(k=0; k<K; k++)
				LL_new = MAX(LL_new,LL_temp[k]);
			sum = exp(LL_K-LL_new);
			for(k=0; k<K; k++)
				sum += exp(LL_temp[k]-LL_new);
			sum_temp = exp(LL_K-LL_new)/sum;
//			printf("%f\t%f\t",num2,sum_temp);
			if(num2 < sum_temp){
//				printf("LL_K%d\n",K);
				if(state_index[t]==K-1){
					printf("changeindex%d\t%d\t%d\n",t,K-1,K);
					for(j=1; j<=K; j++){
						Pimat[t][j] = phmm_old->pi[j];
						mu[t][j] = phmm_old->mu[j];
						sigma2[t][j] = phmm_old->sigma2[j];
					}
					state_index[t] = K;
					if(abs(Pimat[t][1]+Pimat[t][2]+Pimat[t][3]-1.0)>0.01){
						printf("2_to_3_pi_mat_error\n");
						exit(0);
					}
				}
			}
			else{					
				k = 0;
				while(sum_temp <= num2){
					sum_temp += exp(LL_temp[k]-LL_new)/sum;
//					printf("%f\t",sum_temp);
					k ++;
				}
				k--;
//				printf("Ltemp%d\n",k);
				if(state_index[t]==K){
					printf("changeindex%d\t%d\t%d\n",t,K,K-1);
					for(j=0; j<K-1; j++){
						which_two[t][j] = idxmat[k][j];
						Pimat[t][j+1] = Pi_temp[k][j];
						mu[t][j+1] = mu_temp[k][j];
						sigma2[t][j+1] = sig_temp[k][j];
					}
					state_index[t] = K-1;
					if(abs(Pimat[t][1]+Pimat[t][2]-1.0)>0.01){
						printf("3_to_2Pimaterror\n");
						exit(0);
					}
				}
				else if(which_two[t][0]!=idxmat[k][0] || which_two[t][1]!=idxmat[k][1]){
					printf("%d\t change which two states\n",t);
					for(j=0; j<K-1; j++){
						which_two[t][j] = idxmat[k][j];
						Pimat[t][j+1] = Pi_temp[k][j];
						mu[t][j+1] = mu_temp[k][j];
						sigma2[t][j+1] = sig_temp[k][j];
					}
					if(abs(Pimat[t][1]+Pimat[t][2]-1.0)>0.01){
						printf("2_to_2Pimaterror\n");
						exit(0);
					}
					state_index[t] = K-1;
				}
			}
		}

		if (gibbs_count >= burnin){
			if ((gibbs_count-burnin)%thin == 0){
				for(t=0; t<T; t++)
					printf("%d\t",state_index[t]);
				printf("\n");
				for(k=1; k<=K; k++){
				mu0_gibbs[gibbs_save][k] = mu0[k];
				printf("%f\t",mu0[k]);
				sigma20_gibbs[gibbs_save][k] = sigma20[k];
				printf("%f\t",sigma20[k]);
				for(j=1; j<=K;j++){
					P_gibbs[gibbs_save][k][j] = P[k][j];
					printf("%f\t",P[k][j]);
				}
				}
				gibbs_save ++;
				printf("\n");
			}
			if((gibbs_count-burnin)%(thin*10)==0){
					std::string str_name;
					str_name.assign("3_S_paras");
					str_name.append(convertInt(gibbs_count+1));
					str_name.append(".txt");
					cout << str_name << "\n";
					ofstream fout;
					fout.open(str_name.c_str(),ios::app);    // open file for appending
					assert (!fout.fail( ));
					for(t=0; t<T; t++){
						if(state_index[t]==K){
							fout<<t<<"\t";
							for(j=1; j<=K; j++)
								fout<<mu[t][j]<<"\t";
							for(j=1; j<=K; j++)
								fout<<sigma2[t][j]<<"\t";
							for(j=1; j<=K; j++)
								fout<<Pimat[t][j]<<"\t";
							fout<<"\n";
						}
					}
					fout.close( );       //close file
					assert(!fout.fail( ));
					str_name.assign("2_S_paras");
					str_name.append(convertInt(gibbs_count+1));
					str_name.append(".txt");
					cout << str_name << "\n";
					fout.open(str_name.c_str(),ios::app);    // open file for appending
					assert (!fout.fail( ));   
					for(t=0; t<T; t++){
						if(state_index[t]==K-1){
							fout<<t<<"\t";
							fout<<which_two[t][0]<<"\t"<<which_two[t][1]<<"\t";
							for(j=1; j<K; j++)
								fout<<mu[t][j]<<"\t";
							for(j=1; j<K; j++)
								fout<<sigma2[t][j]<<"\t";
							for(j=1; j<K; j++)
								fout<<Pimat[t][j]<<"\t";
							fout<<"\n";
						}
					}
					fout.close( );       //close file
					assert(!fout.fail( ));
					str_name.assign("global_paras");
					str_name.append(convertInt(gibbs_count+1));
					str_name.append(".txt");
					cout << str_name << "\n";
					fout.open(str_name.c_str(),ios::app);    // open file for appending
					assert (!fout.fail( )); 
							for(j=1; j<=K; j++)
								fout<<mu0[j]<<"\t";
							fout<<"\n";
							for(j=1; j<=K; j++)
								fout<<sigma20[j]<<"\t";
							fout<<"\n";
							for(k=1; k<=K; k++){
								for(j=1; j<=K; j++)
								fout<<P[k][j]<<"\t";
								fout<<"\n";
							}
					fout.close( );       //close file
					assert(!fout.fail( ));
			}
		}
		gibbs_count ++;
	}
}



void hier_two_groups_after_select(vector<vector<double> > data, vector<vector<int> > z, int *state_index, int K,
					 double **mu, double **sigma2, double *mu0, double *sigma20, double **mu0_gibbs,
					 double ***P_gibbs, double **sigma20_gibbs, int T, HMM *phmm_temp, int *lengths,
					 double **logalpha, double **logbeta, double **loggamma, double *plogprobinit,
					 double *plogprobfinal, double ***logxi, double *tempA, double nu, double s2_init,
					 double **P, int seed, int num_gibbs, int burnin, int thin, double **Pimat,
					 HMM *phmm_old, double **count)
{	
	gsl_rng * rng = gsl_rng_alloc (gsl_rng_taus);
    gsl_rng_set(rng,seed);

    // initialize the parameters mu0,sigma20 after hierEM_uneqvar() function
	// initialize mu, sigma2 by doing Baum-Welch for each individual trace

	double *s2;
	s2 = new double[K+1];
	int t, k, j, n;
	for(k=1; k<=K; k++)
		s2[k] = s2_init;
	gsl_vector *vv = gsl_vector_alloc(K);
	gsl_permutation * perm = gsl_permutation_alloc(K);
    gsl_permutation * rank = gsl_permutation_alloc(K);
	gsl_vector *vv2 = gsl_vector_alloc(K-1);
	gsl_permutation * perm2 = gsl_permutation_alloc(K-1);
    gsl_permutation * rank2 = gsl_permutation_alloc(K-1);
	int **which_two;
	which_two = new int *[T];
	for(t=0; t<T; t++)
		which_two[t] = new int[2];
	double sum_temp;
	double *LL_temp;
	double LL_new, LL_old;
	LL_temp = new double[K];
	int kk;
	int **idxmat;
	idxmat = new int*[3];
	for(k=0; k<K; k++)
		idxmat[k] = new int[2];
	idxmat[0][0] = 1;
	idxmat[0][1] = 2;
	idxmat[1][0] = 1;
	idxmat[1][1] = 3;
	idxmat[2][0] = 2;
	idxmat[2][1] = 3;
	double **mu_temp;
	double **sig_temp;
	double **Pi_temp;
	mu_temp = new double*[K];
	for(k=0; k<K; k++)
		mu_temp[k] = new double[2];
	sig_temp = new double*[K];
	for(k=0; k<K; k++)
		sig_temp[k] = new double[2];
	Pi_temp = new double*[K];
	for(k=0; k<K; k++)
		Pi_temp[k] = new double[2];
	int num_para_iter = 2;
	int count_para;

	// initialize all the paremeters and fix state indicator of first observation

	for(t=0; t<T; t++){
		initializeHMM_num(phmm_temp, lengths[t], state_index[t]);
		calculateBmatrix(phmm_temp, &data[t][0]);
//		printf("%d\n",state_index[t]);
		if(state_index[t]==K){
		BaumWelch(phmm_temp, &data[t][0], lengths[t], logalpha, logbeta, loggamma, 
			     plogprobinit, plogprobfinal, logxi, tempA, vv, perm, rank, nu, &s2_init);
		which_two[t][0] = 0;
		which_two[t][1] = 0;
		}
		else if(state_index[t]==K-1){
		BaumWelch(phmm_temp, &data[t][0], lengths[t], logalpha, logbeta, loggamma, 
			     plogprobinit, plogprobfinal, logxi, tempA, vv2, perm2, rank2, nu, &s2_init);
		which_two[t][0] = 1;
		which_two[t][1] = 3;
		}
		for(k=1; k<=state_index[t]; k++){
			mu[t][k] = phmm_temp->mu[k];
//			printf("%d\t%d\t%f\n",t,k,mu[t][k]);
			sigma2[t][k] = phmm_temp->sigma2[k];
//			printf("%d\t%d\t%f\n",t,k,sigma2[t][k]);
			// fix z[0] with value of pi from EM
			Pimat[t][k] = phmm_temp->pi[k];
//			printf("%f\t",Pimat[t][k]);
			if(phmm_temp->pi[k]>0.5)
				z[t][0] = k;
		}
//		printf("\n");
//		printf("%d\t",z[t][0]);
	}
	for(t=1; t <= K; t++){
		for(k=1; k <= K; k++)
			P[t][k] = 1.0/K;
		sigma20[t] = 0.001;
	}

	int gibbs_count = 0;
	int gibbs_save = 0;
	double sum;
	int num;
	double num2;
	int Ngibbs = burnin + (num_gibbs - 1) * thin + 1;

	printf("Before Gibbs\n");

	while(gibbs_count < Ngibbs){

		printf("gibbs%d\n",gibbs_count);
		// update mu0

		count_para = 0;
		while(count_para < num_para_iter){
		for(k=1; k <= K; k++){
			sum = 0.0;
			num = 0;
			for(t=0; t< T; t++){					
				if(state_index[t]==K){
					sum += mu[t][k];
					num ++;
				}
				else if(state_index[t]==K-1){
					for(j=0; j<2; j++){
						if(which_two[t][j]==k){
							sum += mu[t][j+1];
							num ++;
						}
					}
				}
				else{
					printf("state index error\n");
					exit(0);
				}
			}
			mu0[k] = sum/num + gsl_ran_gaussian (rng, sqrt(sigma20[k]/num));
//			printf("%f\t%d\tmu0%d\t%f\n",sum,num,k,mu0[k]);
		}
	
		// update sigma20
		for(k=1; k <= K; k++){
			sum = 0.0;
			num = 0;
			for(t=0; t < T; t++){
				if(state_index[t]==K){	
					sum += pow(mu[t][k]-mu0[k],2.0);
					num ++;
				}
				else if(state_index[t]==K-1){
					for(j=0; j<2; j++){
						if(which_two[t][j]==k){
						sum += pow(mu[t][j+1]-mu0[k],2.0);
						num ++;
						}
					}
				}
				else{
					printf("state index error\n");
					exit(0);
				}
			}
			sigma20[k] = sum / gsl_ran_chisq(rng, 1.0 * (num - 2));
//			printf("sigma20%d\t%f\n",k,sigma20[k]);
		}

		// update s2
		for(k=1; k<=K; k++){
			num = 2;
			sum = 0.0;	
			for(t=0; t<T; t++){	
				if(state_index[t]==K){
					num ++;
					sum += 1.0/sigma2[t][k];
//					printf("%f\t",sigma2[t][k]);
				}
				else{
					for(j=0; j<2; j++){
						if(which_two[t][j]==k){
							sum += 1.0/sigma2[t][j+1];
							num ++;
						}
					}
				}
			}
			s2[k] = gsl_ran_chisq(rng,nu*num)/(sum*nu);
//			printf("sum%f\ts2%d\t%f\n",sum,k,s2[k]);
		}

		// update z
			
		for(t=0; t<T; t++){
		phmm_old->N = state_index[t];
		phmm_old->M = lengths[t];
		for(k=1; k<=state_index[t]; k++){
			phmm_old->pi[k] = Pimat[t][k];
			phmm_old->mu[k] = mu[t][k];
			phmm_old->sigma2[k] = sigma2[t][k];
		}
		if(state_index[t]==K){
			for(k=1; k<=state_index[t]; k++)
				for(j=1; j<=state_index[t]; j++)
					phmm_old->A[k][j] = P[k][j];
		}
		else if(state_index[t]==K-1){
			for(k=0; k<state_index[t]; k++){
				sum = 0.0;
				for(j=0; j<state_index[t]; j++)
					sum += P[which_two[t][k]][which_two[t][j]];
				for(j=0; j<state_index[t]; j++)
					phmm_old->A[k+1][j+1] = P[which_two[t][k]][which_two[t][j]]/sum;
			}
		}
		else{
			printf("state index error\n");
			exit(0);
		}
		calculateBmatrix(phmm_old, &data[t][0]);
		sum = 0.0;
		for(k=1; k<=state_index[t]; k++){
			if(phmm_old->pi[k]>=0.5)
				z[t][0] = k;
			if(phmm_old->pi[k]<-0.01 || phmm_old->pi[k]>1.01){
				printf("unreasonablePI\n");
				exit(0);
			}
			sum += phmm_old->pi[k];
		}
		if(sum < 0.99 || sum > 1.01){
			printf("unreasonablePIsum\n");
			exit(0);
		}
//		printf("%d%d%d%d%d\n",z[t][0],z[t][1],z[t][2],z[t][3],z[t][4]);
		samplehiddenstates(phmm_old, &data[t][0], &z[t][0], t+gibbs_count, tempA, logbeta, phmm_temp);
//		printf("%d%d%d%d%d\n",z[t][0],z[t][1],z[t][2],z[t][3],z[t][4]);
		}
		// update parameters for individual traces
		// update mu
		for(t=0; t<T; t++){
			if(state_index[t]==K){
				for(k=1; k<=K; k++){
					sum = 0.0;
					num = 0;
					for(j=0; j<lengths[t]; j++){
						if(z[t][j]==k){
							sum += data[t][j];
							num ++;
						}
					}
					mu[t][k] = (mu0[k]/sigma20[k] + sum/sigma2[t][k])/(1.0/sigma20[k]+1.0*num/sigma2[t][k]) + 
						gsl_ran_gaussian (rng, sqrt(1.0/(1.0/sigma20[k]+1.0*num/sigma2[t][k])));
//					printf("mu%d\t%d\t%f\n",t,k,mu[t][k]);
				}
			}
			else if(state_index[t]==K-1){
				for(k=1; k<K; k++){
					sum = 0.0;
					num = 0;
					for(j=0; j<lengths[t]; j++){
						if(z[t][j]==k){
							sum += data[t][j];
							num ++;
						}
					}
					mu[t][k] = (mu0[which_two[t][k-1]]/sigma20[which_two[t][k-1]] + 
						sum/sigma2[t][k])/(1.0/sigma20[which_two[t][k-1]]+1.0*num/sigma2[t][k]) + 
						gsl_ran_gaussian (rng, sqrt(1.0/(1.0/sigma20[which_two[t][k-1]]+1.0*num/sigma2[t][k])));
//					printf("mu%d\t%d\t%f\n",t,k,mu[t][k]);
				}
			}
			else{
				printf("state index error\n");
				exit(0);
			}
			// sort mu[t]
			if(state_index[t]==K){	
				for(k=0; k<K; k++)
					gsl_vector_set(vv,k,mu[t][k+1]);
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
//				for(k=1; k<=K; k++)
//				printf("mu%d\t%f\t",t,mu[t][k]);
//				printf("\n");
			}
			else if(state_index[t]==K-1){	
				for(k=0; k<K-1; k++)
					gsl_vector_set(vv2,k,mu[t][k+1]);
				gsl_sort_vector_index (perm2, vv2);
				gsl_permutation_inverse (rank2, perm2);
				for(k=0; k<K-1; k++){
					if(rank2->data[k] != k){
						for(n=0; n<lengths[t]; n++){
							if(z[t][n]==k+1){
								z[t][n] = rank2->data[k]+1;
							}
						}
					}
				}
				sort(mu[t]+1,mu[t]+K);
//				for(k=1; k<K; k++)
//				printf("mu%d\t%f\t",t,mu[t][k]);
//				printf("\n");
			}
		}
		// update sigma2
		for(t=0; t<T; t++){
			for(k=1; k<=state_index[t]; k++){
				num2 = nu;
				if(state_index[t]==K)
					sum = s2[k] * nu;
				else
					sum = s2[which_two[t][k-1]] * nu;
				for(j=0; j<lengths[t]; j++){
					if(z[t][j]==k){
						num2 += 1.0;
						sum += pow(data[t][j]-mu[t][k],2.0);
					}
				}
				sigma2[t][k] = sum/gsl_ran_chisq(rng,num2);
//				printf("sig%d\t%d\t%f\n",t,k,sigma2[t][k]);
			}
		}

		// update transition matrix

	for(k=1; k<=K; k++)
		for(j=1; j<=K; j++)
			count[k-1][j-1] = 1.0;
	for(t=0; t<T; t++){	
		if(state_index[t]==K){
		for(n=1; n<lengths[t]; n++){
			count[z[t][n-1]-1][z[t][n]-1] += 1.0;
		}
		}
		else{
//			printf("%d\t%f\t%f\t",t,Pimat[t][1],Pimat[t][2]);
		for(n=1; n<lengths[t]; n++){
//			printf("%d\twhichtwo%d\n",z[t][n-1],which_two[t][z[t][n-1]-1]-1);
			count[which_two[t][z[t][n-1]-1]-1][which_two[t][z[t][n]-1]-1] += 1.0;
		}
		}
	}

	for(k=1; k<=K; k++){
		gsl_ran_dirichlet (rng, K, count[k-1] , tempA);
		for(j=1; j<=K; j++){
			P[k][j] = tempA[j-1];
		}
	}
	count_para ++;
}
	// update state_index
	for(t=0; t<T; t++){
		printf("%d\t%d\n",t,state_index[t]);
		if(state_index[t]==K-1){
			printf("%d\n",t);
			LL_old = log(gsl_ran_gaussian_pdf(data[t][0]-mu[t][z[t][0]],sqrt(sigma2[t][z[t][0]])));
			for(n=1; n<lengths[t]; n++){
				LL_old += log(gsl_ran_gaussian_pdf(data[t][n]-mu[t][z[t][n]],sqrt(sigma2[t][z[t][n]])));
				LL_old += log(P[which_two[t][z[t][n-1]-1]][which_two[t][z[t][n]-1]]/
					(P[which_two[t][z[t][n-1]-1]][which_two[t][0]]+P[which_two[t][z[t][n-1]-1]][which_two[t][1]]));
			}
			for(k=1; k<K; k++){
				LL_old += log(gsl_ran_gaussian_pdf(phmm_old->mu[k]-mu0[which_two[t][k-1]],sqrt(sigma20[which_two[t][k-1]])));
				LL_old +=  log(s2[which_two[t][k-1]]*nu/2.0)*nu/2.0 - gsl_sf_lngamma(nu/2.0) - 
						s2[which_two[t][k-1]]*nu/(2.0*phmm_old->sigma2[k]) -log(phmm_old->sigma2[k])*(nu/2.0+1.0);
			}
			phmm_old->N = K-1;
			phmm_old->M = lengths[t];
			kk = 0;
			while(kk < K){
				if(which_two[t][0]==idxmat[kk][0] && which_two[t][1]==idxmat[kk][1]){
					LL_temp[kk] = LL_old;
					for(k=1; k<K; k++){
						Pi_temp[kk][k-1] = Pimat[t][k];
						mu_temp[kk][k-1] = mu[t][k];
						sig_temp[kk][k-1] = sigma2[t][k];
					}
					if(abs(Pi_temp[kk][0]+Pi_temp[kk][1]-1.0)>0.01){
						printf("2_stay_Pi_temp_sum error\n");
						exit(1);
					}
					kk ++;
				}
				else{
					for(j=1; j<K; j++){
						for(k=1; k<K; k++){
							if(which_two[t][k-1]==idxmat[kk][j-1]){
								phmm_old->pi[j] = Pimat[t][k];
								phmm_old->pi[3-j] = 1.0 - phmm_old->pi[j];
								phmm_old->mu[j] = mu[t][k];
								phmm_old->sigma2[j] = sigma2[t][k];
							}
						}
						if(which_two[t][0]!=idxmat[kk][j-1] && which_two[t][1]!=idxmat[kk][j-1]){
							phmm_old->mu[j] = mu0[idxmat[kk][j-1]] + gsl_ran_gaussian (rng, sqrt(sigma20[idxmat[kk][j-1]]));
							phmm_old->sigma2[j] = nu*s2[idxmat[kk][j-1]]/gsl_ran_chisq(rng,nu);
						}
						for(k=1; k<K; k++)
							phmm_old->A[j][k] = P[idxmat[kk][j-1]][idxmat[kk][k-1]]/
					                    (P[idxmat[kk][j-1]][idxmat[kk][0]]+P[idxmat[kk][j-1]][idxmat[kk][1]]);
					}				
					calculateBmatrix(phmm_old, &data[t][0]);
					samplehiddenstates(phmm_old, &data[t][0], &z[t][0], t+gibbs_count, tempA, logbeta, phmm_temp);
				    LL_new = log(gsl_ran_gaussian_pdf(data[t][0]-phmm_old->mu[z[t][0]],sqrt(phmm_old->sigma2[z[t][0]])));				
					for(n=1; n<lengths[t]; n++){
					    LL_new += log(gsl_ran_gaussian_pdf(data[t][n]-phmm_old->mu[z[t][n]],sqrt(phmm_old->sigma2[z[t][n]])));
					    LL_new += log(phmm_old->A[z[t][n-1]][z[t][n]]);
					}
					for(k=1; k<K; k++){
					    LL_new += log(gsl_ran_gaussian_pdf(phmm_old->mu[k]-mu0[idxmat[kk][k-1]],sqrt(sigma20[idxmat[kk][k-1]])));
					    LL_new +=  log(s2[idxmat[kk][k-1]]*nu/2.0)*nu/2.0 - gsl_sf_lngamma(nu/2.0) - 
						           s2[idxmat[kk][k-1]]*nu/(2.0*phmm_old->sigma2[k]) -log(phmm_old->sigma2[k])*(nu/2.0+1);
					}
					LL_temp[kk] = LL_new;
					for(k=1; k<K; k++){
						Pi_temp[kk][k-1] = phmm_old->pi[k];
						mu_temp[kk][k-1] = phmm_old->mu[k];
						sig_temp[kk][k-1] = phmm_old->sigma2[k];
					}
					if(abs(Pi_temp[kk][0]+Pi_temp[kk][1]-1.0)>0.01){
						printf("pi_temp_sum error\n");
						exit(1);
					}
					kk ++;
				}
			}

			// determine which two states show
            // multinomial distribution propto exp(LL_temp)
            printf("after LLcalculation\n");
			num2 = gsl_ran_flat(rng, 0.0, 1.0);
			LL_new = 0.0;
			for(k=0; k<K; k++)
				LL_new = MAX(LL_new,LL_temp[k]);
			sum = 0.0;
			for(k=0; k<K; k++)
				sum += exp(LL_temp[k]-LL_new);
			sum_temp = 0.0;
//			printf("%f\t%f\t",num2,sum_temp);
			k = 0;
			while(sum_temp <= num2){
				sum_temp += exp(LL_temp[k]-LL_new)/sum;
//				printf("%f\t",sum_temp);
				k ++;
			}
			k--;
//			printf("Ltemp%d\n",k);
			if(state_index[t]==K){
				printf("Error in state index!%d\n",t);
				exit(2);
			}
			else if(which_two[t][0]!=idxmat[k][0] || which_two[t][1]!=idxmat[k][1]){
				printf("%d\t change which two states\n",t);
				for(j=0; j<K-1; j++){
					which_two[t][j] = idxmat[k][j];
					Pimat[t][j+1] = Pi_temp[k][j];
					mu[t][j+1] = mu_temp[k][j];
					sigma2[t][j+1] = sig_temp[k][j];
				}
				if(abs(Pimat[t][1]+Pimat[t][2]-1.0)>0.01){
					printf("2_to_2Pimaterror\n");
					exit(0);
				}
			}
			}
		}

		if (gibbs_count >= burnin){
			if ((gibbs_count-burnin)%thin == 0){
				for(t=0; t<T; t++)
					printf("%d\t",state_index[t]);
				printf("\n");
				for(k=1; k<=K; k++){
				mu0_gibbs[gibbs_save][k] = mu0[k];
				printf("%f\t",mu0[k]);
				sigma20_gibbs[gibbs_save][k] = sigma20[k];
				printf("%f\t",sigma20[k]);
				for(j=1; j<=K;j++){
					P_gibbs[gibbs_save][k][j] = P[k][j];
					printf("%f\t",P[k][j]);
				}
				}
				gibbs_save ++;
				printf("\n");
			}
			if((gibbs_count-burnin)%(thin*10)==0){
					std::string str_name;
					str_name.assign("3_S_paras");
					str_name.append(convertInt(gibbs_count+1));
					str_name.append(".txt");
					cout << str_name << "\n";
					ofstream fout;
					fout.open(str_name.c_str(),ios::app);    // open file for appending
					assert (!fout.fail( ));
					for(t=0; t<T; t++){
						if(state_index[t]==K){
							fout<<t<<"\t";
							for(j=1; j<=K; j++)
								fout<<mu[t][j]<<"\t";
							for(j=1; j<=K; j++)
								fout<<sigma2[t][j]<<"\t";
							for(j=1; j<=K; j++)
								fout<<Pimat[t][j]<<"\t";
							fout<<"\n";
						}
					}
					fout.close( );       //close file
					assert(!fout.fail( ));
					str_name.assign("2_S_paras");
					str_name.append(convertInt(gibbs_count+1));
					str_name.append(".txt");
					cout << str_name << "\n";
					fout.open(str_name.c_str(),ios::app);    // open file for appending
					assert (!fout.fail( ));   
					for(t=0; t<T; t++){
						if(state_index[t]==K-1){
							fout<<t<<"\t";
							fout<<which_two[t][0]<<"\t"<<which_two[t][1]<<"\t";
							for(j=1; j<K; j++)
								fout<<mu[t][j]<<"\t";
							for(j=1; j<K; j++)
								fout<<sigma2[t][j]<<"\t";
							for(j=1; j<K; j++)
								fout<<Pimat[t][j]<<"\t";
							fout<<"\n";
						}
					}
					fout.close( );       //close file
					assert(!fout.fail( ));
					str_name.assign("global_paras");
					str_name.append(convertInt(gibbs_count+1));
					str_name.append(".txt");
					cout << str_name << "\n";
					fout.open(str_name.c_str(),ios::app);    // open file for appending
					assert (!fout.fail( )); 
							for(j=1; j<=K; j++)
								fout<<mu0[j]<<"\t";
							fout<<"\n";
							for(j=1; j<=K; j++)
								fout<<sigma20[j]<<"\t";
							fout<<"\n";
							for(k=1; k<=K; k++){
								for(j=1; j<=K; j++)
								fout<<P[k][j]<<"\t";
								fout<<"\n";
							}
					fout.close( );       //close file
					assert(!fout.fail( ));
			}
		}
		gibbs_count ++;
	}
}


void hier_two_groups(vector<vector<double> > data, vector<vector<int> > z, int *state_index, int K,
					 double **mu, double **sigma2, double *mu0, double *sigma20, double **mu0_gibbs,
					 double ***P_gibbs, double **sigma20_gibbs, int T, HMM *phmm_temp, int *lengths,
					 double **logalpha, double **logbeta, double **loggamma, double *plogprobinit,
					 double *plogprobfinal, double ***logxi, double *tempA, double nu, double s2_init,
					 double **P, int seed, int num_gibbs, int burnin, int thin, double **Pimat,
					 HMM *phmm_old, double **count)
{	
	gsl_rng * rng = gsl_rng_alloc (gsl_rng_taus);
    gsl_rng_set(rng,seed);

    // initialize the parameters mu0,sigma20 after hierEM_uneqvar() function
	// initialize mu, sigma2 by doing Baum-Welch for each individual trace

	double *s2;
	s2 = new double[K+1];
	int t, k, j, n;
	for(k=1; k<=K; k++)
		s2[k] = s2_init;
	gsl_vector *vv = gsl_vector_alloc(K);
	gsl_permutation * perm = gsl_permutation_alloc(K);
    gsl_permutation * rank = gsl_permutation_alloc(K);
	gsl_vector *vv2 = gsl_vector_alloc(K-1);
	gsl_permutation * perm2 = gsl_permutation_alloc(K-1);
    gsl_permutation * rank2 = gsl_permutation_alloc(K-1);
	int **which_two;
	which_two = new int *[T];
	for(t=0; t<T; t++)
		which_two[t] = new int[2];
	double LL_K;
	double sum_temp;
	double *LL_temp;
	double LL_new, LL_old;
	LL_temp = new double[K];
	int kk;
	int **idxmat;
	idxmat = new int*[3];
	for(k=0; k<K; k++)
		idxmat[k] = new int[2];
	idxmat[0][0] = 1;
	idxmat[0][1] = 2;
	idxmat[1][0] = 1;
	idxmat[1][1] = 3;
	idxmat[2][0] = 2;
	idxmat[2][1] = 3;
	double **mu_temp;
	double **sig_temp;
	double **Pi_temp;
	mu_temp = new double*[K];
	for(k=0; k<K; k++)
		mu_temp[k] = new double[2];
	sig_temp = new double*[K];
	for(k=0; k<K; k++)
		sig_temp[k] = new double[2];
	Pi_temp = new double*[K];
	for(k=0; k<K; k++)
		Pi_temp[k] = new double[2];
	int num_para_iter;
	int count_para;
	double *logposterior;

	// initialize all the paremeters and fix state indicator of first observation

	for(t=0; t<T; t++){
		initializeHMM_num(phmm_temp, lengths[t], state_index[t]);
		calculateBmatrix(phmm_temp, &data[t][0]);
//		printf("%d\n",state_index[t]);
		if(state_index[t]==K){
		BaumWelch(phmm_temp, &data[t][0], lengths[t], logalpha, logbeta, loggamma, 
			     plogprobinit, plogprobfinal, logxi, tempA, vv, perm, rank, nu, &s2_init);
		which_two[t][0] = 0;
		which_two[t][1] = 0;
		}
		else if(state_index[t]==K-1){
		BaumWelch(phmm_temp, &data[t][0], lengths[t], logalpha, logbeta, loggamma, 
			     plogprobinit, plogprobfinal, logxi, tempA, vv2, perm2, rank2, nu, &s2_init);
		which_two[t][0] = 1;
		which_two[t][1] = 3;
		}
		for(k=1; k<=state_index[t]; k++){
			mu[t][k] = phmm_temp->mu[k];
//			printf("%d\t%d\t%f\n",t,k,mu[t][k]);
			sigma2[t][k] = phmm_temp->sigma2[k];
//			printf("%d\t%d\t%f\n",t,k,sigma2[t][k]);
			// fix z[0] with value of pi from EM
			Pimat[t][k] = phmm_temp->pi[k];
//			printf("%f\t",Pimat[t][k]);
			if(phmm_temp->pi[k]>0.5)
				z[t][0] = k;
		}
//		printf("\n");
//		printf("%d\t",z[t][0]);
	}
	for(t=1; t <= K; t++){
		for(k=1; k <= K; k++)
			P[t][k] = 1.0/K;
		sigma20[t] = 0.001;
	}

	int gibbs_count = 0;
	int gibbs_save = 0;
	double sum;
	int num;
	double num2;
	int Ngibbs = burnin + (num_gibbs - 1) * thin + 1;
	logposterior = new double[Ngibbs];
	for(k=0; k<Ngibbs; k++)
		logposterior[k] = 0.0;

	printf("Before Gibbs\n");

	while(gibbs_count < Ngibbs){

		printf("gibbs%d\n",gibbs_count);
		// update mu0
		if(gibbs_count < burnin)
			num_para_iter = 20;
		else
			num_para_iter = 2;

		count_para = 0;
		while(count_para < num_para_iter){
		for(k=1; k <= K; k++){
			sum = 0.0;
			num = 0;
			for(t=0; t< T; t++){					
				if(state_index[t]==K){
					sum += mu[t][k];
					num ++;
				}
				else if(state_index[t]==K-1){
					for(j=0; j<2; j++){
						if(which_two[t][j]==k){
							sum += mu[t][j+1];
							num ++;
						}
					}
				}
				else{
					printf("state index error\n");
					exit(0);
				}
			}
			mu0[k] = sum/num + gsl_ran_gaussian (rng, sqrt(sigma20[k]/num));
//			printf("%f\t%d\tmu0%d\t%f\n",sum,num,k,mu0[k]);
		}
	
		// update sigma20
		for(k=1; k <= K; k++){
			sum = 0.0;
			num = 0;
			for(t=0; t < T; t++){
				if(state_index[t]==K){	
					sum += pow(mu[t][k]-mu0[k],2.0);
					num ++;
				}
				else if(state_index[t]==K-1){
					for(j=0; j<2; j++){
						if(which_two[t][j]==k){
						sum += pow(mu[t][j+1]-mu0[k],2.0);
						num ++;
						}
					}
				}
				else{
					printf("state index error\n");
					exit(0);
				}
			}
			sigma20[k] = sum / gsl_ran_chisq(rng, 1.0 * (num - 2));
//			printf("sigma20%d\t%f\n",k,sigma20[k]);
		}

		// update s2
		for(k=1; k<=K; k++){
			num = 2;
			sum = 0.0;	
			for(t=0; t<T; t++){	
				if(state_index[t]==K){
					num ++;
					sum += 1.0/sigma2[t][k];
//					printf("%f\t",sigma2[t][k]);
				}
				else{
					for(j=0; j<2; j++){
						if(which_two[t][j]==k){
							sum += 1.0/sigma2[t][j+1];
							num ++;
						}
					}
				}
			}
			s2[k] = gsl_ran_chisq(rng,nu*num)/(sum*nu);
//			printf("sum%f\ts2%d\t%f\n",sum,k,s2[k]);
		}

		// update z
			
		for(t=0; t<T; t++){
		phmm_old->N = state_index[t];
		phmm_old->M = lengths[t];
		for(k=1; k<=state_index[t]; k++){
			phmm_old->pi[k] = Pimat[t][k];
			phmm_old->mu[k] = mu[t][k];
			phmm_old->sigma2[k] = sigma2[t][k];
		}
		if(state_index[t]==K){
			for(k=1; k<=state_index[t]; k++)
				for(j=1; j<=state_index[t]; j++)
					phmm_old->A[k][j] = P[k][j];
		}
		else if(state_index[t]==K-1){
			for(k=0; k<state_index[t]; k++){
				sum = 0.0;
				for(j=0; j<state_index[t]; j++)
					sum += P[which_two[t][k]][which_two[t][j]];
				for(j=0; j<state_index[t]; j++)
					phmm_old->A[k+1][j+1] = P[which_two[t][k]][which_two[t][j]]/sum;
			}
		}
		else{
			printf("state index error\n");
			exit(0);
		}
		calculateBmatrix(phmm_old, &data[t][0]);
		sum = 0.0;
		for(k=1; k<=state_index[t]; k++){
			if(phmm_old->pi[k]>=0.5)
				z[t][0] = k;
			if(phmm_old->pi[k]<-0.01 || phmm_old->pi[k]>1.01){
				printf("unreasonablePI\n");
				exit(0);
			}
			sum += phmm_old->pi[k];
		}
		if(sum < 0.99 || sum > 1.01){
			printf("unreasonablePIsum\n");
			exit(0);
		}
//		printf("%d%d%d%d%d\n",z[t][0],z[t][1],z[t][2],z[t][3],z[t][4]);
		samplehiddenstates(phmm_old, &data[t][0], &z[t][0], t+gibbs_count, tempA, logbeta, phmm_temp);
//		printf("%d%d%d%d%d\n",z[t][0],z[t][1],z[t][2],z[t][3],z[t][4]);
		}
		// update parameters for individual traces
		// update mu
		for(t=0; t<T; t++){
			if(state_index[t]==K){
				for(k=1; k<=K; k++){
					sum = 0.0;
					num = 0;
					for(j=0; j<lengths[t]; j++){
						if(z[t][j]==k){
							sum += data[t][j];
							num ++;
						}
					}
					mu[t][k] = (mu0[k]/sigma20[k] + sum/sigma2[t][k])/(1.0/sigma20[k]+1.0*num/sigma2[t][k]) + 
						gsl_ran_gaussian (rng, sqrt(1.0/(1.0/sigma20[k]+1.0*num/sigma2[t][k])));
//					printf("mu%d\t%d\t%f\n",t,k,mu[t][k]);
				}
			}
			else if(state_index[t]==K-1){
				for(k=1; k<K; k++){
					sum = 0.0;
					num = 0;
					for(j=0; j<lengths[t]; j++){
						if(z[t][j]==k){
							sum += data[t][j];
							num ++;
						}
					}
					mu[t][k] = (mu0[which_two[t][k-1]]/sigma20[which_two[t][k-1]] + 
						sum/sigma2[t][k])/(1.0/sigma20[which_two[t][k-1]]+1.0*num/sigma2[t][k]) + 
						gsl_ran_gaussian (rng, sqrt(1.0/(1.0/sigma20[which_two[t][k-1]]+1.0*num/sigma2[t][k])));
//					printf("mu%d\t%d\t%f\n",t,k,mu[t][k]);
				}
			}
			else{
				printf("state index error\n");
				exit(0);
			}
			// sort mu[t]
			if(state_index[t]==K){	
				for(k=0; k<K; k++)
					gsl_vector_set(vv,k,mu[t][k+1]);
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
//				for(k=1; k<=K; k++)
//				printf("mu%d\t%f\t",t,mu[t][k]);
//				printf("\n");
			}
			else if(state_index[t]==K-1){	
				for(k=0; k<K-1; k++)
					gsl_vector_set(vv2,k,mu[t][k+1]);
				gsl_sort_vector_index (perm2, vv2);
				gsl_permutation_inverse (rank2, perm2);
				for(k=0; k<K-1; k++){
					if(rank2->data[k] != k){
						for(n=0; n<lengths[t]; n++){
							if(z[t][n]==k+1){
								z[t][n] = rank2->data[k]+1;
							}
						}
					}
				}
				sort(mu[t]+1,mu[t]+K);
//				for(k=1; k<K; k++)
//				printf("mu%d\t%f\t",t,mu[t][k]);
//				printf("\n");
			}
		}
		// update sigma2
		for(t=0; t<T; t++){
			for(k=1; k<=state_index[t]; k++){
				num2 = nu;
				if(state_index[t]==K)
					sum = s2[k] * nu;
				else
					sum = s2[which_two[t][k-1]] * nu;
				for(j=0; j<lengths[t]; j++){
					if(z[t][j]==k){
						num2 += 1.0;
						sum += pow(data[t][j]-mu[t][k],2.0);
					}
				}
				sigma2[t][k] = sum/gsl_ran_chisq(rng,num2);
//				printf("sig%d\t%d\t%f\n",t,k,sigma2[t][k]);
			}
		}

		// update transition matrix

	for(k=1; k<=K; k++)
		for(j=1; j<=K; j++)
			count[k-1][j-1] = 1.0;
	for(t=0; t<T; t++){	
		if(state_index[t]==K){
		for(n=1; n<lengths[t]; n++){
			count[z[t][n-1]-1][z[t][n]-1] += 1.0;
		}
		}
		else{
//			printf("%d\t%f\t%f\t",t,Pimat[t][1],Pimat[t][2]);
		for(n=1; n<lengths[t]; n++){
//			printf("%d\twhichtwo%d\n",z[t][n-1],which_two[t][z[t][n-1]-1]-1);
			count[which_two[t][z[t][n-1]-1]-1][which_two[t][z[t][n]-1]-1] += 1.0;
		}
		}
	}

	for(k=1; k<=K; k++){
		gsl_ran_dirichlet (rng, K, count[k-1] , tempA);
		for(j=1; j<=K; j++){
			P[k][j] = tempA[j-1];
		}
	}
	count_para ++;
}
	// update state_index
	for(t=0; t<T; t++){
//		printf("%d\t%d\n",t,state_index[t]);
		if(state_index[t]==K){
			LL_K = log(gsl_ran_gaussian_pdf(data[t][0]-mu[t][z[t][0]],sqrt(sigma2[t][z[t][0]])));
			for(n=1; n<lengths[t]; n++){
				LL_K += log(gsl_ran_gaussian_pdf(data[t][n]-mu[t][z[t][n]],sqrt(sigma2[t][z[t][n]])));
				LL_K += log(P[z[t][n-1]][z[t][n]]);
			}
			for(k=1; k<=K; k++){
				LL_K += log(gsl_ran_gaussian_pdf(mu[t][k]-mu0[k],sqrt(sigma20[k])));
				LL_K += log(s2[k]*nu/2.0)*nu/2.0 - gsl_sf_lngamma(nu/2.0) - 
					s2[k]*nu/(2.0*sigma2[t][k]) -log(sigma2[t][k])*(nu/2.0+1.0);
			}
			phmm_old->N = K-1;
			phmm_old->M = lengths[t];
			for(kk=1; kk<=K; kk++){
				if(Pimat[t][idxmat[kk-1][0]]+Pimat[t][idxmat[kk-1][1]]!=0){
				    phmm_old->pi[1] = Pimat[t][idxmat[kk-1][0]]/(Pimat[t][idxmat[kk-1][0]]+Pimat[t][idxmat[kk-1][1]]);
					phmm_old->pi[2] = Pimat[t][idxmat[kk-1][1]]/(Pimat[t][idxmat[kk-1][0]]+Pimat[t][idxmat[kk-1][1]]);
				}
				else{
//					printf("warning: pi is not appropriate");
					if(Pimat[t][1]+Pimat[t][2]>0.5){
						phmm_old->pi[1] = 1.0;
						phmm_old->pi[2] = 0.0;
					}
					else{
						phmm_old->pi[1] = 0.0;
						phmm_old->pi[2] = 1.0;
					}
				}
				for(j=1; j<K; j++){
					phmm_old->mu[j] = mu[t][idxmat[kk-1][j-1]];
					phmm_old->sigma2[j] = sigma2[t][idxmat[kk-1][j-1]];
					for(k=1; k<K; k++)
					phmm_old->A[j][k] = P[idxmat[kk-1][j-1]][idxmat[kk-1][k-1]]/
					                    (P[idxmat[kk-1][j-1]][idxmat[kk-1][0]]+P[idxmat[kk-1][j-1]][idxmat[kk-1][1]]);
				}
				calculateBmatrix(phmm_old, &data[t][0]);
				samplehiddenstates(phmm_old, &data[t][0], &z[t][0], t+gibbs_count, tempA, logbeta, phmm_temp);
				LL_new = log(gsl_ran_gaussian_pdf(data[t][0]-phmm_old->mu[z[t][0]],sqrt(phmm_old->sigma2[z[t][0]])));
				for(n=1; n<lengths[t]; n++){
					LL_new += log(gsl_ran_gaussian_pdf(data[t][n]-phmm_old->mu[z[t][n]],sqrt(phmm_old->sigma2[z[t][n]])));
					LL_new += log(phmm_old->A[z[t][n-1]][z[t][n]]);
				}
				for(k=1; k<K; k++){
					LL_new += log(gsl_ran_gaussian_pdf(phmm_old->mu[k]-mu0[idxmat[kk-1][k-1]],sqrt(sigma20[idxmat[kk-1][k-1]])));
					LL_new +=  log(s2[idxmat[kk-1][k-1]]*nu/2.0)*nu/2.0 - gsl_sf_lngamma(nu/2.0) - 
						s2[idxmat[kk-1][k-1]]*nu/(2.0*phmm_old->sigma2[k]) -log(phmm_old->sigma2[k])*(nu/2.0+1);
				}
				LL_temp[kk-1] = LL_new;
				for(k=1; k<K; k++){
					Pi_temp[kk-1][k-1] = phmm_old->pi[k];
					mu_temp[kk-1][k-1] = phmm_old->mu[k];
					sig_temp[kk-1][k-1] = phmm_old->sigma2[k];
				}
				if(abs(Pi_temp[kk-1][0]+Pi_temp[kk-1][1]-1.0)>0.01){
					printf("%f\t%f\n",phmm_old->pi[1],phmm_old->pi[2]);
					printf("Pi_temp_sum_3_to_2_state_error\n");
					exit(0);
				}
			}
		}
		else if(state_index[t]==K-1){
			LL_old = log(gsl_ran_gaussian_pdf(data[t][0]-mu[t][z[t][0]],sqrt(sigma2[t][z[t][0]])));
			for(n=1; n<lengths[t]; n++){
				LL_old += log(gsl_ran_gaussian_pdf(data[t][n]-mu[t][z[t][n]],sqrt(sigma2[t][z[t][n]])));
				LL_old += log(P[which_two[t][z[t][n-1]-1]][which_two[t][z[t][n]-1]]/
					(P[which_two[t][z[t][n-1]-1]][which_two[t][0]]+P[which_two[t][z[t][n-1]-1]][which_two[t][1]]));
			}
			for(k=1; k<K; k++){
				LL_old += log(gsl_ran_gaussian_pdf(phmm_old->mu[k]-mu0[which_two[t][k-1]],sqrt(sigma20[which_two[t][k-1]])));
				LL_old +=  log(s2[which_two[t][k-1]]*nu/2.0)*nu/2.0 - gsl_sf_lngamma(nu/2.0) - 
						s2[which_two[t][k-1]]*nu/(2.0*phmm_old->sigma2[k]) -log(phmm_old->sigma2[k])*(nu/2.0+1.0);
			}
			phmm_old->N = K-1;
			phmm_old->M = lengths[t];
			kk = 0;
			while(kk < K){
				if(which_two[t][0]==idxmat[kk][0] && which_two[t][1]==idxmat[kk][1]){
					LL_temp[kk] = LL_old;
					for(k=1; k<K; k++){
						Pi_temp[kk][k-1] = Pimat[t][k];
						mu_temp[kk][k-1] = mu[t][k];
						sig_temp[kk][k-1] = sigma2[t][k];
					}
					if(abs(Pi_temp[kk][0]+Pi_temp[kk][1]-1.0)>0.01){
						printf("2_stay_Pi_temp_sum error\n");
						exit(1);
					}
					kk ++;
				}
				else{
					for(j=1; j<K; j++){
						for(k=1; k<K; k++){
							if(which_two[t][k-1]==idxmat[kk][j-1]){
								phmm_old->pi[j] = Pimat[t][k];
								phmm_old->pi[3-j] = 1.0 - phmm_old->pi[j];
								phmm_old->mu[j] = mu[t][k];
								phmm_old->sigma2[j] = sigma2[t][k];
							}
						}
						if(which_two[t][0]!=idxmat[kk][j-1] && which_two[t][1]!=idxmat[kk][j-1]){
							phmm_old->mu[j] = mu0[idxmat[kk][j-1]] + gsl_ran_gaussian (rng, sqrt(sigma20[idxmat[kk][j-1]]));
							phmm_old->sigma2[j] = nu*s2[idxmat[kk][j-1]]/gsl_ran_chisq(rng,nu);
						}
						for(k=1; k<K; k++)
							phmm_old->A[j][k] = P[idxmat[kk][j-1]][idxmat[kk][k-1]]/
					                    (P[idxmat[kk][j-1]][idxmat[kk][0]]+P[idxmat[kk][j-1]][idxmat[kk][1]]);
					}				
					calculateBmatrix(phmm_old, &data[t][0]);
					samplehiddenstates(phmm_old, &data[t][0], &z[t][0], t+gibbs_count, tempA, logbeta, phmm_temp);
				    LL_new = log(gsl_ran_gaussian_pdf(data[t][0]-phmm_old->mu[z[t][0]],sqrt(phmm_old->sigma2[z[t][0]])));				
					for(n=1; n<lengths[t]; n++){
					    LL_new += log(gsl_ran_gaussian_pdf(data[t][n]-phmm_old->mu[z[t][n]],sqrt(phmm_old->sigma2[z[t][n]])));
					    LL_new += log(phmm_old->A[z[t][n-1]][z[t][n]]);
					}
					for(k=1; k<K; k++){
					    LL_new += log(gsl_ran_gaussian_pdf(phmm_old->mu[k]-mu0[idxmat[kk][k-1]],sqrt(sigma20[idxmat[kk][k-1]])));
					    LL_new +=  log(s2[idxmat[kk][k-1]]*nu/2.0)*nu/2.0 - gsl_sf_lngamma(nu/2.0) - 
						           s2[idxmat[kk][k-1]]*nu/(2.0*phmm_old->sigma2[k]) -log(phmm_old->sigma2[k])*(nu/2.0+1);
					}
					LL_temp[kk] = LL_new;
					for(k=1; k<K; k++){
						Pi_temp[kk][k-1] = phmm_old->pi[k];
						mu_temp[kk][k-1] = phmm_old->mu[k];
						sig_temp[kk][k-1] = phmm_old->sigma2[k];
					}
					if(abs(Pi_temp[kk][0]+Pi_temp[kk][1]-1.0)>0.01){
						printf("pi_temp_sum error\n");
						exit(1);
					}
					kk ++;
				}
			}
			phmm_old->N = K;
			for(k=1; k<=K; k++){
				for(j=0; j<K-1; j++){
					if(which_two[t][j]==k){
						phmm_old->pi[k] = Pimat[t][j+1];
//						printf("%f\t",Pimat[t][j+1]);
						phmm_old->mu[k] = mu[t][j+1];
						phmm_old->sigma2[k] = sigma2[t][j+1];
					}
				}
				if(which_two[t][0]!=k && which_two[t][1]!=k){
					phmm_old->mu[k] = mu0[k] + gsl_ran_gaussian(rng,sqrt(sigma20[k]));
					phmm_old->sigma2[k] = nu*s2[k]/gsl_ran_chisq(rng,nu);
				}
				for(j=1; j<=K; j++)
					phmm_old->A[k][j] = P[k][j];
			}			
			for(k=1; k<=K; k++){
				if(which_two[t][0]!=k && which_two[t][1]!=k){
					phmm_old->pi[k] = 1.0;
					for(j=1; j<=K; j++)
						if(j!=k)
							phmm_old->pi[k] = phmm_old->pi[k] - phmm_old->pi[j];
				}
			}
			if(abs(phmm_old->pi[1]+phmm_old->pi[2]+phmm_old->pi[3]-1.0)>0.01){
				printf("Phmm_old_pi_sum_error\n");
				exit(0);
			}
			calculateBmatrix(phmm_old, &data[t][0]);
			samplehiddenstates(phmm_old, &data[t][0], &z[t][0], t+gibbs_count, tempA, logbeta, phmm_temp);
			LL_K = log(gsl_ran_gaussian_pdf(data[t][0]-phmm_old->mu[z[t][0]],sqrt(phmm_old->sigma2[z[t][0]])));
			for(n=1; n<lengths[t]; n++){
				LL_K += log(gsl_ran_gaussian_pdf(data[t][n]-phmm_old->mu[z[t][n]],sqrt(phmm_old->sigma2[z[t][n]])));
				LL_K += log(phmm_old->A[z[t][n-1]][z[t][n]]);
			}
			for(k=1; k<=K; k++){
				LL_K += log(gsl_ran_gaussian_pdf(phmm_old->mu[k]-mu0[k],sqrt(sigma20[k])));
				LL_K += log(gsl_ran_chisq_pdf(nu*s2[k]/phmm_old->sigma2[k],nu));
			}
			}
			else{
				printf("state index error\n");
				exit(0);
			}

			// determine whether to change number of states/ change which two states show
            // multinomial distribution propto exp(LL_K,LL_temp)
			num2 = gsl_ran_flat(rng, 0.0, 1.0);
			LL_new = LL_K;
			for(k=0; k<K; k++)
				LL_new = MAX(LL_new,LL_temp[k]);
			sum = exp(LL_K-LL_new);
			for(k=0; k<K; k++)
				sum += exp(LL_temp[k]-LL_new);
			sum_temp = exp(LL_K-LL_new)/sum;
//			printf("%f\t%f\t",num2,sum_temp);
			if(num2 < sum_temp){
//				printf("LL_K%d\n",K);
				if(state_index[t]==K-1){
					printf("changeindex%d\t%d\t%d\n",t,K-1,K);
					for(j=1; j<=K; j++){
						Pimat[t][j] = phmm_old->pi[j];
						mu[t][j] = phmm_old->mu[j];
						sigma2[t][j] = phmm_old->sigma2[j];
					}
					state_index[t] = K;
					if(abs(Pimat[t][1]+Pimat[t][2]+Pimat[t][3]-1.0)>0.01){
						printf("2_to_3_pi_mat_error\n");
						exit(0);
					}
				}
				logposterior[gibbs_count] += LL_K;
			}
			else{					
				k = 0;
				while(sum_temp <= num2){
					sum_temp += exp(LL_temp[k]-LL_new)/sum;
//					printf("%f\t",sum_temp);
					k ++;
				}
				k--;
				logposterior[gibbs_count] += LL_temp[k];
//				printf("Ltemp%d\n",k);
				if(state_index[t]==K){
					printf("changeindex%d\t%d\t%d\n",t,K,K-1);
					for(j=0; j<K-1; j++){
						which_two[t][j] = idxmat[k][j];
						Pimat[t][j+1] = Pi_temp[k][j];
						mu[t][j+1] = mu_temp[k][j];
						sigma2[t][j+1] = sig_temp[k][j];
					}
					state_index[t] = K-1;
					if(abs(Pimat[t][1]+Pimat[t][2]-1.0)>0.01){
						printf("3_to_2Pimaterror\n");
						exit(0);
					}
				}
				else if(which_two[t][0]!=idxmat[k][0] || which_two[t][1]!=idxmat[k][1]){
					printf("%d\t change which two states\n",t);
					for(j=0; j<K-1; j++){
						which_two[t][j] = idxmat[k][j];
						Pimat[t][j+1] = Pi_temp[k][j];
						mu[t][j+1] = mu_temp[k][j];
						sigma2[t][j+1] = sig_temp[k][j];
					}
					if(abs(Pimat[t][1]+Pimat[t][2]-1.0)>0.01){
						printf("2_to_2Pimaterror\n");
						exit(0);
					}
					state_index[t] = K-1;
				}
			}
		}
					
		ofstream fout;

		if (gibbs_count >= burnin){
			if ((gibbs_count-burnin)%thin == 0){
				for(t=0; t<T; t++)
					printf("%d\t",state_index[t]);
				printf("\n");
				for(k=1; k<=K; k++){
				mu0_gibbs[gibbs_save][k] = mu0[k];
				printf("%f\t",mu0[k]);
				sigma20_gibbs[gibbs_save][k] = sigma20[k];
				printf("%f\t",sigma20[k]);
				for(j=1; j<=K;j++){
					P_gibbs[gibbs_save][k][j] = P[k][j];
					printf("%f\t",P[k][j]);
				}
				}
				gibbs_save ++;
				printf("\n");
				fout.open("state_indicators.txt",ios::app);
				assert (!fout.fail( ));
				for(t=0; t<T; t++){
					fout<<t<<"\t";
					fout<<state_index[t]<<"\t";
					fout<<which_two[t][0]<<"\t"<<which_two[t][1]<<"\t";
				}
				fout<<"\n";
				fout.close( );       //close file
				assert(!fout.fail( ));
			}
			if((gibbs_count-burnin)%(thin*100)==0){
					std::string str_name;
					str_name.assign("3_S_paras");
					str_name.append(convertInt(gibbs_count+1));
					str_name.append(".txt");
					cout << str_name << "\n";
					fout.open(str_name.c_str(),ios::app);    // open file for appending
					assert (!fout.fail( ));
					for(t=0; t<T; t++){
						if(state_index[t]==K){
							fout<<t<<"\t";
							for(j=1; j<=K; j++)
								fout<<mu[t][j]<<"\t";
							for(j=1; j<=K; j++)
								fout<<sigma2[t][j]<<"\t";
							for(j=1; j<=K; j++)
								fout<<Pimat[t][j]<<"\t";
							fout<<"\n";
						}
					}
					fout.close( );       //close file
					assert(!fout.fail( ));
					str_name.assign("2_S_paras");
					str_name.append(convertInt(gibbs_count+1));
					str_name.append(".txt");
					cout << str_name << "\n";
					fout.open(str_name.c_str(),ios::app);    // open file for appending
					assert (!fout.fail( ));   
					for(t=0; t<T; t++){
						if(state_index[t]==K-1){
							fout<<t<<"\t";
							fout<<which_two[t][0]<<"\t"<<which_two[t][1]<<"\t";
							for(j=1; j<K; j++)
								fout<<mu[t][j]<<"\t";
							for(j=1; j<K; j++)
								fout<<sigma2[t][j]<<"\t";
							for(j=1; j<K; j++)
								fout<<Pimat[t][j]<<"\t";
							fout<<"\n";
						}
					}
					fout.close( );       //close file
					assert(!fout.fail( ));
					str_name.assign("global_paras");
					str_name.append(convertInt(gibbs_count+1));
					str_name.append(".txt");
					cout << str_name << "\n";
					fout.open(str_name.c_str(),ios::app);    // open file for appending
					assert (!fout.fail( )); 
							for(j=1; j<=K; j++)
								fout<<mu0[j]<<"\t";
							fout<<"\n";
							for(j=1; j<=K; j++)
								fout<<sigma20[j]<<"\t";
							fout<<"\n";
							for(k=1; k<=K; k++){
								for(j=1; j<=K; j++)
								fout<<P[k][j]<<"\t";
								fout<<"\n";
							}
					fout.close( );       //close file
					assert(!fout.fail( ));
			}
		}
		gibbs_count ++;
	}
	ofstream fout;
	fout.open("logposterior.txt",ios::app);    // open file for appending
	assert (!fout.fail( ));
	for(gibbs_count=0; gibbs_count < Ngibbs; gibbs_count++)
		fout<<logposterior[gibbs_count]<<"\n";
	fout.close();
	assert(!fout.fail( ));
}


