

void samplehiddenstates (HMM * phmm, double * O, int *z, int seed, 
						 double *temp, double **logbeta, HMM * phmm_old){
	
	gsl_rng * rng = gsl_rng_alloc (gsl_rng_taus);
    gsl_rng_set(rng,seed);
	double temp_max;
	double r, sum;
	int j, k;

	phmm_old->N = phmm->N;
	phmm_old->M = phmm->M;
	for(k=1; k <= phmm->N; k++){
		phmm_old->pi[k] = phmm->pi[k];
		phmm_old->mu[k] = phmm->mu[k];
		phmm_old->sigma2[k] = phmm->sigma2[k];
		for(j=1; j<=phmm->N; j++)
			phmm_old->A[k][j] = phmm->A[k][j];
		for(j=1; j<=phmm->M; j++)
			phmm_old->B[k][j] = phmm->B[k][j];
		if(phmm_old->pi[k]>=0.5)
			z[0] = k;
	}

	Backwardlog(phmm_old, phmm_old->M, logbeta, temp);
	// keep z[0] fixed
	
	for(k=1; k<phmm_old->M; k++){

			z[k] = 1;
			for(j=1; j<= phmm_old->N; j++){
			    temp[j] = log(gsl_ran_gaussian_pdf(O[k]-phmm_old->mu[j], sqrt(phmm_old->sigma2[j]))) 
				        + logbeta[k+1][j] + log(phmm_old->A[z[k-1]][j]);
			    if(j==1)
				   temp_max = temp[j];
			    else
				   temp_max = MAX (temp_max, temp[j]);
			}
			sum = 0.0;
			for(j=1; j<= phmm_old->N; j++){
				temp[j] = exp(temp[j] - temp_max);
				sum += temp[j];
			}
			for(j=1; j<= phmm_old->N; j++){
				temp[j] = temp[j]/sum;
			}
			r = gsl_rng_uniform_pos(rng);
		    z[k] = 1;
		    sum = temp[1];
		    for(j=2; j<= phmm_old->N; j++){
				if(r > sum){
				z[k]++;
				sum += temp[j];
				}
			}
	}
	gsl_rng_free (rng);
}

void samplenewparameters (HMM * phmm_old, HMM *phmm_new, int *z, double * O, int seed, 
						  double *alpha, double **count, double *tempA, double *sum, 
						  int *num, double *muprior, double *mupriorvar, double nu, 
						  double *s2_ad, gsl_vector *vv, gsl_permutation * perm, gsl_permutation * rank)
{
	gsl_rng * rng = gsl_rng_alloc (gsl_rng_taus);
    gsl_rng_set(rng,seed);

	int i, j, k;
	double s2 = *s2_ad;

	/* Sample pi */
//  keep pi fixed

	/* Sample transition matrix */

	for(j=1; j<= phmm_old->N; j++){
		for(k=1; k<= phmm_old->N; k++){
			count[j][k] = 0.0;
		}
	}
	for(i=1; i< phmm_old->M; i++)
		count[z[i-1]][z[i]] = count[z[i-1]][z[i]]+1.0;
	for(i=0; i< phmm_old->N; i++){
		for(j=0; j< phmm_old->N; j++)
			alpha[j] = 1.0+ 1.0 * count[i+1][j+1];
		gsl_ran_dirichlet (rng, phmm_old->N, alpha , tempA);
		for(j=1;  j<= phmm_new->N; j++){
			phmm_new->A[i+1][j] = tempA[j-1];
		}
	}

	/* Sample mu, sigma2 */

	for(i=0; i<phmm_old->N; i++){
		sum[i] = 0.0;
		num[i] = 0;
	}
	for(i=0; i< phmm_old->M; i++){
		sum[z[i]-1] += O[i];
		num[z[i]-1] ++;
	}
	for(j=1; j<= phmm_old->N; j++){
		phmm_new->mu[j] = (sum[j-1]/phmm_old->sigma2[j] + muprior[j]/mupriorvar[j])/(num[j-1]/phmm_old->sigma2[j] + 1/mupriorvar[j])
		+ gsl_ran_gaussian (rng, 1/sqrt(num[j-1]/phmm_old->sigma2[j] + 1/mupriorvar[j]));
	}
	for(i=0; i<phmm_old->N; i++)
		sum[i] = 0.0;
	for(i=0; i<phmm_old->M; i++)
		sum[z[i]-1] += pow(O[i]-phmm_new->mu[z[i]],2.0);
	for(j=1; j<= phmm_old->N; j++)
		tempA[j] = (sum[j-1] + nu*s2)/gsl_ran_chisq(rng, nu + num[j-1]);
	for(i=0; i<phmm_old->N; i++){
		gsl_vector_set(vv,i,phmm_new->mu[i+1]);	
	}
	gsl_sort_vector_index (perm, vv);
	gsl_permutation_inverse (rank, perm);
	for(i=1; i<=phmm_old->N; i++){
		phmm_new->sigma2[rank->data[i-1]+1] = tempA[i];
	}
	sort(phmm_new->mu+1, phmm_new->mu+phmm_old->N+1);

	/* Sample s2 */

	sum[1] = 0.0;
	for(i=0; i< phmm_old->N; i++)
		sum[1] += 1.0/phmm_new->sigma2[i+1];
	s2 = gsl_ran_gamma(rng, nu*phmm_old->N/2+1,1)*2.0/(nu*sum[1]);
	*s2_ad = s2;
	gsl_rng_free (rng);
	
}

