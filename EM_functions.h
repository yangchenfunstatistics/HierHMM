#define MAX(x,y) ((x)>(y)?(x):(y))
#define MIN(x,y) ((x)<(y)?(x):(y))

void initializeHMM(HMM *phmm, int T, int N)
{
	int i;
	phmm->M = T;
	phmm->N = N;
	phmm->A = new double*[phmm->N+1];
	for(i=0; i<phmm->N+1; i++)
		phmm->A[i] = new double[phmm->N+1];
	phmm->mu = new double[phmm->N+1];
	phmm->pi = new double[phmm->N+1];
	phmm->sigma2 = new double[phmm->N+1];
	phmm->B = new double*[phmm->N+1];
	for(i=0; i<phmm->N+1; i++)
		phmm->B[i] = new double[T+1];

	// for hqs
/*	if(N==2){
		phmm->mu[1] = 0.1;
		phmm->mu[2] = 0.6;
	}
	if(N==3){
		phmm->mu[1] = 0.1;
		phmm->mu[2] = 0.49;
		phmm->mu[3] = 0.65;
	}
	if(N==4){
		phmm->mu[1] = 0.1;
		phmm->mu[2] = 0.45;
		phmm->mu[3] = 0.6;
		phmm->mu[4] = 0.7;
	}
*/	
	// for ffh
	if(N==2){
		phmm->mu[1] = 0.1;
		phmm->mu[2] = 0.6;
	}
	if(N==3){
		phmm->mu[1] = 0.1;
		phmm->mu[2] = 0.34;
		phmm->mu[3] = 0.67;
	}
	if(N==4){
		phmm->mu[1] = 0.08;
		phmm->mu[2] = 0.2;
		phmm->mu[3] = 0.46;
		phmm->mu[4] = 0.7;
	}/**/
	// for ftsy
/*	if(N==2){
		phmm->mu[1] = 0.1;
		phmm->mu[2] = 0.7;
	}
	if(N==3){
		phmm->mu[1] = 0.1;
		phmm->mu[2] = 0.39;
		phmm->mu[3] = 0.8;
	}
	if(N==4){
		phmm->mu[1] = 0.067;
		phmm->mu[2] = 0.2;
		phmm->mu[3] = 0.53;
		phmm->mu[4] = 0.83;
	}
*/	int j, k;
	for(j=1; j<=N; j++){
		phmm->sigma2[j] = 0.01;
		phmm->pi[j] = 1.0/N;
	}
	for(j=1; j<=N; j++){
		for(k=1; k<=N; k++)
			phmm->A[j][k] = 1.0/N;
	}
	for(j=1; j<=N; j++)
		for(k=1; k<=T; k++)
			phmm->B[j][k] = -0.5;
	phmm->logprobf = new double[1];
}

void initializeHMM_num(HMM *phmm, int T, int N)
{
	phmm->M = T;
	phmm->N = N;

	// for hqs
/*	if(N==2){
		phmm->mu[1] = 0.1;
		phmm->mu[2] = 0.6;
	}
	if(N==3){
		phmm->mu[1] = 0.1;
		phmm->mu[2] = 0.49;
		phmm->mu[3] = 0.65;
	}
	if(N==4){
		phmm->mu[1] = 0.1;
		phmm->mu[2] = 0.45;
		phmm->mu[3] = 0.6;
		phmm->mu[4] = 0.7;
	}
*/	
	// for ffh
	if(N==2){
		phmm->mu[1] = 0.1;
		phmm->mu[2] = 0.6;
	}
	if(N==3){
		phmm->mu[1] = 0.1;
		phmm->mu[2] = 0.34;
		phmm->mu[3] = 0.67;
	}
	if(N==4){
		phmm->mu[1] = 0.08;
		phmm->mu[2] = 0.2;
		phmm->mu[3] = 0.46;
		phmm->mu[4] = 0.7;
	}/**/
	// for ftsy
/*	if(N==2){
		phmm->mu[1] = 0.1;
		phmm->mu[2] = 0.7;
	}
	if(N==3){
		phmm->mu[1] = 0.1;
		phmm->mu[2] = 0.39;
		phmm->mu[3] = 0.8;
	}
	if(N==4){
		phmm->mu[1] = 0.067;
		phmm->mu[2] = 0.2;
		phmm->mu[3] = 0.53;
		phmm->mu[4] = 0.83;
	}
*/	int j, k;
	for(j=1; j<=N; j++){
		phmm->sigma2[j] = 0.01;
		phmm->pi[j] = 1.0/N;
	}
	for(j=1; j<=N; j++){
		for(k=1; k<=N; k++)
			phmm->A[j][k] = 1.0/N;
	}
}

void calculateBmatrix(HMM *phmm, double *O)
{
	int i, t;
	for(i=1; i<=phmm->N; i++)
		for(t=1; t<=phmm->M; t++)
			phmm->B[i][t] = -1.0/2.0*log(2.0*PI*phmm->sigma2[i])- 
			           pow(O[t-1] - phmm->mu[i],2.0)/(2.0*phmm->sigma2[i]);
}
void Backwardlog(HMM *phmm, int T, double **logbeta, double *temp)
{
	int i, j, t;
	double sum, maxtemp;
	for(i=1; i<=phmm->N; i++)
		logbeta[T][i] = 0.0;
	for(t=T-1; t>=1; t--){
		for(i=1; i<=phmm->N; i++){
			sum = 0.0;
			for(j=1; j<=phmm->N; j++){
				temp[j] = log(phmm->A[i][j]) + phmm->B[j][t+1] + logbeta[t+1][j];
				if(j==1) maxtemp = temp[j];
				else maxtemp = MAX(maxtemp,temp[j]);
			}
			for(j=1; j<=phmm->N; j++)
				sum += exp(temp[j] - maxtemp);
			logbeta[t][i] = log(sum) + maxtemp;
		}
	}
}
void Forwardlog(HMM *phmm, int T, double **logalpha, double *temp)
{
	int i, j, k, t;
	double sum, maxtemp;
	for(i=1; i<=phmm->N; i++)
		logalpha[1][i] = log(phmm->pi[i]) + phmm->B[i][1];
	for(t=2; t<=T; t++){
		for(i=1; i<=phmm->N; i++){
			sum = 0.0;
			for(j=1; j<=phmm->N; j++){
				temp[j] = log(phmm->A[j][i]) + phmm->B[i][t] + logalpha[t-1][j];
				if(j==1) maxtemp = temp[j];
				else maxtemp = MAX(maxtemp,temp[j]);
			}
			for(j=1; j<=phmm->N; j++)
				sum += exp(temp[j] - maxtemp);
			logalpha[t][i] = log(sum) + maxtemp;
		}
	}
	double logprobftemp = 0.0;
	for(j=1; j<=phmm->N; j++){
		if(j==1) maxtemp = logalpha[T][j];
		else maxtemp = MAX(maxtemp, logalpha[T][j]);
	}
	for(k=1; k<=phmm->N; k++)
		logprobftemp += exp(logalpha[T][k] - maxtemp);
	phmm->logprobf[1] = log(logprobftemp) + maxtemp;
}
void ComputeGamma(HMM *phmm, int T, double **logalpha, double **logbeta, double **loggamma)
{
	int i, j, k, t;
	double denominator, maxtemp;
	for(t=1; t<=T; t++){
		denominator = 0.0;
		for(j=1; j<=phmm->N; j++){
			loggamma[t][j] = logalpha[t][j] + logbeta[t][j];
			if(j==1) maxtemp = loggamma[t][j];
			else maxtemp = MAX(maxtemp, loggamma[t][j]);
		}
		for(k=1; k<=phmm->N; k++)
			denominator += exp(loggamma[t][k] - maxtemp);
		for(i=1; i<=phmm->N; i++)
			loggamma[t][i] = loggamma[t][i] - log(denominator) - maxtemp;
	}
}
void ComputeXi(HMM *phmm, int T, double **logalpha, double **logbeta, double ***logxi)
{
	int i, j, t;
	double sum, maxtemp;
	maxtemp = logalpha[1][1] + logbeta[2][1] + log(phmm->A[1][1]) + phmm->B[1][2];
	for(t=1; t<T; t++){
		sum = 0.0;
		for(i=1; i<=phmm->N; i++){
			for(j=1; j<=phmm->N; j++){
				logxi[t][i][j] = logalpha[t][i] + logbeta[t+1][j] +
					            log(phmm->A[i][j]) + phmm->B[j][t+1];
				maxtemp = MAX(maxtemp, logxi[t][i][j]);
			}
		}
		for(i=1; i<=phmm->N; i++)
			for(j=1; j<=phmm->N; j++)
				sum += exp(logxi[t][i][j] - maxtemp);
		for(i=1; i<=phmm->N; i++)
			for(j=1; j<=phmm->N; j++)
				logxi[t][i][j] = logxi[t][i][j] - log(sum) - maxtemp;
	}
}
#define DELTA 0.000001
void BaumWelch(HMM *phmm, double *O, int T, double **logalpha, double **logbeta,
			   double **loggamma, double *plogprobinit, double *plogprobfinal,
			   double ***logxi, double *temp, gsl_vector *vv, gsl_permutation * perm, 
			   gsl_permutation * rank, double nu, double *s22)
{
	int i, j, t;
	double s2 = *s22;
	double tempsum, numeratorA, denominatorA, numeratorB, denominatorB, numeratorC;
	double delta, deltaprev, logprobprev;
	deltaprev = 10e-10;
	calculateBmatrix(phmm,O);
	Forwardlog(phmm, T, logalpha, temp);
	*plogprobinit = phmm->logprobf[1];
	Backwardlog(phmm, T, logbeta, temp);
	ComputeGamma(phmm, T, logalpha, logbeta, loggamma);
	ComputeXi(phmm, T, logalpha, logbeta, logxi);
	logprobprev = phmm->logprobf[1];

	do{
		tempsum = 0.0;
		for(i=1; i<=phmm->N; i++){
			phmm->pi[i] = exp(loggamma[1][i]);
			tempsum += phmm->pi[i];
		}
		for(i=1; i<=phmm->N; i++)
			phmm->pi[i] = phmm->pi[i]/tempsum;
		for(i=1; i<=phmm->N; i++)
		{
			denominatorA = 0.0;
			for(t=1; t<T; t++)
				denominatorA += exp(loggamma[t][i]);
			for(j=1; j<=phmm->N; j++){
				numeratorA = 0.0;
				for(t=1; t<T; t++)
					numeratorA += exp(logxi[t][i][j]);
				// for ffh 0.00001 0.99999
				// for hqs ftsy
				if(numeratorA/denominatorA<10e-7){
					phmm->A[i][j] = 0.00001 + 0.99999 * numeratorA/denominatorA;
				}
				else{
					phmm->A[i][j] = numeratorA/denominatorA;
				}
			}
			denominatorB = denominatorA + exp(loggamma[T][i]);
			numeratorB = 0.0;
			numeratorC = 0.0;
			for(t=1; t<=T; t++){
				numeratorB += O[t-1] * exp(loggamma[t][i]);
				numeratorC += pow(O[t-1] - phmm->mu[i],2.0) * exp(loggamma[t][i]);
			}
			phmm->mu[i] = numeratorB/denominatorB;
			temp[i] = (numeratorC + nu*s2)/(denominatorB+nu+2.0);
		}
		for(i=0; i<phmm->N; i++){
		gsl_vector_set(vv,i,phmm->mu[i+1]);
		}
		gsl_sort_vector_index (perm, vv);
		gsl_permutation_inverse (rank, perm);
		sort(phmm->mu+1, phmm->mu+phmm->N+1);
		for(i=1; i<=phmm->N; i++){
			phmm->sigma2[rank->data[i-1]+1] = temp[i];
			tempsum += 1.0/temp[i];
		}
		s2 = phmm->N/tempsum;
		calculateBmatrix(phmm, O);
		Forwardlog(phmm, T, logalpha, temp);
		Backwardlog(phmm, T, logbeta, temp);
		ComputeGamma(phmm, T, logalpha, logbeta, loggamma);
		ComputeXi(phmm, T, logalpha, logbeta, logxi);
		delta = phmm->logprobf[1] - logprobprev;
		logprobprev = phmm->logprobf[1];
	}
	while (delta > DELTA);
	*plogprobfinal = phmm->logprobf[1];
	*s22= s2;
}


