#include <stdio.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <cmath>
#include <assert.h>

//stantard containers
#include <string>
#include <deque>
#include <queue>
#include <list>
#include <vector>
#include <set>
#include <map>
#include <iterator>
#include <stack>
//algorithm
#include <algorithm>
//stream
#include <iostream>
#include <fstream>
#include <sstream>
//limits
#include <limits>

#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_heapsort.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_sf_gamma.h>

using namespace std;
#include "datainputoutput.h"

#define Thqs 142
#define Lhqs 238651

//#define Thqs 208
//#define Lhqs 259623

//#define Thqs 138
//#define Lhqs 147770

#define PI 3.1415926

int HQSlength[Thqs];
double HQSfret[Lhqs];
vector<vector<double> > HQS;

class HMM {
public:
	int N;		    /* number of states */
	int M; 		    /* number of observations */
	double	**A;	/* transition matrix */
	double	**B;	/* LOG probability density of observed values at different states*/
	double  *mu;
	double  *sigma2;
	double	*pi;	/* starting point distribution */
	double * logprobf;   /* marginal log likelihood p(observation|parameters), from forward probability */
};  

#include "EM_functions.h"
#include "Gibbs_functions.h"
#include "hier_functions.h"
#include "hier_two_groups.h"

int main(){

	int x;
	x=readdata(Thqs,Lhqs,HQSlength, HQSfret,
		"Ffhlengths.txt",
		"Ffhfret.txt");
	if(x==1){
		printf("Sucessfully read data Ffh.\n");
	}
    HQS = putinvector(HQSlength, HQSfret, Thqs);

	int K = 3;
	int state_index[Thqs];
	readindex(Thqs, state_index, "num_state_ffh.txt");
	int Thqstemp = 0;
	int i, j, k;
	for(x=0; x<Thqs; x++){
		if(state_index[x]>K)
			state_index[x] = K;
	}
	for(x=0; x<Thqs; x++){
		if(state_index[x]==K){
			Thqstemp++;
		}
	}

	int *HQSlengthtemp;
	HQSlengthtemp =  new int[Thqstemp];
	vector<vector<double> > HQStemp;
	i = 0;
	for(x=0; x<Thqs; x++){
		if(state_index[x]==K){
			HQSlengthtemp[i] = HQSlength[x];
			HQStemp.push_back(vector<double>());
			for(j=0; j<HQSlengthtemp[i]; j++)
				HQStemp[i].push_back(HQS[x][j]);
//			printf("%f\t",HQStemp[i][0]);
			i++;
		}
	}


	int num_gibbs = 2000;
	int burnin = 2000;
	int thin = 10;
	double nu = 1.0;
	double s2 = 0.01;
	int mlength;
	int max_num_state = 4;
	HMM **hmm;
	double * temp;
	temp = new double [max_num_state+1];
	mlength = HQSlength[0];
	for(i=1; i<Thqs; i++)
		mlength = MAX(mlength,HQSlength[i]);
	double **count;
	double *sum;
	double *row_P;
	double *mu0;
	double *sigma20;

	double ***logbeta_vec;
	double ***logalpha_vec;
	double ***loggamma_vec;
	double ****logxi_vec;


   /******** Allocate Space ******/

	logbeta_vec = new double **[Thqs];
	for(i=0; i<Thqs; i++){
		logbeta_vec[i] = new double*[mlength+1];
		for(j=0; j<mlength+1; j++)
			logbeta_vec[i][j] = new double[max_num_state+1];
	}
	logalpha_vec = new double **[Thqs];
	for(i=0; i<Thqs; i++){
		logalpha_vec[i] = new double*[mlength+1];
		for(j=0; j<mlength+1; j++)
			logalpha_vec[i][j] = new double[max_num_state+1];
	}
	loggamma_vec = new double **[Thqs];
	for(i=0; i<Thqs; i++){
		loggamma_vec[i] = new double*[mlength+1];
		for(j=0; j<mlength+1; j++)
			loggamma_vec[i][j] = new double[max_num_state+1];
	}
	logxi_vec = new double ***[Thqs];
	for(i=0; i<Thqs; i++){
		logxi_vec[i] = new double**[mlength+1];
		for(j=0; j<mlength+1; j++){
			logxi_vec[i][j] = new double*[max_num_state+1];
			for(k=0; k<max_num_state+1; k++)
				logxi_vec[i][j][k] = new double[max_num_state+1];
		}
	}
	sum = new double[max_num_state+1];
	row_P = new double[max_num_state+1];
	mu0 = new double[max_num_state+1];
	sigma20 = new double[max_num_state+1];
	count = new double*[max_num_state+1];
	for(i=0; i<max_num_state+1; i++)
		count[i] = new double[max_num_state+1];
	hmm = new HMM*[Thqs];
	for(i=0; i<Thqs; i++){
		hmm[i] = new HMM[1];
		initializeHMM(hmm[i],mlength,max_num_state);
	}
/*
	printf("beforeEM\n");
	hierEM_uneqvar(HQSlengthtemp, Thqstemp, HQStemp, mu0, sigma20, logbeta_vec, logalpha_vec, 
		   loggamma_vec, logxi_vec, sum, count, row_P, nu, s2, K, hmm, 
		   "EM_Pmu0sigam20.txt", "EM_indtraces.txt");
*/
	delete[] row_P;

	for(i=0; i<Thqs; i++){
		delete[] hmm[i];
	}
	delete[] hmm;

	for(i=0; i<Thqs; i++){
		for(j=0; j<mlength+1; j++)
			delete[] logbeta_vec[i][j];
		delete[] logbeta_vec[i];
	}
	delete[] logbeta_vec;

	for(i=0; i<Thqs; i++){
		for(j=0; j<mlength+1; j++)
			delete[] logalpha_vec[i][j];
		delete[] logalpha_vec[i];
	}
	delete[] logalpha_vec;

	for(i=0; i<Thqs; i++){
		for(j=0; j<mlength+1; j++)
			delete[] loggamma_vec[i][j];
		delete[] loggamma_vec[i];
	}
	delete[] loggamma_vec;

	for(i=0; i<Thqs; i++){
		for(j=0; j<mlength+1; j++){
			for(k=0; k<max_num_state+1; k++)
				delete[] logxi_vec[i][j][k];
			delete[] logxi_vec[i][j];
		}
		delete[] logxi_vec[i];
	}
	delete[] logxi_vec;

	printf("Deleteall EM parameters\n");

	double ** logbeta;
	double ** logalpha;
	double ** loggamma;
	double *** logxi;
	int *num;
	double **P;
	double **mu;
	double **sigma2_unequalvar;
	vector<vector<int> > Z;
	vector<vector<int> >z;
	double **Pimat;
	double ***P_gibbs;
	double ***mu_gibbs;
	double ***sigma2_gibbs_unequalvar;
	double **mu0_gibbs;
	double **sigma20_gibbs;

	num = new int[max_num_state+1];

	int seed = K;
	for(i=0; i<Thqstemp; i++)
	{
     Z.push_back(vector<int>());
	 for(j=0; j<HQSlengthtemp[i]; j++) Z[i].push_back(0);
	}
	for(i=0; i<Thqs; i++)
	{
     z.push_back(vector<int>());
	 for(j=0; j<HQSlength[i]; j++) z[i].push_back(0);
	}

	HMM * phmm_old;
	HMM * phmm_temp;
	phmm_old = new HMM[1];
	phmm_temp = new HMM[1];
	initializeHMM(phmm_old,mlength,max_num_state);
	initializeHMM(phmm_temp,mlength,max_num_state);
	double * plogprobfinal;  // final marginal log likelihood in Baum Welch
	double * plogprobinit;   // initial marginal log likelihood in Baum Welch
	plogprobfinal = new double[1];
	plogprobinit = new double[1];

	P_gibbs = new double **[num_gibbs];
	for(i=0; i<num_gibbs; i++){
		P_gibbs[i] = new double *[max_num_state+1];
		for(j=0; j<max_num_state+1; j++)
			P_gibbs[i][j] = new double [max_num_state+1];
	}
	mu_gibbs = new double **[num_gibbs];
	for(i=0; i<num_gibbs; i++){
		mu_gibbs[i] = new double *[Thqs];
		for(j=0; j<Thqs; j++)
			mu_gibbs[i][j] = new double [max_num_state+1];
	}
	sigma2_gibbs_unequalvar = new double **[num_gibbs];
	for(i=0; i<num_gibbs; i++){
		sigma2_gibbs_unequalvar[i] = new double *[Thqs];
		for(j=0; j<Thqs; j++)
			sigma2_gibbs_unequalvar[i][j] = new double [max_num_state+1];
	}
	mu0_gibbs = new double *[num_gibbs];
	for(i=0; i<num_gibbs; i++){
		mu0_gibbs[i] = new double [max_num_state+1];
	}
	sigma20_gibbs = new double *[num_gibbs];
	for(i=0; i<num_gibbs; i++){
		sigma20_gibbs[i] = new double [max_num_state+1];
	}
	P = new double*[max_num_state+1];
	for(i=0; i<max_num_state+1; i++)
		P[i] = new double[max_num_state+1];
	mu = new double*[Thqs];
	for(i=0; i<Thqs; i++)
		mu[i] = new double[max_num_state+1];
	sigma2_unequalvar = new double*[Thqs];
	for(i=0; i<Thqs; i++)
		sigma2_unequalvar[i] = new double[max_num_state+1];
	Pimat = new double*[Thqs];
	for(i=0; i<Thqs; i++)
		Pimat[i] = new double[max_num_state+1];
	logbeta = new double*[mlength+1];
	for(i=0; i<mlength+1; i++)
		logbeta[i] = new double[max_num_state+1];
	logalpha = new double*[mlength+1];
	for(i=0; i<mlength+1; i++)
		logalpha[i] = new double[max_num_state+1];
	loggamma = new double*[mlength+1];
	for(i=0; i<mlength+1; i++)
		loggamma[i] = new double[max_num_state+1];
	logxi = new double**[mlength+1];
	for(i=0; i<mlength+1; i++){
		logxi[i] = new double*[max_num_state+1];
		for(j=0; j<max_num_state+1; j++)
			logxi[i][j] = new double[max_num_state+1];
	}
	hier_two_groups(HQS, z, state_index, K, mu, sigma2_unequalvar, mu0, sigma20, mu0_gibbs,
		P_gibbs, sigma20_gibbs, Thqs, phmm_temp, HQSlength, logalpha, logbeta, 
		loggamma, plogprobinit, plogprobfinal, logxi, temp, nu, s2, P, K, num_gibbs, 
		burnin, thin, Pimat, phmm_old, count);
/*	
	hiergibbsunequalvar(P, mu, sigma2_unequalvar, K, Thqstemp,
			   HQStemp, Z, seed, count, HQSlengthtemp, mu0, sigma20,
			   sum, num, temp, s2, nu,
			   Pimat, logbeta, phmm_temp, phmm_old,
			   logalpha, loggamma, logxi,
			   plogprobinit, plogprobfinal, num_gibbs,
			   burnin, thin, P_gibbs, mu_gibbs,
			   sigma2_gibbs_unequalvar, mu0_gibbs, sigma20_gibbs);
*/
	// print gibbs results to file
	     
	ofstream fout;
    fout.open("P_gibbs.txt",ios::app);    // open file for appending
    assert (!fout.fail( ));   
	for(i=0; i<num_gibbs; i++){
		for(j=1; j<=K; j++)
			for(k=1; k<=K; k++){
				if(k != j)
				fout<<P_gibbs[i][j][k]<<"\t";
				if(k==K-1 && j==K)
				fout<<"\n";
			}
	}
    fout.close( );       //close file
    assert(!fout.fail( ));

    fout.open("mu0_gibbs.txt",ios::app);    // open file for appending
    assert (!fout.fail( ));   
	for(i=0; i<num_gibbs; i++){
		for(j=1; j<=K; j++)
			fout<<mu0_gibbs[i][j]<<"\t";
		fout<<"\n";
	}
    fout.close( );       //close file
    assert(!fout.fail( ));

	fout.open("sigma20_gibbs.txt",ios::app);    // open file for appending
    assert (!fout.fail( ));   
	for(i=0; i<num_gibbs; i++){
		for(j=1; j<=K; j++)
			fout<<sigma20_gibbs[i][j]<<"\t";
		fout<<"\n";
	}
    fout.close( );       //close file
    assert(!fout.fail( ));
/*
	fout.open("musigma2_gibbs.txt",ios::app);    // open file for appending
    assert (!fout.fail( ));   
	for(i=0; i<num_gibbs; i++){
		for(j=0; j<Thqstemp; j++){
			for(k=1; k<=K; k++)
			fout<<mu_gibbs[i][j][k]<<"\t";
			for(k=1; k<=K; k++)
			fout<<sigma2_gibbs_unequalvar[i][j][k]<<"\t";
		}
		fout<<"\n";
	}
    fout.close( );       //close file
    assert(!fout.fail( ));
*/
	return(0);
}


