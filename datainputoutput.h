
#define max_length 1000

/* Function to read data
T_length: number of traces
L_fret: total number of observations (all traces)
Length: vector for lengths of individual traces
fret: vector of all the FRET values from all traces
dir_length: directory of the txt file storing Length
dir_fret: directory of the txt file storing fret
*/

int readdata(int T_length, int L_fret, int *Length, double *fret, const char* dir_length, const char* dir_fret){

	 int i=0;
	 int j=0;
     char str[max_length];

     // read lengths

     fstream file_op(dir_length,ios::in);
     while(!file_op.eof() && i<T_length) {
     file_op.getline(str,max_length); 
     if(str!=NULL){
		 Length[i]=atoi(str);
//	     printf("%d\n",Length[i]);
	 }
	 i++;
	 }

	file_op.close();
    

    // read Fret data

	i=0;
	fstream file_opf(dir_fret,ios::in);
     while(!file_opf.eof() && i<L_fret) {
     file_opf.getline(str,max_length); 
     if(str!=NULL){
		 fret[i]=atof(str);
	 }
	 i++;
	 }
	file_opf.close();
	return(1);
}


/* Funtion for reading a vector (index) from directory */

int readindex(int T, int * index, const char* directory){
 
	 int i=0;
     char str[max_length];

     fstream file_op(directory,ios::in);
     while(!file_op.eof() && i<T) {
     file_op.getline(str,max_length); 
     if(str!=NULL){
		 index[i]=atoi(str);
	 }
	 i++;
	 }

	file_op.close();
	return(i-1);
}

/* Put the data in the form of vector DATA[i][j] means j+1 observation from trace i+1 */
/* T is number of traces, length is the vector for individual lengths */

vector<vector<double> > putinvector(int *length, double *fret, int T){

	vector<vector<double> > DATA;
	int i, j;
	for(i=0; i<T; i++)
	{
     DATA.push_back(vector<double>());
	 for(j=0; j<length[i]; j++) DATA[i].push_back(0);
	}
	int k=0;
	i=0;
	while(i<T){
		j=0;
		while(j<length[i]){
			DATA[i][j]=fret[k];
			k++;
			j++;
		}
		i++;
	}
    return(DATA);
}

/* Function for output into file */

#include <assert.h>
int appendtofile(int length, double * mu, double * sigma2, const char * filename)
{
     ofstream fout;
     fout.open(filename,ios::app);    // open file for appending
     assert (!fout.fail( ));     

	 if(length==3){
         fout<<mu[1]<<" "<<mu[2]<<" "<<mu[3]<<" "<<sigma2[1]<<" "<<sigma2[2]<<" "<<sigma2[3]<<endl;
     }
	 if(length==2){
         fout<<mu[1]<<" "<<mu[2]<<" "<<sigma2[1]<<" "<<sigma2[2]<<endl;
     }
 	 if(length==4){
         fout<<mu[1]<<" "<<mu[2]<<" "<<mu[3]<<" "<<mu[4]<<" "<<sigma2[1]<<" "<<sigma2[2]<<" "<<sigma2[3]<<" "<<sigma2[4]<<endl;
     }
     fout.close( );       //close file
     assert(!fout.fail( ));
     return 0;
}
string convertInt(int number)
{
   stringstream ss;//create a stringstream
   ss << number;//add number to the stream
   return ss.str();//return a string with the contents of the stream
}

