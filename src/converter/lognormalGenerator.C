#include <cstring>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <math.h>
#include <iomanip>
#include <assert.h>
#include <vector>
#include <time.h> 
#include <stdlib.h>
#include <ctime>

using namespace std;
  
#define IA 16807
#define IM 2147483647
#define AM ( 1.0 / IM )
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV ( 1 + ( IM - 1 ) / NTAB )
#define EPS 1.2e-7
#define RNMX ( 1.0 - EPS )

int main(int argc, char* argv[]);
double ran1(long *idum);
double normalCdfInverse(double cdf, double a, double b);
double normal01CdfInverse(double p);
double dpolyValue(int n, double a[], double x);
void indexx(int n, double arr[], int indx[]);
int *ivector(int nl, int nh);
void free_ivector(int *v, int nl, int nh);

int main(int argc, char* argv[])
{
  if(argc != 3){
    printf("1 output file is needed\n");
    printf("gausGenerator gauss.in cumu.dat\n");
    std::exit(1);
  }

  cout << "Start Gaussian generator \n";  
  
  std::ofstream gaussFile(argv[1]);
  if (!gaussFile) {
    cout << "Can't open output file " << argv[1] << "\n";
    exit(1);
  } 
  else{
    cout << "Open output file "<<argv[1] << "\n";                                                           
  }

   std::ofstream cumuFile(argv[2]);
  if (!cumuFile) {
    cout << "Can't open input file " << argv[2] << "\n";
  } 
  else{
    cout << "Open input file "<<argv[2] << "\n";                                                           
    }
  
  int nRandom = 200000;  
  double mean = 2.139;
  double coefficientOfVariation = 0.099751345;

  double *randomNumbers = new double [nRandom];
  //  double *sumArray = new double [nRandom];

  sleep(2);

  long randomIntegerOne = -time(NULL);
  
  for (int i=0;i<nRandom;i++){
    randomNumbers[i] = normalCdfInverse(ran1(& randomIntegerOne), mean, coefficientOfVariation * mean)*0.000001;
    if(randomNumbers[i] < 0.0000001 ){
      randomNumbers[i] = 0.0000001;
    } 
    else if
( randomNumbers[i] > 0.00002){
      randomNumbers[i] = 0.00002;
      
	}
  }    
  
  for (int i =0;i<nRandom;i++){
     gaussFile << randomNumbers[i]*0.5 <<"\n";    
   }



  //Test the indexx-function
  int* conversion = new int[nRandom];
  double *sortedRandomNumbers = new double [nRandom];

  //Sort the field using indexx 
  conversion--;
  randomNumbers--;
  sortedRandomNumbers--;
  indexx(nRandom, randomNumbers, conversion);
  for (int i =1;i<=nRandom;i++){
    sortedRandomNumbers[i] = randomNumbers[conversion[i]];    
  }
  randomNumbers++;
  sortedRandomNumbers++;
  conversion++;


  printf("Sorting done\n");

  
  double ratio;

  ratio = 0;
  for (int i =0;i<nRandom;i++){
    ratio += 1./nRandom;
    cumuFile << sortedRandomNumbers[i] <<" "<< ratio <<"\n";    
  }
 
}




double ran1(long *idum)
{
    long k;
    static long iy = 0;
    static long iv [ NTAB ];
    float temp;

    if ( * idum <= 0 || !iy ) {
        if ( -( * idum ) < 1 ) {
            * idum = 1;
        } else {
            * idum = -( * idum );
        }

        for ( int j = NTAB + 7; j >= 0; j-- ) {
            k = ( * idum ) / IQ;
            * idum = IA * ( * idum - k * IQ ) - IR * k;
            if ( * idum < 0 ) {
                * idum += IM;
            }

            if ( j < NTAB ) {
                iv [ j ] = * idum;
            }
        }

        iy = iv [ 0 ];
    }

    k = ( * idum ) / IQ;
    * idum = IA * ( * idum - k * IQ ) - IR * k;
    if ( * idum < 0 ) {
        * idum += IM;
    }

    int j = iy / NDIV;
    iy = iv [ j ];
    iv [ j ] = * idum;
    if ( ( temp = AM * iy ) > RNMX ) {
        return RNMX;
    } else {
        return temp;
    }
}

double normalCdfInverse(double cdf, double a, double b)
{
  
    double x;
    double x2;
    if ( cdf < 0.0 || 1.0 < cdf ) {
      printf("Error: in normalCdfInverse. Values outside range 0-1.\n");
      exit;
    }

    x2 = normal01CdfInverse(cdf);
    x = exp(a + b * x2);
    
    return x;
}

double  normal01CdfInverse(double p)
{
    double a [ 8 ] = {
        3.3871328727963666080,     1.3314166789178437745e+2,
        1.9715909503065514427e+3,  1.3731693765509461125e+4,
        4.5921953931549871457e+4,  6.7265770927008700853e+4,
        3.3430575583588128105e+4,  2.5090809287301226727e+3
    };
    double b [ 8 ] = {
        1.0,                       4.2313330701600911252e+1,
        6.8718700749205790830e+2,  5.3941960214247511077e+3,
        2.1213794301586595867e+4,  3.9307895800092710610e+4,
        2.8729085735721942674e+4,  5.2264952788528545610e+3
    };
    double c [ 8 ] = {
        1.42343711074968357734,     4.63033784615654529590,
        5.76949722146069140550,     3.64784832476320460504,
        1.27045825245236838258,     2.41780725177450611770e-1,
        2.27238449892691845833e-2,  7.74545014278341407640e-4
    };
    double const1 = 0.180625;
    double const2 = 1.6;
    double d [ 8 ] = {
        1.0,                        2.05319162663775882187,
        1.67638483018380384940,     6.89767334985100004550e-1,
        1.48103976427480074590e-1,  1.51986665636164571966e-2,
        5.47593808499534494600e-4,  1.05075007164441684324e-9
    };
    double e [ 8 ] = {
        6.65790464350110377720,     5.46378491116411436990,
        1.78482653991729133580,     2.96560571828504891230e-1,
        2.65321895265761230930e-2,  1.24266094738807843860e-3,
        2.71155556874348757815e-5,  2.01033439929228813265e-7
    };
    double f [ 8 ] = {
        1.0,                        5.99832206555887937690e-1,
        1.36929880922735805310e-1,  1.48753612908506148525e-2,
        7.86869131145613259100e-4,  1.84631831751005468180e-5,
        1.42151175831644588870e-7,  2.04426310338993978564e-15
    };
    double q;
    double r;
    double split1 = 0.425;
    double split2 = 5.0;
    double value;

    if ( p <= 0.0 ) {
        value = -HUGE_VAL;
        return value;
    }

    if ( 1.0 <= p ) {
        value = HUGE_VAL;
        return value;
    }

    q = p - 0.5;

    if ( fabs(q) <= split1 ) {
        r = const1 - q * q;
        value = q * dpolyValue(8, a, r) / dpolyValue(8, b, r);
    } else {
        if ( q < 0.0 ) {
            r = p;
        } else {
            r = 1.0 - p;
        }

        if ( r <= 0.0 ) {
            value = -1.0;
            printf("LocalGaussianRandomGenerator :: normal01CdfInverse - r < 0.0!");
	    exit;
        }

        r = sqrt( -log(r) );
        if ( r <= split2 ) {
            r = r - const2;
            value = dpolyValue(8, c, r) / dpolyValue(8, d, r);
        } else {
            r = r - split2;
            value = dpolyValue(8, e, r) / dpolyValue(8, f, r);
        }

        if ( q < 0.0 ) {
            value = -value;
        }
    }

    return value;
}


double  dpolyValue(int n, double a[], double x)
{
    int i;
    double value;
    value = 0.0;
    for ( i = n - 1; 0 <= i; i-- ) {
        value = value * x + a [ i ];
    }

    return value;
}



#define SWAP(a,b) itemp=(a); (a)=(b);(b)=itemp;
#define M 7
#define NSTACK 50

void indexx(int n, double arr[], int indx[])
{
  int i, indxt, ir=n, itemp, j, k, l=1;
  int jstack=0,*istack;
  double a;
  
  istack=ivector(1,NSTACK);

  for(j=1;j<=n;j++){
    indx[j] = j;
  }
  for(;;) {
    if (ir-l <M){
      for (j=l+1;j<=ir;j++){
        indxt=indx[j];
        a=arr[indxt];
        for(i=j-1;i>=l;i--){
          if (arr[indx[i]] <= a) break;
          indx[i+1]=indx[i];
        }
        indx[i+1]=indxt;
      }
      if (jstack == 0) break;
      ir=istack[jstack--];
      l=istack[jstack--];
    } else{
      k=(l+ir) >> 1;
      SWAP(indx[k],indx[l+1]);
      if(arr[indx[l]] > arr[indx[ir]]){
        SWAP(indx[l],indx[ir])
          }
      if (arr[indx[l+1]] >arr[indx[ir]]){
        SWAP(indx[l+1], indx[ir])
          }
      if(arr[indx[l]] >arr[indx[l+1]]){
        SWAP(indx[l],indx[l+1])
          }
      i=l+1;
      j=ir;
      indxt=indx[l+1];
      a=arr[indxt];
      for(;;){
        do i++; while (arr[indx[i]] <a );
        do j--; while (arr[indx[j]] >a );
        if (j < i) break;
        SWAP(indx[i],indx[j])
          }
      indx[l+1]=indx[j];
      indx[j]=indxt;
      jstack +=2;
      if (jstack >NSTACK){
        printf("Error in indexx: NSTACK too small in indexx.");
        std::exit(1);
      }
      if (ir-i+1 >= j-l){
        istack[jstack]=ir;
        istack[jstack-1]=i;
        ir=j-1;
      }
      else{
        istack[jstack]=j-1;
        istack[jstack-1]=l;
        l=i;
      }
    }
  }
  free_ivector(istack,1,NSTACK);
}


#define NR_END 1
#define FREE_ARG char*
int *ivector(int nl, int nh)
{
  int *v;
  
  v=(int *) malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
  if (!v){
    printf("allocation failure in ivector()");
    std::exit(1);
  }
  return v-nl+NR_END;
}
 

void free_ivector(int *v, int nl, int nh)
{
  free((FREE_ARG) (v+nl-NR_END));
}
