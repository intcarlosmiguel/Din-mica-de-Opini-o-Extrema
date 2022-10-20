#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include "kleinberg.h"
#include "ziff_newman.h"
#include "mtwister.h"

double retorno(double* q, double* media,int* n,double* a){
	int m = (*media > 0) ? 1 : -1;
	int u = (*q > 0) ? 1 : -1;
	if((m==u) && (fabs(*q)-fabs(*media)<0.00001)) return *q;
    if((m**media>u**q)&&(u==m)){
        *n = *n+1;
        return *media;
    }
    if(*q>0 && *media<(1-*a)**q){
	    *n = *n+1;
	    return *media + *a**q;
    }
    if(*q<0 && *media>(1-*a)**q){
	    *n = *n+1;
	    return *media + *a**q;
    }
    return *q;
}
int interactions(double* antes,double* depois,double* a,int* L,int* linkagem,int* N ){
    int i = 0,j;
    int n = 0;
    double media = 0;
	int V = 5;
    for(i = 0;i<*N;i++){
        for(j=0;j<V;j++) media += antes[vizinho(&i,&j,L,linkagem)];
        media = media/V;
        depois[i] = retorno(&antes[i],&media,&n,a);
        media = 0;
    }
    return n;
}
double *process(double* antes,double*a,int*L,int*linkagem,int* N){
    double* depois = (double*) malloc(*N*sizeof(double));
    int n = 1,t = 0;
    while(n!=0){
	    if((t+1)%2!=0) n = interactions(antes,depois,a,L,linkagem,N); // para i+1 ímpar
		else n = interactions(depois,antes,a,L,linkagem,N); // para i+1 par
		t++;
		if(((double)n/ *N >0.9) &&(t>=25)) break;
		//else printf("%d\n",n);
		if(n<(3**N)/100) break;
		//else J = 0;
		//if(J==100) break;
	}
	return depois;
}
void inicia(double* alpha,int*L,int* seed,int* linkagem,double *a,double** retorno){
	int N = *L**L;
	double m = 0,f0,No;
	int q_e = 0,q = 0,j,big = 0; 
	for(int i=0;i<=1000;i+=10){
		int* dados = (int*) malloc(0*sizeof(int));
		double* antes = (double*) malloc(N*sizeof(double));
		double* copy = (double*) malloc(N*sizeof(double));
		f0 = (double) i/1000;
		No = (double) f0*N;
		int N0 = (int) No;
		for(j=0;j<N;j++){
			init_genrand64(*seed + i**L + j);
			double rands = genrand64_real1();
			antes[j] = -rands;
		}
		while(N0!=0){
			int u = genrand64_real1()*N;
			if(antes[u]<0){
				antes[u] = -antes[u];
				N0--;
			}
		}
		double* final = process(antes,a,L,linkagem,&N);
		for(j=0;j<N;j++){
		    if(final[j]>0){
		        q++;
		        if(final[j]>0.5){
		            q_e++;
		            dados = (int*)realloc(dados,q_e*sizeof(int));
					dados[q_e-1] = j;
				}
		    }
		    m += final[j];
		}
		m = fabs(m)/N;
		int cluster = Ziff(dados,&q_e,L,seed,linkagem,&big);
		//printf("%f\t%f\t%f\t%f\t%f\t%f\n",f0,(double)q/N ,(double)q_e/N ,(double)big/N,(double)cluster/N ,m);
		retorno[0][(int) (i/10)] = (double)f0;
		retorno[1][(int) (i/10)] += (double)q/N;
		retorno[2][(int) (i/10)] += (double)q_e/N;
		retorno[3][(int) (i/10)] += (double)big/N;
		retorno[4][(int) (i/10)] += (double)cluster/N;
		retorno[5][(int) (i/10)] += (double)m;
		free(antes);
		free(final);
		free(dados);
		free(copy);
		q_e = 0,q = 0,m = 0,big = 0;
	}
}
int main(int argc,char *argv[ ]){

    int L = atoi(argv[1]);
	int N = L*L;
	double alpha = (atof(argv[2]))/100;
	double a = atof(argv[3]);
	int seed = atoi(argv[4]);
	int tam = atoi(argv[5]);

	int i;
	int* linkagem = malloc(N*sizeof(int));
	FILE *fp = NULL;
	char file[50];

	double** retorno = (double**) malloc(6*sizeof(double*));
	for (i = 0; i < 6; i++) retorno[i] = (double*) malloc(101*sizeof(double));
	for (int j = 0; j < tam; j++){
		init_genrand64(seed*(i+1));
		for(i=0;i<N;i++)linkagem[i] = long_link(&i,&alpha,&L);
		inicia(&alpha,&L,&seed,linkagem,&a,retorno);
		if(j%10==0) printf("Rodou: %d vezes - (%f,%d)\n",j+1,alpha,L);
	}

	sprintf(file, "./dados/dados_%d_%d.txt",L,(int) alpha);
	fp = fopen(file ,"w");
	for (int i = 0; i < 101; i++){
		for (int j = 0; j < 6; j++){
			if(j!=5){
				if(j==0) fprintf(fp,"%.5f\t",retorno[j][i]);
				else fprintf(fp,"%.5f\t",retorno[j][i]/tam);
			}
			else fprintf(fp,"%.5f\n",retorno[j][i]/tam);
		}
	} 
    return 0;
}
