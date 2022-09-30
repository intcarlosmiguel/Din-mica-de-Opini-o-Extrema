#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
//#include "kleinberg.h"
//#include "ziff_newman.h"
#include "mtwister.h"

int X(int* i,int* L){
    return *i%*L;
}
int Y(int* i,int* L){
    return *i/ *L;
}

int vizinho(int* i, int* y,int* L,int*linkagem){
    switch(*y){
        case 0: //UP//
            if (Y(i,L)!=*L-1)return *i+*L;
            else return X(i,L);
        case 1: //Down//
            if(Y(i,L)!=0) return *i-*L;
            else return (X(i,L) + *L*(*L-1));
        case 2: //Right//
            if(X(i,L)!=*L-1) return *i+1;
            else return abs(*i-(*L-1));
        case 3: //Left//
            if (X(i,L)!=0) return *i-1;
            else return *i+*L-1;
        case 4:
            return linkagem[*i];
    }
}
double sinal(double* i){
    switch(genrand64_int63()%2){
        case 0:
            return *i;
        case 1:
            return -*i;
    }
}
int power_law(int* R,double* alpha){
	double rand = genrand64_real1();
	double CN = pow(*R,1 -*alpha) - pow(2,1 -*alpha);
	if(*alpha==1)return 2*pow(*R/2,rand);
    else return pow(rand*CN + pow(2,1 -*alpha),1/(1-*alpha));
}
int long_link(int* site,double* alpha,int* L){
    double ALPHA = *alpha - 1;
	int L_ = *L-1;
    int r = power_law(&L_,&ALPHA),x,y;
	double r_ = genrand64_real1()*r;
    x = sinal(&r_);
    if(X(site,L)+x<0) x +=*L;
    if(X(site,L)+x>*L-1) x-=*L;
	double r__ = r - abs(x);
    y = sinal(&r__);
    if(Y(site,L)+y<0) y +=*L;
    if(Y(site,L)+y>*L-1) y -= *L;
    return *site + x + y**L;
}

int findroot(int* i,int* ptr){
    if (ptr[*i]<0) return *i;
    return ptr[*i] = findroot(&ptr[*i],ptr);
}

void percolate(int x,int* y,int* t,int*ptr,int*big ,int*r_maior){
	int r1,r2;
	r1 = findroot(y,ptr);
	r2 = findroot(&x,ptr);
	if(r2!=r1){
		if (ptr[r1]>ptr[r2]){
		    ptr[r2] += ptr[r1];
		    ptr[r1] = r2;
		    r1 = r2;
		} 
		else{
		    ptr[r1] += ptr[r2];
		    ptr[r2] = r1;
		}
		if(-ptr[r1]>*big){
			*big = -ptr[r1];
			*r_maior = r1;
		}
	}
}
void neighbors(int* dados,int* tam,int*ptr,int* EMPTY,int* big,int* r_maior,int* L,int* linkagem){
    int i,j,site,V = 5;
    for(i=0;i<*tam;i++){
		site = dados[i];
        ptr[site] = -1;
        for(j=0;j<V-1;j++) if(ptr[vizinho(&site,&j,L,linkagem)]!= *EMPTY) percolate((vizinho(&site,&j,L,linkagem)),&site,&i,ptr,big,r_maior);
    }
}
void permutation(int* tam,int* dados,int* seed,int* big){
    int i,j,x;
    for(i=0; i<*tam; i++){
		init_genrand64(*seed + i);
        j = i + (*tam-i)*genrand64_real2();
        x = dados[i];
        dados[i] = dados[j];
        dados[j] = x;
    }
}
int Ziff(int* dados,int* tam,int* L,int* seed,int* linkagem,int* big){
    int N = *L**L;
	int EMPTY = -(N+1);
	int r_maior = EMPTY;
	int* ptr = (int*) malloc(N*sizeof(int));
	int i = 0;
    for(i=0;i<N;i++) ptr[i] = EMPTY;
    permutation(tam,dados,seed,big);
    neighbors(dados,tam,ptr,&EMPTY,big,&r_maior,L,linkagem);
	int cluster=0;
	for(i=0;i<N;i++) if((ptr[i]!=EMPTY)&&(ptr[i]<0)&&(-ptr[i]>cluster)&&(i!=r_maior)) cluster = -ptr[i];
	free(ptr);
	return cluster;
}
double retorno(double* q, double* media,int* n,double* a){
	int m,u;
	if(*media>0) m = 1;
	else m = -1;
	if(q>0) u = 1;
	else u = -1;
	if(fabs(*q-*media)<0.000001) return *q;
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
    int n = 0, N0 = *N - 1;
    double media = 0;
	double m,u;
	int V = 5;
    while(*N/2!=i){
		//int vizinho(int* i, int* y,int* L,int*linkagem);
        for(j=0;j<V;j++) media += antes[vizinho(&i,&j,L,linkagem)];
        media = media/V;
        depois[i] = retorno(&antes[i],&media,&n,a);
		//printf("%f\n",depois[i]);
        media = 0;
        for(j=0;j<V;j++) media += antes[vizinho(&N0,&j,L,linkagem)];
        media = media/V;
        depois[N0] = retorno(&antes[N0],&media,&n,a);
        media = 0;
        i++;
        N0--;
    }
    return n;
}
double *process(double* antes,double*a,int*L,int*linkagem,int* N){
    double* depois = (double*) malloc(*N*sizeof(double));
    int n = 1,J = 0,t = 0;
    while(n!=0){
	    if((t+1)%2!=0) n = interactions(antes,depois,a,L,linkagem,N); // para i+1 ímpar
		else n = interactions(depois,antes,a,L,linkagem,N); // para i+1 par
		t++;
		if(n<0.5**N) break;
		else J = 0;
		if(J==100) break;
	}
	return depois;
}
void inicia(double* alpha,int*L,int* seed,int* linkagem,double *a,double** retorno){
	int N = *L**L;
	float r1 = 0,r2 = 0,r3 = 0,r4 = 0,r5 = 0,r6 = 0;
	double m = 0,f0,No;
	int q_e = 0,q = 0,j,t=0,pos = 0,check1,check2,big = 0; 
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
		q_e = 0,q = 0,j,t=0,pos = 0,m = 0,big = 0;
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
		printf("Rodou: %d vezes\n",j+1);
	}

	sprintf(file, "dados_%d_%d.txt",L,(int) alpha);
	fp = fopen(file ,"w");
	for (int i = 0; i < 101; i++){
		for (int j = 0; j < 6; j++){
			if(j!=5){
				if(j==0) fprintf(fp,"%f\t",retorno[j][i]);
				else fprintf(fp,"%f\t",retorno[j][i]/tam);
			}
			else fprintf(fp,"%f\n",retorno[j][i]/tam);
		}
	} 
    return 0;
}