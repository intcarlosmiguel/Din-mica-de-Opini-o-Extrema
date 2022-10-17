#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "mtwister.h"
#include <unistd.h>


#define d 2
#define V 5
double time_spent = 0.0;
int L,N;
double a;
int seed;
int EMPTY;
int r_maior;
int* linkagem;
int* dados;
int* ptr;
//int* order;
int big = 0;
int X(int i){
    return i%L;
}
int Y(int i){
    return i/L;
}
int vizinho(int i, int y){
    switch(y){
        case 0: //UP//
            if (Y(i)!=L-1)return i+L;
            else return X(i);
        case 1: //Down//
            if(Y(i)!=0) return i-L;
            else return (X(i)+L*(L-1));
        case 2: //Right//
            if(X(i)!=L-1) return i+1;
            else return abs(i-(L-1));
        case 3: //Left//
            if (X(i)!=0) return i-1;
            else return i+L-1;
        case 4:
            return linkagem[i];
    }
    printf("Erro!");
    return 0;
}
double sinal(double i){
    switch(genrand64_int63()%2){
        case 0:
            return i;
        case 1:
            return -i;
    }
    printf("Erro!");
    return 0;
}
int power_law(int R,double alpha){
	double rand = genrand64_real1();
	double CN = pow(R,1 -alpha) - pow(2,1 -alpha);
	if(alpha==1)return 2*pow(R/2,rand);
    else return pow(rand*CN + pow(2,1 -alpha),1/(1-alpha));
}
int long_link(int site,double alpha){
    double ALPHA = alpha - 1;
    int r = power_law(L-1,ALPHA),x,y;
    x = sinal(genrand64_real1()*r);
    if(X(site)+x<0) x +=L;
    if(X(site)+x>L-1) x-=L;
    y = sinal(r - abs(x));
    if(Y(site)+y<0) y +=L;
    if(Y(site)+y>L-1) y -= L;
    return site + x + y*L;
}
int findroot(int i){
    if (ptr[i]<0) return i;
    return ptr[i] = findroot(ptr[i]);
}
void percolate(int x,int y,int t){
	int r1,r2;
	r1 = findroot(y);
	r2 = findroot(x);
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
		if(-ptr[r1]>big){
			big = -ptr[r1];
			r_maior = r1;
		}
	}
}
void neighbors(int* dados,int tam){
    int i,j,site;
    for(i=0;i<tam;i++){
		site = dados[i];
        ptr[site] = -1;
        for(j=0;j<V-1;j++) if(ptr[vizinho(site,j)]!= EMPTY) percolate(vizinho(site,j),site,i);
    }
}
void orderandpointer(){
    int i;
    for(i=0;i<N;i++) ptr[i] = EMPTY;
}
void permutation(int tam){
    int i,j,x;
    for(i=0; i<tam; i++){
		init_genrand64(seed + i);
        j = i + (tam-i)*genrand64_real2();
        x = dados[i];
        dados[i] = dados[j];
        dados[j] = x;
    }
}
int Ziff(int* dados,int tam){
	EMPTY = -(N+1);
	r_maior = EMPTY;
	ptr = (int*) malloc(N*sizeof(int));
	big = 0;
    orderandpointer();
    permutation(tam);
    neighbors(dados,tam);
	int i,cluster=0;
	for(i=0;i<N;i++) if((ptr[i]!=EMPTY)&&(ptr[i]<0)&&(-ptr[i]>cluster)&&(i!=r_maior)) cluster = -ptr[i];
	free(ptr);
	return cluster;
}
double retorno(double q, double media,double *n){
	int m = (media > 0) ? 1 : -1;
	int u = (q > 0) ? 1 : -1;
	//if((m==u) && (fabs(q)-fabs(media)<0.00001) && (fabs(q)>fabs(media))) return q;
    if((m*media>u*q)&&(u==m)){
    	*n= *n + 1;
        return media;
    }
    if(q>0 && media<(1-a)*q){
    	*n= *n + 1;
	    return media + a*q;
    }
    if(q<0 && media>(1-a)*q){
    	*n= *n + 1;
	    return media + a*q;
    }
    return q;
}
double interactions(double* antes,double* depois){
    int i = 0,j;
    double media = 0;
	double n = 0.0;

	for(i=0;i<N;i++){
        media = 0;
        for(j=0;j<V;j++) media += antes[vizinho(i,j)];
        media = media/V;
        depois[i] = retorno(antes[i],media,&n);
    }

    return n/N;
}
double *process(double* antes){
    double* depois = (double*) malloc(N*sizeof(double));
    int t = 0, contador = 0;
	double n = 1.0,cont = 0;
    while(n!=0){
	    if((t+1)%2!=0) n = interactions(antes,depois); // para i+1 ímpar
		else n = interactions(depois,antes); // para i+1 par
		t++;
		if(fabs(cont-n)<0.0001) contador++;
		if(contador==50) break;
		cont = n;
	}
	return depois;
}
void inicia(double alpha,double** retorno){

	double m = 0,f0,No;
	int q_e = 0,q = 0, j;
	
	for(int i=0;i<=1000;i+=10){
	
		dados = (int*) malloc(0*sizeof(int));
		
		double* antes = (double*) malloc(N*sizeof(double));
		double* copy = (double*) malloc(N*sizeof(double));
		
		f0 = (double) i/1000;
		No = (double) f0*N;
		int N0 = (int) No;
		
		for(j=0;j<N;j++){
			init_genrand64(seed + i*L + j);
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
		
		double* final = process(antes);
		
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
		int cluster = Ziff(dados,q_e);
		
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
		
		q_e = 0,q = 0,m = 0;
	}
}

void create_file(double alpha,int L,int tam,double** retorno){
	FILE *fp = NULL;
	char file[20];

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
}

int main(int argc,char *argv[ ]){
	// omp_set_num_threads(4);
	L = atoi(argv[1]);
	double alpha = (atof(argv[2]))/100;
	a = (atof(argv[3]))/100;
	seed = atoi(argv[4]);
	int tam = atoi(argv[5]);
	
	N = L*L;

	double** retorno = (double**) malloc(6*sizeof(double*));
	for (int i = 0; i < 6; i++) retorno[i] = (double*) malloc(101*sizeof(double));
	
	for (int j = 0; j < tam; j++){
		init_genrand64(seed*(j+1));
		linkagem = (int*) malloc(N*sizeof(int));
		for(int i=0;i<N;i++) linkagem[i] = long_link(i,alpha);
		inicia(alpha,retorno);
		free(linkagem);
		printf("Rodou: %d\n",j+1);
	}

	create_file(alpha,L,tam,retorno);
	
    return 0;
}
