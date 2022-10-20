#include "mtwister.h"
#include <stdlib.h>
#include <stdio.h>
#include "ziff_newman.h"
#include "kleinberg.h"
#include "mtwister.h"

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
	for(i=0;i<N;i++) if((ptr[i]!=EMPTY)&&(ptr[i]<0)&&(-ptr[i]>cluster)&&(findroot(&i,ptr)!=r_maior)) cluster = -ptr[i];
	free(ptr);
	return cluster;
}

