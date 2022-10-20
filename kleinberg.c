#include "kleinberg.h"
#include "ziff_newman.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mtwister.h"


int X(int* i,int* L){
    return *i%*L;
}
int Y(int* i,int* L){
    return *i/ *L;
}

int vizinho(int* i, int* y,int* L,int* linkagem){
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
        default:
            return linkagem[*i];

    }
}
double sinal(double* i){
    int retorno = (genrand64_int63()%2) ? *i : -*i;
    return retorno;
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
