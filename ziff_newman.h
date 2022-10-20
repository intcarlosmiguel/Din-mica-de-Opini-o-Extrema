#ifndef ZIFF_NEWMAN_H
#define ZIFF_NEWMAN_H

int Ziff(int* dados,int* tam,int* L,int* seed,int* linkagem,int* big);
void permutation(int* tam,int* dados,int* seed,int* big);
void neighbors(int* dados,int* tam,int*ptr,int* EMPTY,int*big,int*r_maior,int*L,int*linkagem);
void percolate(int x,int* y,int* t,int*ptr,int*big ,int*r_maior);
int findroot(int* i,int* ptr);

#endif
