#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <istream>
#include <vector>
#include <cstring>
namespace smallfilter {

typedef struct { int q,k,mod; uint8_t add[81][81],mul[81][81],neg[81],cube[81],square[81]; int chi[81]; } Field;
typedef struct { Field *f; int b,c,d,mb,bbc,bc,ccbd,md,bd,mcd,dd; uint8_t left[729][3],right[729][3]; } Test;
static Field fields[4];
static Test tests[40];
static int pow3(int n) { int v=1; while(n--)v*=3; return v; }
static int rawadd(int a,int b,int k) { int v=0,p=1; for(int i=0;i<k;i++,p*=3,a/=3,b/=3)v+=((a+b)%3)*p;return v; }
static int rawmul(int a,int b,int k,int mod) {
 int aa[4]={0},bb[4]={0},m[4]={0},c[8]={0};
 for(int i=0;i<k;i++,a/=3,b/=3,mod/=3){aa[i]=a%3;bb[i]=b%3;m[i]=mod%3;}
 for(int i=0;i<k;i++)for(int j=0;j<k;j++)c[i+j]=(c[i+j]+aa[i]*bb[j])%3;
 for(int i=2*k-2;i>=k;i--)for(int j=0;j<k;j++)c[i-k+j]=(c[i-k+j]+6-c[i]*m[j])%3;
 int v=0;for(int i=k-1;i>=0;i--)v=3*v+c[i];return v;
}
static Field *field(int q,int mod) {
 int k=0;while(pow3(k)<q)k++;
 Field *f=&fields[k-1];
 if(f->q){if(f->mod!=mod)abort();return f;}
 f->q=q;f->k=k;f->mod=mod;
 for(int a=0;a<q;a++)for(int b=0;b<q;b++){f->add[a][b]=rawadd(a,b,k);f->mul[a][b]=rawmul(a,b,k,mod);}
 for(int a=0;a<q;a++){f->neg[a]=f->mul[a][2];f->square[a]=f->mul[a][a];f->cube[a]=f->mul[f->square[a]][a];f->chi[a]=a? -1:0;}
 for(int a=1;a<q;a++)f->chi[f->square[a]]=1;
 return f;
}
static int fieldnorm(Test *p,int u,int v,int w) {
 Field *f=p->f;
 #define A(a,b) f->add[(a)][(b)]
 #define M(a,b) f->mul[(a)][(b)]
 int u2=f->square[u],v2=f->square[v],w2=f->square[w];
 int z=A(f->cube[u],M(u2,A(M(p->mb,v),M(p->bbc,w))));
 z=A(z,M(u,A(A(M(p->c,v2),M(p->bc,M(v,w))),M(p->ccbd,w2))));
 z=A(z,M(p->md,f->cube[v]));z=A(z,M(p->bd,M(v2,w)));
 z=A(z,M(p->mcd,M(v,w2)));z=A(z,M(p->dd,f->cube[w]));
 return z;
 #undef A
 #undef M
}
std::vector<int> run(std::istream &in,int dim) {
 int nt;in>>nt;if(dim<1||dim>12||nt>40)throw std::runtime_error("filter dimensions");
 int leftdim=dim/2,rightdim=dim-leftdim,ls=pow3(leftdim),rs=pow3(rightdim);
 for(int j=0;j<nt;j++){
  Test *p=&tests[j];std::memset(p,0,sizeof(*p));
  int q,mod,b,c,d;in>>q>>mod>>b>>c>>d;Field *f=p->f=field(q,mod);p->b=b;p->c=c;p->d=d;
  p->mb=f->neg[b];p->bbc=f->add[f->square[b]][c];p->bc=f->neg[f->mul[b][c]];
  p->bd=f->mul[b][d];p->ccbd=f->add[f->square[c]][p->bd];p->md=f->neg[d];p->mcd=f->neg[f->mul[c][d]];p->dd=f->square[d];
  int basis[12][3];for(int i=0;i<dim;i++)for(int z=0;z<3;z++)in>>basis[i][z];
  for(int side=0;side<2;side++){
   uint8_t (*tab)[3]=side?p->right:p->left;int off=side?leftdim:0,sz=side?rs:ls;
   for(int n=1;n<sz;n++){int i=0,pow=1;while(3*pow<=n){pow*=3;i++;}for(int z=0;z<3;z++)tab[n][z]=f->add[tab[n-pow][z]][basis[off+i][z]];}
  }
 }
 std::vector<int> keep;
 for(int top=1;top<pow3(dim);top*=3)for(int code=top;code<2*top;code++){
  int lo=code%ls,hi=code/ls,oddchi=0,ok=1;
  for(int j=0;j<nt;j++){
   Test *p=&tests[j];Field *f=p->f;int xyz[3];
   for(int z=0;z<3;z++)xyz[z]=f->add[p->left[lo][z]][p->right[hi][z]];
   int ch=f->chi[fieldnorm(p,xyz[0],xyz[1],xyz[2])];
   if(!ch)continue;
   if(f->k%2==0){if(ch<0){ok=0;break;}}
   else {if(oddchi&&ch!=oddchi){ok=0;break;}oddchi=ch;}
  }
  if(ok)keep.push_back(code);
 }
 return keep;
}
}
