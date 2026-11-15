#include <NTL/lzz_pXFactoring.h>
#include <NTL/lzz_pEXFactoring.h>
#include <NTL/mat_lzz_p.h>
#include <array>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <time.h>
#include "filter.hpp"
using namespace NTL;
using Poly=zz_pX;
using Base=zz_pE;
using Ext=std::vector<Base>;
using Series=std::vector<Ext>;
using Rows=std::vector<std::vector<int>>;
static void require(bool b,const char *s) {if(!b)throw std::runtime_error(s);}
static double cpu() {timespec t{};clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&t);return t.tv_sec+t.tv_nsec*1e-9;}
static Poly readpoly(std::istream &in) {
 int n,c;in>>n;require(n>=0&&n<1000,"polynomial length");Poly p;
 for(int i=0;i<n;i++){in>>c;SetCoeff(p,i,c);}return p;
}
static Poly cube(const Poly &p) {Poly q;for(long i=0;i<=deg(p);i++)SetCoeff(q,3*i,coeff(p,i));return q;}
static Poly normpoly(const std::array<Poly,3>& a,const std::array<Poly,3>& p) {
 const auto &u=a[0],&v=a[1],&w=a[2],&b=p[0],&c=p[1],&d=p[2];
 Poly u2=sqr(u),v2=sqr(v),w2=sqr(w);
 return cube(u)-b*u2*v+(sqr(b)+c)*u2*w+c*u*v2-b*c*u*v*w+(sqr(c)+b*d)*u*w2
        -d*cube(v)+b*d*v2*w-c*d*v*w2+sqr(d)*cube(w);
}
static bool square_constant(Poly f) {
 if(IsZero(f))return false;
 MakeMonic(f);if(deg(f)%2)return false;
 auto fs=SquareFreeDecomp(f);
 for(long i=0;i<fs.length();i++)if(fs[i].b%2)return false;
 return true;
}
struct Tower {
 int r,bd;zz_pEX mod;
 explicit Tower(const zz_pEX &h):r(deg(h)),bd(deg(zz_pE::modulus())),mod(h){}
 Ext zero()const{return Ext(r);}
 Ext integer(int a)const{Ext x(r);x[0]=a;return x;}
 Ext base(const Base &a)const{Ext x(r);x[0]=a;return x;}
 Ext root()const{Ext x(r);if(r==1)x[0]=-coeff(mod,0);else x[1]=1;return x;}
 bool iszero(const Ext &a)const{for(const auto &x:a)if(!IsZero(x))return false;return true;}
 Ext add(const Ext&a,const Ext&b)const{Ext c(r);for(int i=0;i<r;i++)c[i]=a[i]+b[i];return c;}
 Ext neg(const Ext&a)const{Ext c(r);for(int i=0;i<r;i++)c[i]=-a[i];return c;}
 Ext mul(const Ext&a,const Ext&b)const{
  Ext c(2*r-1);for(int i=0;i<r;i++)for(int j=0;j<r;j++)c[i+j]+=a[i]*b[j];
  for(int i=2*r-2;i>=r;i--)for(int j=0;j<r;j++)c[i-r+j]-=c[i]*coeff(mod,j);
  c.resize(r);return c;
 }
 Ext inverse(const Ext&a)const{
  zz_pEX p;for(int i=0;i<r;i++)SetCoeff(p,i,a[i]);
  auto q=InvMod(p,mod);Ext b(r);for(int i=0;i<r;i++)b[i]=coeff(q,i);return b;
 }
 std::vector<int> flatten(const Ext&a)const{
  std::vector<int> v;for(const auto &x:a)for(int j=0;j<bd;j++)v.push_back(rep(coeff(rep(x),j)));return v;
 }
 Series series(int n)const{return Series(n,zero());}
 Series constant(const Ext&a,int n)const{auto v=series(n);v[0]=a;return v;}
 Series add(const Series&a,const Series&b)const{auto c=series(a.size());for(size_t i=0;i<a.size();i++)c[i]=add(a[i],b[i]);return c;}
 Series neg(const Series&a)const{auto c=series(a.size());for(size_t i=0;i<a.size();i++)c[i]=neg(a[i]);return c;}
 Series mul(const Series&a,const Series&b)const{
  auto c=series(a.size());for(size_t i=0;i<a.size();i++)for(size_t j=0;i+j<a.size();j++)c[i+j]=add(c[i+j],mul(a[i],b[j]));return c;
 }
 Series inverse(const Series&a)const{
  auto b=series(a.size());b[0]=inverse(a[0]);
  for(size_t i=1;i<a.size();i++){auto c=zero();for(size_t j=1;j<=i;j++)c=add(c,mul(a[j],b[i-j]));b[i]=neg(mul(c,b[0]));}return b;
 }
 Series polynomial(const Poly&p,const Series&t)const{
  auto a=series(t.size());for(long i=deg(p);i>=0;i--){a=mul(a,t);a[0]=add(a[0],integer(rep(coeff(p,i))));}return a;
 }
 Series cubic(const Series&x,const Series&b,const Series&c,const Series&d)const{
  auto xc=series(x.size());for(size_t i=0;3*i<x.size();i++)xc[3*i]=mul(mul(x[i],x[i]),x[i]);
  return add(add(xc,mul(b,mul(x,x))),add(mul(c,x),d));
 }
};
struct PlaceData {int degree,m,steps;Rows rows;};
struct Input {
 std::string key;int g,dim;std::array<Poly,3> pars;Poly den;std::vector<std::array<Poly,3>> basis;
 std::vector<int> keep;std::set<int> squares;std::vector<PlaceData> exceptional;
};
static Input readinput(std::istream &in,bool exhaustive) {
 Input a;std::string magic;in>>magic>>a.key>>a.g>>a.dim;require(magic=="D2CPP1","input version");
 require(a.g>=0&&a.g<=10&&a.dim==a.g+1,"genus and dimension");
 for(auto &p:a.pars)p=readpoly(in);
 a.den=readpoly(in);
 a.basis.resize(a.dim);for(auto &v:a.basis)for(auto &p:v)p=readpoly(in);
 a.keep=smallfilter::run(in,a.dim);
 if(exhaustive){a.keep.clear();for(int top=1;top<smallfilter::pow3(a.dim);top*=3)for(int n=top;n<2*top;n++)a.keep.push_back(n);}
 int n,c;in>>n;for(int i=0;i<n;i++){in>>c;a.squares.insert(c);}
 in>>n;for(int i=0;i<n;i++){
  PlaceData p;in>>p.degree>>p.m>>p.steps;p.rows.assign(p.degree*p.steps,std::vector<int>(a.dim));
  for(auto &row:p.rows)for(auto &v:row)in>>v;
  a.exceptional.push_back(p);
 }
 require(bool(in),"truncated input");return a;
}
static std::vector<int> digits(int code,int dim){std::vector<int> v(dim);for(auto &x:v){x=code%3;code/=3;}return v;}
static int dot(const std::vector<int>&a,const std::vector<int>&b){int s=0;for(size_t i=0;i<a.size();i++)s+=a[i]*b[i];return s%3;}
static int rank3(const Rows &rows,int dim) {
 mat_zz_p m;m.SetDims(rows.size(),dim);for(size_t i=0;i<rows.size();i++)for(int j=0;j<dim;j++)m[i][j]=rows[i][j];return gauss(m);
}
static bool exceptional_test(const Input&a,const std::vector<int>&cs,Rows &rows,int &zero_degree) {
 for(const auto &p:a.exceptional){
  int k=0;
  for(;k<p.steps;k++){bool nz=false;for(int j=0;j<p.degree;j++)if(dot(cs,p.rows[k*p.degree+j])){nz=true;break;}if(nz)break;}
  require(k<p.steps,"exceptional expansion precision");if(k%2)return false;
  zero_degree+=k*p.degree;
  rows.insert(rows.end(),p.rows.begin(),p.rows.begin()+(k/2)*p.degree);
 }
 return true;
}
static bool ordinary_test(const Input&a,const Poly&prime,int exponent,const std::vector<int>&cs,Rows&rows,int&zero_degree) {
 zz_pE::init(prime);Base alpha=to_zz_pE(Poly(INIT_MONO,1));
 zz_pEX cub;SetCoeff(cub,3,1);for(int j=0;j<3;j++)SetCoeff(cub,2-j,eval(conv<zz_pEX>(a.pars[j]),alpha));
 auto factors=CanZass(cub);
 for(long j=0;j<factors.length();j++){
  require(factors[j].b==1,"ordinary cubic ramification");Tower field(factors[j].a);
  int n=exponent/field.r+1;auto t=field.constant(field.base(alpha),n);if(n>1)t[1]=field.integer(1);
  auto b=field.polynomial(a.pars[0],t),c=field.polynomial(a.pars[1],t),d=field.polynomial(a.pars[2],t);
  auto x=field.constant(field.root(),n);
  auto invder=field.inverse(field.add(field.mul(field.integer(2),field.mul(b[0],x[0])),c[0]));
  for(int k=1;k<n;k++)x[k]=field.neg(field.mul(field.cubic(x,b,c,d)[k],invder));
  auto x2=field.mul(x,x);auto iden=field.inverse(field.polynomial(a.den,t));
  std::vector<Series> bs;
  for(const auto &v:a.basis){
   auto e=field.add(field.add(field.polynomial(v[0],t),field.mul(field.polynomial(v[1],t),x)),field.mul(field.polynomial(v[2],t),x2));
   bs.push_back(field.mul(e,iden));
  }
  int val=0;
  for(;val<n;val++){
   auto z=field.zero();for(int i=0;i<a.dim;i++)if(cs[i])z=field.add(z,field.mul(field.integer(cs[i]),bs[i][val]));
   if(!field.iszero(z))break;
  }
  require(val<n,"ordinary expansion precision");if(val%2)return false;
  zero_degree+=val*field.r*field.bd;
  for(int k=0;k<val/2;k++){
   std::vector<std::vector<int>> flat;for(const auto &v:bs)flat.push_back(field.flatten(v[k]));
   for(int z=0;z<field.r*field.bd;z++){std::vector<int> row(a.dim);for(int i=0;i<a.dim;i++)row[i]=flat[i][z];rows.push_back(row);}
  }
 }
 return true;
}
struct Answer {int d2,filtered,normpasses,evenpasses;double exacttime;};
static Answer calculate(const Input&a,std::ostream *details) {
 double start=cpu();Poly bad=a.den*(sqr(a.pars[0])*sqr(a.pars[1])-cube(a.pars[1])-cube(a.pars[0])*a.pars[2]);
 // In characteristic 3 the cubic discriminant is b^2 c^2 - c^3 - b^3 d.
 long num=1,den=1;int normpasses=0,evenpasses=0;
 for(int code:a.keep){
  if(a.squares.count(code))continue;
  auto cs=digits(code,a.dim);std::array<Poly,3> f;
  for(int i=0;i<a.dim;i++)if(cs[i])for(int j=0;j<3;j++)f[j]+=a.basis[i][j]*cs[i];
  Poly nrm=normpoly(f,a.pars);if(!square_constant(nrm*a.den))continue;normpasses++;
  Rows rows;int zeros=0;
  if(!exceptional_test(a,cs,rows,zeros))continue;
  MakeMonic(nrm);auto factors=CanZass(nrm);bool ok=true;
  for(long i=0;i<factors.length();i++){
   if(IsZero(bad%factors[i].a))continue;
   if(!(ordinary_test(a,factors[i].a,factors[i].b,cs,rows,zeros))){ok=false;break;}
  }
  if(!ok)continue;
  // Omitted exceptional places have degree greater than g and cannot occur in an even zero divisor.
  if(zeros!=2*a.g)continue;
  int ell=a.dim-rank3(rows,a.dim);require(ell>=1,"effective half divisor");
  long weightden=smallfilter::pow3(ell)-1;long newnum=num*weightden+2*den,newden=den*weightden;
  long g=std::gcd(newnum,newden);num=newnum/g;den=newden/g;evenpasses++;
  if(details)*details<<a.key<<','<<code<<','<<ell<<'\n';
 }
 require(den==1&&num>0&&(num&(num-1))==0,"weighted count must be a power of 2");
 int d2=0;while((1L<<d2)<num)d2++;
 return {d2,int(a.keep.size()),normpasses,evenpasses,cpu()-start};
}
#include "prepare.hpp"
int main(int argc,char **argv) {
 if(argc==2&&std::string(argv[1])=="--stream"){
  zz_p::init(3);SetSeed(ZZ(20261114));std::string line;
  std::cout<<"key,genus,d2,status,total_cpu\n"<<std::setprecision(9);
  while(std::getline(std::cin,line)){
   if(line.empty())continue;
   std::istringstream in(line),label(line);std::string key;label>>key;double tm=cpu();
   try{auto a=native_prepare(in);auto filters=filter_text(a);std::istringstream fi(filters);
    a.keep=smallfilter::run(fi,a.dim);auto v=calculate(a,nullptr);
    std::cout<<a.key<<','<<a.g<<','<<v.d2<<",ok,"<<cpu()-tm<<std::endl;
   }catch(const std::exception&e){std::cerr<<key<<": "<<e.what()<<'\n';
    std::cout<<key<<",,,fallback,"<<cpu()-tm<<std::endl;}
  }
  return 0;
 }

 try{
  zz_p::init(3);SetSeed(ZZ(20261114));bool exhaustive=false,native=false;std::string exportdir;int repeat=1;std::ofstream detail;
  std::vector<std::string> files;
  for(int i=1;i<argc;i++){
   std::string arg=argv[i];
   if(arg=="--native")native=true;
   else if(arg=="--export"&&i+1<argc)exportdir=argv[++i];
   else if(arg=="--exhaustive")exhaustive=true;
   else if(arg=="--repeat"&&i+1<argc)repeat=std::stoi(argv[++i]);
   else if(arg=="--details"&&i+1<argc){detail.open(argv[++i]);detail<<"key,code,ell\n";}
   else if(arg=="--manifest"&&i+1<argc){std::ifstream in(argv[++i]);require(bool(in),"manifest open");std::string file;while(in>>file)files.push_back(file);}
   else files.push_back(arg);
  }
  require(!files.empty()&&repeat>=1,"usage: d2 [--manifest file] [--repeat n] [--exhaustive] [--details file] inputs...");
  std::cout<<"key,genus,d2,filter_candidates,norm_passes,even_passes,preparation_cpu,filter_cpu,exact_cpu,total_cpu\n"<<std::setprecision(9);
  for(int j=0;j<repeat;j++)for(const auto &file:files){
   double tm=cpu();std::ifstream in(file);require(bool(in),"input open");Input a;double prep=0;
   if(native){a=native_prepare(in);auto filters=filter_text(a);
    if(!exportdir.empty()){std::ofstream out(exportdir+"/"+a.key+".txt");require(bool(out),"export open");writeinput(out,a,filters);}
    prep=cpu()-tm;std::istringstream fi(filters);a.keep=smallfilter::run(fi,a.dim);
    if(exhaustive){a.keep.clear();for(int top=1;top<smallfilter::pow3(a.dim);top*=3)for(int code=top;code<2*top;code++)a.keep.push_back(code);}
   }else a=readinput(in,exhaustive);
   double setup=cpu()-tm-prep;
   auto v=calculate(a,detail.is_open()?&detail:nullptr);
   std::cout<<a.key<<','<<a.g<<','<<v.d2<<','<<v.filtered<<','<<v.normpasses<<','<<v.evenpasses<<','<<prep<<','<<setup<<','<<v.exacttime<<','<<prep+setup+v.exacttime<<std::endl;
  }
 }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}
}
