#include "core.hpp"
#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>
using namespace NTL;
P decode(std::string s,int q) {
 Z n(s);if(n<0)throw std::runtime_error("negative polynomial code");P a;
 for(long i=0;n>0;++i){unsigned long c=mpz_fdiv_q_ui(n.get_mpz_t(),n.get_mpz_t(),q);SetCoeff(a,i,c);}return a;
}
std::string encode(const P& a,int q) {Z n=0;for(long i=deg(a);i>=0;--i){n*=q;n+=rep(coeff(a,i));}return n.get_str();}
R::R(P a,P b):n(a),d(b) {
 if(IsZero(d))throw std::runtime_error("zero rational denominator");
 if(IsZero(n)){d=P(1);return;}P g=GCD(n,d);n/=g;d/=g;zz_p c=LeadCoeff(d);n/=c;d/=c;
}
R operator+(const R&a,const R&b){return R(a.n*b.d+b.n*a.d,a.d*b.d);}
R operator-(const R&a,const R&b){return R(a.n*b.d-b.n*a.d,a.d*b.d);}
R operator-(const R&a){return R(-a.n,a.d);}
R operator*(const R&a,const R&b){return R(a.n*b.n,a.d*b.d);}
R operator/(const R&a,const R&b){return R(a.n*b.d,a.d*b.n);}
bool operator==(const R&a,const R&b){return a.n==b.n&&a.d==b.d;}
R power(R a,int n){R b(1);while(n){if(n&1)b=b*a;a=a*a;n>>=1;}return b;}
bool operator==(const Point&a,const Point&b){if(a.identity||b.identity)return a.identity==b.identity;return a.x==b.x&&a.y==b.y;}
Curve::Curve(int q0,const std::array<P,5>& a0):q(q0),a(a0) {
 const auto&a1=a[0],&a2=a[1],&a3=a[2],&a4=a[3],&a6=a[4];
 b2=a1*a1+4*a2;b4=a1*a3+2*a4;b6=a3*a3+4*a6;
 b8=a1*a1*a6+4*a2*a6-a1*a3*a4+a2*a3*a3-a4*a4;
 c4=b2*b2-24*b4;delta=-b2*b2*b8-8*b4*b4*b4-27*b6*b6+9*b2*b4*b6;
 if(IsZero(delta))throw std::runtime_error("singular generic fibre");
 int w[]={2,4,6,8,12};for(int j=0;j<5;++j){if(deg(a[j])>w[j])throw std::runtime_error("coefficient outside K3 chart");reverse(inf[j],a[j],w[j]);}
}
Point Curve::neg(const Point& p)const {return p.identity?p:Point(p.x,-p.y-R(a[0])*p.x-R(a[2]));}
Point Curve::add(const Point& p,const Point& s)const {
 if(p.identity)return s;
 if(s.identity)return p;
 R l;
 if(p.x==s.x){if(s==neg(p))return {};l=(R(3)*p.x*p.x+R(2)*R(a[1])*p.x+R(a[3])-R(a[0])*p.y)/(R(2)*p.y+R(a[0])*p.x+R(a[2]));}
 else l=(s.y-p.y)/(s.x-p.x);
 R nu=p.y-l*p.x;R x=l*l+R(a[0])*l-R(a[1])-p.x-s.x;
 return Point(x,-(l+R(a[0]))*x-nu-R(a[2]));
}
Point Curve::mul(Point p,Z n)const {if(n<0){p=neg(p);n=-n;}Point s;while(n>0){if(mpz_odd_p(n.get_mpz_t()))s=add(s,p);n/=2;if(n>0)p=add(p,p);}return s;}
bool Curve::contains(const Point&p)const {return p.identity || p.y*p.y+R(a[0])*p.x*p.y+R(a[2])*p.y==p.x*p.x*p.x+R(a[1])*p.x*p.x+R(a[3])*p.x+R(a[4]);}
constexpr int infinity=1000000;
static int valuation(P a,const P& pi){if(IsZero(a))return infinity;int v=0;while(IsZero(a%pi)){a/=pi;++v;}return v;}
static int valuation(const R&a,const P& pi){return a.zero()?infinity:valuation(a.n,pi)-valuation(a.d,pi);}
static R chart(const R& a,int weight){if(a.zero())return a;P n,d;reverse(n,a.n);reverse(d,a.d);int e=weight+deg(a.d)-deg(a.n);if(e>=0)LeftShift(n,n,e);else LeftShift(d,d,-e);return R(n,d);}
Q Curve::local(const Point& p,const Place& pl)const {
 if(p.identity)return 0;
 const auto& co=pl.infinity?inf:a;
 R x=pl.infinity?chart(p.x,4):p.x,y=pl.infinity?chart(p.y,6):p.y;
 R a1(co[0]),a2(co[1]),a3(co[2]),a4(co[3]),a6(co[4]);
 R bb2=a1*a1+R(4)*a2,bb4=a1*a3+R(2)*a4,bb6=a3*a3+R(4)*a6;
 R bb8=a1*a1*a6+R(4)*a2*a6-a1*a3*a4+a2*a3*a3-a4*a4;
 auto val=[&](const R& f){return valuation(f,pl.pi);};
 int f2=val(R(2)*y+a1*x+a3);
 if(val(R(3)*x*x+R(2)*a2*x+a4-a1*y)<=0||f2<=0)return Q(std::max(0,-val(x)))+rat(pl.vd,6);
 if(val(bb2*bb2-R(24)*bb4)==0){Q alpha=std::min(rat(1,2),rat(f2,pl.vd));return pl.vd*alpha*(alpha-1)+rat(pl.vd,6);}
 int f3=val(R(3)*power(x,4)+bb2*power(x,3)+R(3)*bb4*x*x+R(3)*bb6*x+bb8);
 Q lambda=f3>=3*f2 ? rat(-2*f2,3):rat(-f3,4);
 return lambda+rat(pl.vd,6);
}
Q Curve::height(const Point& p)const {
 if(p.identity)return 0;
 Q h=4+deg(p.x.d)+std::max(0,int(deg(p.x.n)-deg(p.x.d)-4));
 for(auto pl:places){R x=pl.infinity?chart(p.x,4):p.x;Q pole=std::max(0,-valuation(x,pl.pi));h+=pl.degree*(local(p,pl)-rat(pl.vd,6)-pole);}
 if(h<0)throw std::runtime_error("negative height");
 return h;
}
Q Curve::pairing(const Point& p,const Point&s)const {return (height(add(p,s))-height(p)-height(s))/2;}
Q determinant(Matrix a){Q d=1;int n=a.size();for(int i=0;i<n;++i){int j=i;while(j<n&&a[j][i]==0)++j;if(j==n)return 0;if(j!=i){std::swap(a[j],a[i]);d=-d;}Q x=a[i][i];d*=x;for(int k=i+1;k<n;++k){Q c=a[k][i]/x;for(int l=i;l<n;++l)a[k][l]-=c*a[i][l];}}return d;}
std::vector<Q> solve(Matrix a,std::vector<Q>b){int n=a.size();for(int i=0;i<n;++i){int j=i;while(j<n&&a[j][i]==0)++j;if(j==n)throw std::runtime_error("singular height matrix");std::swap(a[j],a[i]);std::swap(b[j],b[i]);Q x=a[i][i];for(int k=i;k<n;++k)a[i][k]/=x;b[i]/=x;for(int k=0;k<n;++k)if(k!=i){Q c=a[k][i];for(int l=i;l<n;++l)a[k][l]-=c*a[i][l];b[k]-=c*b[i];}}return b;}
static Matrix product(const Matrix& c,const Matrix& g){int m=c.size(),n=g.size();Matrix h(m,std::vector<Q>(m));for(int i=0;i<m;++i)for(int j=0;j<m;++j)for(int k=0;k<n;++k)for(int l=0;l<n;++l)h[i][j]+=c[i][k]*g[k][l]*c[j][l];return h;}
static Matrix identity(int n){Matrix a(n,std::vector<Q>(n));for(int i=0;i<n;++i)a[i][i]=1;return a;}
Lattice build(const Curve& e,const std::vector<Point>&pts){
 Lattice out;std::vector<Q> h;std::vector<int> idx(pts.size());std::iota(idx.begin(),idx.end(),0);
 for(const auto&p:pts)h.push_back(e.height(p));
 std::stable_sort(idx.begin(),idx.end(),[&](int i,int j){return h[i]<h[j];});
 for(int i:idx){if(h[i]==0)continue;int r=out.basis.size();std::vector<Q> v;for(const auto&p:out.basis)v.push_back(e.pairing(pts[i],p));auto c=solve(out.gram,v);Q residual=h[i];for(int j=0;j<r;++j)residual-=v[j]*c[j];
  if(residual!=0){if(residual<0)throw std::runtime_error("negative residual");for(int j=0;j<r;++j)out.gram[j].push_back(v[j]);v.push_back(h[i]);out.gram.push_back(v);out.basis.push_back(pts[i]);continue;}
  Z den=1;for(const auto&x:c)mpz_lcm(den.get_mpz_t(),den.get_mpz_t(),x.get_den_mpz_t());if(den==1)continue;
  fmpz_mat_t a,u,hnf;fmpz_mat_init(a,r+1,r);fmpz_mat_init(u,r+1,r+1);fmpz_mat_init(hnf,r+1,r);
  for(int j=0;j<r;++j){fmpz_set_mpz(fmpz_mat_entry(a,j,j),den.get_mpz_t());Q scaled=c[j]*den;fmpz_set_mpz(fmpz_mat_entry(a,r,j),scaled.get_num_mpz_t());}
  fmpz_mat_hnf_transform(hnf,u,a);auto gens=out.basis;gens.push_back(pts[i]);std::vector<Point> nb;Matrix cm(r,std::vector<Q>(r));
  for(int j=0;j<r;++j){Point p;for(int k=0;k<=r;++k){Z z;fmpz_get_mpz(z.get_mpz_t(),fmpz_mat_entry(u,j,k));if(z!=0)p=e.add(p,e.mul(gens[k],z));}nb.push_back(p);for(int k=0;k<r;++k){Z z;fmpz_get_mpz(z.get_mpz_t(),fmpz_mat_entry(hnf,j,k));cm[j][k]=Q(z)/den;}}
  out.basis=nb;out.gram=product(cm,out.gram);fmpz_mat_clear(a);fmpz_mat_clear(u);fmpz_mat_clear(hnf);
 }return out;
}
static std::pair<std::vector<Q>,Matrix> gs(const Matrix& g){int n=g.size();std::vector<Q>b(n);Matrix mu(n,std::vector<Q>(n));for(int i=0;i<n;++i){b[i]=g[i][i];for(int j=0;j<i;++j){Q x=g[i][j];for(int k=0;k<j;++k)x-=mu[i][k]*mu[j][k]*b[k];mu[i][j]=x/b[j];b[i]-=mu[i][j]*mu[i][j]*b[j];}}return {b,mu};}
Q babai(Matrix g){int n=g.size(),k=1;while(k<n){for(int j=k-1;j>=0;--j){auto [b,m]=gs(g);long z=floorq(m[k][j]+rat(1,2));if(z){auto t=identity(n);t[k][j]=-z;g=product(t,g);}}
 auto [b,m]=gs(g);if(b[k]>=(rat(3,4)-m[k][k-1]*m[k][k-1])*b[k-1])++k;else{auto t=identity(n);std::swap(t[k],t[k-1]);g=product(t,g);k=std::max(1,k-1);}}
 Q sum=0;for(auto x:gs(g).first)sum+=x;return sum/4;
}
