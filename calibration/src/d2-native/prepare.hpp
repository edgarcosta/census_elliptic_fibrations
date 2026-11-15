#include <sstream>
using Bivar = std::map<std::pair<int,int>, Ext>;

static Series smul(const Tower& f,const Series& a,const Series& b) {
 if(f.r!=1)return f.mul(a,b);
 zz_pEX aa,bb;
 for(size_t i=0;i<a.size();i++){SetCoeff(aa,i,a[i][0]);SetCoeff(bb,i,b[i][0]);}
 auto cc=MulTrunc(aa,bb,a.size());auto c=f.series(a.size());
 for(size_t i=0;i<a.size();i++)c[i][0]=coeff(cc,i);
 return c;
}
static Series sinv(const Tower& f,const Series& a) {
 if(f.r!=1)return f.inverse(a);
 zz_pEX aa;for(size_t i=0;i<a.size();i++)SetCoeff(aa,i,a[i][0]);
 auto bb=InvTrunc(aa,a.size());auto b=f.series(a.size());
 for(size_t i=0;i<a.size();i++)b[i][0]=coeff(bb,i);
 return b;
}
static void addterm(Bivar& p,int i,int j,const Ext& v,const Tower& f) {
 if(f.iszero(v))return;
 auto key=std::make_pair(i,j);auto it=p.find(key);
 if(it==p.end())p[key]=v;
 else {it->second=f.add(it->second,v);if(f.iszero(it->second))p.erase(it);}
}
static int binom3(int n,int k) {
 int v=1;
 while(n||k){int a=n%3,b=k%3;if(b>a)return 0;if(a==2&&b==1)v=v*2%3;n/=3;k/=3;}
 return v;
}
static Bivar centered(const std::array<Poly,3>& pars,const Ext& alpha,const Ext& beta,const Tower& f) {
 Bivar p;std::vector<Ext> ap(13,f.integer(1)),bp(4,f.integer(1));
 for(int i=1;i<=12;i++)ap[i]=f.mul(ap[i-1],alpha);
 for(int i=1;i<=3;i++)bp[i]=f.mul(bp[i-1],beta);
 for(int j=0;j<=3;j++){
  Poly c=j==3?Poly(1):pars[2-j];
  for(int k=0;k<=deg(c);k++)if(!IsZero(coeff(c,k)))
   for(int u=0;u<=k;u++)for(int v=0;v<=j;v++){
    int b=rep(coeff(c,k))*binom3(k,u)*binom3(j,v)%3;
    if(b)addterm(p,u,v,f.mul(f.integer(b),f.mul(ap[k-u],bp[j-v])),f);
   }
 }
 return p;
}
static Series evalbi(const Bivar& p,const Series& t,const Series& x,const Tower& f) {
 int mt=0,mx=0;for(const auto& [ij,c]:p){mt=std::max(mt,ij.first);mx=std::max(mx,ij.second);}
 std::vector<Series> tp{f.constant(f.integer(1),t.size())},xp=tp;
 for(int k=1;k<=mt;k++)tp.push_back(smul(f,tp.back(),t));
 for(int k=1;k<=mx;k++)xp.push_back(smul(f,xp.back(),x));
 auto out=f.series(t.size());
 for(const auto& [ij,c]:p){auto s=smul(f,tp[ij.first],xp[ij.second]);for(size_t k=0;k<t.size();k++)out[k]=f.add(out[k],f.mul(c,s[k]));}
 return out;
}
static Bivar derivbi(const Bivar& p,int variable,const Tower& f) {
 Bivar q;for(const auto& [ij,c]:p){int n=variable?ij.second:ij.first;if(n%3)addterm(q,ij.first-(variable==0),ij.second-(variable==1),f.mul(c,f.integer(n%3)),f);}return q;
}
struct Transform {bool vertical;Ext slope;};
struct Branch {
 Poly prime;zz_pEX h;Ext alpha,beta;Bivar equation;std::vector<Transform> path;
 bool infinity=false;int degree=0,ramification=0;
};
static Bivar liftbi(Bivar p,int r) {for(auto& [ij,c]:p)c.resize(r);return p;}
static Bivar blowup(const Bivar& p,int multiplicity,bool vertical,const Ext& slope,const Tower& f) {
 Bivar q;
 for(const auto& [ij,c]:p){
  int n=vertical?ij.first:ij.second;
  auto pow=f.integer(1);std::vector<Ext> powers{pow};for(int k=1;k<=n;k++){pow=f.mul(pow,slope);powers.push_back(pow);}
  for(int j=0;j<=n;j++){int b=binom3(n,j);if(b)addterm(q,ij.first+ij.second-multiplicity,j,f.mul(c,f.mul(f.integer(b),powers[n-j])),f);}
 }
 return q;
}
static void resolve(Branch branch,std::vector<Branch>& out,int& delta,int depth=0) {
 require(depth<30,"normalization blowup limit");Tower f(branch.h);
 int mult=1000;for(const auto& [ij,c]:branch.equation)mult=std::min(mult,ij.first+ij.second);
 require(mult>=1&&mult<=3,"plane multiplicity");
 if(mult==1){if(!branch.path.empty())branch.ramification=0;branch.degree=f.bd*f.r;out.push_back(std::move(branch));return;}
 delta+=mult*(mult-1)/2*f.bd*f.r;
 require(f.r==1,"singular point over an extended tangent field");
 zz_pEX tangent;
 for(const auto& [ij,c]:branch.equation)if(ij.first+ij.second==mult)SetCoeff(tangent,ij.second,c[0]);
 require(!IsZero(tangent),"nonzero tangent cone");
 int td=deg(tangent);MakeMonic(tangent);auto factors=CanZass(tangent);
 for(long i=0;i<factors.length();i++){
  Branch child=branch;child.h=factors[i].a;Tower ff(child.h);
  child.alpha.resize(ff.r);child.beta.resize(ff.r);child.equation=liftbi(child.equation,ff.r);
  for(auto& tr:child.path)tr.slope.resize(ff.r);
  auto slope=ff.root();child.equation=blowup(child.equation,mult,false,slope,ff);
  child.path.push_back({false,slope});resolve(std::move(child),out,delta,depth+1);
 }
 if(td<mult){
  branch.equation=blowup(branch.equation,mult,true,f.zero(),f);
  branch.path.push_back({true,f.zero()});resolve(std::move(branch),out,delta,depth+1);
 }
}
static std::pair<Series,Series> expand_branch(const Branch& b,int n) {
 zz_pE::init(b.prime);Tower f(b.h);auto t=f.series(n),x=t;
 auto it=b.equation.find({0,1});bool solve_x=it!=b.equation.end()&&!f.iszero(it->second);
 auto derivative=derivbi(b.equation,solve_x?1:0,f);
 if(n>1)(solve_x?t:x)[1]=f.integer(1);
 for(int size=2;size/2<n;size*=2){
  int m=std::min(size,n);auto tt=t,xx=x;tt.resize(m);xx.resize(m);
  auto error=evalbi(b.equation,tt,xx,f),der=evalbi(derivative,tt,xx,f);
  auto correction=smul(f,error,sinv(f,der));auto& target=solve_x?x:t;
  for(int k=0;k<m;k++)target[k]=f.add(target[k],f.neg(correction[k]));
 }
 for(auto it2=b.path.rbegin();it2!=b.path.rend();++it2){
  auto shifted=x;shifted[0]=f.add(shifted[0],it2->slope);auto prod=smul(f,t,shifted);
  if(it2->vertical){x=t;t=std::move(prod);}else x=std::move(prod);
 }
 t[0]=f.add(t[0],b.alpha);x[0]=f.add(x[0],b.beta);
 return {t,x};
}
static zz_pEX specialize(const std::array<Poly,3>& pars,const Base& alpha) {
 zz_pEX p;SetCoeff(p,3,1);for(int j=0;j<3;j++)SetCoeff(p,2-j,eval(conv<zz_pEX>(pars[j]),alpha));return p;
}
static std::vector<Branch> branches(const std::array<Poly,3>& pars,const Poly& prime,bool infinity,int& delta) {
 zz_pE::init(prime);Base alpha=to_zz_pE(Poly(INIT_MONO,1));auto fs=CanZass(specialize(pars,alpha));
 std::vector<Branch> out;
 for(long i=0;i<fs.length();i++){
  Tower f(fs[i].a);Branch b;b.prime=prime;b.h=fs[i].a;b.alpha=f.base(alpha);b.beta=f.root();b.infinity=infinity;
  b.ramification=fs[i].b;b.equation=centered(pars,b.alpha,b.beta,f);resolve(std::move(b),out,delta);
 }
 int total=0;
 for(auto& b:out){
  if(!b.ramification){auto tx=expand_branch(b,5);Tower f(b.h);int e=1;while(e<5&&f.iszero(tx.first[e]))e++;require(e<=3,"local ramification degree");b.ramification=e;}
  total+=b.ramification*b.degree;
 }
 require(total==3*deg(prime),"normalization branch degree");return out;
}
static Poly reverse_weight(const Poly& p,int weight) {
 require(deg(p)<=weight,"weighted degree bound");Poly q;for(int i=0;i<=deg(p);i++)SetCoeff(q,weight-i,coeff(p,i));return q;
}
static std::vector<Poly> irreducibles(int degree) {
 static std::map<int,std::vector<Poly>> cache;
 if(cache.count(degree))return cache[degree];
 std::vector<Poly> ps;
 for(int code=0;code<smallfilter::pow3(degree);code++){
  Poly p;SetCoeff(p,degree);int z=code;for(int i=0;i<degree;i++,z/=3)SetCoeff(p,i,z%3);
  if(IterIrredTest(p))ps.push_back(p);
 }
 return cache[degree]=ps;
}
struct Normalization {
 int genus=10,delta_inf=0;Poly index,prime;int m=0;zz_pEX chosen;
 std::vector<Branch> places;
};
static Normalization normalize(const std::array<Poly,3>& pars) {
 Normalization n;n.index=1;
 Poly disc=sqr(pars[0])*sqr(pars[1])-cube(pars[1])-cube(pars[0])*pars[2];require(!IsZero(disc),"inseparable cubic");
 MakeMonic(disc);auto fs=CanZass(disc);
 for(long i=0;i<fs.length();i++){
  if(fs[i].b==1&&deg(fs[i].a)>10)continue;
  int delta=0;auto bs=branches(pars,fs[i].a,false,delta);require(delta%deg(fs[i].a)==0,"local index degree");
  n.genus-=delta;n.index*=power(fs[i].a,delta/deg(fs[i].a));
  n.places.insert(n.places.end(),bs.begin(),bs.end());
 }
 std::array<Poly,3> ip;for(int j=0;j<3;j++)ip[j]=reverse_weight(pars[j],4*(j+1));
 auto inf=branches(ip,Poly(INIT_MONO,1),true,n.delta_inf);n.genus-=n.delta_inf;
 n.places.insert(n.places.end(),inf.begin(),inf.end());require(n.genus>=0&&n.genus<=10,"normalization genus");
 if(n.genus==0)return n;
 std::vector<int> degrees;
 for(int d=std::min(5,n.genus);d>=1;d--)if(n.genus%d==0)degrees.push_back(d);
 if(n.genus>5)degrees.push_back(n.genus);
 for(int d:degrees)for(const auto& p:irreducibles(d)){
  if(IsZero(disc%p))continue;
  zz_pE::init(p);auto ff=CanZass(specialize(pars,to_zz_pE(Poly(INIT_MONO,1))));
  for(long i=0;i<ff.length();i++)if(deg(ff[i].a)==1){
   n.prime=p;n.chosen=ff[i].a;n.m=n.genus/d;int delta=0;auto bs=branches(pars,p,false,delta);
   n.places.insert(n.places.end(),bs.begin(),bs.end());return n;
  }
 }
 throw std::runtime_error("no unramified divisor of degree genus");
}
static Series evalpoly(const Poly& p,const Series& t,const Tower& f) {
 auto v=f.series(t.size());for(int k=deg(p);k>=0;k--){v=smul(f,v,t);v[0]=f.add(v[0],f.integer(rep(coeff(p,k))));}return v;
}
static int valuation(const Poly& p,const Poly& prime) {
 if(IsZero(p))return 10000;
 Poly q=p;int n=0;while(IsZero(q%prime)){q/=prime;n++;}return n;
}
struct RRData {
 std::vector<std::array<Poly,3>> basis;mat_zz_p vectors;
 std::vector<std::pair<int,int>> monomials;Poly den;int bound;
};
static RRData riemann_roch(const std::array<Poly,3>& pars,const Normalization& n,int multiple) {
 (void)pars;RRData rr;rr.den=n.index*power(n.prime,multiple*n.m);rr.bound=deg(rr.den)+n.delta_inf;
 for(int j=0;j<3;j++)for(int k=0;k<=rr.bound-4*j;k++)rr.monomials.emplace_back(j,k);
 Rows rows;
 for(const auto& b:n.places){
  bool chosen=!b.infinity&&b.prime==n.prime&&b.h==n.chosen;
  int order=b.infinity?n.delta_inf*b.ramification:valuation(rr.den,b.prime)*b.ramification;
  if(chosen)order-=multiple*n.m;
  if(order==0)continue;
  auto [t,x]=expand_branch(b,order);Tower f(b.h);
  std::vector<Series> tp{f.constant(f.integer(1),order)};
  for(int k=1;k<=rr.bound;k++)tp.push_back(smul(f,tp.back(),t));
  std::array<Series,3> xp{f.constant(f.integer(1),order),x,smul(f,x,x)};
  Rows block(order*b.degree,std::vector<int>(rr.monomials.size()));
  for(size_t i=0;i<rr.monomials.size();i++){
   auto [j,k]=rr.monomials[i];if(b.infinity)k=rr.bound-4*j-k;
   auto e=smul(f,tp[k],xp[j]);
   for(int v=0;v<order;v++){auto flat=f.flatten(e[v]);for(int z=0;z<b.degree;z++)block[v*b.degree+z][i]=flat[z];}
  }
  rows.insert(rows.end(),block.begin(),block.end());
 }
 mat_zz_p mat;mat.SetDims(rr.monomials.size(),rows.size());
 for(size_t i=0;i<rows.size();i++)for(size_t j=0;j<rr.monomials.size();j++)mat[j][i]=rows[i][j];
 kernel(rr.vectors,mat);gauss(rr.vectors);
 for(long i=0;i<rr.vectors.NumRows();i++){
  std::array<Poly,3> v;for(size_t j=0;j<rr.monomials.size();j++)SetCoeff(v[rr.monomials[j].first],rr.monomials[j].second,rr.vectors[i][j]);rr.basis.push_back(v);
 }
 return rr;
}
static int square_code(const RRData& rr,const std::array<Poly,3>& f) {
 vec_zz_p target;target.SetLength(rr.monomials.size());
 for(size_t i=0;i<rr.monomials.size();i++)target[i]=coeff(f[rr.monomials[i].first],rr.monomials[i].second);
 std::vector<int> cs;
 for(long i=0;i<rr.vectors.NumRows();i++){
  long pivot=0;while(pivot<rr.vectors.NumCols()&&IsZero(rr.vectors[i][pivot]))pivot++;
  require(pivot<rr.vectors.NumCols(),"Riemann-Roch basis independence");auto c=target[pivot]/rr.vectors[i][pivot];cs.push_back(rep(c));
  for(long j=pivot;j<rr.vectors.NumCols();j++)target[j]-=c*rr.vectors[i][j];
 }
 require(IsZero(target),"square in L(2D0)");int last=int(cs.size())-1;while(last>=0&&!cs[last])last--;require(last>=0,"nonzero square");
 int code=0;for(int i=0;i<=last;i++)code+=(cs[i]*cs[last]%3)*smallfilter::pow3(i);return code;
}
static std::array<Poly,3> square_function(const std::array<Poly,3>& f,const std::array<Poly,3>& p) {
 std::array<Poly,5> c;for(int i=0;i<3;i++)for(int j=0;j<3;j++)c[i+j]+=f[i]*f[j];
 for(int i=4;i>=3;i--)for(int j=0;j<3;j++)c[i-1-j]-=c[i]*p[j];
 return {c[0],c[1],c[2]};
}
static PlaceData expansion(const Branch& b,const RRData& rr,const Normalization& n) {
 int m=!b.infinity&&b.prime==n.prime&&b.h==n.chosen?n.m:0;
 int dv=b.infinity?n.delta_inf*b.ramification:valuation(rr.den,b.prime)*b.ramification;
 int steps=2*n.genus/b.degree+1;auto [t,x]=expand_branch(b,steps+dv);Tower f(b.h);
 Poly denominator=b.infinity?reverse_weight(rr.den,deg(rr.den)):rr.den;
 auto ds=evalpoly(denominator,t,f);
 if(b.infinity){auto up=f.constant(f.integer(1),t.size());for(int k=0;k<n.delta_inf;k++)up=smul(f,up,t);ds=smul(f,ds,up);}
 for(int k=0;k<dv;k++)require(f.iszero(ds[k]),"denominator valuation");
 require(!f.iszero(ds[dv]),"denominator leading term");
 auto unit=f.series(steps);for(int k=0;k<steps;k++)unit[k]=ds[k+dv];auto inverse=sinv(f,unit),x2=smul(f,x,x);
 PlaceData p{b.degree,m,steps,Rows(steps*b.degree,std::vector<int>(rr.basis.size()))};
 for(size_t j=0;j<rr.basis.size();j++){
  auto v=rr.basis[j];if(b.infinity)for(int i=0;i<3;i++)v[i]=reverse_weight(v[i],rr.bound-4*i);
  auto ns=f.add(f.add(evalpoly(v[0],t,f),smul(f,evalpoly(v[1],t,f),x)),smul(f,evalpoly(v[2],t,f),x2));
  int shift=dv-2*m;require(shift>=0,"expansion shift");
  for(int k=0;k<shift;k++)require(f.iszero(ns[k]),"Riemann-Roch pole constraint");
  auto ss=f.series(steps);for(int k=0;k<steps;k++)ss[k]=ns[k+shift];ss=smul(f,ss,inverse);
  for(int k=0;k<steps;k++){auto flat=f.flatten(ss[k]);for(int z=0;z<b.degree;z++)p.rows[k*b.degree+z][j]=flat[z];}
 }
 return p;
}
static int small_eval(const Poly& p,int z,const smallfilter::Field* f) {
 int n=0;for(long i=deg(p);i>=0;i--)n=f->add[f->mul[n][z]][rep(coeff(p,i))];return n;
}
static int encode_poly(const Poly& p) {int code=0;for(long i=deg(p);i>=0;i--)code=code*3+rep(coeff(p,i));return code;}
static std::string filter_text(const Input& a) {
 std::vector<std::string> tests;
 for(int degree=1;degree<=4;degree++){
  auto primes=irreducibles(degree);int q=smallfilter::pow3(degree);
  int mod=degree==1?0:encode_poly(primes[0]-Poly(INIT_MONO,degree));auto f=smallfilter::field(q,mod);
  for(const auto& prime:primes){
   int alpha=-1;for(int z=0;z<q;z++)if(small_eval(prime,z,f)==0){alpha=z;break;}require(alpha>=0,"residue root");
   int den=small_eval(a.den,alpha,f);if(!den)continue;int inverse=1;while(f->mul[den][inverse]!=1)inverse++;
   std::ostringstream s;s<<q<<' '<<mod;for(const auto& p:a.pars)s<<' '<<small_eval(p,alpha,f);
   for(const auto& v:a.basis)for(const auto& p:v)s<<' '<<int(f->mul[small_eval(p,alpha,f)][inverse]);
   tests.push_back(s.str());
  }
 }
 std::ostringstream s;s<<tests.size()<<'\n';for(const auto& test:tests)s<<test<<'\n';return s.str();
}
static bool irreducible_cubic(const std::array<Poly,3>& p) {
 for(const auto& prime:irreducibles(1)){
  zz_pE::init(prime);auto cub=specialize(p,to_zz_pE(Poly(INIT_MONO,1)));
  if(IterIrredTest(cub))return true;
 }
 for(int code=0;code<243;code++){
  Poly root;int v=code;for(int i=0;i<5;i++,v/=3)SetCoeff(root,i,v%3);
  if(IsZero(cube(root)+p[0]*sqr(root)+p[1]*root+p[2]))return false;
 }
 return true;
}
static Input native_prepare(std::istream& in) {
 Input a;in>>a.key;for(auto& p:a.pars)p=readpoly(in);require(bool(in),"native cubic input");
 for(int i=0;i<3;i++)require(deg(a.pars[i])<=4*(i+1),"K3 coefficient degrees");
 require(irreducible_cubic(a.pars),"fallback required: reducible cubic");
 auto n=normalize(a.pars);a.g=n.genus;a.dim=a.g+1;
 if(!a.g){a.den=1;a.basis={std::array<Poly,3>{Poly(1),Poly(),Poly()}};a.squares.insert(1);return a;}
 auto rr=riemann_roch(a.pars,n,2);require(int(rr.basis.size())==a.dim,"L(2D0) dimension");a.basis=rr.basis;a.den=rr.den;
 auto small=riemann_roch(a.pars,n,1);require(!small.basis.empty()&&small.basis.size()<=6,"L(D0) dimension");
 for(int top=1;top<smallfilter::pow3(small.basis.size());top*=3)for(int code=top;code<2*top;code++){
  auto cs=digits(code,small.basis.size());std::array<Poly,3> v;
  for(size_t i=0;i<small.basis.size();i++)for(int j=0;j<3;j++)v[j]+=small.basis[i][j]*cs[i];
  v=square_function(v,a.pars);for(auto& p:v){require(IsZero(p%n.index),"integral square denominator");p/=n.index;}
  a.squares.insert(square_code(rr,v));
 }
 for(const auto& b:n.places)if(b.degree<=a.g||(!b.infinity&&b.prime==n.prime&&b.h==n.chosen))a.exceptional.push_back(expansion(b,rr,n));
 return a;
}
static void writepoly(std::ostream& out,const Poly& p) {out<<deg(p)+1;for(int i=0;i<=deg(p);i++)out<<' '<<rep(coeff(p,i));out<<'\n';}
static void writeinput(std::ostream& out,const Input& a,const std::string& filters) {
 out<<"D2CPP1 "<<a.key<<' '<<a.g<<' '<<a.dim<<'\n';for(const auto& p:a.pars)writepoly(out,p);writepoly(out,a.den);
 for(const auto& v:a.basis)for(const auto& p:v)writepoly(out,p);
 out<<filters;
 out<<a.squares.size();for(int code:a.squares)out<<' '<<code;out<<'\n'<<a.exceptional.size()<<'\n';
 for(const auto& p:a.exceptional){out<<p.degree<<' '<<p.m<<' '<<p.steps<<'\n';for(const auto& row:p.rows){for(int v:row)out<<v<<' ';out<<'\n';}}
}
