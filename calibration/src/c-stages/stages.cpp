#include "core.hpp"
#include "../lfun-tate/tate.hpp"
#include "../search-c/search.hpp"
#include <boost/json/src.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <ctime>
#include <iostream>
#include <map>
#include <sstream>
namespace j=boost::json;
std::string value_string(const j::value& x){return x.is_string()?std::string(x.as_string()):j::serialize(x);}
j::object parse_input(const std::string& line){
 boost::property_tree::ptree tree;std::istringstream stream(line);boost::property_tree::read_json(stream,tree);
 j::object out;
 for(const char* name:{"q","key","b"})out[name]=tree.get<std::string>(name);
 for(const char* name:{"codes","L"}){j::array a;for(const auto& v:tree.get_child(name))a.emplace_back(v.second.data());out[name]=a;}
 if(auto points=tree.get_child_optional("points")){j::array a;for(const auto& pt:*points){j::array c;for(const auto& v:pt.second)c.emplace_back(v.second.data());a.push_back(c);}out["points"]=a;}
 return out;
}
long integer(const j::value& x){return std::stol(value_string(x));}
j::array matrix_json(const Matrix& a){j::array r;for(const auto& row:a){j::array s;for(const auto& x:row)s.emplace_back(x.get_str());r.push_back(s);}return r;}
Point point(const j::value&v,int q){auto c=v.as_array();if(c.empty())return {};if(c.size()!=4)throw std::runtime_error("expected four coordinate codes");return Point(R(decode(value_string(c[0]),q),decode(value_string(c[1]),q)),R(decode(value_string(c[2]),q),decode(value_string(c[3]),q)));}
j::array point_json(const Point&p,int q){if(p.identity)return {};return {encode(p.x.n,q),encode(p.x.d,q),encode(p.y.n,q),encode(p.y.d,q)};}
static double cpu(){timespec t;clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&t);return t.tv_sec+1e-9*t.tv_nsec;}
struct Analytic {int rank;Q product;};
Analytic analytic(const j::array& coefficients,int q,const std::vector<LocalData>& local,int torsion){
 std::vector<Q> f;for(const auto&v:coefficients)f.emplace_back(value_string(v));
 if(f.empty()||f.front()!=1)throw std::runtime_error("expected signed L-polynomial with constant term one");
 int rank=0;
 for(;;){std::vector<Q> quo(f.size()-1);for(size_t i=0;i<quo.size();++i)quo[i]=f[i]+(i?q*quo[i-1]:Q(0));if(quo.empty()||f.back()+q*quo.back()!=0)break;f=quo;++rank;}
 Q z=0;for(auto i=f.rbegin();i!=f.rend();++i)z=z/q+*i;
 Q product=q*z*torsion*torsion;for(const auto& l:local)product/=l.tamagawa;
 return {rank,product};
}
std::pair<Q,int> contribution(const LocalData& l){
 const auto& s=l.kodaira;
 if(s=="II"||s=="II*"||s=="I0"||s=="I1")return {0,1};
 if(s=="III")return {rat(1,2),2};
 if(s=="III*")return {rat(3,2),2};
 if(s=="IV")return {rat(2,3),3};
 if(s=="IV*")return {rat(4,3),3};
 int n=std::stoi(s.substr(1));if(s.back()=='*')return {1+rat(n,4),4};return {rat((n/2)*(n-n/2),n),n};
}
std::vector<Point> search(const Curve& e,const std::array<uint64_t,5>& codes,int b){
 if(b<0)return {};
 if(b>(e.q==2?7:6))throw std::runtime_error("native section bound exceeded; Magma required");
 Stats stats;std::ostringstream out;if(e.q==2)search2(codes,b,out,stats);else search3(codes,b,out,stats);
 auto rows=j::parse("["+out.str()+"]").as_array();std::vector<Point> pts;
 for(const auto& v:rows){auto a=v.as_array();P den=decode(value_string(a[1]),e.q);pts.emplace_back(R(decode(value_string(a[0]),e.q),NTL::power(den,2)),R(decode(value_string(a[2]),e.q),NTL::power(den,3)));}
 return pts;
}
std::vector<Point> torsion(const Curve&e,const std::array<uint64_t,5>&codes,const Q&cmax){
 std::vector<Point> result(1);if(cmax<4)return result;
 for(const auto&p:search(e,codes,floorq((cmax-4)/2)))if(e.height(p)==0)result.push_back(p);
 return result;
}
j::array torsion_invariants(const Curve&e,const std::vector<Point>&pts){
 int exponent=1;
 for(const auto&p:pts){Point s;int order=0;do{s=e.add(s,p);++order;if(order>int(pts.size()))throw std::runtime_error("torsion enumeration is not closed");}while(!s.identity);exponent=std::lcm(exponent,order);}
 if(pts.size()==1)return {};
 int n=pts.size()/exponent;return n==1?j::array{exponent}:j::array{n,exponent};
}
struct Saturation {bool resolved=true,ok=true;j::array checks,primes,unresolved_primes;};
Saturation saturated(const Curve&e,const std::array<uint64_t,5>&codes,const Lattice&lat,const Q&cmax,int d,const std::vector<Point>&tors,std::map<int,std::vector<Point>>& cache,int cap=100){
 Saturation out;std::map<std::string,Q> heights;int r=lat.basis.size();if(r==0)return out;
 Q n=determinant(lat.gram);for(int i=0;i<r;++i)n*=d;
 if(n.get_den()!=1||n<=0)throw std::runtime_error("nonintegral lattice discriminant");
 if(!n.get_num().fits_ulong_p()){out.resolved=false;return out;}
 unsigned long nn=n.get_num().get_ui();std::vector<int> primes;
 for(unsigned long p=2;p<=nn/p;++p)if(nn%p==0){int power=0;do{nn/=p;++power;}while(nn%p==0);if(power>=2)primes.push_back(p);}
 struct DivisionClass {std::vector<int> v;Q h;long bound;};
 struct Prime {int ell;bool class_limit=false;long bound=0;std::vector<DivisionClass> classes;};
 std::vector<Prime> plans;
 for(int ell:primes){
  Prime plan;plan.ell=ell;
  unsigned long classes=1;for(int i=0;i<r;++i){if(classes>20000UL/ell){plan.class_limit=true;break;}classes*=ell;}
  if(plan.class_limit){plans.push_back(plan);continue;}
  for(unsigned long code=1;code<classes;++code){unsigned long z=code;std::vector<int> v(r);for(int i=0;i<r;++i){v[i]=z%ell;z/=ell;}int first=0;while(v[first]==0)++first;if(v[first]!=1)continue;
   Q h=0;for(int i=0;i<r;++i)for(int k=0;k<r;++k)h+=v[i]*v[k]*lat.gram[i][k];
   long bound=floorq((h/ell/ell-4+cmax)/2);
   if(plan.classes.empty()||bound>plan.bound)plan.bound=bound;
   plan.classes.push_back({v,h,bound});
  }
  plans.push_back(std::move(plan));
 }
 std::stable_sort(plans.begin(),plans.end(),[](const Prime&a,const Prime&b){if(a.class_limit!=b.class_limit)return !a.class_limit;return a.bound<b.bound;});
 for(auto& plan:plans){
  int ell=plan.ell;bool complete=!plan.class_limit,divisible=false;
  std::stable_sort(plan.classes.begin(),plan.classes.end(),[](const DivisionClass&a,const DivisionClass&b){return a.bound<b.bound;});
  for(const auto& cls:plan.classes){
   const auto&v=cls.v;const auto&h=cls.h;long bound=cls.bound;
   if(bound>std::min(cap,e.q==2?7:6)){complete=false;continue;}
   Point p;std::map<std::string,int> multiples;
   if(bound>=0){
    auto found=cache.lower_bound(bound);
    if(found==cache.end())found=cache.emplace(bound,search(e,codes,bound)).first;
    for(const auto&pt:found->second){
     auto label=j::serialize(point_json(pt,e.q));auto it=heights.find(label);
     if(it==heights.end())it=heights.emplace(label,e.height(pt)).first;
     if(it->second==h/ell/ell)++multiples[j::serialize(point_json(e.mul(pt,ell),e.q))];
    }
    if(!multiples.empty())for(int i=0;i<r;++i)if(v[i])p=e.add(p,e.mul(lat.basis[i],v[i]));
   }
   for(const auto&t:tors){int count=multiples.empty()?0:multiples[j::serialize(point_json(e.add(p,t),e.q))];if(count)divisible=true;j::array vv;for(int x:v)vv.push_back(x);out.checks.push_back({ell,vv,point_json(t,e.q),count,bound});}
  }
  j::object status={{"ell",ell},{"bound",plan.class_limit?j::value():j::value(plan.bound)},{"status",divisible?"not saturated":complete?"saturated":"unresolved"}};
  if(!complete)status["reason"]=plan.class_limit?"coefficient vector limit":"native section bound exceeded";
  out.primes.push_back(status);
  if(divisible)out.ok=false;
  else if(!complete)out.unresolved_primes.push_back(status);
 }
 out.resolved=!out.ok||out.unresolved_primes.empty();
 return out;
}
struct Options {bool supplied=false,audit=false,bundle=false;};
j::object run(const j::object& row,const Options& opts){
 double start=cpu();int q=integer(row.at("q"));if(q!=2&&q!=3)throw std::runtime_error("expected field 2 or 3");NTL::zz_p::init(q);
 std::array<P,5>a;std::array<uint64_t,5>codes;Surface surface(q);surface.key=std::stoull(value_string(row.at("key")));
 for(int i=0;i<5;++i){codes[i]=std::stoull(value_string(row.at("codes").as_array()[i]));a[i]=decode(std::to_string(codes[i]),q);for(int k=0;k<=24;++k)surface.a[i][k]=NTL::rep(NTL::coeff(a[i],k));}
 Curve e(q,a);double setup=1000*(cpu()-start);
 start=cpu();bool minimal;auto local=tate(surface,&minimal);if(!minimal)throw std::runtime_error("requires everywhere minimal K3 model");Q cmax=0;int d=1;j::array local_json;
 for(const auto& l:local){auto [c,dd]=contribution(l);cmax+=l.degree*c;d=std::lcm(d,dd);bool inf=l.place=="inf";P pi=decode(inf?std::to_string(q):l.place,q);e.places.push_back({pi,l.vd,l.degree,inf});local_json.push_back({l.place,l.degree,l.kodaira,l.tamagawa,l.vd,l.conductor});}
 double local_ms=1000*(cpu()-start);
 start=cpu();auto tors=torsion(e,codes,cmax);auto tinv=torsion_invariants(e,tors);double torsion_ms=1000*(cpu()-start);
 start=cpu();auto ai=analytic(row.at("L").as_array(),q,local,tors.size());double bsd_ms=1000*(cpu()-start);
 int b=integer(row.at("b"));std::vector<Point> pts;j::array points_json;
 double search_ms=0;
 if(ai.rank>0&&!opts.supplied){start=cpu();pts=search(e,codes,b);search_ms=1000*(cpu()-start);}
 start=cpu();
 if(row.contains("points")){
  std::vector<Point> given;for(const auto&v:row.at("points").as_array()){auto p=point(v,q);if(!e.contains(p))throw std::runtime_error("point outside curve");given.push_back(p);}
  if(ai.rank>0&&!opts.supplied){auto key=[&](const Point&p){return j::serialize(point_json(p,q));};std::vector<std::string>x,y;for(const auto&p:pts)x.push_back(key(p));for(const auto&p:given)y.push_back(key(p));std::sort(x.begin(),x.end());std::sort(y.begin(),y.end());if(x!=y)throw std::runtime_error("provided point set differs from exhaustive native search");}
  pts=given;
 }else if(opts.supplied)throw std::runtime_error("--supplied requires complete points");
 auto lat=build(e,ai.rank?pts:std::vector<Point>{});double lattice_ms=1000*(cpu()-start);
 start=cpu();Q mu=babai(lat.gram),det=determinant(lat.gram),h0=6+2*b-cmax;
 bool certified=ai.rank==0||(int(lat.basis.size())==ai.rank&&mu<h0);
 Q sha=ai.product/det;
 if(certified&&(sha.get_den()!=1||sha<1||!mpz_perfect_square_p(sha.get_num_mpz_t())))throw std::runtime_error("BSD square check failed");
 std::map<int,std::vector<Point>> cache;
 if(ai.rank>0)cache[b]=pts;
 Saturation sat;double saturation_ms=0;
 if(certified&&!lat.basis.empty()){double ts=cpu();sat=saturated(e,codes,lat,cmax,d,tors,cache,opts.bundle?b:100);saturation_ms=1000*(cpu()-ts);if(sat.resolved&&!sat.ok)throw std::runtime_error("division contradicts covering certificate");}
 double cert_ms=1000*(cpu()-start);
 j::array basis;for(const auto&p:lat.basis)basis.push_back(point_json(p,q));
 j::object result={{"q",q},{"key",surface.key},{"local",local_json},{"torsion",tinv},{"rank",ai.rank},{"sha_reg",ai.product.get_str()},{"Cmax",cmax.get_str()},{"D",d},{"gram",matrix_json(lat.gram)},{"basis",basis},{"mu",mu.get_str()},{"certified",certified},{"ms",j::array{setup,local_ms,torsion_ms,bsd_ms,lattice_ms,cert_ms}},{"search_ms",search_ms},{"saturation_ms",saturation_ms}};
 if(opts.bundle){
  j::array sections,heights,torsion_points;
  for(const auto&p:pts){sections.push_back(point_json(p,q));heights.emplace_back(e.height(p).get_str());}
  for(const auto&p:tors)torsion_points.push_back(point_json(p,q));
  result["sections"]=sections;result["heights"]=heights;result["torsion_points"]=torsion_points;
  result["bound"]=b;result["missing_height_bound"]=h0.get_str();result["determinant"]=det.get_str();
  if(!certified){double ts=cpu();sat=saturated(e,codes,lat,cmax,d,tors,cache,b);saturation_ms=1000*(cpu()-ts);}
  result["saturation_ms"]=saturation_ms;result["saturation_resolved"]=sat.resolved;
  result["saturated"]=sat.resolved?j::value(sat.ok):j::value();result["division_checks"]=sat.checks;
  result["saturation_primes"]=sat.primes;result["unresolved_primes"]=sat.unresolved_primes;
 }
 if(opts.audit){Matrix pg;for(const auto&p:pts){std::vector<Q> v;for(const auto&s:pts)v.push_back(e.pairing(p,s));pg.push_back(v);}result["point_gram"]=matrix_json(pg);
  if(!certified){double ts=cpu();sat=saturated(e,codes,lat,cmax,d,tors,cache,opts.bundle?b:100);saturation_ms=1000*(cpu()-ts);}result["saturation_ms"]=saturation_ms;result["saturation_resolved"]=sat.resolved;result["saturated"]=sat.resolved?j::value(sat.ok):j::value();result["division_checks"]=sat.checks;result["saturation_primes"]=sat.primes;result["unresolved_primes"]=sat.unresolved_primes;
  result["negative_passed"]=true;
  if(!lat.basis.empty()&&sat.resolved&&sat.ok){auto negative=lat;negative.basis[0]=e.mul(negative.basis[0],2);for(size_t i=0;i<negative.gram.size();++i){negative.gram[0][i]*=2;negative.gram[i][0]*=2;}auto ns=saturated(e,codes,negative,cmax,d,tors,cache);result["negative_passed"]=ns.resolved?j::value(ns.ok):j::value();result["negative_checks"]=ns.checks;result["negative_primes"]=ns.primes;result["negative_unresolved_primes"]=ns.unresolved_primes;}
 }
 return result;
}
int main(int argc,char**argv){try{Options opts;for(int i=1;i<argc;++i){std::string arg=argv[i];if(arg=="--supplied")opts.supplied=true;else if(arg=="--audit")opts.audit=true;else if(arg=="--bundle")opts.bundle=true;else throw std::runtime_error("usage: stages [--supplied] [--audit]");}std::string line;while(std::getline(std::cin,line)){if(line.empty())continue;auto row=parse_input(line);std::cout<<j::serialize(run(row,opts))<<std::endl;}}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
