#pragma once
#include <NTL/lzz_pX.h>
#include <gmpxx.h>
#include <algorithm>
#include <array>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>
using P=NTL::zz_pX;
using Q=mpq_class;
using Z=mpz_class;
using Matrix=std::vector<std::vector<Q>>;
inline Q rat(long n,long d=1) {Q x(n,d);x.canonicalize();return x;}
inline long floorq(const Q& x) {Z z;mpz_fdiv_q(z.get_mpz_t(),x.get_num_mpz_t(),x.get_den_mpz_t());if(!z.fits_slong_p())throw std::runtime_error("bound too large");return z.get_si();}
P decode(std::string s,int q);
std::string encode(const P& a,int q);
struct R {
 P n,d;
 R(long a=0):n(a),d(1){}
 R(const P& a):n(a),d(1){}
 R(P a,P b);
 bool zero()const{return IsZero(n);}
};
R operator+(const R&,const R&);R operator-(const R&,const R&);R operator-(const R&);
R operator*(const R&,const R&);R operator/(const R&,const R&);
bool operator==(const R&,const R&);R power(R,int);
struct Point {R x,y;bool identity=true;Point()=default;Point(R x0,R y0):x(x0),y(y0),identity(false){};};
bool operator==(const Point&,const Point&);
struct Place {P pi;int vd,degree;bool infinity;};
struct Curve {
 int q;
 std::array<P,5> a,inf;
 P b2,b4,b6,b8,c4,delta;
 std::vector<Place> places;
 Curve(int q0,const std::array<P,5>& a0);
 Point neg(const Point&)const;
 Point add(const Point&,const Point&)const;
 Point mul(Point,Z)const;
 bool contains(const Point&)const;
 Q height(const Point&)const;
 Q pairing(const Point&,const Point&)const;
 Q local(const Point&,const Place&)const;
};
Q determinant(Matrix);
std::vector<Q> solve(Matrix,std::vector<Q>);
struct Lattice {std::vector<Point> basis;Matrix gram;};
Lattice build(const Curve&,const std::vector<Point>&);
Q babai(Matrix);
