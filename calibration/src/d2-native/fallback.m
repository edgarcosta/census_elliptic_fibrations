SetQuitOnError(true); SetNthreads(1); SetColumns(0);
record:=eval Read(IN);key,pars:=Explode(record);
k:=GF(3);R:=PolynomialRing(k);F<t>:=FunctionField(k);P<x>:=PolynomialRing(F);
cub:=x^3+(F!(R!pars[1]))*x^2+(F!(R!pars[2]))*x+F!(R!pars[3]);
assert IsSeparable(cub);SetSeed(key mod 2147483647);tm:=Cputime();
d2:=0;genus:=0;
for pair in Factorization(cub) do
 assert pair[2] eq 1;
 if Degree(pair[1]) eq 1 then continue;end if;
 B:=FunctionField(pair[1]);genus+:=Genus(B);
 inv:=ClassGroupAbelianInvariants(B:Proof:=true);
 d2+:=#[v:v in inv|v ne 0 and IsEven(v)];
 printf "FACTORS|%o|%o\n",Degree(pair[1]),inv;
end for;
printf "FALLBACK|%o|%o|%o|%o\n",key,genus,d2,Cputime(tm);
print "COMPLETE";quit;
