SetQuitOnError(true);SetNthreads(1);SetColumns(0);
AttachSpec("src/census/spec");
q:=StringToInteger(q);
R<s,u>:=PolynomialRing(GF(q),2);P<t>:=PolynomialRing(GF(q));
lines:=Split(Read(input),"\n");
for line in lines do
 fields:=[StringToInteger(v):v in Split(line," ")];
 W:=aInvariantsFromKey(R,fields[1]);
 a:=[P|Evaluate(w,[t,P!1]):w in W];
 enc:=[&+[Integers()|Integers()!Coefficient(w,j)*q^j:j in [0..Max(0,Degree(w))]]:w in a];
 assert enc eq fields[2..6];
end for;
printf "DECODE_OK|%o\n",#lines;
quit;
