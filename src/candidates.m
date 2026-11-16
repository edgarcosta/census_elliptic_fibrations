function rational_corrections(row)
    c := row`Tamagawa;
    values := {Rationals() | 0};
    case row`Type:
        when "I":
            n := row`Index;
            if c eq n then
                values := {Rationals() | i*(n-i)/n : i in [0..n-1]};
            elif IsEven(n) then
                Include(~values, n/4);
            end if;
        when "I*":
            if c in [2, 4] then Include(~values, 1); end if;
            if c eq 4 then Include(~values, 1+row`Index/4); end if;
        when "III": Include(~values, 1/2);
        when "III*": Include(~values, 3/2);
        when "IV":
            if c eq 3 then Include(~values, 2/3); end if;
        when "IV*":
            if c eq 3 then Include(~values, 4/3); end if;
    end case;
    return {Rationals() | row`Degree*x : x in values};
end function;

function lattice_parity(discriminant, rho)
    if IsOdd(rho) then return IsEven(discriminant); end if;
    return discriminant mod 4 eq 0
        or ((-1)^(rho div 2+1)*discriminant) mod 4 eq 1;
end function;

intrinsic FibrationCandidates(data::Rec, fibres::SeqEnum) -> SeqEnum, SeqEnum
{The sorted admissible square Sha orders and matching rational regulators from
a FibrationInvariants record and complete FibrationLocalData records of the same
elliptic K3 surface. The supplied data are assumed correct.
See also FibrationInvariants}
    rank := data`Rank;
    rho := data`PicardNumber;
    N := data`ArtinTateProduct;
    product := Rationals()!data`ShaRegulatorProduct;
    require Type(rank) eq RngIntElt and rank in [0..20]:
        "rank must be between 0 and 20";
    require Type(rho) eq RngIntElt and rho in [2..22]:
        "Picard number must be between 2 and 22";
    require IsCoercible(Integers(), N) and N gt 0 and product gt 0:
        "Artin-Tate product must be a positive integer and Sha-regulator product positive";
    N := Integers()!N;
    squares := [Integers() | 1];
    for factor in Factorization(N) do
        squares := [s*factor[1]^(2*i) : s in squares, i in [0..factor[2] div 2]];
    end for;
    Sort(~squares);
    squares := [s : s in squares | lattice_parity(N div s, rho)];
    if rank eq 0 then
        squares := [s : s in squares | s eq product];
    else
        corrections := {Rationals() | 0};
        for row in fibres do
            corrections := {a+b : a in corrections, b in rational_corrections(row)};
        end for;
        if rank eq 1 then
            squares := [s : s in squares | &or[product/s+c ge 4
                and IsIntegral((product/s+c-4)/2) : c in corrections]];
        else
            minimum := Min([4+2*Max(0, Floor((c-4)/2)+1)-c : c in corrections]);
            if rank le 5 then
                gamma_powers := [Rationals() | 1, 4/3, 2, 4, 8];
                squares := [s : s in squares | product/s ge minimum^rank/gamma_powers[rank]];
            else
                // Hermite's bound gamma_r <= (4/3)^((r-1)/2), raised to rank r.
                squares := [s : s in squares | (product/s)*(4/3)^(rank*(rank-1) div 2)
                    ge minimum^rank];
            end if;
        end if;
    end if;
    return squares, [Rationals() | product/s : s in squares];
end intrinsic;

intrinsic FibrationCandidates(E::CrvEll[FldFunRat]) -> SeqEnum, SeqEnum
{The sorted admissible square Sha orders and matching rational regulators of an
everywhere minimal elliptic K3 polynomial model E over GF(2)(t) or GF(3)(t).
A singleton determines both values. See also FibrationInvariants}
    return FibrationCandidates(FibrationInvariants(E), FibrationLocalData(E));
end intrinsic;
