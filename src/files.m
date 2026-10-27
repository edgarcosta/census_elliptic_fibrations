import "models.m": valid_ring, prime_field_key, prefix_degrees;
import "forms.m": coefficients;
import "census.m": prefixes, selection_error, selections, selected_representative_keys;

function output_paths(stem, indices)
    return [ stem cat ".csv" ] cat &cat[ [ Sprintf("%o_%o.bin",stem,IntegerToString(10000+i)[2..5]),
        Sprintf("%o_%o.stats",stem,IntegerToString(10000+i)[2..5]) ] : i in indices ];
end function;

intrinsic WriteFibrationRepresentativeKeys(R::RngMPol, stem::MonStgElt :
    PrefixIndices := [], GenericFibreKinds := [], Minimalities := []) -> RngIntElt
{The number of representative keys written over R to new files with stem.
The parameter PrefixIndices (default []) selects all prefixes when empty.
The parameter GenericFibreKinds (default []) selects "smooth", "regular cusp",
and "nonregular" when empty, or a subset of these case-sensitive values.
The parameter Minimalities (default []) selects "everywhere minimal" and
"nonminimal somewhere" when empty, or a subset of these case-sensitive values.
See also FibrationRepresentativeKeys, aInvariantsFromKey}
    require valid_ring(R): "R must be a binary polynomial ring over GF(2) or GF(3)";
    P := prefixes(R);
    message := selection_error(PrefixIndices,GenericFibreKinds,Minimalities,#P);
    require message eq "": message;
    indices,kinds,minimalities := selections(PrefixIndices,GenericFibreKinds,Minimalities,#P);
    require &and[ not FileExists(path) : path in output_paths(stem,indices) ]:
        "An output path already exists";
    summary := Open(stem cat ".csv","w");
    Write(summary,"file_format_version,key_encoding_version,characteristic,prefix_index,prefix_form_1_degree,prefix_form_1_key,prefix_form_2_degree,prefix_form_2_key,generic_fibre_kind,minimality,number_of_keys\n");
    Flush(summary);
    total := 0;
    p := Characteristic(R);
    degrees := prefix_degrees(R);
    for i in indices do
        keys,cells := selected_representative_keys(R,P[i],kinds,minimalities);
        binary := Open(Sprintf("%o_%o.bin",stem,IntegerToString(10000+i)[2..5]),"wb");
        for lo in [1..#keys by 4096] do
            WriteBytes(binary,&cat[ [ (key div 256^j) mod 256 : j in [0..7] ]
                : key in keys[lo..Min(lo+4095,#keys)] ]);
        end for;
        Flush(binary);
        delete binary;
        form_keys := [ prime_field_key(coefficients(P[i][j],degrees[j]),p) : j in [1,2] ];
        metadata := [ "file_format_version=1","key_encoding_version=1",
            Sprintf("characteristic=%o",p),Sprintf("prefix_index=%o",i),
            Sprintf("prefix_form_1_degree=%o",degrees[1]),Sprintf("prefix_form_1_key=%o",form_keys[1]),
            Sprintf("prefix_form_2_degree=%o",degrees[2]),Sprintf("prefix_form_2_key=%o",form_keys[2]),
            "generic_fibre_kinds=" cat Join(kinds,"|"),
            "minimalities=" cat Join(minimalities,"|"),Sprintf("number_of_keys=%o",#keys) ];
        for kind in kinds do
            for minimality in minimalities do
                n := cells[<kind,minimality>];
                Write(summary,Sprintf("1,1,%o,%o,%o,%o,%o,%o,%o,%o,%o\n",
                    p,i,degrees[1],form_keys[1],degrees[2],form_keys[2],kind,minimality,n));
                Append(~metadata,Sprintf("cell=%o|%o|%o",kind,minimality,n));
            end for;
        end for;
        Flush(summary);
        sidecar := Open(Sprintf("%o_%o.stats",stem,IntegerToString(10000+i)[2..5]),"w");
        Write(sidecar,Join(metadata,"\n") cat "\ncomplete=1\n");
        Flush(sidecar);
        delete sidecar;
        total +:= #keys;
        delete keys;
        vprintf FibrationCensus,1: "Wrote prefix %o\n",i;
    end for;
    Flush(summary);
    delete summary;
    return total;
end intrinsic;
