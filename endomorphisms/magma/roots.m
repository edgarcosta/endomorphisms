intrinsic RootsPari(f::RngUPolElt, K::Fld) -> .
{Given a polynomial f and a field K, finds all roots of f in K.}

assert BaseRing(f) eq Rationals();
assert BaseRing(K) eq Rationals();
g := DefiningPolynomial(K);
cmd := Sprintf(
"{f = Pol(Vecrev(%o),'x); g = Pol(Vecrev(%o),'y); K = nfinit(g); print1(apply(z->vector(poldegree(g),i, polcoeff(z,i-1)),lift(nfroots(g,f))))}",
Coefficients(f), Coefficients(g));
s := Pipe("gp -q", cmd);
rts := [ K ! rt : rt in eval(s) ];
return [ rt : rt in rts | Evaluate(f, rt) eq 0 ];

end intrinsic;


intrinsic HasRootPari(f::RngUPolElt, K::Fld) -> .
{Given a polynomial f and a field K, determines if f has a root in K.}

return #RootsPari(f, K) ne 0;

end intrinsic;


intrinsic AutomorphismGroupPari(K::Fld) -> GrpPerm, RngIntElt, Map
{Similar to usual function, but outsources to Pari for better performance.}

assert BaseRing(K) eq Rationals();
rts := RootsPari(DefiningPolynomial(K), K); S := Sym(#rts);

sigmas := [ ];
for rt in rts do
    h := hom< K -> K | rt >;
    if h(K`base_gen) eq K`base_gen then
        L := [ ];
        for i in [1..#rts] do
            for j in [1..#rts] do
                if h(rts[i]) eq rts[j] then
                    Append(~L, j);
                end if;
            end for;
        end for;
        Append(~sigmas, S ! L);
    end if;
end for;
Gp := sub< S | sigmas >;

for i in [1..#rts] do
    if rts[i] eq K.1 then
        n0 := i;
    end if;
end for;

tups := [ ]; hs := [ ];
for sigma in Gp do
    j := Eltseq(sigma)[n0];
    h := hom< K -> K | rts[j] >;
    Append(~tups, <sigma, h>);
    Append(~hs, h);
end for;
hs := Set(hs);
Gphi := map< Gp -> hs | tups >;

return Gp, 0, Gphi;

end intrinsic;


intrinsic FactorizationPari(f::RngUPolElt, K::Fld) -> .
{Given a polynomial f and a field K, finds the factorization of f over K.}

assert BaseRing(f) eq Rationals();
assert BaseRing(K) eq Rationals();
g := DefiningPolynomial(K);
cmd := Sprintf(
"{f = Pol(Vecrev(%o),'x); g = Pol(Vecrev(%o),'y); K = nfinit(g); apply(h->apply(c->vector(poldegree(g),i,polcoeff(c,i-1)),lift(Vecrev(h))),nffactor(K,f)[,1]~)",
Coefficients(f), Coefficients(g));
s := Pipe("gp -q", cmd);

R := PolynomialRing(K);
seqs := eval(s);
facs := [ &+[ (K ! seq[i])*R.1^(i - 1) : i in [1..#seq] ] : seq in seqs ];
return facs;

end intrinsic;


intrinsic SplittingFieldPari(f::RngUPolElt) -> .
{Splitting field of f calculated using Pari.}

assert BaseRing(f) eq Rationals();
cmd := Sprintf(
"{f = Pol(Vecrev(%o),'x); K = nfinit(f); nfsplitting(K)",
Coefficients(f), Coefficients(f));
s := Pipe("gp -q", cmd);
R<x> := PolynomialRing(BaseRing(f));
return Polredbestabs(NumberField(eval(s)));

end intrinsic;
