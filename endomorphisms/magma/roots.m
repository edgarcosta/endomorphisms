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


intrinsic AutomorphismGroupPari(K::Fld) -> .
{Like the usual function, but outsources to Pari for better performance, and the second return value is usually 0 because we do not use it. Stores its results.}

if assigned K`aut then
    Gp, Gf, Gphi := Explode(K`aut);
    return Gp, Gf, Gphi;
end if;

assert BaseRing(K) eq Rationals();
/* Special case of QQ  */
if IsQQ(K) then
    Gp, Gf, Gphi := AutomorphismGroup(K);
    K`aut := [* Gp, 0, Gphi *];
    return Gp, 0, Gphi ;
end if;

/* Determine all roots of the defining polynomial */
rts := RootsPari(DefiningPolynomial(K), K); S := Sym(#rts);

sigmas := [ ]; tups := [ ]; hs := [ ];
for rt in rts do
    /* Automorphism as map */
    h := hom< K -> K | rt >;
    /* Test if it fixes the base, and if so determine the associated
     * permutation representation */
    if h(K`base_gen) eq K`base_gen then
        L := [ ];
        for i in [1..#rts] do
            for j in [1..#rts] do
                if h(rts[i]) eq rts[j] then
                    Append(~L, j);
                end if;
            end for;
        end for;
        sigma := S ! L;
        Append(~sigmas, sigma);
        Append(~hs, h);
        Append(~tups, <sigma, h>);
    end if;
end for;
Gp := sub< S | sigmas >;
hs := Set(hs);
Gphi := map< Gp -> hs | tups >;

/* Assign and return */
K`aut := [* Gp, 0, Gphi *];
return Gp, 0, Gphi;

end intrinsic;
