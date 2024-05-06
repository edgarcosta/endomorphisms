intrinsic RootsPari(f::RngUPolElt, K::Fld) -> .
{Given a polynomial f and a field K, finds all roots of f in K.}

  assert Type(BaseRing(f)) in [FldRat, RngInt];
  assert Type(BaseRing(K)) in [FldRat, RngInt];
  rts := [];
  try
    g := DefiningPolynomial(K);
    cmd := Sprintf(
    "{f = Pol(Vecrev(%o),'x); g = Pol(Vecrev(%o),'y); K = nfinit(g); print1(apply(z->vector(poldegree(g),i, polcoeff(z,i-1)),lift(nfroots(g,f))))}",
    Coefficients(f), Coefficients(g));
    s := Pipe("gp -q -D timer=0", cmd);
    rts := [ K ! rt : rt in eval(s) ];
    rts := [ rt : rt in rts | Evaluate(f, rt) eq 0 ];
  catch e
    vprintf EndoFind : "WARNING: Need gp at command line for RootsPari!\n";
    for pair in Roots(f, K) do
      for i in [1..pair[2]] do
        Append(~rts, pair[1]);
      end for;
    end for;
  end try;
  return rts;

end intrinsic;


intrinsic HasRootPari(f::RngUPolElt, K::Fld) -> .
{Given a polynomial f and a field K, determines if f has a root in K.}

return #RootsPari(f, K) ne 0;

end intrinsic;


intrinsic FactorizationPari(f::RngUPolElt, K::Fld) -> .
{Given a polynomial f and a field K, finds the factorization of f over K.}

  assert Type(BaseRing(f)) in [FldRat, RngInt];
  assert Type(BaseRing(K)) in [FldRat, RngInt];
  facs := [ ];
  try
  g := DefiningPolynomial(K);
  cmd := Sprintf(
  "{f = Pol(Vecrev(%o),'x); g = Pol(Vecrev(%o),'y); K = nfinit(g); print1(apply(h->apply(c->vector(poldegree(g),i,polcoeff(c,i-1)),lift(Vecrev(h))),nffactor(K,f)[,1]~))",
  Coefficients(f), Coefficients(g));
  s := Pipe("gp -q -D timer=0", cmd);

  R := PolynomialRing(K);
  seqs := eval(s);
  facs := [ &+[ (K ! seq[i])*R.1^(i - 1) : i in [1..#seq] ] : seq in seqs ];
  catch e
    vprintf EndoFind : "WARNING: Need gp at command line for FactorizationPari!\n";
    for pair in Factorization(f, K) do
      for i in [1..pair[2]] do
        Append(~facs, pair[1]);
      end for;
    end for;
  end try;
  return facs;
end intrinsic;


intrinsic SplittingFieldPari(f::RngUPolElt : roots := false) -> .
{Splitting field of f calculated using Pari.}

  //return SplittingField(f);
  f0 := f;
  F := BaseRing(f);
  assert F eq Rationals();
  try
    cmd := Sprintf(
    "{f = Pol(Vecrev(%o),'x); print1(nfsplitting(f))",
    Coefficients(f), Coefficients(f));
    s := Pipe("gp -q -D timer=0", cmd);
    R<x> := PolynomialRing(F);
    L := Polredbestabs(NumberField(eval(s)));
    f := R ! DefiningPolynomial(L);
  catch e
    vprintf EndoFind : "WARNING: Need gp at command line for SplittingFieldPari!\n";
    f := DefiningPolynomial(SplittingField(f));
  end try;
  L := NumberFieldExtra(f);
  if not roots then
    return L;
  end if;
  return L, RootsPari(f0, L);
end intrinsic;


intrinsic SplittingFieldPari(fs::SeqEnum[RngUPolElt] : roots := true) -> .
{Splitting field of f calculated using Pari.}

  L := SplittingField(&*fs);
  if not roots then
    return L;
  end if;
  rss := [ RootsPari(f, L) : f in fs ];
  return L, rss;

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
