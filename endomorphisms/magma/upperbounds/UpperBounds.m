// Section references throughout this file are to Costa, Mascot, Sijsling,
// Voight, arXiv:1705.09248, Math. Comp. 88 (2019) 1303-1339.
// exposes some of Honda--Tate theory necessary to produce tight upper bounds
// For more details see Section 7.2

intrinsic EndomorphismAlgebra(f::RngUPolElt) -> Tup
{The geometric endomorphism data of the abelian variety with Frobenius
 polynomial f over a finite field: the dimension dim_Q End(Abar), the degree of
 the least field over which all endomorphisms of Abar are defined, and the
 sorted isogeny decomposition [<m_i, m_i * deg(c_i), c_i>]}
    if IsMonic(f) then
        f := Reverse(f);
    end if;
    require Coefficient(f, 0) eq 1: "f is not a Weil polynomial";
    d := Degree(f);
    genus := Integers()!(Degree(f)/2);
    // The leading coefficient is q^genus, so we need the exact genus-th root.
    // IsPower/1 would instead demand a nontrivial perfect power, rejecting
    // every genus-1 Weil polynomial over a prime field.
    b, q := IsPower(Coefficient(f, d), genus);
    require b: "f is not a Weil polynomial";
    T := Parent(f).1;
    fof := TensorCharacteristicPolynomial(f,f);
    g := Evaluate(fof, ChangeRing(T, Rationals())/q);

    dimtotal := 0;
    fieldext := 1;

    for factorpower in Factorization(g) do
        factor, power := Explode(factorpower);
        b, k := IsCyclotomicPolynomial(factor);
        if b then
            dimtotal +:= power * Degree(factor);
            fieldext := LCM(fieldext, k);
        end if;
    end for;

    // fieldext is the least extension degree defining all geometric endomorphisms
    // (CMSV, Lemma 7.2.7(b)); factoring fext therefore gives the geometric
    // isogeny factors (Remark 7.2.13).
    fext := PowerCharacteristicPolynomial(f, fieldext);

    endo := Sort([
        <power, power * Degree(factor), factor>
        where factor, power := Explode(factorpower)
        : factorpower in Factorization(fext)]);

    return dimtotal, fieldext, endo;
end intrinsic;

// exposes some of the functionality mentioned in Section 7.3 and Section 7.4
intrinsic EndomorphismAlgebraEtaBound(frob_list::SeqEnum[RngUPolElt] : eta_char0 := false)
    -> BoolElt, MonStgElt, RngIntElt, RngIntElt, SeqEnum
{The eta and t narrowing step of Section 7.3, as
 <success, message, eta_char0, t, eta_lower> over one Frobenius polynomial per
 prime. eta_char0 bounds eta = sum_i e_i n_i^2 dim(A_i) of (7.3.16), NOT the
 algebra dimension sum_i n_i^2 e_i^2 [L_i:Q]; a quartic-CM surface has 2 and 4}
    require #frob_list ne 0: "frob_list must not be empty";
    g := Degree(frob_list[1]) div 2;
    if eta_char0 cmpeq false then
        eta := 4 * g * g;
    else
        eta := 2 * eta_char0;
    end if;

    t := g;
    eta_lower := [];
    for f in frob_list do
        dimtotal, _, endo := EndomorphismAlgebra(f);
        if dimtotal lt eta then
            eta := dimtotal;
            t := #endo;
            eta_lower := [];
        end if;
        if dimtotal eq eta then
            if #endo lt t then
                t := #endo;
                eta_lower := [];
            end if;
            if #endo eq t then
                Append(~eta_lower, endo);
            end if;
        end if;
    end for;

    if #eta_lower eq 0 then
        return false,
               "We did not manage to find any prime where eta(A_p) = 2 * eta(A)",
               0, 0, [];
    end if;

    return true, "", eta div 2, t, eta_lower;
end intrinsic;

intrinsic EndomorphismAlgebraCenterBounds(eta::RngIntElt, t::RngIntElt, eta_lower::SeqEnum)
    -> BoolElt, MonStgElt, SeqEnum, RngIntElt
{The center-bounding step of Sections 7.3-7.4, as
 <success, message, output, total_dim> with output a sequence of
 <ejnj, njdimAj, Lj, RRj>. Fails when a factor has no greatest common subfield,
 since no single field then bounds its center. See also RealRepresentationBound}
    require #eta_lower gt 0: "eta_lower must not be empty";
    QQT := PolynomialRing(Rationals());

    // Step 1: build the canonical multiset of (m, m*deg(h)) shapes from the
    // first prime; assert every other prime gives the same multiset.
    multiset0 := Sort([<x, y> where x, y, _ := Explode(elt) : elt in eta_lower[1]]);
    for endo in eta_lower do
        ms := Sort([<x, y> where x, y, _ := Explode(elt) : elt in endo]);
        if ms ne multiset0 then
            return false,
                   "We only managed to find an upper bound for eta. If the upper bound for eta indeed is eta, then the number of factors is a strict upper bound",
                   [], 0;
        end if;
    end for;

    // Step 2: group factor polynomials by (m, m*deg(h)) pair, indexed by prime.
    // Coerce to Q[T] but skip Polredabs at this stage: SubfieldsPolynomials
    // applies it internally where it matters, and the Polredabs fallback (used
    // when PARI/gp is unavailable) destroys integer coefficients here.
    pair_set := SequenceToSet(multiset0);
    frob_factors := AssociativeArray();
    for pair in pair_set do
        frob_factors[pair] := [[QQT |] : i in [1..#eta_lower]];
    end for;
    for i in [1..#eta_lower] do
        for elt in eta_lower[i] do
            x, y, hpj := Explode(elt);
            Append(~frob_factors[<x, y>][i], QQT ! hpj);
        end for;
    end for;

    // Step 3: for each pair, bound the center via FieldIntersectionMatrix and
    // assemble the output tuple.
    output := [];
    total_dim := 0;
    for pair in pair_set do
        ejnj := pair[1];
        njdimAj := pair[2] div 2;
        L := FieldIntersectionMatrix(frob_factors[pair]);
        for Lj in L do
            // A center-bound candidate must contain every other candidate. With
            // eta and t correct, CMSV Corollary 7.4.4 puts the true center among
            // these candidates. Incomparable maximal candidates do not certify a
            // bound; additional primes may resolve the ambiguity.
            Ljmax_poly := Lj[1];
            if IsZero(Ljmax_poly) then
                return false,
                       "The common subfields of a factor have no unique maximal element, so its center is not bounded by any single field",
                       [], 0;
            end if;
            if Degree(Ljmax_poly) eq 1 then
                Ljmax := RationalsAsNumberField();
            else
                Ljmax := NumberField(Ljmax_poly);
            end if;
            RRj := RealRepresentationString(njdimAj, Ljmax, ejnj);
            Append(~output, <ejnj, njdimAj, Lj, RRj>);
            total_dim +:= ejnj^2 * Degree(Ljmax);
        end for;
    end for;

    return true,
           "We have putatively computed eta and t. Under this assumption, we bounded the corresponding centers.",
           output, total_dim;
end intrinsic;

intrinsic EndomorphismAlgebraUpperBound(frob_list::SeqEnum[RngUPolElt] : eta_char0 := false)
    -> BoolElt, MonStgElt, RngIntElt, RngIntElt, SeqEnum, RngIntElt
{The Section 7 upper bound, as
 <success, message, eta_char0, t, output, total_dim>. The eta bound is
 unconditional once the eta step succeeds (Corollary 7.3.19(a)); the factor and
 center bounds assume eta and t are correct, and exact centers need Mumford-Tate}
    ok, msg, eta_c, t, eta_lower := EndomorphismAlgebraEtaBound(
        frob_list : eta_char0 := eta_char0);
    if not ok then
        return false, msg, 0, 0, [], 0;
    end if;
    ok2, msg2, output, total_dim := EndomorphismAlgebraCenterBounds(2 * eta_c, t, eta_lower);
    if not ok2 then
        return false, msg2, eta_c, t, [], 0;
    end if;
    return true, msg2, eta_c, t, output, total_dim;
end intrinsic;

intrinsic RealRepresentationBound(frob_list::SeqEnum[RngUPolElt]) -> SeqEnum
{Convenience wrapper extracting just the RR-representation strings for each
 simple factor. Returns an empty sequence when the upper bound cannot be
 established. Mirrors Sage's RR_upper_bound.}
    ok, _, _, _, output, _ := EndomorphismAlgebraUpperBound(frob_list);
    if not ok then
        return [];
    end if;
    return [tup[4] : tup in output];
end intrinsic;

intrinsic EndomorphismAlgebraUpperBound(C::Crv, B::RngIntElt : eta_char0 := false)
    -> BoolElt, MonStgElt, RngIntElt, RngIntElt, SeqEnum, RngIntElt
{The Section 7 upper bound for the curve C, from the L-polynomials at good
 primes below B. Any B ge 1 is allowed: with no usable prime below B the bound
 fails rather than raising, per CLV (arXiv:1906.02803) Algorithm 5.1}
    frob_list := [pair[2] : pair in LPolynomials(C, B)];
    if #frob_list eq 0 then
        return false,
               "No prime below the bound yielded a usable L-polynomial",
               0, 0, [], 0;
    end if;
    return EndomorphismAlgebraUpperBound(frob_list : eta_char0 := eta_char0);
end intrinsic;

intrinsic RealRepresentationBound(C::Crv, B::RngIntElt) -> SeqEnum
{Curve-level overload of RealRepresentationBound. Returns an empty sequence
 when no prime below B has good reduction, matching the failure convention of
 the SeqEnum form}
    frob_list := [pair[2] : pair in LPolynomials(C, B)];
    if #frob_list eq 0 then
        return [];
    end if;
    return RealRepresentationBound(frob_list);
end intrinsic;

/* Needs to be fixed
// exposes some of the functionality mentioned in Section 7.3 and Section 7.4
intrinsic EndomorphismAlgebraUpperBound(frob_list::SeqEnum[RngUPolElt] : eta_char0 := false) -> Tup
  {given a list of Frobenius polynomials return a Tuple ...}
  degrees := {Degree(f) : f in frob_list};
  require #degrees eq 1: "the Frobenius should have all degree 2*genus";
  g := Integers()!(Degree(frob_list[1])/2);

  if eta_char0 cmpeq false then
    eta := 4*g*g; // max value for eta in finite characteristic
  else
    eta := 2*eta_char0;
  end if;


  t := g; // max value for t
  eta_lower := [];
  for f in frob_list do
    dimtotal, _, endo :=EndomorphismAlgebra(f);
    if dimtotal lt eta then // erase all previous upper bounds
      eta := dimtotal;
      t := #endo;
      eta_lower := [];
    end if;
    if dimtotal eq eta then // Corollary 7.3.19
      if #endo lt t then
        t := #endo;
        eta_lower := [];
      end if;
      if #endo eq t then
        Append(~eta_lower, endo);
      end if;
    end if;
  end for;

  if #eta_lower eq 0 then
        return false, "We did not manage to find any prime where eta(A_p) = 2 * eta(A)";
  end if;

  eta_char0 := Integers()!(eta/2);
  multiset_char0 := Multiset([<x, y> where x,y,_ := Explode(elt) : elt in eta_lower[1]]);
  frob_factors := AssociativeArray();
  for pair in MultisetToSet(multiset_char0) do
    frob_factors[pair] := [];
  end for;


  for endo in eta_lower do
    multiset := Multiset([<x, y> where x,y,_ := Explode(elt) : elt in endo]);
    if multiset_char0 ne multiset then
      // we only managed to bound eta
      message := "We only managed to find an upper bound for eta.";
      message cat:= " If the upper bound for eta indeed is eta, then the number of factors is a strict upper bound";
      return false, message, eta_char0, t;
    end if;
    for elt in endo do
      // endo[j] = mpj, mpj*deg(hpj), hpj
      x, y, hpj := Explode(elt);
      Append(~frob_factors[<x,y>], hpj);
    end for;
  end for;

  // it looks like we have a consistent upper bound for eta and t
  message := "We have putatively computed eta and t.";
  message cat:= " Under this assumption, we bounded the corresponding centers.";

  // We can now try to bound the center of each factor


  // FIXME
  return "a", "b";

end intrinsic;
*/
