// exposes some of Honda--Tate theory necessary to produce tight upper bounds
// For more details see Section 7.2

intrinsic EndomorphismAlgebra(f::RngUPolElt) -> Tup
{Given a Frobenius polynomial  f of an Abelian variety f.
Return the triple. The first item is dim_Q ( End(Abar ).
The second item the degree of field ext where all endomorphism are defined.
The third item, a sorted list [ (m_i, m_i * deg(c_i), c_i) : i in [1..t] ] that represents the geometric isogeny decomposition over an algebraic closure, where
Abar = (A_1)^n_1 x ... x (A_k)^n_t
and
det(1 - T Frob^k | H^1(Ai^n_i)) = c_i (T)^m_i
}
    if IsMonic(f) then
        f := Reverse(f);
    end if;
    require Coefficient(f, 0) eq 1: "f is not a Weil polynomial";
    d := Degree(f);
    genus := Integers()!(Degree(f)/2);
    b, p, ga := IsPower(Coefficient(f, d));
    require b: "f is not a Weil polynomial";
    q := p^(Integers()!(ga/genus));
    T := Parent(f).1;
    fof := TensorCharacteristicPolynomial(f,f);
    g := Evaluate(fof, T/q);

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

    fext := PowerCharacteristicPolynomial(f, fieldext);

    endo := Sort([
        <power, power * factor.degree(), factor>
        where factor, power := Explode(factorpower)
        : factorpower in Factorization(fext)]);

    return dimtotal, fieldext, endo;
end intrinsic;


// exposes some of the functionality mentioned in Section 7.3 and Section 7.4
intrinsic EndomorphismAlgebraUpperBound(frob_list::SeqEnum[RngUPolElt] : eta_char0 := false) -> Tup
  {given a list of Frobenius polynomials return a Tuple ...}
  degrees := {Degree(f) : f in frob_list*};
  require #degrees eq 1, "the Frobenius should have all degree 2*genus";
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
  return "a", "b"

end intrinsic;
