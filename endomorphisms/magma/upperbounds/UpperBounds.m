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

function MonicReciprocal(f)
    g := Parent(f) ! Reverse(f);
    return g / LeadingCoefficient(g);
end function;

// Literal subfields of M are compared as Q-vector subspaces of M. Field-object
// equality is unsuitable: two separately constructed copies of the same
// literal subfield compare unequal, while distinct conjugate subfields must stay
// distinct.
function FieldSpace(M, K)
    V := VectorSpace(Rationals(), Degree(M));
    if Degree(K) eq 1 then
        return sub<V | V ! Eltseq(M ! 1)>;
    end if;
    return sub<V | [V ! Eltseq(M ! b) : b in Basis(K)]>;
end function;

// Candidate coefficient fields for hp inside the fixed anchor M. Products of
// irreducible factors are essential: a normic divisor can be reducible even
// when none of its irreducible factors is normic. The normic-divisor criterion
// of arXiv:1906.02803, Proposition 2.3, makes the complete divisor family
// downward closed: every subfield occurs as another coefficient field.
function NormicCandidates(M, hp)
    gp := MonicReciprocal(hp);
    fac := Factorization(PolynomialRing(M) ! gp);
    assert forall{f : f in fac | f[2] eq 1};
    irr := [f[1] : f in fac];
    out := [* *];
    for mask in [1 .. 2^#irr - 1] do
        h := Parent(irr[1]) ! 1;
        for i in [1 .. #irr] do
            if BitwiseAnd(mask, 2^(i - 1)) ne 0 then
                h *:= irr[i];
            end if;
        end for;
        if Degree(gp) mod Degree(h) ne 0 then
            continue;
        end if;
        L := sub<M | Coefficients(h)>;
        if Degree(L) * Degree(h) ne Degree(gp) then
            continue;
        end if;
        W := FieldSpace(M, L);
        K := Degree(L) eq 1 select RationalsAsNumberField() else L;
        if not exists{c : c in out | c[1] eq W} then
            Append(~out, <W, K>);
        end if;
    end for;
    return out;
end function;

function PossibilityToken(token_lists)
    encoded := [Join(tokens, "+") : tokens in token_lists];
    Sort(~encoded);
    return "oneof{" cat Join(encoded, "|") cat "}";
end function;

function CenterBounds(eta, t, eta_lower)
    QQT := PolynomialRing(Rationals());

    multiset0 := Sort([<x, y> where x, y, _ := Explode(elt) : elt in eta_lower[1]]);
    for endo in eta_lower do
        ms := Sort([<x, y> where x, y, _ := Explode(elt) : elt in endo]);
        if ms ne multiset0 then
            return false,
                   "We only managed to find an upper bound for eta. If the upper bound for eta indeed is eta, then the number of factors is a strict upper bound",
                   [], 0;
        end if;
    end for;

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

    output := [];
    total_dim := 0;
    for pair in Sort(SetToSequence(pair_set)) do
        ejnj := pair[1];
        njdimAj := pair[2] div 2;
        rows := frob_factors[pair];
        for hq in rows[1] do
            if Degree(hq) eq 1 then
                maximal := [* <1, RationalsAsNumberField()> *];
            else
                M := NumberField(MonicReciprocal(hq));
                if #rows eq 1 then
                    maximal := [* <FieldSpace(M, M), M> *];
                else
                    survivors := [* *];
                    for i in [2 .. #rows] do
                        at_prime := [* *];
                        for hp in rows[i] do
                            for candidate in NormicCandidates(M, hp) do
                                if not exists{c : c in at_prime |
                                              c[1] eq candidate[1]} then
                                    Append(~at_prime, candidate);
                                end if;
                            end for;
                        end for;
                        if i eq 2 then
                            survivors := at_prime;
                        else
                            survivors := [* c : c in survivors |
                                exists{d : d in at_prime | d[1] eq c[1]} *];
                        end if;
                    end for;

                    // Fix an embedding of the true centre into M. The
                    // normic-divisor criterion of arXiv:1906.02803,
                    // Proposition 2.3, puts that literal copy in every
                    // per-prime set, hence under a maximal survivor here.
                    maximal := [* c : c in survivors |
                        not exists{d : d in survivors |
                            c[1] ne d[1] and c[1] subset d[1]} *];
                end if;
            end if;

            d := Max([Degree(c[2]) : c in maximal]);
            token_lists := [];
            polynomials := [];
            for candidate in maximal do
                tokens := RealRepresentationString(
                    njdimAj, candidate[2], ejnj);
                if tokens notin token_lists then
                    Append(~token_lists, tokens);
                end if;
                if Degree(candidate[2]) eq 1 then
                    Append(~polynomials, QQT.1);
                else
                    Append(~polynomials, QQT ! DefiningPolynomial(
                        AbsoluteField(candidate[2])));
                end if;
            end for;
            RRj := #token_lists eq 1 select token_lists[1]
                   else [PossibilityToken(token_lists)];
            Append(~output, <ejnj, njdimAj, <d, polynomials>, RRj>);
            total_dim +:= ejnj^2 * d;
        end for;
    end for;

    return true,
           "We have putatively computed eta and t. Under this assumption, we bounded the corresponding centers.",
           output, total_dim;
end function;

intrinsic EndomorphismAlgebraCenterBounds(eta::RngIntElt, t::RngIntElt, eta_lower::SeqEnum)
    -> BoolElt, MonStgElt, SeqEnum, RngIntElt
{The center-bounding step of Sections 7.3-7.4, as
 <success, message, output, total_dim> with output a sequence of
 <ejnj, njdimAj, Lj, RRj>, where Lj gives a degree bound and the surviving
 candidate fields. See also RealRepresentationBound}
    require #eta_lower gt 0: "eta_lower must not be empty";
    return CenterBounds(eta, t, eta_lower);
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
    ok2, msg2, output, total_dim := CenterBounds(2 * eta_c, t, eta_lower);
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
