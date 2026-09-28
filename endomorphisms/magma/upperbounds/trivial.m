function normalize_reduction_polynomial(f, expected_degree)
    coefficients := [];
    for coefficient in Coefficients(f) do
        integral, value := IsCoercible(Integers(), coefficient);
        if not integral then
            error "polynomials must have integral coefficients";
        end if;
        Append(~coefficients, value);
    end for;

    QT := PolynomialRing(Rationals());
    f := QT ! coefficients;
    degree := Degree(f);
    if not (degree ge 2 and IsEven(degree) and
            (expected_degree eq 0 or degree eq expected_degree)) then
        error "polynomials must have common positive even degree";
    end if;

    if IsMonic(f) then
        P := f;
        L := QT ! Reverse(f);
    else
        if Coefficient(f, 0) ne 1 then
            error "each polynomial must have constant term 1 or be monic";
        end if;
        L := f;
        P := QT ! Reverse(f);
    end if;

    genus := degree div 2;
    is_power, q := IsPower(Integers() ! Coefficient(P, 0), genus);
    if not (is_power and q gt 1) then
        error "polynomial does not determine a prime-power field size";
    end if;
    is_prime_power, _, _ := IsPrimePower(q);
    if not is_prime_power then
        error "polynomial does not determine a prime-power field size";
    end if;
    if not forall{i : i in [0..genus] |
            Coefficient(L, degree - i) eq
                q^(genus - i) * Coefficient(L, i)} then
        error "polynomial does not satisfy the Frobenius functional equation";
    end if;

    return L, P, degree;
end function;

function reduction_record(f, expected_degree, geometric)
    L, P, degree := normalize_reduction_polynomial(f, expected_degree);
    squarefree := IsSquarefree(P);
    irreducible := squarefree and IsIrreducible(P);
    dimension := 0;
    if geometric and squarefree then
        dimension, _, _ := EndomorphismAlgebra(L);
    end if;
    return <L, P, squarefree, irreducible, dimension>, degree;
end function;

function reduction_is_eligible(record, degree, geometric)
    return record[3] and (not geometric or record[5] eq degree);
end function;

function witness_factor(anchor_field, record, degree, geometric)
    if not reduction_is_eligible(record, degree, geometric) then
        return false, 0;
    end if;
    MT := PolynomialRing(anchor_field);
    for factor_power in Factorization(record[2]) do
        h := factor_power[1];
        if IsIrreducible(MT ! h) then
            return true, h;
        end if;
    end for;
    return false, 0;
end function;

function trivial_endomorphism_certificate(records, sources, degree, geometric)
    if #records lt 2 then
        vprint UpperBounds: "At least two reductions are needed";
        return false, _;
    end if;

    anchor_index := 0;
    for i in [1..#records] do
        if records[i][4] and
                (not geometric or records[i][5] eq degree) then
            anchor_index := i;
            break;
        end if;
    end for;
    if anchor_index eq 0 then
        vprint UpperBounds: "No eligible irreducible anchor was found";
        return false, _;
    end if;

    anchor := records[anchor_index];
    anchor_field := NumberField(anchor[2]);
    for witness_index in [1..#records] do
        if witness_index eq anchor_index then
            continue;
        end if;
        found, h := witness_factor(
            anchor_field, records[witness_index], degree, geometric);
        if found then
            return true, <geometric, sources[anchor_index],
                sources[witness_index], anchor[2],
                records[witness_index][2], h>;
        end if;
    end for;

    vprint UpperBounds: "No eligible witness factor stayed irreducible";
    return false, _;
end function;

intrinsic CanProveTrivialEndomorphismRing(
    frob_list::SeqEnum[RngUPolElt] : Geometric := false) -> BoolElt, Tup
{A sufficient test and witness proving that the endomorphism ring is Z over
the base field from the exact good-reduction polynomials in frob_list. True
proves the claim; false means it was not established. The parameter Geometric
(default false) instead tests over an algebraic closure}
    if IsEmpty(frob_list) then
        vprint UpperBounds: "No reduction polynomials were supplied";
        return false, _;
    end if;

    records := [* *];
    degree := 0;
    for f in frob_list do
        record, degree := reduction_record(f, degree, Geometric);
        Append(~records, record);
    end for;
    return trivial_endomorphism_certificate(
        records, [1..#records], degree, Geometric);
end intrinsic;

intrinsic CanProveTrivialEndomorphismRing(
    C::Crv, B::RngIntElt : Geometric := false) -> BoolElt, Tup
{A sufficient test and witness proving that the Jacobian of C has endomorphism
ring Z over the base field, using exact good-reduction polynomials at primes
below B. True proves the claim; false means it was not established. The
parameter Geometric (default false) instead tests over an algebraic closure}
    require BaseRing(C) cmpeq Rationals(): "C must be defined over Q";
    require B ge 1: "B must be at least 1";
    genus := Genus(C);
    require genus ge 1: "C must have positive genus";

    records := [* *];
    sources := [];
    degree := 2 * genus;
    p := 2;
    while p lt B do
        if p^genus le 10^6 then
            try
                L := LPolynomial(ChangeRing(C, GF(p)) : Al := "Naive");
                record, degree := reduction_record(L, degree, Geometric);
                Append(~records, record);
                Append(~sources, p);
                found, certificate := trivial_endomorphism_certificate(
                    records, sources, degree, Geometric);
                if found then
                    return true, certificate;
                end if;
            catch e
                vprintf UpperBounds: "p = %o skipped: exact reduction unavailable\n", p;
            end try;
        else
            vprintf UpperBounds: "p = %o skipped: p^g exceeds 10^6\n", p;
        end if;
        p := NextPrime(p);
    end while;
    return false, _;
end intrinsic;
