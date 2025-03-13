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
    genus := d div 2;
    b, p, ga := IsPower(Coefficient(f, d));
    require b: "f is not a Weil polynomial";
    q := p^(ga div genus);
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
        <power, power * Degree(factor), factor>
        where factor, power := Explode(factorpower)
        : factorpower in Factorization(fext)]);

    return dimtotal, fieldext, endo;
end intrinsic;
