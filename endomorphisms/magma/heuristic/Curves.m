/***
 *  Properties of curves
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa       (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo   (davide.lombardo@unipi.it)
 *            Jeroen Sijsling   (jeroen.sijsling@uni-ulm.de)
 *            Andrew Sutherland (drew@math.mit.edu)
 *
 *  See LICENSE.txt for license details.
 */

declare attributes Crv : plane_model, period_matrix, riesrf, geo_endo_rep_CC, geo_endo_rep, base_endo_rep, base_point, ghpols;
forward Homogenize;


intrinsic PlaneCurve(F::RngMPolElt) -> Crv
{Returns the plane curve defined by F, which can be given affinely or
projectively.}

R := Parent(F); gens := GeneratorsSequence(R);
if #gens eq 3 then
    return Curve(Scheme(ProjectiveSpace(R), F));
end if;
S := PolynomialRing(BaseRing(R), 3);
x := S.1; y := S.2; z := S.3;
Fhom := S ! (z^(Degree(F)) * Evaluate(F, [ x/z, y/z ]));
return Curve(Scheme(ProjectiveSpace(S), Fhom));

end intrinsic;


function Homogenize(fs)
// Given polys f(x0,x1,...x_(n-1)) in n (unweighted) variables returns homogeneous polynomial g(x_0,...,x_n) = x_n^deg(f)*f(x0/xn,...,x_(n-1)/ x_n) in n+1 variables.
    R := Universe(fs);  F:= BaseRing(R);
    assert Set(VariableWeights(R)) eq {1};
    n := Rank(R)+1;
    S := PolynomialRing(F,n);
    K := FieldOfFractions(S);
    return [S!(K.n^(Degree(f))*Evaluate(f,[K.i/K.n:i in [1..n-1]])): f in fs];
end function;


intrinsic GHCurve(q::RngMPolElt,f::RngMPolElt) -> Crv
{ Given conic q(x,y,z) and homogeneous poynomial f(x,y,z) of degree >= 3, returns curve [q(x,y,z)=0,w^2=f(z,y,x)] in [2,1,1,1] weighted        projective space. }
// Written by Andrew Sutherland
    require Parent(q) eq Parent(f): "Polynomials must lie in the same polynomial ring.";
    if Rank(Parent(q)) eq 2 then 
      R0 := PolynomialRing(BaseRing(Parent(q)),3);
      q := Homogenization(Evaluate(q,[R0.1,R0.2]),R0.3);
      f := Homogenization(Evaluate(f,[R0.1,R0.2]),R0.3);
    end if;
    R := Parent(q); F := BaseRing(R);
    require VariableWeights(R) eq [1,1,1]: "Polynomials must lie in an unweighted polynomial ring of rank 2 or 3.";
    require IsHomogeneous(q) and Degree(q) eq 2: "First input should be a conic.";
    require IsHomogeneous(f) and IsEven(Degree(f)) and Degree(f) ge 4: "Second input should be bivariate or homogeneous trivariate polynomial  whose homogenization is of even degree at least 4.";
    P<w,x,y,z> := ProjectiveSpace(F,[2,1,1,1]);
    return Curve(P,[Evaluate(q,[x,y,z]),w^2-Evaluate(f,[x,y,z])]);
end intrinsic;


intrinsic IsGHCurve(C::Crv) -> BoolElt, RngMPolElt, RngMPolElt
{ Determines whether the curve C is (exactly, not just isomorphic to a curve) of the form [q(x,y,z)=0,w^2=f(x,y,z)] in [2,1,1,1]-weighted      projective space.  If so, also returns q an f. }
// Written by Andrew Sutherland
    if VariableWeights(CoordinateRing(Ambient(C))) ne [2,1,1,1] then return false,_,_; end if;
    p := Sort(DefiningPolynomials(C),func<a,b|Degree(a)-Degree(b)>);
    if Degree(p[1]) ne 2 or Degree(p[2]) ne 4 then return false,_,_; end if;
    q,h := Explode(p);
    R := Parent(h);  F := BaseRing(R);
    f := R.1^2 - h;
    for g in [q,f] do if Evaluate(g,[0,R.2,R.3,R.4]) ne g then return false,_,_; end if; end for;
    R<x,y,z> := PolynomialRing(F,3);
    return true,Evaluate(q,[0,x,y,z]), Evaluate(f,[0,x,y,z]);
end intrinsic;


intrinsic CurveType(X::Crv) -> MonStgElt
{Returns a string that describes the type of curve that X belongs to, which is
one of "hyp", "genhyp", "plane" and "gen".}

if Type(X) eq CrvHyp then
    return "hyp";
elif IsGHCurve(X) then
    _, q, f := IsGHCurve(X);
    X`ghpols := [q, f];
    return "genhyp";
elif Type(X) eq CrvPln then
    return "plane";
else
    return "gen";
end if;

end intrinsic;


intrinsic ChangeRingCurve(X::Crv, h::.) -> Crv
{Returns X base changed by a morphism h whose domain is the base ring of X.}

if Type(X) eq CrvHyp then
    fK, hK := HyperellipticPolynomials(X);
    L := Codomain(h); S := PolynomialRing(L);
    if fK eq 0 then fL := S ! 0; else fL := &+[ h(Coefficient(fK, d))*S.1^d : d in [0..Degree(fK)] ]; end if;
    if hK eq 0 then hL := S ! 0; else hL := &+[ h(Coefficient(hK, d))*S.1^d : d in [0..Degree(hK)] ]; end if;
    return HyperellipticCurve(fL, hL);

elif Type(X) eq CrvPln then
    FK := DefiningPolynomials(X)[1];
    L := Codomain(h); S := PolynomialRing(L, 3);
    FL := &+[ h(MonomialCoefficient(FK, mon))*Monomial(S, Exponents(mon)) : mon in Monomials(FK) ];
    return PlaneCurve(FL);

elif IsGHCurve(X) then
    _, qK, fK := IsGHCurve(X);
    L := Codomain(h); S := PolynomialRing(L, 3);
    qL := &+[ h(MonomialCoefficient(qK, mon))*Monomial(S, Exponents(mon)) : mon in Monomials(qK) ];
    fL := &+[ h(MonomialCoefficient(fK, mon))*Monomial(S, Exponents(mon)) : mon in Monomials(fK) ];
    XL := GHCurve(qL, fL);
    XL`ghpols := [qL, fL];
    return XL;
else
	eqns := DefiningPolynomials(X);
	L := Codomain(h);
	S := PolynomialRing(L, Rank(Parent(eqns[1])));
	new_eqns := [];
	for FK in eqns do
	    FL := &+[ h(MonomialCoefficient(FK, mon)) * Monomial(S, Exponents(mon))
		      : mon in Monomials(FK)];
	    Append(~new_eqns, FL);
	end for;
	return Curve(Proj(S), new_eqns);
end if;
end intrinsic;


function CurveDescriptionHyperelliptic(X)

if Genus(X) eq 1 then
    desc := "ell";
else
    desc := "hyp";
end if;

K := BaseRing(X);
K_seq := FieldDescriptionExtra(K);
field := K_seq;

f, h := HyperellipticPolynomials(X);
f_seq := Eltseq(f); h_seq := Eltseq(h);
f_seq_seq := [ ElementDescriptionExtra(coeff) : coeff in f_seq ];
h_seq_seq := [ ElementDescriptionExtra(coeff) : coeff in h_seq ];
if #h_seq_seq eq 0 then
    coeffs := f_seq_seq;
else
    coeffs := [ f_seq_seq, h_seq_seq ];
end if;

return [* desc, field, coeffs *];

end function;


function CurveDescriptionPlane(X)

desc := "pln";

K := BaseRing(X);
K_seq := FieldDescriptionExtra(K);
field := K_seq;

f := DefiningPolynomials(X)[1];
mons := Monomials(f);
coeffsexps := [ ];
for mon in mons do
    coeff := MonomialCoefficient(f, mon);
    coeff_seq := ElementDescriptionExtra(coeff);
    exp := Exponents(mon);
    Append(~coeffsexps, [* coeff_seq, exp *]);
end for;

return [* desc, field, coeffsexps *];

end function;


intrinsic CurveDescription(X::Crv) -> List
{Returns a string description of the curve X.}

if Type(X) eq CrvHyp then
    return CurveDescriptionHyperelliptic(X);
elif Type(X) eq CrvPln then
    return CurveDescriptionPlane(X);
else
    error "No description for general curves yet";
end if;

end intrinsic;


intrinsic PlaneModel(X::Crv) -> .
{Returns a plane model for more complicated curves.}

if assigned X`plane_model then return X`plane_model; end if;
if Type(X) eq CrvHyp or Type(X) eq CrvPln then X`plane_model := X; return X`plane_model; end if;

F := DefiningPolynomial(AlgorithmicFunctionField(FunctionField(X)));
L := Parent(F); K := BaseRing(L);

R := PolynomialRing(BaseRing(K));
S := PolynomialRing(BaseRing(K), 2);
h := hom< R -> S | S.1 >;

F0 := S ! 0; coeffs := Coefficients(F);
for i in [1..#coeffs] do
    F0 +:= h(R ! coeffs[i])*S.2^(i - 1);
end for;

X`plane_model := PlaneCurve(F0);
return X`plane_model;

end intrinsic;


/* These functions are for the SageMath interface */
intrinsic HyperellipticCurveExtra(f::RngUPolElt, h::., prec::RngIntElt) -> Crv
{Returns the hyperelliptic curve over the rationals with precision prec defined
by f and h. Only relevant in the Sage interface.}

QQ := RationalsExtra(prec); RQQ := PolynomialRing(QQ);
return HyperellipticCurve(RQQ ! f, RQQ ! h);

end intrinsic;


intrinsic PlaneCurveExtra(F::RngMPolElt, prec::RngIntElt) -> Crv
{Returns the plane curve over the rationals with precision prec defined by F,
which can be given affinely or projectively. Only relevant in the Sage interface.}

QQ := RationalsExtra(prec); RQQ := PolynomialRing(QQ, #GeneratorsSequence(Parent(F)));
return PlaneCurve(RQQ ! F);

end intrinsic;
