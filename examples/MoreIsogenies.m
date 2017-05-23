// load "Initialize.m";
AttachSpec("../spec");
//AttachSpec("oldenburg/periodmatrices/spec");

function TestIsogeny()
Cpx := ComplexFieldExtra(500);

    R<x> := PolynomialRing(Cpx);


    f := x^8+3*x^2+2;
    h := 0;
    g := 4*f+h^2;
    Crv := HyperellipticCurve(g);

    J := AnalyticJacobian(Crv);
    a,m := EndomorphismRing(J);
    P := BigPeriodMatrix(J);


    f2 := x^4+3*x+2;
    f2 := x^5+3*x^2+2*x;
    h2 := 0;
    g2 := 4*f2+h2^2;
    Crv2 := HyperellipticCurve(g2);

    J2 := AnalyticJacobian(Crv2);
    a2,m2 := EndomorphismRing(J2);
    P2 := BigPeriodMatrix(J2);


GeometricIsogenyBasisApproximations(Transpose(P), Transpose(P2));

return true;
end function;


function MatrixRatio(M1,M2);
    for i in [1..#Rows(M1)] do
        for j in [1..#Rows(Transpose(M1))] do
            if M2[i][j] ne 0 then
                return M1[i][j]/M2[i][j];
            end if;
        end for;
    end for;
end function;


function RationalizeMatrix(M);
	N := Matrix(Rationals(), 0, #Rows(Transpose(M)), []);
	for r in Rows(M) do
		N := VerticalJoin(N, Matrix(Rationals(), 1, #Rows(Transpose(M)), [ BestApproximation(Re(r[i]),10^10) : i in [1..#Rows(Transpose(M))] ]) );
	end for;
	return N;
end function;

function PullbackPolarization()
    Cpx := ComplexFieldExtra(300);
    C<i> := ComplexFieldExtra(4);
    R<x> := PolynomialRing(Cpx);

    PrincipalPolarization := Matrix(Integers(),6,6,[0,0,0,1,0,0,  0,0,0,0,1,0,  0,0,0,0,0,1,  -1,0,0,0,0,0,  0,-1,0,0,0,0,  0,0,-1,0,0,0]);




    // examples in which the map C \to E has degree 2. They really sort of work, except that TwistingFactor() gives all sort of errors

    f := x^8-x^2-1; // Runtime error in 'Characteristic': Bad argument types
		    // Argument types given: PowSeqEnum

    f := x^8+3*x^2-1; // same as previous one

    f := x^8+x^6+1;	// Runtime error in sequence construction: No valid universe containing all elements

 


    f := x^4+2*x^3+2*x^2+x;
    h := x^4+x^3+1; // No valid universe...

    f := -x^4;
    h := x^4+x^3+x+1; // No valid universe...



    f := x^7-2*x^6+2*x^5-x^4-x^3+x^2-x;
    h := x^2+1; // IsNormal complains

    f := x^7 + x^6 + x^5 + x^3 + x^2 + x;
    h := x^4 + x^2 + 1; // PolynomializeMatrix(As[1]); Runtime error in '[]': Illegal null sequence




    f := x^8+3*x^2+2;
    h := 0;




    f := -x^6-x^5-2*x^4-x^3-x^2;
    h := x^4+x^2+1; // PolynomializeMatrix(As[1]); Runtime error in '[]': Illegal null sequence

   h := 0;
    f := x^8+2*x^4-3*x^2-1; // Runtime error in '[]': Illegal null sequence




    // degree 3 map towards EC
    f := -2*x^7 - 4*x^6 + 3*x^4 + x^3 - 2*x^2 - x;
    h := x^2 + x + 1;

    f := x^7 - x^6 + 2*x^4 - 3*x^3 + 2*x^2 - x;
    h := x^4 + x^2 + 1;


    // covering of degree 5! The smallest field of definition of the genus 2 curve appears to be BIG - do we really believe that? Yes, I guess we do.
    f := -x^7+2*x^6-2*x^5+3*x^4+x^3+2*x^2+x;
    h := x^4+x^3+x+1;




    g := 4*f+h^2;
    Crv := HyperellipticCurve(g);


    J := AnalyticJacobian(Crv);
    a,m := EndomorphismRing(J);
    P := Transpose(BigPeriodMatrix(J));


    "== Genus 3 curve ==";
    print "Discriminant = ", BestApproximation(Discriminant(Crv), 10^10);


    // PSmall := Submatrix(P,[1,2,3],[4,5,6])^(-1)*Submatrix(P,[1,2,3],[1,2,3]);
    // print "Original small matrix", ChangeRing(PSmall,C);

    mRat := RationalizeMatrix(Transpose(m[1]));

    if mRat eq IdentityMatrix(Rationals(),3) then
	    mRat := RationalizeMatrix(Transpose(m[2]));
    end if;


    ev := Eigenvalues(mRat);

    for v in ev do
	if v[2] eq 2 then
		ECEigenvalue := v[1];
	end if;
    end for;



    idem1 := mRat - ECEigenvalue;
    d := MatrixRatio(idem1^2, idem1);
    idem2 := d-idem1;
    d := Abs(d);
    /*
    "== idem1 ==";
    print idem1;
    "== idem2 ==";
    print idem2;
    */
    "== Degree of the map C to E ==";
    print d;

    idem1 := ChangeRing(idem1, Cpx);
    idem2 := ChangeRing(idem2, Cpx);



    LatticeEC := P*idem1;
    LatticeSurface := P*idem2;
    // print ChangeRing(SplitMatrix(LatticeEC),C);

    "==== Surface lattice ====";
    print ChangeRing(LatticeSurface, C);

    "==== Elliptic curve lattice ====";
    print ChangeRing(LatticeEC, C);
    
    "==== Preliminary saturation ====";

    /*LatticeEC := CombineMatrix( SubmatrixOfRank(SplitMatrix(LatticeEC), 2, "Rows"));
    LatticeSurface := CombineMatrix( SubmatrixOfRank(SplitMatrix(LatticeSurface), 4, "Rows"));*/

    LatticeEC := SubmatrixOfRank(SplitMatrix(LatticeEC), 2 : ColumnsOrRows := "Rows");
    LatticeSurface := SubmatrixOfRank(SplitMatrix(LatticeSurface), 4 : ColumnsOrRows := "Rows");

    LatticeEC := CombineMatrix( SaturateLattice( SplitMatrix(P*idem1) , LatticeEC) , Cpx ) ;
    LatticeSurface := CombineMatrix( SaturateLattice( SplitMatrix(P*idem2) , LatticeSurface), Cpx );

    PSplit := SplitMatrix(P);

    BaseChangeMatrixEC := Matrix(Integers(), 0, 6, []);
    BaseChangeMatrixSurface := Matrix(Integers(), 0, 6, []);

    BaseChangeMatrixComplete := Matrix(Integers(), 0, 6, []);

    "==== Elliptic Curve lattice in terms of the original lattice ====";
    for i in Rows(LatticeEC) do
        iSplit := SplitMatrix(i);
        vec := Transpose(PSplit)^(-1) * Transpose(iSplit);
	BaseChangeMatrixEC := VerticalJoin(BaseChangeMatrixEC, Matrix(Integers(),1,6,[BestApproximation(vec[j][1],10^10) : j in [1..6]]));
	BaseChangeMatrixComplete := VerticalJoin(BaseChangeMatrixComplete, Matrix(Integers(),1,6,[BestApproximation(vec[j][1],10^10) : j in [1..6]]));  
    end for;
    print BaseChangeMatrixEC;

    "==== Pullback of canonical polarization to elliptic curve ====";
    print BaseChangeMatrixEC * PrincipalPolarization * Transpose(BaseChangeMatrixEC);
    
    "==== Abelian surface lattice in terms of the original lattice ====";
    for i in Rows(LatticeSurface) do
        iSplit := SplitMatrix(i);
        vec := Transpose(PSplit)^(-1) * Transpose(iSplit);
	BaseChangeMatrixSurface := VerticalJoin(BaseChangeMatrixSurface, Matrix(Integers(),1,6,[BestApproximation(vec[j][1],10^10) : j in [1..6]]));
	BaseChangeMatrixComplete := VerticalJoin(BaseChangeMatrixComplete, Matrix(Integers(),1,6,[BestApproximation(vec[j][1],10^10) : j in [1..6]]));  
    end for;
    print BaseChangeMatrixSurface;

    "==== Pullback of canonical polarization to abelian surface ====";
    PullbackPolarizationSurface := BaseChangeMatrixSurface * PrincipalPolarization * Transpose(BaseChangeMatrixSurface);
    print PullbackPolarizationSurface ;

    "==== Frobenius form ====";
    TypePolarization, baseChange := FrobeniusFormAlternating(Matrix(Integers(), PullbackPolarizationSurface));
    print TypePolarization;

    "==== Canonizing Abelian Surface Lattice ====";
    LatticeSurface := Matrix(Cpx, baseChange) * LatticeSurface ;

    "==== Repeating computations with canonized lattice ====";
    BaseChangeMatrixSurface := Matrix(Integers(), 0, 6, []);
    for i in Rows(LatticeSurface) do
        iSplit := SplitMatrix(i);
        vec := Transpose(PSplit)^(-1) * Transpose(iSplit);
	BaseChangeMatrixSurface := VerticalJoin(BaseChangeMatrixSurface, Matrix(Integers(),1,6,[BestApproximation(vec[j][1],10^10) : j in [1..6]]));
    end for;
    PullbackPolarizationSurface := BaseChangeMatrixSurface * PrincipalPolarization * Transpose(BaseChangeMatrixSurface);
    
    CanonizedBaseChangeMatrix := BaseChangeMatrixSurface;

    "==== Canonized lattice ====";
    print BaseChangeMatrixSurface;

    "==== Canonized polarization ====";
    print PullbackPolarizationSurface ;

    
    MinimalExtensionFound := d^10;


    "==== Trying all isogenies that factor [d] ====";
    dTorsionCart := CartesianPower([0..d-1], 4);
    dTorsion := {t : t in dTorsionCart};
    TorsionSubschemes := Subsets(dTorsion, 2);

    PSplit := SplitMatrix(P);
    LatticeSurfaceSplit := SplitMatrix(LatticeSurface);
    R := Parent(LatticeSurfaceSplit[1][1]);

    AlreadyEncountered := [];

    for tS in TorsionSubschemes do

        // print "Isogeny corresponding to", tS;
        
	for t in tS do
		iSplit := Matrix( R, 1, 6, [ &+ [ t[i]/d * LatticeSurfaceSplit[i][j] : i in [1..4] ] : j in [1..6] ] );
		EnlargedLattice := VerticalJoin ( SplitMatrix(LatticeSurface), iSplit );
	end for;

	IsogenousLattice := CombineMatrix( SaturateLattice (EnlargedLattice , SplitMatrix(LatticeSurface) ), Cpx );
	vec := Transpose(PSplit)^(-1) * Transpose(iSplit);

	// "==== Computations with isogenous lattice ====";

	BaseChangeMatrixSurface := Matrix(Rationals(), 0, 6, []);
	for i in Rows(IsogenousLattice) do
		iSplit := SplitMatrix(i);
		vec := Transpose(PSplit)^(-1) * Transpose(iSplit);
		BaseChangeMatrixSurface := VerticalJoin(BaseChangeMatrixSurface, Matrix(Rationals(),1,6,[BestApproximation(vec[j][1],10^10) : j in [1..6]]));
        end for;
	PullbackPolarizationSurface := BaseChangeMatrixSurface * ChangeRing(PrincipalPolarization, Rationals()) * Transpose(BaseChangeMatrixSurface);



	AlreadyEncounteredBool := false;
	for ae in AlreadyEncountered do
		if ae eq BaseChangeMatrixSurface then
			AlreadyEncounteredBool := true;
		end if;
	end for;

	if not AlreadyEncounteredBool then
	AlreadyEncountered := Append(AlreadyEncountered, BaseChangeMatrixSurface);
	
	/*
	print "== Lattice ==";
	BaseChangeMatrixSurface ;

	print "== Polarization ==";
	TypePolarization, baseChange := FrobeniusFormAlternating(Matrix(Integers(), PullbackPolarizationSurface));
	print TypePolarization;
	*/
	
	
	
	TypePolarization, baseChange := FrobeniusFormAlternating(Matrix(Integers(), PullbackPolarizationSurface));
	if TypePolarization eq d*Matrix(Integers(), 4,4, [0,0,1,0,  0,0,0,1,  -1,0,0,0,  0,-1,0,0] ) then
		/* print "== Lattice ==";
		BaseChangeMatrixSurface ;

		print "== Polarization ==";
		print TypePolarization; */

		PrincipallyPolarizedLattice := Matrix(Cpx, baseChange) * IsogenousLattice ;


		PEllHugeSplit := SplitMatrix(PrincipallyPolarizedLattice);
    
		PreliminaryLatticeMatrix := SubmatrixOfRank(PrincipallyPolarizedLattice, 2 : ColumnsOrRows := "Columns"); // extract g columns (i.e. decide which projection to use)
		    
		PreliminaryLatticeMatrix := Transpose(PreliminaryLatticeMatrix); // necessary before calling SaturateLattice
		PEllBig := Transpose(PrincipallyPolarizedLattice);

		PrincipallyPolarizedLatticeC2 := SaturateLattice(PEllBig, PreliminaryLatticeMatrix);

    		SmallPeriodMatrix := Submatrix(PrincipallyPolarizedLatticeC2, [1,2], [3,4])^(-1) * Submatrix(PrincipallyPolarizedLatticeC2, [1,2], [1,2]);
		// print "== Computing Rosenhain invariants ==";
		ri := RosenhainInvariants(SmallPeriodMatrix);

		RCpx<xCC> := PolynomialRing(Cpx);
		fCC := xCC * (xCC - 1) * &*[ xCC - rosenCC : rosenCC in ri ];
		g2sCC := G2Invariants(HyperellipticCurve(fCC));



		"==== Check: is there a geometric isogeny? ====";

		// print #Rows(P), #Rows(Transpose(P));
		PeriodMatrixQ := Transpose(BigPeriodMatrix(AnalyticJacobian(HyperellipticCurve(fCC))));

		GeometricIsogenyBasisApproximations(PeriodMatrixQ, P)[2];

		MaximalDegree := 0;
		for r in g2sCC do
			// print r;
			f := MinimalPolynomial(r, Integers()!25);
			MaximalDegree := Max(MaximalDegree, Degree(f));
			if Degree(f) eq 0 then
				MaximalDegree := d^10;
			end if;
			if MaximalDegree ge MinimalExtensionFound then
				break;
			end if;
		end for;
		if MaximalDegree lt MinimalExtensionFound then
			MinimalExtensionFound := MaximalDegree;
			BestPPLattice := PrincipallyPolarizedLatticeC2;
			BestSmallPeriodMatrix := SmallPeriodMatrix;
		end if;
		print "Degree of G2 invariants:", MaximalDegree;
		print "Current minimal degree of extension:", MinimalExtensionFound;
		if MinimalExtensionFound eq 1 then
			break;
		end if;


	end if;
	end if;

    end for;

    if MinimalExtensionFound ne d^10 then
	    "Found an isogeny to a Jacobian over an extension of degree", MinimalExtensionFound;
	    "Computing corresponding genus 2 curve";

            BestPPLattice := Transpose(BestPPLattice);
	    PPLattice := Matrix(Cpx, 4, 2, [BestPPLattice[3], BestPPLattice[4], BestPPLattice[1], BestPPLattice[2] ] );

	    K := Rationals();
	    SetInfinitePlace(K, InfinitePlaces(K)[1]);

	    return FactorReconstructG2(PPLattice, K);
    else
	    "Could not find an isogeny to a Jacobian over any small degree extension";
	    return false;
    end if;

end function;


// TestIsogeny();
PullbackPolarization();
