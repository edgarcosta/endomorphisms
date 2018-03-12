AttachSpec("../../endomorphisms/magma/spec");

CC := ComplexFieldExtra(100);
RR := RealField(CC);

L := Matrix(RR, [[1,0,0],[0,0,1]]);
print "";
print InvertibleSubmatrix(L);
print "";
print InvertibleSubmatrix(Transpose(L));
print "";

print "";
print SubmatrixOfRank(L, 2 : ColumnsOrRows := "Columns");
print "";
print SubmatrixOfRank(L, 2 : ColumnsOrRows := "Rows");
print "";

L := Matrix(RR, [[1,0,0],[0,1,0]]);
M := Matrix(RR, [[1/4,0,0],[-1/2,1/5,0]]);
B, T, U := SaturateLattice(L, M : ColumnsOrRows := "Rows");
print "";
print B;
print "";
print T;
print "";
print U;

L := Matrix(RR, [[1,0],[0,1]]);
M := Matrix(RR, [[1/4,0,0],[-1/2,1/5,0]]);
B, T, U := SaturateLattice(L, M : ColumnsOrRows := "Columns");
print "";
print B;
print "";
print T;
print "";
print U;
print "";

exit;
