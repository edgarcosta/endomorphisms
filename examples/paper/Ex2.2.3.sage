"""
 *  Some examples of endomorphism rings of genus 2 hyperelliptic curves
 *
 *  Copyright (C) 2016  J.R. Sijsling (sijsling@gmail.com)
 *
 *  Distributed under the terms of the GNU General License (GPL)
 *                  http://www.gnu.org/licenses/
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License along
 *  with this program; if not, write to the Free Software Foundation, Inc., 51
 *  Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
"""

# Add if no initilization script set:
load('../../Initialize.sage')

# Ambient ring:
F = QQ
CCSmall = magma.ComplexField(5)
R.<x> = PolynomialRing(F)

# Curve:
f = x^5 - x^4 + 4*x^3 - 8*x^2 + 5*x - 1
h = R(0)
X = mHyperellipticCurve(f, h)

print X
Endo = EndomorphismData(X, prec = 100, have_oldenburg = False)

print "Period matrix:"
P = magma.Transpose(Endo._P_)
print magma.ChangeRing(P, CCSmall)

print "Endomorphism over QQ (sqrt (2)):"
R.<t> = PolynomialRing(F)
K.<s> = NumberField(t^2 - 2)
overK = Endo.over_field(K)
endodict = overK.full()
M = endodict['representation'][1]['tangent']
MCC = endodict['representation'][1]['approx']
R = endodict['representation'][1]['homology']
print M
print R

#CC = magma.BaseRing(P)
#MCC = magma.ChangeRing(MCC, CC)
#RCC = magma.ChangeRing(R, CC)
#print magma.ChangeRing(MCC * P - P * RCC, CCSmall)
