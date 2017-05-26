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
R.<x> = PolynomialRing(F)

# Curve:
f = 5*x^6 + 10*x^3 - 4*x + 1
h = R(0)
X = mHyperellipticCurve(f, h)

print X
Endo = EndomorphismData(X, prec = 300, have_oldenburg = False)

print "Endomorphism over QQ (sqrt (5)):"
R.<t> = PolynomialRing(F)
K.<s> = NumberField(t^2 - 5)
overK = Endo.over_field(K)
endodict = overK.full()
print endodict['representation'][0]['tangent'] + endodict['representation'][1]['tangent']
print endodict['representation'][0]['homology'] + endodict['representation'][1]['homology']
