/***
 *  Sage version of raw Magma data
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@math.u-psud.fr)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */


intrinsic SagifyDescription(obj::.) -> MonStgElt
{Returns a conversion of the list description obj into a list that, when
evaluated in Sage, corresponds with the old list.}

case Type(obj):
    when List: if obj eq [* *] then return "[ ]"; else return "[" cat &cat[ SagifyDescription(x) cat ",": x in obj ] cat "]"; end if;
    when SeqEnum: if obj eq [ ] then return "[ ]"; else return "[" cat &cat[ SagifyDescription(x) cat ",": x in obj ] cat "]"; end if;
    when RngIntElt: return Sprint(obj);
    when FldRatElt: return Sprint(obj);
    when MonStgElt: return "'" cat Sprint(obj) cat "'";
end case;

end intrinsic;
