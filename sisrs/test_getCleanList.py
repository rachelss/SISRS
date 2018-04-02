#!/usr/bin/env python2
from specific_genome import getCleanList

def test_answer():
    #Test getCleanList (Shared function)
    assert getCleanList('A','***.......,..^7.^7.^7.^7.^7.^7,') == list('***AAAAAAAAAAAAAAAA')
    assert getCleanList('A','t,.*,.-2TT..,,......,,,.,,,') == list('TAA*AAAAAAAAAAAAAAAAAAA')
    assert getCleanList('G','.,,,,*.$..$.,,,,,,,.') == list('GGGGG*GGGGGGGGGGGG')
    assert getCleanList('T',',+1a,+1a,+1a.+1A.+1A,+1a.+1A') == list('TTTTTTT')
    assert getCleanList('T','AAaA*A*aA**aaAa') == list('AAAA*A*AA**AAAA')
    assert getCleanList('T','..,....,,,,.,') == list('TTTTTTTTTTTTT')
    assert getCleanList('A',',,....+18AGTTAACCCTAAGGGACC,+18agttaaccctaagggacc,+18agttaaccctaagggacc') == list('AAAAAAAA')
    assert getCleanList('T','***.*,*.***,*') == list('***T*T*T***T*')
    assert getCleanList('T','*************') == list('*************')
