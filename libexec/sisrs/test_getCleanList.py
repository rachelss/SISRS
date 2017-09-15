#!/usr/bin/env python2
from specific_genome import getCleanList

def test_answer():
    assert getCleanList('A','***.......,..^7.^7.^7.^7.^7.^7,') == list('DDDAAAAAAAAAAAAAAAA')
    assert getCleanList('A','t,.*,.-2TT..,,......,,,.,,,') == list('TAADAAAAAAAAAAAAAAAAAAA')
    assert getCleanList('G','.,,,,*.$..$.,,,,,,,.') == list('GGGGGDGGGGGGGGGGGG')
    assert getCleanList('T',',+1a,+1a,+1a.+1A.+1A,+1a.+1A') == list('TTTTTTT')
    assert getCleanList('T','AAaA*A*aA**aaAa') == list('AAAADADAADDAAAA')
    assert getCleanList('T','..,....,,,,.,') == list('TTTTTTTTTTTTT')
    assert getCleanList('A',',,....+18AGTTAACCCTAAGGGACC,+18agttaaccctaagggacc,+18agttaaccctaagggacc') == list('AAAAAAAA')
