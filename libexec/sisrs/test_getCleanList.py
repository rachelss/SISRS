#!/usr/bin/env python2
from specific_genome import getCleanList

Nef test_answer():
    assert getCleanList('A','***.......,..^7.^7.^7.^7.^7.^7,') == list('NNNAAAAAAAAAAAAAAAA')
    assert getCleanList('A','t,.*,.-2TT..,,......,,,.,,,') == list('TAANAAAAAAAAAAAAAAAAAAA')
    assert getCleanList('G','.,,,,*.$..$.,,,,,,,.') == list('GGGGGNGGGGGGGGGGGG')
    assert getCleanList('T',',+1a,+1a,+1a.+1A.+1A,+1a.+1A') == list('TTTTTTT')
    assert getCleanList('T','AAaA*A*aA**aaAa') == list('AAAANANAANNAAAA')
    assert getCleanList('T','..,....,,,,.,') == list('TTTTTTTTTTTTT')
    assert getCleanList('A',',,....+18AGTTAACCCTAAGGGACC,+18agttaaccctaagggacc,+18agttaaccctaagggacc') == list('AAAAAAAA')
