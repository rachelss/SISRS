#!/usr/bin/env python2
from specific_genome import getCleanList
from specific_genome import getFinalBase_Specific
from get_pruned_dict import getFinalBase_Pruned

def test_answer():
    assert getCleanList('A','***.......,..^7.^7.^7.^7.^7.^7,') == list('***AAAAAAAAAAAAAAAA')
    assert getCleanList('A','t,.*,.-2TT..,,......,,,.,,,') == list('TAA*AAAAAAAAAAAAAAAAAAA')
    assert getCleanList('G','.,,,,*.$..$.,,,,,,,.') == list('GGGGG*GGGGGGGGGGGG')
    assert getCleanList('T',',+1a,+1a,+1a.+1A.+1A,+1a.+1A') == list('TTTTTTT')
    assert getCleanList('T','AAaA*A*aA**aaAa') == list('AAAA*A*AA**AAAA')
    assert getCleanList('T','..,....,,,,.,') == list('TTTTTTTTTTTTT')
    assert getCleanList('A',',,....+18AGTTAACCCTAAGGGACC,+18agttaaccctaagggacc,+18agttaaccctaagggacc') == list('AAAAAAAA')
    assert getCleanList('T','***.*,*.***,*') == list('***T*T*T***T*')
    assert getCleanList('T','*************') == list('*************')

    assert getFinalBase_Specific(getCleanList('A','***.......,..^7.^7.^7.^7.^7.^7,'))=='A'
    assert getFinalBase_Pruned(getCleanList('A','***.......,..^7.^7.^7.^7.^7.^7,'),3,1)=='N'

    assert getFinalBase_Specific(getCleanList('A','t,.*,.-2TT..,,......,,,.,,,'))=='A'
    assert getFinalBase_Pruned(getCleanList('A','t,.*,.-2TT..,,......,,,.,,,'),3,1)=='N'

    assert getFinalBase_Specific(getCleanList('G','.,,,,*.$..$.,,,,,,,.'))=='G'
    assert getFinalBase_Pruned(getCleanList('G','.,,,,*.$..$.,,,,,,,.'),3,1)=='N'
    assert getFinalBase_Pruned(getCleanList('G','.,,,,*.$..$.,,,,,,,.'),3,0.93)=='G'

    assert getFinalBase_Specific(getCleanList('T',',+1a,+1a,+1a.+1A.+1A,+1a.+1A'))=='T'
    assert getFinalBase_Pruned(getCleanList('T',',+1a,+1a,+1a.+1A.+1A,+1a.+1A'),3,1)=='T'
    assert getFinalBase_Pruned(getCleanList('T',',+1a,+1a,+1a.+1A.+1A,+1a.+1A'),8,1)=='N'

    assert getFinalBase_Specific(getCleanList('T','AAaA*A*aA**aaAa'))=='A'
    assert getFinalBase_Pruned(getCleanList('T','AAaA*A*aA**aaAa'),3,1)=='N'

    assert getFinalBase_Specific(getCleanList('T','..,....,,,,.,'))=='T'
    assert getFinalBase_Pruned(getCleanList('T','..,....,,,,.,'),3,1)=='T'

    assert getFinalBase_Specific(getCleanList('A',',,....+18AGTTAACCCTAAGGGACC,+18agttaaccctaagggacc,+18agttaaccctaagggacc'))=='A'
    assert getFinalBase_Pruned(getCleanList('A',',,....+18AGTTAACCCTAAGGGACC,+18agttaaccctaagggacc,+18agttaaccctaagggacc'),3,1)=='A'

    assert getFinalBase_Specific(getCleanList('T','***.*,*.***,*'))=='N'
    assert getFinalBase_Pruned(getCleanList('T','***.*,*.***,*'),3,1)=='N'

    assert getFinalBase_Specific(getCleanList('T','*************'))=='N'
    assert getFinalBase_Pruned(getCleanList('T','*************'),3,1)=='-'
