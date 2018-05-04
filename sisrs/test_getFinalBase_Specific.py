#!/usr/bin/env python3
from .specific_genome import getCleanList
from .specific_genome import getFinalBase_Specific

def test_answer():

    #Test getFinalBase for specific_genome
    assert getFinalBase_Specific(getCleanList('A','***.......,..^7.^7.^7.^7.^7.^7,'))=='A'
    assert getFinalBase_Specific(getCleanList('A','t,.*,.-2TT..,,......,,,.,,,'))=='A'
    assert getFinalBase_Specific(getCleanList('G','.,,,,*.$..$.,,,,,,,.'))=='G'
    assert getFinalBase_Specific(getCleanList('T',',+1a,+1a,+1a.+1A.+1A,+1a.+1A'))=='T'
    assert getFinalBase_Specific(getCleanList('T','AAaA*A*aA**aaAa'))=='A'
    assert getFinalBase_Specific(getCleanList('T','..,....,,,,.,'))=='T'
    assert getFinalBase_Specific(getCleanList('A',',,....+18AGTTAACCCTAAGGGACC,+18agttaaccctaagggacc,+18agttaaccctaagggacc'))=='A'
    assert getFinalBase_Specific(getCleanList('T','***.*,*.***,*'))=='N'

    #Test how to handle ties
    assert getFinalBase_Specific(getCleanList('T','TtTtTAAaAa'))=='A'
    assert getFinalBase_Specific(getCleanList('T','GgGggAAaAa'))=='A'
    assert getFinalBase_Specific(getCleanList('T','TtTtTgGggG'))=='G'

    #Test how to handle deletions
    assert getFinalBase_Specific(getCleanList('T','*************'))=='N'
    assert getFinalBase_Specific(getCleanList('T','*****AAaAa'))=='A'
