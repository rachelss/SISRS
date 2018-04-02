#!/usr/bin/env python2
from specific_genome import getCleanList
from get_pruned_dict_Memory import getFinalBase_Pruned

def test_answer():

    assert getFinalBase_Pruned(getCleanList('A','***.......,..^7.^7.^7.^7.^7.^7,'),3,1,0,0,0)[0]=='N'
    assert getFinalBase_Pruned(getCleanList('A','t,.*,.-2TT..,,......,,,.,,,'),3,1,0,0,0)[0]=='N'
    assert getFinalBase_Pruned(getCleanList('G','.,,,,*.$..$.,,,,,,,.'),3,1,0,0,0)[0]=='N'
    assert getFinalBase_Pruned(getCleanList('G','.,,,,*.$..$.,,,,,,,.'),3,0.93,0,0,0)[0]=='G' #Test threshold
    assert getFinalBase_Pruned(getCleanList('T',',+1a,+1a,+1a.+1A.+1A,+1a.+1A'),3,1,0,0,0)[0]=='T'
    assert getFinalBase_Pruned(getCleanList('T',',+1a,+1a,+1a.+1A.+1A,+1a.+1A'),8,1,0,0,0)[0]=='N'   #Test minread
    assert getFinalBase_Pruned(getCleanList('T','AAaA*A*aA**aaAa'),3,1,0,0,0)[0]=='N'
    assert getFinalBase_Pruned(getCleanList('T','..,....,,,,.,'),3,1,0,0,0)[0]=='T'
    assert getFinalBase_Pruned(getCleanList('A',',,....+18AGTTAACCCTAAGGGACC,+18agttaaccctaagggacc,+18agttaaccctaagggacc'),3,1,0,0,0)[0]=='A'
    assert getFinalBase_Pruned(getCleanList('T','***.*,*.***,*'),3,1,0,0,0)[0]=='N'

    #Test how to handle deletions
    assert getFinalBase_Pruned(getCleanList('T','*************'),3,1,0,0,0)[0]=='-'
