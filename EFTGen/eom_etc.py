from index_contractions import encode_lorentz_grouping

def EOM_remove(field_multiset, derivative_assignment, contraction_multiset):
    # loop through contractions in contraction_multiset; for each, check if it is a
    # contraction between a D and a V, T, or V_p. Do this by getting the symbol grouping
    # for field_multiset and derivative_assignment using encode_lorentz_groupings().

    # do not remove \gamma^{\mu} \bar{\psi} D_{\mu} \psi term that generates Dirac field EoMs
    if field_multiset == {'D':1, 'V':1} and derivative_assignment == [[(0,1)]] and contraction_multiset == {(0,1):1}:
        return False
    symbol_grouping, _ = encode_lorentz_grouping(field_multiset, derivative_assignment)
    symbol_grouping_flattened = [x for group in symbol_grouping for x in group]
    #print(symbol_grouping_flattened)
    for contraction in contraction_multiset.keys():
        index_0 = contraction[0]
        index_1 = contraction[1]
        #print(index_0)
        #print(index_1)

        symbol_0 = symbol_grouping_flattened[index_0]
        symbol_1 = symbol_grouping_flattened[index_1]

        # if both symbols are in the same group and one is a D and another a V, T, or V_p,
        # return True.

        # first, check if both symbols and indices are in the same group of basic operators - i.e.,
        # belong to the same Dirac bilinear Lorentz tensor. determine whether two basic operators
        # belong to the same group by setting boundary indices for each group and seeing
        # whether index_0 and index_1 lie between the same pair of boundary indices.
        #index_group_boundaries = [sum([len(symbol_grouping[i]) for j in range(i+1)]) - 1 for i in range(len(symbol_grouping))]


        index_group_boundaries = [sum([len(symbol_grouping[i]) for i in range(j+1)])-1 for j in range(len(symbol_grouping))]
        #print("index_group_boundaries: " + str(index_group_boundaries))
        if index_0 <= index_group_boundaries[0]:
            index_0_group = 0
        elif index_0 > index_group_boundaries[-1]:
            print('ERROR in EOM_remove: index_0 value greater than length of list of basic operators.')
            return

        if index_1 <= index_group_boundaries[0]:
            index_1_group = 0
        elif index_1 > index_group_boundaries[-1]:
            print('ERROR in EOM_remove: index_1 value greater than length of list of basic operators.')
            return

        else:
            for i in range(len(index_group_boundaries)-1):
                if index_0 > index_group_boundaries[i] and index_0 <= index_group_boundaries[i+1]:
                    index_0_group = i+1
                if index_1 > index_group_boundaries[i] and index_1 <= index_group_boundaries[i+1]:
                    index_1_group = i+1

        same_group = index_0_group==index_1_group

        # next
        one_is_D = 'D' in symbol_grouping_flattened[index_0] or 'D' in symbol_grouping_flattened[index_1]
        one_is_V = symbol_grouping_flattened[index_0]=='V' or symbol_grouping_flattened[index_1]=='V'
        one_is_T = symbol_grouping_flattened[index_0]=='T' or symbol_grouping_flattened[index_1]=='T'
        one_is_Vp = symbol_grouping_flattened[index_0]=='Vp' or symbol_grouping_flattened[index_1]=='Vp'
        one_is_F = symbol_grouping_flattened[index_0]=='F' or symbol_grouping_flattened[index_1]=='F'

        if same_group and one_is_D and (one_is_V or one_is_T or one_is_Vp):
            return True
        if same_group and one_is_D and one_is_F:
            return True
    # if no individual contraction meets the conditions for removal, return False.
    return False



def DD_comm_remove(field_multiset, derivative_assignment, contraction_multiset):
    # determine whether the same contraction between a group of D's and any F or T
    # (antisymmetric)- not necessarily an F or T that the D's act on - occurs more than
    # once. this yields a commutator of D's which yields an F. if so, remove the operator
    # - i.e., return True. otherwise, return False.

    # determine which contractions occur more than once.
    contractions_list = [item[0] for item in contraction_multiset.items() if item[1] > 1]

    # determine symbols associated with each index
    symbol_grouping, _ = encode_lorentz_grouping(field_multiset, derivative_assignment)
    symbol_grouping_flattened = [x for group in symbol_grouping for x in group]

    # for each contraction in relevant contractions, determine whether it is between a
    # group of D's and an F or T. if so, return True.
    for contraction in contractions_list:
        index_0 = contraction[0]
        index_1 = contraction[1]

        one_is_D = 'D' in symbol_grouping_flattened[index_0] or 'D' in symbol_grouping_flattened[index_1]
        one_is_T = symbol_grouping_flattened[index_0]=='T' or symbol_grouping_flattened[index_1]=='T'
        one_is_F = symbol_grouping_flattened[index_0]=='F' or symbol_grouping_flattened[index_1]=='F'

        if one_is_D and (one_is_F or one_is_T):
            return True

    # if no contraction that occurs more than once is between a group of D's and an F or T,
    # return False.
    return False



def FT_selfcontracted_remove(field_multiset, derivative_assignment, contraction_multiset):
    # determine whether an F or a T is contracted with itself. if so, remove operator
    # i.e., return True. otherwise, return False.

    # determine symbols associated with each index
    symbol_grouping, _ = encode_lorentz_grouping(field_multiset, derivative_assignment)
    symbol_grouping_flattened = [x for group in symbol_grouping for x in group]

    # find indices i associated with F or T.
    F_or_T_indices = [i for i in range(len(symbol_grouping_flattened)) if (symbol_grouping_flattened[i]=='F' or symbol_grouping_flattened[i]=='T')]

    # determine if there is a contraction (i,i) for any of these indices. if so, return True
    for i in F_or_T_indices:
        if (i,i) in contraction_multiset.keys():
            return True
    # if no contraction is found of an F with itself or a T with itself, return False.
    return False
