import itertools

def encode_lorentz_grouping(field_multiset, derivative_assignment):
    # check compatibility of field_multiset and derivative_assignment.
    # first, check that number of non-derivative elements in field_multiset is equal to the number of sub-lists in
    # derivative_assignment.
    #print('IN ENCODE_LORENTZ_GROUPINGS()')
    #print('field_multiset: ' + str(field_multiset))
    nonderivative_fields = list(field_multiset.keys())
    if 'D' in nonderivative_fields:
        nonderivative_fields.remove('D')
    if len(nonderivative_fields) != len(derivative_assignment):
        print('In encode_lorentz_grouping(). Arguments are incompatible.')

    lorentz_grouping = []
    num_field_types = len(derivative_assignment)
    for i in range(num_field_types):
        field_type = nonderivative_fields[i]
        derivative_dist = derivative_assignment[i]
        num_copies = len(derivative_dist)
        for j in range(num_copies):
            # for F-type fields
            if field_type == 'F':
                nD = derivative_dist[j]
                derivative_group = [nD*'D' for i in range(1) if nD*'D'] #a bit hacky. leaves list empty if nD is zero
                group = derivative_group + [field_type]
            # for bilinear fields
            else:
                two_tuple = derivative_dist[j]
                nD_left = two_tuple[0]
                nD_right = two_tuple[1]
                derivative_group_left = [nD_left*'D' for i in range(1) if nD_left*'D'] #a bit hacky. leaves list empty if nD_left is zero
                #print(derivative_group_left)
                derivative_group_right = [nD_right*'D' for i in range(1) if nD_right*'D']
                #print(derivative_group_right)
                #print([field_type])
                group = derivative_group_left + [field_type] + derivative_group_right

            lorentz_grouping.append(group)

    # group indices in nested list by group type, group copy
    index_grouping_by_type = []
    index = 0
    num_groups = len(lorentz_grouping)
    group_type = []
    for i in range(num_groups):
        #print('i: ' + str(i))
        group = lorentz_grouping[i]
        #print(group)
        num_items = len(group)

        group_indices = []
        for j in range(num_items):
            group_indices.append(index)
            #print(index)
            index += 1
        #print(group_indices)

        if i == 0:
            #print('i==0')
            group_type = [group_indices]
            #print('group_type: ' + str(group_type))
            if i == num_groups - 1:
                    #print('i == num_groups - 1')
                    index_grouping_by_type.append(group_type)
        else:
            if group == lorentz_grouping[i-1]:
                #print('group == lorentz_grouping[i-1]')
                group_type.append(group_indices)
                #print('group_type: ' + str(group_type))
                if i == num_groups - 1:
                    #print('i == num_groups - 1')
                    index_grouping_by_type.append(group_type)
            else:
                index_grouping_by_type.append(group_type)
                group_type = [group_indices]
                if i == num_groups - 1:
                    #print('i == num_groups - 1')
                    index_grouping_by_type.append(group_type)
        #print(group_type)


    return lorentz_grouping, index_grouping_by_type



def list_lorentz_ranks(lorentz_grouping_flattened):
    ranks_list = []
    lorRanks = {'D': 1, 'F': 2, 'S': 0, 'V': 1, 'T': 2, 'Vp': 1, 'Sp': 0}
    n_items = len(lorentz_grouping_flattened)
    for i in range(n_items):
        group = lorentz_grouping_flattened[i]
        group_lorRank = 0
        if group[0]=='D':
            for j in range(len(group)):
                group_lorRank += lorRanks[group[j]]
        else:
            group_lorRank += lorRanks[group]
        ranks_list.append(group_lorRank)
    return ranks_list



def generate_lorentz_contractions_from_ranks(lorentz_ranks_list, i0):
    # include lorentz_grouping_flattened as argument to exclude self contractions of F and T fields
    # base case: if the total number of free Lorentz indices is 0 or 1, return empty list of contraction sets
    #print(lorentz_grouping_flattened)
    #print(lorentz_ranks_list)
    if sum(lorentz_ranks_list) < 2:
        contraction_multisets_list = [{}]
        return contraction_multisets_list

    # find first field of non-zero rank for first index of contraction
    N_groups = len(lorentz_ranks_list)
    for i in range(N_groups):
        if lorentz_ranks_list[i] >= 1:
            i_start = i
            lorentz_ranks_list[i_start] -= 1
            #print('i_start: ' + str(i_start))
            break

    # if lower limit for second contraction index is less than lower limit for first index, set
    # former equal to latter.
    if i0 <= i_start:
        i0 = i_start

    #print('i0: ' + str(i0))

    contraction_multisets_list = []
    # find all ways of contracting first non-zero lorentz rank field with other fields or itself.
    for i in range(i0, N_groups):
        #print(i)

        if lorentz_ranks_list[i] >= 1:
            contraction = (i_start, i)

            lorentz_ranks_list_old = lorentz_ranks_list.copy()

            # decrement lorentzRanks_list_old
            lorentz_ranks_list_old[i] -= 1

            # if i_start of decremented list is the same as the non-decremented list, only include contractions
            # where second index is greater than or equal to i. otherwise, i_start will be larger and we can
            # include all values of the second index.
            if lorentz_ranks_list_old[i_start] > 0:
                contraction_multisets_list_old = generate_lorentz_contractions_from_ranks(lorentz_ranks_list_old, i)
            else:
                contraction_multisets_list_old = generate_lorentz_contractions_from_ranks(lorentz_ranks_list_old, 0)

            for contraction_multiset_old in contraction_multisets_list_old:
                contraction_multiset = contraction_multiset_old.copy()
                if contraction in contraction_multiset.keys():
                    contraction_multiset[contraction] += 1
                else:
                    contraction_multiset[contraction] = 1
                contraction_multisets_list.append(contraction_multiset)

    return contraction_multisets_list



def generate_lorentz_contractions_with_repeats(field_multiset, derivative_assignment):
    #print('IN GENERATE_LORENTZ_CONTRACTIONS()')
    lorentz_grouping,_ = encode_lorentz_grouping(field_multiset, derivative_assignment)
    lorentz_grouping_flattened = [x for group in lorentz_grouping for x in group]
    #print('lorentz_grouping_flattened: ' + str(lorentz_grouping_flattened))
    lorentz_ranks_list = list_lorentz_ranks(lorentz_grouping_flattened)
    #print('lorentz_ranks_list: ' + str(lorentz_ranks_list))
    contraction_multisets_list = generate_lorentz_contractions_from_ranks(lorentz_ranks_list, 0)
    return contraction_multisets_list



def flatten(index_grouping_by_type):
    return [x for group_type in index_grouping_by_type for copy in group_type for x in copy]

def replace_indices(contraction_tuple, index_map):
    #print('IN REPLACE_INDICES()')
    #print('index_map: ' + str(index_map))
    #print('contraction_tuple: ' + str(contraction_tuple))
    contraction_tuple_new = tuple(sorted((index_map[contraction_tuple[0]], index_map[contraction_tuple[1]])))
    return contraction_tuple_new

def substitute_permuted_indices(lorentz_contraction, index_map):
    #print('IN SUBSTITUTE_PERMUTED_INDICES()')
    lorentz_contraction_new = {}
    for contraction_tuple in lorentz_contraction.keys():
        #print('contraction_tuple: ' + str(contraction_tuple))
        contraction_tuple_new = replace_indices(contraction_tuple, index_map)
        #print('contraction_tuple_new: ' + str(contraction_tuple_new))
        lorentz_contraction_new[contraction_tuple_new] = lorentz_contraction[contraction_tuple]
    return lorentz_contraction_new

def permuted_contractions(lorentz_contraction, index_grouping_by_type):
    #print('IN PERMUTED_CONTRACTIONS()')
    # generates all permutations of lorentz_contraction equivalent under index_grouping_by_type. for each group
    # type with more than one copy - i.e., for each sublist with more than one subsublist -

    '''
    # for each group type, generate a list of all permutations of indices associated with copies of that group type.
    # assemble a list of these lists, called permutations_by_type
    permutations_by_type = []
    for group_type in index_grouping_by_type:
        permuted_contractions = []
        permutations_list = list(itertools.permutations(group_type))
        permutations_by_type.append(permutations_list)
    '''

    # generate all possible sequences of permutations, where each sequence consists of one permutation for each group
    # type.
    #print('index_grouping_by_type: ' + str(index_grouping_by_type))
    N_group_types = len(index_grouping_by_type)
    #print('N_group_types')
    full_permutations_list = [[]]
    counter = 0
    while counter < N_group_types:
        #print('counter: ' + str(counter))
        full_permutations_list_old = full_permutations_list.copy()
        full_permutations_list = []
        group_type = index_grouping_by_type[counter]
        permutations_list = list(itertools.permutations(group_type))
        #print(permutations_list)
        for full_permutation in full_permutations_list_old:
            #print('full_permutation: ' + str(full_permutation))
            for permutation in permutations_list:
                #print('permutation: ' + str(permutation))
                full_permutation_new = full_permutation + [list(permutation)]
                full_permutations_list.append(full_permutation_new)
        #print(full_permutations_list)
        counter += 1
    #print('full_permutations_list: ' + str(full_permutations_list))

    # for each possible sequence of permutations, generate an index map from original indices to permuted ones
    permuted_contractions_list = []
    index_grouping_by_type_flattened = flatten(index_grouping_by_type)
    for full_permutation in full_permutations_list:
        #print('full_permutation: ' + str(full_permutation))
        index_grouping_by_type_permuted_flattened = flatten(full_permutation)
        #print('index_grouping_by_type_permuted_flattened: ' + str(index_grouping_by_type_permuted_flattened))
        index_map = {i: index_grouping_by_type_permuted_flattened[i] for i in index_grouping_by_type_flattened}
        #print('index_map: ' + str(index_map))
        lorentz_contraction_permuted = substitute_permuted_indices(lorentz_contraction, index_map)
        permuted_contractions_list.append(lorentz_contraction_permuted)

    return permuted_contractions_list



def lorentz_contractions_equiv(lorentz_contraction1, lorentz_contraction2, field_multiset, derivative_assignment):
    #print('IN LORENTZ_CONTRACTIONS_EQUIV()')
    #print('field_multiset: ' + str(field_multiset))
    #print('derivative_assignment: ' + str(derivative_assignment))
    _, index_grouping_by_type = encode_lorentz_grouping(field_multiset, derivative_assignment)
    #print('index_grouping_by_type: ' + str(index_grouping_by_type))
    if lorentz_contraction2 in permuted_contractions(lorentz_contraction1, index_grouping_by_type):
        return True
    else:
        return False

def generate_lorentz_contractions(field_multiset, derivative_assignment):
    # THIS IS THE MAIN FUNCTION OF THIS SECTION
    #print('IN GENERATE_LORENTZ_CONTRACTIONS()')
    contraction_multisets_list = generate_lorentz_contractions_with_repeats(field_multiset, derivative_assignment)
    N_contractions = len(contraction_multisets_list)
    #print('N_contractions: ' + str(N_contractions))
    # start list of unique elements of contraction_multisets_list with first element
    contraction_multisets_list_unique = [contraction_multisets_list[0]]

    for i in range(1,N_contractions):
        #print('i: ' + str(i))
        contraction_multiset = contraction_multisets_list[i]
        N_unique = len(contraction_multisets_list_unique)
        #print('N_unique: ' + str(N_unique))

        def already_in_list():
            for j in range(N_unique):
                contraction_multiset_unique = contraction_multisets_list_unique[j]
                equiv = lorentz_contractions_equiv(contraction_multiset, contraction_multiset_unique, field_multiset, derivative_assignment)
                if equiv:
                    return True
            return False

        if already_in_list():
            continue
        else:
            contraction_multisets_list_unique.append(contraction_multiset)

    return contraction_multisets_list_unique
            
