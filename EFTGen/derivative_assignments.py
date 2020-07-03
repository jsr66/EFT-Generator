from itertools import permutations
from copy import deepcopy



def tuples_sum(nbval,total,order=True) :
    """
        Generate all the tuples L of nbval positive or nul integer
        such that sum(L)=total.
        The tuples may be ordered (decreasing order) or not
    """
    if nbval == 0 and total == 0 : yield tuple() ; return #raise StopIteration
    if nbval == 1 : yield (total,) ; return #raise StopIteration
    if total==0 : yield (0,)*nbval ; return #raise StopIteration
    for start in range(total,0,-1) :
        for qu in tuples_sum(nbval-1,total-start) :
            if qu[0]<=start :
                sol=(start,)+qu
                if order : yield sol
                else :
                    l=set()
                    for p in permutations(sol,len(sol)) :
                        if p not in l :
                            l.add(p)
                            yield p



def generate_pair_partitions(inner_partition, i_start):
    # INPUT:
    # inner_partition: sorted list (non-increasing order) of numbers that sum to some n_Di, an element of the outer
    # partition of n_D.
    # OUTPUT:
    # pair_partitions_list: a list of lists of tuples, where each tuple has two elements that sum to the corresponding
    # element of the inner_partition. The number of tuples in each sublist matches the number of elements in
    # inner_partition.
    # EXPLANATION:
    # recursive function that returns a list of all distinct ways to partition elements of inner partition into two
    # non-negative integers (indicating number of derivatives acting on psi_bar and psi, respectively).

    # base case
    #print('RECURSION DEPTH: ' + str(len(inspect.stack(0))-29))
    if len(inner_partition) == 1:
        pair_partitions = list(tuple(reversed(x)) for x in tuples_sum(2, inner_partition[0], order=True))[i_start:]
        pair_partitions_list = [[item] for item in pair_partitions]
        return pair_partitions_list

    # if two successive elements of inner_partition are equal, then different orderings of
    # corresponding pair partitions [..., (a,b), (c,d), ...] and [..., (c,d), (a,b), ...] are equal. To avoid
    # repeats, need to modify recursion step for case when successive elements are equal.
    if inner_partition[1] == inner_partition[0]:
        pair_partitions_list = []
        pair_partitions = list(tuple(reversed(x)) for x in tuples_sum(2, inner_partition[0], order=True))
        for i in range(i_start, len(pair_partitions)):
            pair_partitions_list_old = generate_pair_partitions(inner_partition[1:], i)
            #print('RECURSION DEPTH: ' + str(len(inspect.stack(0))-29))
            extension = pair_partitions[i]
            for j in range(len(pair_partitions_list_old)):
                pair_partition_old = pair_partitions_list_old[j]
                pair_partition_new = [extension] + pair_partition_old
                pair_partitions_list.append(pair_partition_new)
        return pair_partitions_list

    # recursion step is simpler in case where successive elements are not equal
    else:
        pair_partitions_list = []
        pair_partitions = list(tuple(reversed(x)) for x in tuples_sum(2, inner_partition[0], order=True))[i_start:]
        pair_partitions_list_old = generate_pair_partitions(inner_partition[1:], 0)
        for i in range(len(pair_partitions)):
            extension = pair_partitions[i]
            for j in range(len(pair_partitions_list_old)):
                pair_partition_old = pair_partitions_list_old[j]
                pair_partition_new = [extension] + pair_partition_old
                pair_partitions_list.append(pair_partition_new)
        return pair_partitions_list



def generate_extended_pair_partitions(inner_partitions_list, first_is_F):
    # INPUT:
    # inner_partitions_list: list of sorted tuples (each in non-increasing order) where each tuple has m_i positive
        # integers that sum to n_i, where the n_i sum to n, the total number of derivatives.
    # first_is_F: indicates whether first sublist of inner_partitions_list corresponds to F-type fields and so
        # does not need to be further partitioned, by contrast with bilinear fields.
    # OUTPUT:
    # extended_pair_partitions_list: list of lists of m_i-tuples of 2-tuples, each sublist an extended pair
    # partition, only represents Dirac bilinear fields, not F fields.
    # EXPLANATION:
    # recursively generates a list of all sequences (lists) of pair partitions, where each pair partition is an m_i-tuple of
    # 2-tuples, one m_i-tuple for each fielnd type. each sequence contains one pair partition for each
    # of the tuples in inner_partitions_list (one such tuple serves as an argument for generate_pair_partitions()).
    # if first inner partition in inner_partition_list corresponds to distribution of derivatives among F-type fields, omit
    # first partition. NOTE: output lists do not include numbers of derivatives acting on F-type fields. Full
    # derivative assignment is output by generate_full_partitions(inner_partitions_list, first_is_F).

    # remove any inner partition for F-type fields to focus on pair partitions of bilinear fields
    if first_is_F == True:
        inner_partitions_list = inner_partitions_list[1:]
        # if remaining list of inner partitions is empty, return empty list
        if len(inner_partitions_list)==0:
            return []

    # base case
    if len(inner_partitions_list) == 1:
        partition = inner_partitions_list[0]
        pair_partitions_list = generate_pair_partitions(partition, 0)
        extended_pair_partitions_list = [[pair_partitions_list[i]] for i in range(len(pair_partitions_list))]
        return extended_pair_partitions_list

    # main body
    extended_pair_partitions_list = []
    first_partition = inner_partitions_list[0]
    first_pair_partitions_list = generate_pair_partitions(first_partition, 0) # list of all pair partitions of first inner partition
    extended_pair_partitions_list_old = generate_extended_pair_partitions(inner_partitions_list[1:], False)
    for i in range(len(first_pair_partitions_list)):
        extension = first_pair_partitions_list[i]
        for j in range(len(extended_pair_partitions_list_old)):
            extended_pair_partition_old = extended_pair_partitions_list_old[j]
            extended_pair_partition_new = [extension] + extended_pair_partition_old
            extended_pair_partitions_list.append(extended_pair_partition_new)

    return extended_pair_partitions_list



def generate_full_partitions(inner_partitions_list, first_is_F):
    # INPUT:
    # inner_partitions_list: list of sorted tuples (each in non-increasing order) where each tuple has m_i positive
        # integers that sum to n_i, where the n_i sum to n, the total number of derivatives.
    # first_is_F: indicates whether first sublist of inner_partitions_list corresponds to F-type fields and so
        # does not need to be further partitioned, by contrast with bilinear fields.
    # OUTPUT:
    # full_partitions_list: list of lists of m_i-tuples. if first_is_F is True, first sublist of full_partitions_list
    # is equal to first sublist of inner_partitions_list. remaining sublists of full_partitions_list are lists of
    # 2-tuples, each sublist a pair partition for a Dirac bilinear.
    if first_is_F:
        F_inner_partition = inner_partitions_list[0]
        full_partitions_list = []
        if len(inner_partitions_list) > 1: # if inner_partitions_list contains more than just F fields
            extended_pair_partitions_list = generate_extended_pair_partitions(inner_partitions_list, first_is_F)
            for pair_partitions_list in extended_pair_partitions_list:
                full_partition = [list(F_inner_partition)] + pair_partitions_list
                full_partitions_list.append(full_partition)
        elif len(inner_partitions_list)==1:
            full_partition = [list(F_inner_partition)]
            full_partitions_list.append(full_partition)
        else:
            return []
        return full_partitions_list
    else:
        # if first sublist inner_partitions_list is not for F-type fields (because these are absent), simply
        # return the pair partitions output by generate_extended_pair_partitions()
        full_partitions_list = generate_extended_pair_partitions(inner_partitions_list, first_is_F)
        return full_partitions_list



def generate_inner_partitions(outer_partition, field_multiplicities):
    # INPUT:
    # outer_partition: list of numbers that sum to number of derivatives n_D with n_fields elements.
    # field_multiplicities: list of multiplicities of fields, same length as outer_partition
    # OUTPUT:
    # subpartitions_list: list of lists of tuples, with each tuple a sorted partition of the corresponding element
    # of outer partition
    # EXPLANATION:
    # recursive function that returns a list of all distinct ways to sub-partition elements of outer partition.
    # if there are k_1 distinct ways to partition n1, k_2 ways to partition n2, ..., k_d ways to partition nd,
    # then there are k_1 x k_2 x ... x k_d distinct ways to sub-partition elements of outer partition.

    if len(outer_partition) != len(field_multiplicities):
        print('Error: Lengths of outer_partition and field_multiplicities arguments must be equal.')
        return None

    # base case
    if len(outer_partition)==1:
        n = outer_partition[0]
        m = field_multiplicities[0]
        return [[subpartition_tuple] for subpartition_tuple in tuples_sum(m, n, order=True)]

    # main body
    subpartitions_list = [] # stores all distinct subpartitions
    n = outer_partition[0]
    m = field_multiplicities[0]
    inner_partitions = list(tuples_sum(m, n, order=True))
    subpartitions = generate_inner_partitions(outer_partition[1:], field_multiplicities[1:])
    for i in range(len(inner_partitions)):
        for j in range(len(subpartitions)):
            subpartition_extended = [inner_partitions[i]] + subpartitions[j]
            subpartitions_list.append(subpartition_extended)

    return subpartitions_list



def generate_derivative_assignments(field_multiset):
    # INPUT:
    # field_multiset: encodes the set of fields and derivatives making up the term
    #     e.g., field_multiset = {D:3, F:2, S:1, T:2};
    #     if a field is absent, simply omit it - e.g., write {D:3, F:2, S:1, T:2}, not {D:3, F:2, S:1, T:2, V:0}
    # OUTPUT:
    # derivative_assignments_list: A dictionary of lists of lists of tuples, where each sublist of each list is a
    # distinct assignment of derivatives to fields.
    # EXPLANATION: generates all distinct ways of assigning derivatives to fields in field_multiset. The fields
    # S, V, T, Vp, Sp all contain two fields, a psi and a psi_bar. A derivative may be placed on either psi or psi_bar.
    # More specifically, generate all 'outer partitions' among the field types - in the above example, F, psi_bar_S, psi_S, psi_bar_T, psi_T.
    # For each outer partition, generate all 'inner_partitions' or 'subpartitions.' Given the
    # outer partition of n, (n1, n2, n3, n4, n5), where each of the five fields fi occurs with multiplicity mi, a
    # subpartition is a list of tuples [(), (), (), (), ()], where the ith tuple is a sorted partition of ni with
    # mi elements. The function generates all such lists of tuples

    #print(field_multiset)
    # extract number of derivatives
    if 'D' in field_multiset.keys():
        n_D = field_multiset['D']
    else:
        n_D = 0

    # find length of outer partition
    #bilinear_keys = [k for k in list(field_multiset.keys()) if k != 'D' and k != 'F'] # keys for bilinear elements
    #bilinear_multiplicities = [field_multiset[m] for m in list(field_multiset.keys()) if m != 'D' and m != 'F']
    #n_fields = sum([1 if 'F' in field_multiset.keys() else 0]) \ # use list comprehension on list of one element
    #    + sum([2 for key in bilinears]) # group psi and psi_bar together as one field for now

    #if 'F' in field_multiset.keys():
    #    F_multiplicity = [field_multiset['F']]
    #else:
    #    F_multiplicity = []
    #bilinear_multiplicities = [x for pair in zip(bilinear_multiplicities,bilinear_multiplicities) for x in pair]
    # for each bilinear, multiplicities are the same for psi and psi_bar
    field_multiplicities = [field_multiset[x] for x in field_multiset.keys() if x != 'D']

    n_fields = len(field_multiplicities)
    #print('n_fields: ' + str(n_fields))
    #print('bilinear_multiplicities: ' + str(bilinear_multiplicities))
    #print('field_multiplicities: ' + str(field_multiplicities))


    #list of lists, each sublist an outer partition, where different orderings are distinct
    outer_partitions_list = list(tuples_sum(n_fields, n_D, order=False))

    # for each outer partition, generate a list of sorted inner partitions (ordering does not matter, so one may
    # as well sort them) where the sum of each inner partition is the corresponding element of the outer partition.
    # this is a list of lists of lists.
    #inner_partitions_dict = {}
    derivative_assignments_dict = {}
    for i in range(len(outer_partitions_list)):
        outer_partition = outer_partitions_list[i]
        #print('outer_partition: ' + str(outer_partition))
        inner_partitions_list = generate_inner_partitions(outer_partition, field_multiplicities)
        #inner_partitions_dict[outer_partition] = inner_partitions_list
        # for each outer partition and inner partition, generate a list of pair partitions
        inner_partitions_dict = {}
        for j in range(len(inner_partitions_list)):
            inner_partition = tuple(inner_partitions_list[j])
            #print('inner_partition: ' + str(inner_partition))
            full_partitions_list = generate_full_partitions(inner_partition, 'F' in field_multiset.keys())
            inner_partitions_dict[inner_partition] = full_partitions_list
        derivative_assignments_dict[outer_partition] = inner_partitions_dict

    return derivative_assignments_dict



def IBP_remove(derivative_assignment):
    # INPUT:
    # derivative_assignment: list of lists, each sublist corresponding to a different field type. if F-type
    # fields are represented, first sublist should be a list of integers. any other lists, representing Dirac
    # bilinears, should be lists of 2-tuples.
    # OUTPUT:
    # remove_bool: True if derivative_assignment meets requirements for removal via IBP, False otherwise
    # EXPLANATION:
    # checks whether a derivative assignment is to be removed under the specified prescription - i.e.,
    # whether the maximum number of derivatives acting on any field is unique, and whether any field that comes
    # earlier in the ordering of field types has exactly one fewer derivative acting on it.

    # flatten derivative assignment list to facilitate checks below
    first_is_F = type(derivative_assignment[0][0])==int
    if first_is_F:
        flattened_list_F = [x for x in derivative_assignment[0]]
        if len(derivative_assignment) > 1: # if both Dirac bilinears and F fields are included
            flattened_list_bilinears = [x for inner_partition in derivative_assignment[1:] for pair in inner_partition for x in pair]
            flattened_list = flattened_list_F + flattened_list_bilinears
        else: # if only F-type fields are included
            flattened_list = flattened_list_F
    else: # if only Dirac bilinear fields are included
        flattened_list = [x for inner_partition in derivative_assignment for pair in inner_partition for x in pair]


    # check that maximal number of derivatives only acts on one field
    max_derivs = max(flattened_list)
    unique_max = sum([x == max_derivs for x in flattened_list]) == 1
    #print('unique_max: ' + str(unique_max ))

    # check that no field before this field has one fewer derivative acting on it
    max_index = flattened_list.index(max_derivs)
    #print('max_index: ' + str(max_index))
    if max_index > 0:
        one_less_before = max(flattened_list[:max_index]) == max_derivs - 1
    else:
        one_less_before = False
    #print('one_less_before: ' + str(one_less_before))

    if unique_max and not one_less_before:
        return True

    return False



def CPT_transformation_factor(field_multiset, derivative_assignment, symmetry):
    # under C, if D acts on a Dirac bilinear, then it gets a minus sign from IBP in switching it from psibar to psi.
    # under T, if D acts on a Dirac bilinear, then it comes with an i, which transforms with an additional minus sign because of antilinearity.
    # for each derivative in derivative assignment:
    # - if symmetry is C, if that derivative acts on a Dirac bilinear,
    #   include an additional factor of minus one for each derivative to account for IBP switch.
    # - if symmetry is T, if that derivative acts on a Dirac bilinear, include an additional factor of -1 to
    #   account for the additional factor of i attached to derivative and antilinearity.
    # - if symmetry is P, simply multiply relevant factors in transformation_factor_dict
    if symmetry not in ['C', 'P', 'T']:
        print('"symmetry" argument must be "C", "P", or "T".')
        return None

    transformation_factor_dict = {
    'C': {'D': 1, 'F': -1, 'S': 1, 'V': -1, 'T': -1, 'Vp': 1, 'Sp': 1}, #change D to -1 to account for ibp manipulation
    'P': {'D': 1, 'F': 1, 'S': 1, 'V': 1, 'T': 1, 'Vp': -1, 'Sp': -1},
    'T': {'D': -1, 'F': -1, 'S': 1, 'V': 1, 'T': -1, 'Vp': 1, 'Sp': -1} #note D's acting on psi come with i while D's acting on F don't. this interacts with antilinearity of T
    }

    if symmetry == 'C' or symmetry == 'T':
        nD = field_multiset['D'] if 'D' in field_multiset.keys() else 0
        if 'F' in field_multiset.keys(): #subtract number of derivatives acting on F since these aren't affected by IBP or anti-linearity
            derivative_assignment_F = derivative_assignment[0]
            nD_on_F = sum(derivative_assignment_F)
            nD_on_bilinear = nD - nD_on_F
            if nD_on_bilinear < 0:
                print('In CPT_transformation_factor: nD_on_bilinear cannot be negative.')
            factor = (-1)**nD_on_bilinear
        else: # all derivatives acting on Dirac bilinears
            factor = (-1)**nD

        for field_label in field_multiset.keys():
            num_fields = field_multiset[field_label]
            field_factor = transformation_factor_dict[symmetry][field_label]
            factor *= field_factor**num_fields

    if symmetry == 'P':
        factor = 1
        for field_label in field_multiset.keys():
            num_fields = field_multiset[field_label]
            field_factor = transformation_factor_dict[symmetry][field_label]
            factor *= field_factor**num_fields

    return factor



def generate_derivative_assignments_IBP_CPT_filtered(field_multiset, IBP=True, C=True, P=True, T=True):
    # INPUT:
    # field_multiset: dictionary containing multiplicities o f derivative and field types
    # OUTPUT:
    # reduced_derivative_assignments_dict: a dictionary containing an IBP-reduced set of derivative assignments
    # deleted: a list containing all removed derivative assignments
    # EXPLANATION:
    # scans through all derivative assignments and performs IBP reduction. it removes one operator for each
    # IBP relation, which appears in one and only one IBP relation. in practice, it does this by removing
    # any derivative assignments for which the maximum number of derivatives acting on any field is unique, and
    # for which there is not one fewer derivative acting on any field that occurs earlier in the specified ordering
    # of field types. it is worth noting that there is an element of conventionality to this prescription, and that
    # other IBP reductions are possible that remove different sets of derivative assignments.

    # generate derivative assignments and make copy of dictionary to use for deletion
    derivative_assignments_dict = generate_derivative_assignments(field_multiset)
    derivative_assignments_dict_reduced = deepcopy(derivative_assignments_dict)
    #print(derivative_assignments_dict_reduced)
    # declare empty list to store removed derivative assignments
    IBP_deleted = []
    C_deleted = []
    P_deleted = []
    T_deleted = []
    #first_is_F = 'F' in field_multiset.keys()

    # loop through derivative assignments and remove from copy of dictionary
    # those derivative assignments that meet the criteria for removal, appending them to
    # the list of deleted items
    for outer_partition in derivative_assignments_dict.keys():
        #print('outer partition: ' + str(outer_partition))
        inner_partitions_dict = derivative_assignments_dict[outer_partition]
        for inner_partition in inner_partitions_dict.keys():
            #print('inner partition: ' + str(inner_partition))
            full_partitions_list = inner_partitions_dict[inner_partition]
            full_partitions_list_reduced = derivative_assignments_dict_reduced[outer_partition][inner_partition]
            for derivative_assignment in full_partitions_list:
                #print(derivative_assignment)
                IBP_delete = IBP_remove(derivative_assignment)
                C_delete = C and CPT_transformation_factor(field_multiset, derivative_assignment, 'C')!=1
                P_delete = P and CPT_transformation_factor(field_multiset, derivative_assignment, 'P')!=1
                T_delete = T and CPT_transformation_factor(field_multiset, derivative_assignment, 'T')!=1
                delete = False
                if IBP_delete:
                    #print('IBP_delete')
                    IBP_deleted.append(derivative_assignment)
                    delete = True
                if C_delete:
                    #print('C_delete')
                    C_deleted.append(derivative_assignment)
                    delete = True
                if P_delete:
                    #print('P_delete')
                    P_deleted.append(derivative_assignment)
                    delete= True
                if T_delete:
                    #print('T_delete')
                    T_deleted.append(derivative_assignment)
                    delete = True
                    #print(len(derivative_assignments_dict[outer_partition][inner_partition]))
                if delete:
                    full_partitions_list_reduced.remove(derivative_assignment)
                    #print(len(derivative_assignments_dict[outer_partition][inner_partition]))

    return derivative_assignments_dict_reduced, IBP_deleted, C_deleted, P_deleted, T_deleted
