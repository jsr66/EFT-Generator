from frozendict import frozendict
from field_multisets import generate_field_multisets
from derivative_assignments import generate_derivative_assignments_IBP_CPT_filtered
from index_contractions import generate_lorentz_contractions, encode_lorentz_grouping
from eom_etc import EOM_remove, DD_comm_remove, FT_selfcontracted_remove



def generate_operators(massDim, C=True, P=True, T=True, IBP=True, EOM=True, DD_comm=True, FT_selfcontracted=True, verbose = False):
    # EXPLANATION: generate dictionary of operators in reduced basis, and dictionaries of operators removed according to various
    # criteria.
    if verbose:
        print('IN GENERATE_OPERATORS()')
    if massDim < 2 or type(massDim) != int:
        print('massDim argument must be an int greater than or equal to 2.')
        return

    field_symbols = ['D', 'F', 'S', 'V', 'T', 'Vp', 'Sp']
    field_massDims = [1, 2, 3, 3, 3 ,3 ,3]
    field_lorentzRanks = [1, 2, 0, 1, 2, 1, 0]

    field_multisets_dict = generate_field_multisets(massDim, field_symbols, field_massDims, field_lorentzRanks)
    if verbose:
        print('field/derivative multisets generated.')
    field_multiset_list = field_multisets_dict[massDim]

    # declare empty dictionary to store final list of generated operators
    operator_dict = {}

    # declare dictionaries to store deleted derivative assignments
    IBP_deleted_dict = {}
    C_deleted_dict = {}
    P_deleted_dict = {}
    T_deleted_dict = {}

    # declare dictionaries to store deleted contractions
    EOM_deleted_dict = {}
    DDcomm_deleted_dict = {}
    FTselfcontracted_deleted_dict = {}

    for field_multiset in field_multiset_list:
        # for each field multiset, generate all derivative assignments
        if verbose:
            print('')
            print('')
            print('')
            print('field multiset: ' + str(field_multiset))
        # declare empty dictionary for each multiset, with derivative assignments as keys and list of full contractions as the value for each key
        operator_dict_field_multiset = {}

        # generate derivative assignments, filtered as specified by IBP, C, P, and T
        derivative_assignments_dict, IBP_deleted, C_deleted, P_deleted, T_deleted = generate_derivative_assignments_IBP_CPT_filtered(field_multiset, IBP, C, P, T)

        # populate dictionaries of deleted derivative assignments
        IBP_deleted_dict[frozendict(field_multiset)] = IBP_deleted
        C_deleted_dict[frozendict(field_multiset)] = C_deleted
        P_deleted_dict[frozendict(field_multiset)] = P_deleted
        T_deleted_dict[frozendict(field_multiset)] = T_deleted

        # for each field multiset, declare empty dictionary with derivative assignments as keys and a list of
        # deleted full contractions for each key
        EOM_deleted_dict[frozendict(field_multiset)] = {}
        DDcomm_deleted_dict[frozendict(field_multiset)] = {}
        FTselfcontracted_deleted_dict[frozendict(field_multiset)] = {}

        for outer_partition in derivative_assignments_dict.keys():
            for inner_partition in derivative_assignments_dict[outer_partition].keys():
                for derivative_assignment in derivative_assignments_dict[outer_partition][inner_partition]:
                    def make_hashable(derivative_assignment):
                        return tuple(tuple(sublist) for sublist in derivative_assignment)

                    if verbose:
                        print('')
                        print('derivative_assignment: ' + str(derivative_assignment))

                    lorentz_contractions_unique = generate_lorentz_contractions(field_multiset, derivative_assignment)

                    ############### filter out eom-redundant operators
                    lorentz_contractions_unique_filtered = []

                    # for each of eom, ddcomm, and ft self contracted deletions, declare empty list of deleted items for
                    # each field multiset and derivative assignment
                    EOM_deleted_dict[frozendict(field_multiset)][make_hashable(derivative_assignment)] = []
                    DDcomm_deleted_dict[frozendict(field_multiset)][make_hashable(derivative_assignment)] = []
                    FTselfcontracted_deleted_dict[frozendict(field_multiset)][make_hashable(derivative_assignment)] = []

                    # determine whether to remove full contractions on the basis of eom, DD commutator, or self
                    # contraction of antisymmetric operators
                    for i in range(len(lorentz_contractions_unique)):
                        contraction_multiset = lorentz_contractions_unique[i]
                        if verbose:
                            print('contraction_multiset: ' + str(contraction_multiset))
                        eom_remove = EOM and EOM_remove(field_multiset, derivative_assignment, contraction_multiset)
                        dd_comm_remove = DD_comm and DD_comm_remove(field_multiset, derivative_assignment, contraction_multiset)
                        ft_selfcontracted_remove = FT_selfcontracted and FT_selfcontracted_remove(field_multiset, derivative_assignment, contraction_multiset)

                        if eom_remove or dd_comm_remove or ft_selfcontracted_remove:
                            if eom_remove:
                                EOM_deleted_dict[frozendict(field_multiset)][make_hashable(derivative_assignment)].append(contraction_multiset)
                                if verbose:
                                    print('removed based on EoM')
                            if dd_comm_remove:
                                DDcomm_deleted_dict[frozendict(field_multiset)][make_hashable(derivative_assignment)].append(contraction_multiset)
                                if verbose:
                                    print('removed based on [D,D]~F')
                            if ft_selfcontracted_remove:
                                FTselfcontracted_deleted_dict[frozendict(field_multiset)][make_hashable(derivative_assignment)].append(contraction_multiset)
                                if verbose:
                                    print('removed based on self contraction of antisymmetric operator')
                            continue
                        else:
                            lorentz_contractions_unique_filtered.append(contraction_multiset)

                    ###################
                    operator_dict_field_multiset[make_hashable(derivative_assignment)] = lorentz_contractions_unique_filtered

        operator_dict[frozendict(field_multiset)] = operator_dict_field_multiset

    return operator_dict, IBP_deleted_dict, C_deleted_dict, P_deleted_dict, T_deleted_dict, EOM_deleted_dict, DDcomm_deleted_dict, FTselfcontracted_deleted_dict



def display_operators(massDim, C=True, P=True, T=True, IBP=True, EOM=True, DD_comm=True, FT_selfcontracted=True, verbose=False, show_removed=False):
    # displays contents of operator_dict produced in generate_operators()
    if massDim < 2 or type(massDim) != int:
        print('massDim argument must be an int greater than or equal to 2.')
        return

    operator_dict, IBP_deleted_dict, C_deleted_dict, P_deleted_dict, T_deleted_dict, EOM_deleted_dict, DDcomm_deleted_dict, FTselfcontracted_deleted_dict  = generate_operators(massDim, C, P, T, IBP, EOM, DD_comm, FT_selfcontracted, verbose)

    file = open('operators_' + 'massDim_' + str(massDim) + '_C_' + str(C) + '_P_' + str(P) + '_T_' + str(T) + '_IBP_' + str(IBP) + '_EOM_' + str(EOM) + '_DDcomm_' + str(DD_comm) + '_FTself_' + str(FT_selfcontracted) + '.txt', 'w')
    for field_multiset in operator_dict.keys():
        if not operator_dict[field_multiset]:
            continue
        for derivative_assignment in operator_dict[field_multiset].keys():
            symbol_grouping, _ = encode_lorentz_grouping(field_multiset, derivative_assignment)
            if not operator_dict[field_multiset][derivative_assignment]:
                    continue
            for lorentz_contraction_multiset in operator_dict[field_multiset][derivative_assignment]:
                #print to console and file
                print('')
                print('', file=file)
                print('')
                print('', file=file)

                print('CONSTITUENT FIELDS: ')
                print('CONSTITUENT FIELDS: ', file=file)
                print(dict(field_multiset))
                print(dict(field_multiset), file=file)
                print('DERIVATIVE ASSIGNMENT: ')
                print('DERIVATIVE ASSIGNMENT: ', file=file)
                print(derivative_assignment)
                print(derivative_assignment, file=file)
                print('GROUPING OF LORENTZ OBJECTS: ')
                print('GROUPING OF LORENTZ OBJECTS: ', file=file)
                print(symbol_grouping)
                print(symbol_grouping, file=file)
                print('LORENTZ INDEX CONTRACTIONS: ')
                print('LORENTZ INDEX CONTRACTIONS: ', file=file)
                print(lorentz_contraction_multiset)
                print(lorentz_contraction_multiset, file=file)

    if show_removed:
        print('')
        print('', file=file)
        print('')
        print('', file=file)
        print('')
        print('', file=file)
        print('')
        print('', file=file)
        print('####################################### OPERATORS REMOVED IN BASIS REDUCTION #####################################')
        print('####################################### OPERATORS REMOVED IN BASIS REDUCTION #####################################', file=file)
        basis_reduction_variables = {'C': C, 'P': P, 'T': T, 'IBP': IBP, 'EOM': EOM, 'DD_comm': DD_comm, 'FT_selfcontracted': FT_selfcontracted}
        basis_reduction_dictionaries = {'C': C_deleted_dict, 'P': P_deleted_dict, 'T': T_deleted_dict, 'IBP': IBP_deleted_dict, 'EOM': EOM_deleted_dict, 'DD_comm': DDcomm_deleted_dict, 'FT_selfcontracted': FTselfcontracted_deleted_dict}
        basis_reduction_operator_type = {'C': 'DERIVATIVE ASSIGNMENT', 'P': 'DERIVATIVE ASSIGNMENT', 'T': 'DERIVATIVE ASSIGNMENT', 'IBP': 'DERIVATIVE ASSIGNMENT', 'EOM': 'FULL LORENTZ INDEX CONTRACTION', 'DD_comm': 'FULL LORENTZ INDEX CONTRACTION', 'FT_selfcontracted': 'FULL LORENTZ INDEX CONTRACTION'}
        def make_hashable(derivative_assignment):
            return tuple(tuple(sublist) for sublist in derivative_assignment)
        for key in basis_reduction_variables.keys():
            if basis_reduction_variables[key]:
                print('')
                print('', file=file)
                print('')
                print('', file=file)
                print(str(basis_reduction_operator_type[key]) + 'S REMOVED BASED ON ' + str(key) + ':')
                print(str(basis_reduction_operator_type[key]) + 'S REMOVED BASED ON ' + str(key) + ':', file=file)
                for field_multiset in basis_reduction_dictionaries[key].keys():
                    if basis_reduction_dictionaries[key][field_multiset]:
                        if basis_reduction_operator_type[key] == 'DERIVATIVE ASSIGNMENT':
                            print('')
                            print('', file=file)
                            print(dict(field_multiset))
                            print(dict(field_multiset), file=file)
                            for derivative_assignment in basis_reduction_dictionaries[key][field_multiset]:
                                print(derivative_assignment)
                                print(derivative_assignment, file=file)

                        elif basis_reduction_operator_type[key] == 'FULL LORENTZ INDEX CONTRACTION':
                            for derivative_assignment in basis_reduction_dictionaries[key][field_multiset]:
                                if basis_reduction_dictionaries[key][field_multiset][make_hashable(derivative_assignment)]:
                                    for full_contraction in basis_reduction_dictionaries[key][field_multiset][make_hashable(derivative_assignment)]:
                                        symbol_grouping, _ = encode_lorentz_grouping(field_multiset, derivative_assignment)
                                        print('')
                                        print('', file=file)
                                        print(dict(field_multiset))
                                        print(dict(field_multiset), file=file)
                                        print(derivative_assignment)
                                        print(derivative_assignment, file=file)
                                        print('Lorentz index grouping: ' + str(symbol_grouping))
                                        print('Lorentz index grouping: ' + str(symbol_grouping), file=file)
                                        print(full_contraction)
                                        print(full_contraction, file=file)
    file.close()
    return




if __name__ == "__main__":

    # Have user input mass dimension
    print('Input mass dimension:')
    massDim = input()
    try:
        massDim = int(massDim)
    except:
        pass

    while type(massDim) != int or (type(massDim) == int and int(massDim) <= 0):
        print('ERROR: Input must be a positive integer')
        print('Input mass dimension:')
        try:
            massDim = int(input())
        except:
            massDim = input()

    # Have user input whether to impose C symmetry
    print('Impose charge conjugation invariance (C)? (y/n)')
    C_string = input()
    while C_string != 'y' and C_string != 'n':
        print('ERROR: Input must be "y" or "n".')
        C_string = input()
    if C_string == 'y':
        C = True
    if C_string == 'n':
        C = False

    # Have user input whether to impose P symmetry
    print('Impose parity invariance (P)? (y/n)')
    P_string = input()
    while P_string != 'y' and P_string != 'n':
        print('ERROR: Input must be "y" or "n".')
        P_string = input()
    if P_string == 'y':
        P = True
    if P_string == 'n':
        P = False

    # Have user input whether to impose T symmetry
    print('Impose time reversal invariance (T)? (y/n)')
    T_string = input()
    while T_string != 'y' and T_string != 'n':
        print('ERROR: Input must be "y" or "n".')
        T_string = input()
    if T_string == 'y':
        T = True
    if T_string == 'n':
        T = False

    # Have user input whether to impose IBP reduction
    print('Perform integration by parts basis reduction? (y/n)')
    IBP_string = input()
    while IBP_string != 'y' and IBP_string != 'n':
        print('ERROR: Input must be "y" or "n".')
        IBP_string = input()
    if IBP_string == 'y':
        IBP = True
    if IBP_string == 'n':
        IBP = False

    # Have user input whether to impose EoM reduction
    print('Perform equations of motion basis reduction? (y/n)')
    EOM_string = input()
    while EOM_string != 'y' and EOM_string != 'n':
        print('ERROR: Input must be "y" or "n".')
        EOM_string = input()
    if EOM_string == 'y':
        EOM = True
    if EOM_string == 'n':
        EOM = False

    # Have user input whether to impose reduction by covariant derivative commutator
    print('Perform basis reduction using covariant derivative commutator? (y/n)')
    DDcomm_string = input()
    while DDcomm_string != 'y' and DDcomm_string != 'n':
        print('ERROR: Input must be "y" or "n".')
        DDcomm_string = input()
    if DDcomm_string == 'y':
        DD_comm = True
    if DDcomm_string == 'n':
        DD_comm = False

    # Have user input whether to impose reduction by covariant derivative commutator
    print('Remove operators with self contraction of antisymmetric tensors? (y/n)')
    FTselfcontracted_string = input()
    while FTselfcontracted_string != 'y' and FTselfcontracted_string != 'n':
        print('ERROR: Input must be "y" or "n".')
        FTselfcontracted_string = input()
    if FTselfcontracted_string == 'y':
        FT_selfcontracted = True
    if FTselfcontracted_string == 'n':
        FT_selfcontracted = False

    # Have user input whether to show progress of operator generation
    print('Print progress of operator generation? (y/n)')
    verbose_string = input()
    while verbose_string != 'y' and verbose_string != 'n':
        print('ERROR: Input must be "y" or "n".')
        verbose_string = input()
    if verbose_string == 'y':
        verbose = True
    if verbose_string == 'n':
        verbose = False

    # Have user input whether to show operators that were removed
    print('Print removed operators? (y/n)')
    removed_string = input()
    while removed_string != 'y' and removed_string != 'n':
        print('ERROR: Input must be "y" or "n".')
        removed_string = input()
    if removed_string == 'y':
        show_removed = True
    if removed_string == 'n':
        show_removed = False



    display_operators(massDim, C, P, T, IBP, EOM, DD_comm, FT_selfcontracted, verbose, show_removed)


    '''
    ################################ TEST OF generate_field_multisets() ################################
    massDim = 50
    field_symbols = ['D', 'F', 'S', 'V', 'T', 'Vp', 'Sp']
    field_massDims = [1, 2, 3, 3, 3 ,3 ,3]
    field_lorentzRanks = [1, 2, 0, 1, 2, 1, 0]
    field_multisets_dict = generate_field_multisets(massDim, field_symbols, field_massDims, field_lorentzRanks)
    for (key, value) in field_multisets_dict.items():
        print(key)
        print(value)
    '''

    '''
    ################################ TEST OF generate_derivative_assignments_IBP_CPT_filtered() ################################
    from frozendict import frozendict
    
    massDim = 30
    field_symbols = ['D', 'F', 'S', 'V', 'T', 'Vp', 'Sp']
    field_massDims = [1, 2, 3, 3, 3 ,3 ,3]
    field_lorentzRanks = [1, 2, 0, 1, 2, 1, 0]
    field_multisets_dict = generate_field_multisets(massDim, field_symbols, field_massDims, field_lorentzRanks)
    field_multiset_list = field_multisets_dict[massDim]
    derivative_assignments = {}
    for i in range(len(field_multiset_list)):
        print('')
        print('')
        print(str(i) + '/' + str(len(field_multiset_list)))
        field_multiset = field_multiset_list[i]
        print(field_multiset)
        derivative_assignments_dict, IBP_deleted, C_deleted, P_deleted, T_deleted = generate_derivative_assignments_IBP_CPT_filtered(field_multiset, IBP=False, C=False, P=False, T=False)
        derivative_assignments[frozendict(field_multiset)] = derivative_assignments_dict
    '''
