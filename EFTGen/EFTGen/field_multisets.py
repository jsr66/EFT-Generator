def totalMassDim(field_multiplicities, field_massDims):
    if len(field_multiplicities) != len(field_massDims):
        print('totalMassDim(): Input lists must have same length')
        return None
    n_types = len(field_multiplicities)
    total_massDim = 0
    for i in range(n_types):
        total_massDim += field_multiplicities[i]*field_massDims[i]
    return total_massDim

def totalLorentzRank(field_multiplicities, field_lorentzRanks):
    if len(field_multiplicities) != len(field_lorentzRanks):
        print(field_multiplicities)
        print(field_lorentzRanks)
        print('totalLorentzRank(): Input lists must have same length')
        return None
    n_types = len(field_multiplicities)
    total_lorentzRank = 0
    for i in range(n_types):
        total_lorentzRank += field_multiplicities[i]*field_lorentzRanks[i]
    return total_lorentzRank

def convert_to_dict(field_multiplicities, field_types, remove_zeros):
    # here, field_types = ['D', 'F', 'S', 'V', 'T', 'Vp', 'Sp'], but may be adjusted to include more
    # field types (for example, in theories with multiple Fermi fields or non-Abelian gauge fields) if necessary.
    if len(field_multiplicities) != len(field_types):
        print('Length of field_multiplicities input to convert_to_dict() must be equal to length of field_types.')
        return None
    n_types = len(field_types)
    field_multiset = {}
    if remove_zeros:
        for i in range(n_types):
            if field_multiplicities[i] != 0:
                field_multiset[field_types[i]] = field_multiplicities[i]
        return field_multiset
    else:
        for i in range(n_types):
            field_multiset[field_types[i]] = field_multiplicities[i]
        return field_multiset

def generate_field_multisets(massDim, field_symbols, field_massDims, field_lorentzRanks):
    n_types = len(field_symbols)
    zero_list = [0 for i in range(len(field_symbols))]
    front = [zero_list]
    field_multiplicities_dict = {}

    counter = 0
    while front:
        #print(counter)
        front_new = front.copy()
        for i in range(len(front)):
            #print('i: ' + str(i))
            multiplicities_list = front[i]
            #print(multiplicities_list)
            # find index last non-zero element
            if multiplicities_list == len(field_symbols)*[0]:
                index_last_nonzero = 0
            else:
                index_last_nonzero = [i for (i, m) in enumerate(multiplicities_list) if m != 0][-1]

            for j in range(index_last_nonzero, n_types):
                multiplicities_list_new = multiplicities_list.copy()
                multiplicities_list_new[j] += 1
                M = totalMassDim(multiplicities_list_new, field_massDims)
                L = totalLorentzRank(multiplicities_list_new, field_lorentzRanks)
                if M < massDim:
                    #print('front_new before append:' + str(front_new))
                    front_new.append(multiplicities_list_new)
                    #print('front_new after append:' + str(front_new))
                if M <= massDim and L%2 == 0:
                    try:
                        field_multiplicities_dict[M].append(multiplicities_list_new)
                    except:
                        field_multiplicities_dict[M] = [multiplicities_list_new]
            #print('len(front): ' + str(len(front)))
            #print('len(front_new): ' + str(len(front_new)))
            #del front_new[i]
            front_new.remove(multiplicities_list)
        front = front_new
        counter += 1

    # create field multisets dict, where each value is a dict rather than a list
    field_multisets_dict = {}

    for key in field_multiplicities_dict.keys():
        field_multisets_dict[key] = []
        field_multiplicities_list = field_multiplicities_dict[key]
        for field_multiplicities in field_multiplicities_list:
            # if number of derivatives (field_multiplicities[0]) equal to mass dimension (key), do not include
            # since multiset is all derivatives
            if field_multiplicities[0] == key:
                continue
            field_multiset = convert_to_dict(field_multiplicities, field_symbols, True)
            field_multisets_dict[key].append(field_multiset)
    return field_multisets_dict
