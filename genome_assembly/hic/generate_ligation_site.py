# generate all possible combinations for Ns in the supplied ligation sites to be included in juicer.sh.
# for convenience if you have a restriction enzyme cocktail. 

# for example: DpnII-HinFI-MseI-DdeI.
# input ligation sites (set in re parameter)
# re = ['GATCGATC', 'GANTANTC', 'TTATAA', 'CTNATNAG']
# will print output as:
# 'GATCGATC|GAATAATC|GAATATTC|GAATACTC|GAATAGTC|GATTAATC|GATTATTC|GATTACTC|GATTAGTC|GACTAATC|GACTATTC|GACTACTC|GACTAGTC|GAGTAATC|GAGTATTC|GAGTACTC|GAGTAGTC|TTATAA|CTAATAAG|CTAATTAG|CTAATCAG|CTAATGAG|CTTATAAG|CTTATTAG|CTTATCAG|CTTATGAG|CTCATAAG|CTCATTAG|CTCATCAG|CTCATGAG|CTGATAAG|CTGATTAG|CTGATCAG|CTGATGAG'

import itertools
re = ['GATCGATC', 'GANTANTC', 'TTATAA', 'CTNATNAG']
result = []
for sequence in re:
    if 'N' in sequence:
        n_indices = [i for i, char in enumerate(sequence) if char == 'N']
        n_count = len(n_indices)
        combinations = list(itertools.product('ATCG', repeat=n_count))
        for comb in combinations:
            new_sequence = sequence
            for i, char in zip(n_indices, comb):
                new_sequence = new_sequence[:i] + char + new_sequence[i+1:]
            result.append(new_sequence)
    else:
        result.append(sequence)
print(result)
print("|".join(result))