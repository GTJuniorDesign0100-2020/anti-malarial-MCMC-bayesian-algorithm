import numpy as np

def recodeallele(alleles_definitions_subset,proposed):

    # verify shapes.
    if alleles_definitions_subset.shape != (1,2):
        raise("Improper alleles_definition_subset.")
    if len(proposed.shape) > 1 or proposed.shape[2] != 1:
        raise("Improper proposed vector.")

    condition = (proposed > alleles_definitions_subset[0] & proposed <= alleles_definitions_subset[1])
    result = proposed[condition]
    return result