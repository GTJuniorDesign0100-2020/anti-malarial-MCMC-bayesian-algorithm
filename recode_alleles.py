import numpy as np

def recodeallele(alleles_definitions_subset,proposed):

    # verify shapes.
    if len(alleles_definitions_subset.shape) != 2 or alleles_definitions_subset.shape[1] != 2:
        raise ValueError(f"Improper alleles_definition_subset shape {alleles_definitions_subset.shape} (expected (:,2))")
    if len(proposed.shape) > 1:
        raise ValueError("Improper proposed vector shape {proposed.shape} (expected 1D vector)")

    # alleles_definitions_subset ranges guaranteed to be non-overlapping, so it will always fall within either 0 or exactly 1 of the ranges (i.e. rows)
    result = np.argwhere(np.logical_and(proposed > alleles_definitions_subset[:,0], proposed <= alleles_definitions_subset[:,1]))
    result = result.reshape(-1)
    if (result.size == 0):
        result = np.nan
    else:
        return result[0]
    return result
