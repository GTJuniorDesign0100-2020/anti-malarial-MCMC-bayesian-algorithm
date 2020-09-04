from random import uniform
import numpy as np
import pandas as pd

def switch_hidden(x):
    z = random.uniform((10,))
    if sum() > 0: #if hidden alleles exist TODO: Figure out comma notation in arrays
        if len(np.where(x == 1, np.concatonate())) > 1: #TODO same as above
            chosen = np.random.choice(np.where(x == 1, np.concatonate()), 1, False) #TODO same as above w/
        else:
            chosen = np.where(x == 1, np.concatonate()) #TODO same as above
    if classification[x] == 0: #reinfection
        if chosen <= nloci * maxMoi: #day0 hidden allele
            chosenLocus = math.ceil(chosen / maxMoi)
            old = recoded0[x][chosen]
