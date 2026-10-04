
import numpy as np

def approximately_equal(new_data, old_data):
    """
    (doc to be written)
    """
    absolute_difference = 0
    for i in [0, 1] if len(new_data) == 3 else [0, 1, 2]:
        new = np.array(new_data[i])
        old = np.array(old_data[i])
        absolute_difference += np.sum(np.abs(new-old))/np.sum(np.abs(new)+np.abs(old))
    print(absolute_difference)
    return absolute_difference