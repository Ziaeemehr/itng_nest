import numpy as np

def adj_list_to_matrix(adj_list):
    A = np.asarray(adj_list)
    
    n = np.max(np.asarray(adj_list))
    adj_matrix = np.zeros((n,n))

    for i in range(A.shape[0]):
        for j, k in adj_list[i]:
            adj_matrix[j-1,k-1] = 1
    return adj_matrix


adj_list = [((1, 4), (1, 5), (1, 6), (1, 7)),
            ((2, 4), (2, 5), (2, 6), (2, 7)),
            ((3, 4), (3, 5), (3, 6), (3, 7))]

print adj_list_to_matrix(adj_list)
