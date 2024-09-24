import numpy as np
import pandas as pd
import math
import sys
import mykmeanssp as km


def ReadData(fileName):
    # Read the file, splitting by lines
    try:
        f = open(fileName, 'r')
    except:
        print_error()
    lines = f.read().splitlines()
    f.close()
    items = []
    for i in range(0, len(lines)):
        line = lines[i].split(',')
        itemFeatures = []
        for j in range(len(line)):
            # Convert feature value to float
            if line[j] != '':
                v = float(line[j])
                # Add feature value to dict
                itemFeatures.append(v)
        items.append(itemFeatures)
    for i in range(len(items)):
        try:
            items.remove([])
        except:
            pass
    return items


def getfirstkeigen(eigen, k):  # sort the eigenvectors
    eigen = eigen[:, eigen[len(eigen) - 1, :].argsort()]
    return eigen[:, :k]


def print_matrix(matrix):  # printing matrix by 4f format
    for vec in matrix:
        for i in range(len(vec) - 1):
            print(str("%.4f" % vec[i]) + ",", end="")
        print(str("%.4f" % vec[len(vec) - 1]))


### test ###
def normal_matrix(mat):
    n = len(mat) - 1
    indexes = [4, 8, 2, 7, 5, 0, 3, 9, 6, 1]
    temp = [[0 for i in range(n)] for i in range(n + 1)]
    for i in range(n + 1):
        for j in range(n):
            temp[i][j] = mat[i][indexes[j]]
    print_matrix(temp)
    res = temp[:n][:]
    for i in range(n):
        mult = 1 / res[n - 1][i]
        for j in range(n):
            res[j][i] *= mult
    return res


### test ###

def initial_kmeans_pp(items, k):  # initializing the first centroids indexes according to HW2
    centroids_indexes = []
    index_lst = [i for i in range(len(items))]
    items_dst = np.full(len(items), float('inf'))
    p_lst = np.full(len(items), (1 / len(items)))
    for i in range(k):
        if i == 0:
            centroids_indexes.append(np.random.choice(index_lst))
        else:
            centroids_indexes.append(np.random.choice(index_lst, p=p_lst))
        total = 0
        for j in range(len(items)):
            new_dst = euclidean_distance(items[(centroids_indexes[i])], items[j])
            if new_dst < items_dst[j]:
                items_dst[j] = new_dst
            total += items_dst[j]
        for j in range(len(items)):
            p_lst[j] = items_dst[j] / total
    return centroids_indexes


def euclidean_distance(x, y):  # calculate euclidean distance
    sum = 0
    for i in range(len(x)):
        sum += math.pow(x[i] - y[i], 2)
    return math.sqrt(sum)


def eigen_gap_heuristic(e_values, n):  # calculate the optimal k according to eigengap_heuristic
    e_values_sorted = np.sort(e_values)
    e_values_diff = np.diff(e_values_sorted)[0]
    e_values_diff = e_values_diff[:math.ceil(int(n / 2))]
    k = np.argmax(e_values_diff) + 1
    return k


def calc_spkmeans(goal, file_name, k=None):  # main algorithem
    items = ReadData(file_name)
    if not (len(items) > 0):
        print_error()
    d = len(items[0])
    n = len(items)
    if goal == 'spk':  # spk case
        lm = km.gl(items)
        eigen = km.jacobi(lm)
        for i in range(len(items)):  # change -0.0000 eigenvalues to 0 and multiply eigen vector by -1
            if ("%.4f" % eigen[len(items)][i]) == "-0.0000":
                eigen[len(items)][i] = 0
                for j in range(len(eigen) - 1):
                    eigen[j][i] = eigen[j][i] * -1
        eigen = np.array(eigen)
        eigenvalues = eigen[len(eigen) - 1:, :]
        if k == None:
            k = eigen_gap_heuristic(eigenvalues, len(eigen) - 1)
        try:
            k = float(k)
        except:
            print_error()
        if not (1 < k < n and (float(k)).is_integer()):
            print_error()
        k = int(k)
        u_mat = getfirstkeigen(eigen, k)
        u_mat = u_mat[:len(u_mat) - 1, :]
        u_mat = np.ndarray.tolist(u_mat)
        centroids_indexes = initial_kmeans_pp(u_mat, k)
        for i in range(len(centroids_indexes) - 1):
            print(str(centroids_indexes[i]) + ",", end="")
        print(str(centroids_indexes[len(centroids_indexes) - 1]))
        lst_centroids = km.spk(u_mat, centroids_indexes, k)  # passing k and not dim
        print_matrix(lst_centroids)
        return 1
    elif goal == 'wam':  # wam case
        w = km.wam(items)
        print_matrix(w)
        return 1
    elif goal == 'ddg':  # ddg case
        d = km.ddg(items)
        print_matrix(d)
        return 1
    elif goal == 'gl':  # gl case
        l = km.gl(items)
        print_matrix(l)
        return 1
    elif goal == 'jacobi':  # jacobi case
        if d != n:
            print_error()
        eigen = km.jacobi(items)
        for i in range(len(items)):
            if ("%.4f" % eigen[len(items)][i]) == "-0.0000":
                eigen[len(items)][i] = 0
                for j in range(len(eigen) - 1):
                    eigen[j][i] = eigen[j][i] * -1
        for i in range(len(items) - 1):
            print(str("%.4f" % eigen[len(items)][i]) + ",", end="")
        print(str("%.4f" % eigen[len(items)][len(items) - 1]))
        eigen = np.array(eigen)
        eigen = eigen[:len(items), :]
        print_matrix(eigen)
        return 1
    else:
        print_error()


def print_error():
    print("An Error Has Occurred")
    exit(1)


np.random.seed(0)
if len(sys.argv) == 4:
    calc_spkmeans(sys.argv[2], sys.argv[3], sys.argv[1])
elif len(sys.argv) == 3:
    calc_spkmeans(sys.argv[1], sys.argv[2])
else:
    print_error()

#########