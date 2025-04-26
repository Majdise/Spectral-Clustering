import sys
import numpy as np
import pandas as pd
import mykmeanssp


def validate_arguments():
    # read the CMD arguments and check if valid
    cmd_arguments_list = sys.argv
    arguments_cnt = len(cmd_arguments_list)
    assert arguments_cnt == 4 or arguments_cnt == 5, "4 or 5 arguments needed"
    assert cmd_arguments_list[1].isnumeric() and int(cmd_arguments_list[1]) > 0, "k should be an integer > 0"
    if arguments_cnt == 5:
        assert cmd_arguments_list[2].isnumeric() and int(cmd_arguments_list[2]) > 0, "max_iter should be an integer > 0"
    return cmd_arguments_list, arguments_cnt


def read_files_and_merge(file_name_1, file_name_2):
    file1_pd = pd.read_csv(file_name_1, header=None)
    file2_pd = pd.read_csv(file_name_2, header=None)
    # inner join according to the first column
    result_pd = pd.merge(file1_pd, file2_pd, on=0)
    # sort the file using the 1st column
    result_pd = result_pd.sort_values(result_pd.columns[0])
    result_pd = result_pd.reset_index(drop=True)
    # remove the 1st column that we joined and sorted according to it.
    result_pd = result_pd.drop(result_pd.columns[0], axis=1)
    return result_pd


def get_variables():
    # initialize variables with the given arguments
    max_iteration = 300
    K = int(cmd_arguments[1])
    if arguments_count == 5:
        max_iteration = int(cmd_arguments[2])
        file1 = cmd_arguments[3]
        file2 = cmd_arguments[4]
    else:
        file1 = cmd_arguments[2]
        file2 = cmd_arguments[3]
    return K, max_iteration, file1, file2


def distance_between_2_vectors(vec1, vec2):
    return (np.linalg.norm(vec1 - vec2)) ** 2


def kmeans_pp():
    z = 1
    while z < k:
        i = 0
        j = 0
        while i < points_count:
            index = centroids_indices[z - 1]
            current_dist = (np.linalg.norm(np_vector_list[i] - np_vector_list[index])) ** 2
            if z == 1:
                extended_np_vector_list[i][vector_len] = current_dist
            else:
                extended_np_vector_list[i][vector_len] = min(current_dist, extended_np_vector_list[i][vector_len])
            i += 1
        dists_sum = extended_np_vector_list.sum(axis=0)[vector_len]
        while j < points_count:
            extended_np_vector_list[j][vector_len + 1] = extended_np_vector_list[j][vector_len] / dists_sum
            j += 1

        n_array = np.arange(points_count)
        centroids_indices.append(np.random.choice(n_array, p=extended_np_vector_list[:, vector_len + 1]))
        z += 1


def print_centroids(centroids):
    centroids = np.asarray(centroids).round(4)
    centroids = np.array_split(centroids, k)
    centroids = [centroids[i].tolist() for i in range(len(centroids))]

    for cent in centroids:
        print(','.join(str(val) for val in cent))


if __name__ == "__main__":
    cmd_arguments, arguments_count = validate_arguments()
    k, max_iter, file_name1, file_name2 = get_variables()

    vector_list = read_files_and_merge(file_name1, file_name2)
    vector_len = vector_list.shape[1]  # num of element in every vector

    points_to_cluster_init = vector_list.to_numpy().flatten().tolist()

    np_vector_list = vector_list.to_numpy()
    points_count = len(np_vector_list)  # num of vectors

    # validate that k < n
    assert k < points_count, "k should be < n (number of datapoints)"

    # extend the vector matrix with 2 columns for the distance and the property
    vector_list['distCol'] = np.nan
    vector_list['propCol'] = np.nan
    extended_np_vector_list = vector_list.to_numpy()

    # get initial random index to start with
    np.random.seed(0)
    rand_index = np.random.randint(0, (points_count - 1))
    centroids_indices = [rand_index]

    # run the kmeans_pp algorithm and update the 'centroids_indices'
    kmeans_pp()

    # create
    centroids_lists = [np_vector_list[ind].tolist() for ind in centroids_indices]
    centroids_init = [num for sublist in centroids_lists for num in sublist]

    num_of_lines = len(vector_list)
    centroid_indexes_string = ','.join([str(num) for num in centroids_indices])
    print(centroid_indexes_string)

    # run the C code of the kmeans (sending to it initial centroids that we calculate above using kmeans_pp algorithem)
    centroids_result = mykmeanssp.fit(k, max_iter, num_of_lines, vector_len, centroids_init, points_to_cluster_init)

    print_centroids(centroids_result)