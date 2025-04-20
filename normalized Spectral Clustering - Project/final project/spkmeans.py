import numpy as np
import sys
import math
import spkmeans


# Class Cluster 
class Cluster :
    def __init__(self, datapoint) :
        self.centroid = datapoint
        self.dataPointsSum = [0] * len(datapoint)
        self.size = 0

    
    def addVectorToCluster(self, vector) :
        self.size += 1
        for i in range(len(vector)) :
            self.dataPointsSum[i] += vector[i]

    
    def euclideanDistance(self, vector) :
        squareSum = 0
        for i in range(len(vector)) :
            squareSum += (self.centroid[i] - vector[i]) ** 2
        return squareSum ** 0.5


    def getNewCentroid(self) :
        if (self.size == 0) :
            return self.centroid
        centroid =[]
        for i in range(len(self.dataPointsSum)) :
            centroid.append(self.dataPointsSum[i] / self.size)
        return centroid

    
    def updateCentroid(self, centroid) :
        self.centroid = centroid
        self.dataPointsSum = [0] * len(self.centroid)
        self.size = 0
 

def parse_input():
    if len(sys.argv) == 4:
        try:
            k_int = int(sys.argv[1])
            goal_str = str(sys.argv[2])
            file_name_str = str(sys.argv[3])
        except ValueError:
            assert False, "An Error Has Occurred"
    else:
        try:
            k_int = 0
            goal_str = str(sys.argv[1])
            file_name_str = str(sys.argv[2])
        except ValueError:
            assert False, "An Error Has Occurred"
    return k_int, goal_str, file_name_str



def read_data_points_from_file(file_to_open):
    try:
        file = open(file_to_open, 'r')
    except IOError:
        assert False, "An Error Has Occurred"
    lines_arr = file.readlines()
    data_points_array = []

    for line in lines_arr:
        dot = [float(word) for word in line.split(sep=",")]
        data_points_array.append(dot)

    return data_points_array


def minusZero(val_to_check):
    if -0.00005 < val_to_check < 0:
        val_to_check = -0.00001
    return val_to_check


def print_list_of_lists(list_of_lists):
    tmp_lst = []
    for row in list_of_lists:
        for val in row:
            tmp_val = minusZero(val)
            tmp_lst.append(format(tmp_val, '.4f'))
        print(','.join(str(tmpVal) for tmpVal in tmp_lst))
        tmp_lst = []


def random_nearest_cluster_index_by_weights(weights_arr):
    weights_len = len(weights_arr)
    return int(np.random.choice(range(weights_len), 1, replace=False, p=weights_arr))


def get_cluster_distance(data_point, cluster):
    cluster_as_vector = np.array(cluster)
    data_point_as_vector = np.array(data_point)
    distance = np.linalg.norm(cluster_as_vector - data_point_as_vector)
    return distance


def find_initial_clusters(data_points_mat, k_val):
    n = len(data_points_mat)
    random_index = np.random.choice(n)
    first_cluster = data_points_mat[random_index]
    clusters_list = [first_cluster]
    init_indexes = [random_index]
    distances_list = [float("inf")] * n
    for i in range(1, k_val):
        for j, data_point in enumerate(data_points_mat):
            d = get_cluster_distance(data_point, clusters_list[i - 1])
            distances_list[j] = min(distances_list[j], d)
        sum_distances = sum(distances_list)
        weights = [dist / sum_distances for dist in distances_list]
        index = random_nearest_cluster_index_by_weights(weights)
        init_indexes.append(index)
        clusters_list.append(data_points_mat[index])
    assert (k_val == len(init_indexes))
    return init_indexes


def initializeClusters(inputAsMat, initial_indices) :
    initialClusters =[]
    for i in initial_indices :
        initialClusters.append(Cluster(inputAsMat[i]))
    return initialClusters


def updateClusters(clusters) :
    for cluster in clusters :
        newCentroid = cluster.getNewCentroid()
        cluster.updateCentroid(newCentroid)

    
def printCentroids(clusters) :
    for cluster in clusters :
        centroidAsStr = ",".join(str('%.4f' % point) for point in cluster.centroid)
        print(centroidAsStr)


def kmeans(inputAsMat,k, initial_indices, iter) :
    clusters = initializeClusters(inputAsMat, initial_indices)
    iterations = 0
    while (iterations < iter) :
        iterations += 1
        for vector in inputAsMat :
            minDistance = float("inf")
            for i in range(len(clusters)) :
                currDst = clusters[i].euclideanDistance(vector)
                if (currDst < minDistance) :
                    minDistance = currDst
                    minClusterIndex = i
            clusters[minClusterIndex].addVectorToCluster(vector)
        updateClusters(clusters)
    return clusters
    
    
if __name__ == "__main__":
    np.random.seed(0)
    k, goal, file_name = parse_input()

    data_points = read_data_points_from_file(file_name)

    if goal == 'wam':
        wam_matrix_result = spkmeans.wam(data_points)
        print_list_of_lists(wam_matrix_result)

    elif goal == 'ddg':
        ddg_matrix_result = spkmeans.ddg(data_points)
        print_list_of_lists(ddg_matrix_result)

    elif goal == 'gl':
        gl_matrix_result = spkmeans.gl(data_points)
        print_list_of_lists(gl_matrix_result)

    elif goal == 'jacobi':
        jacobi_eigen_values_and_vectors = spkmeans.jacobi(data_points)
        print_list_of_lists(jacobi_eigen_values_and_vectors)

    elif goal == 'spk':
        t_matrix, k_updated = spkmeans.createU(data_points, k)
        initial_indexes = find_initial_clusters(t_matrix, k_updated)
        iter = 300
        clusters = kmeans(t_matrix, k_updated,initial_indexes, iter)

        # print the initial clusters indexes
        print(','.join(str(tmpVal) for tmpVal in initial_indexes))
        # print the clusters
        printCentroids(clusters)

