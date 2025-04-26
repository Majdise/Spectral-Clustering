import sys

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

    
def validInput(args) :
    if (len(args) != 2 and len(args) != 3) :
        print("An Error Has Occurred")
        return False
    if ((not args[0].isnumeric()) or int(args[0]) <= 1) :
        print("Invalid number of clusters!")
        return False
    if (len(args) == 3) :
        if ((not args[1].isnumeric()) or int(args[1]) <= 1 or int(args[1]) >= 1000) :
            print("Invalid maximum iteration!")
            return False
    return True


def readDataPoints(filePath) :
    try :
        dataPoints =[]
        file = open(filePath, "r")
        while (True) :
            line = file.readline()
            if (not line) :
                break
            line = line.rstrip('\n').split(",")
            for i in range(len(line)) :
                line[i] = float(line[i])
            dataPoints.append(line)
    except :
        sys.exit("An Error Has Occurred")
    file.close()
    return dataPoints


def initializeClusters(inputAsMat, k) :
    initialClusters =[]
    for i in range(k) :
        initialClusters.append(Cluster(inputAsMat[i]))
    return initialClusters


def updateClustersAndCheckConvergence(clusters) :
    eps = 10**(-3)
    converged = True
    for cluster in clusters :
        newCentroid = cluster.getNewCentroid()
        if ((converged) and (cluster.euclideanDistance(newCentroid) >= eps)) :
            converged = False
        cluster.updateCentroid(newCentroid)
    return converged


def kmeans(inputAsMat, k, iter) :
    clusters = initializeClusters(inputAsMat, k)
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
        if (updateClustersAndCheckConvergence(clusters)) :
            break
    return clusters


def printCentroids(clusters) :
    for cluster in clusters :
        centroidAsStr = ",".join(str('%.4f' % point) for point in cluster.centroid)
        print(centroidAsStr)


def main(args) :
    assert(validInput(args))
    k = int(args[0])
    if (len(args) == 2) :
        iter = 200
        filePath = args[1]
    else :
        iter = int(args[1])
        filePath = args[2]
    inputAsMat = readDataPoints(filePath)
    if (k >= len(inputAsMat)) :
        sys.exit("Invalid number of clusters!")
    clusters = kmeans(inputAsMat, k, iter)
    printCentroids(clusters)


if (__name__ == "__main__") :
    main(sys.argv[1:])