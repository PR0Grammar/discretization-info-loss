import math;
import statistics;
import copy;

# H_n
def getEntropy(probDist):
    print("H("+ str(probDist) + ") = ")
    print("-1 * (")
    ans = 0
    
    for i in range(len(probDist)):
        p = probDist[i]
        ans += (p * math.log(p, 2))
        
        if(i == len(probDist) - 1):
            print(str(p) + " * log_2(" + str(p) +")")
            print(")")
        else:
            print(str(p) + " * log_2(" + str(p) +") +")

    ans = -(ans) 
    print("= " + str(ans))
    return ans
    
def createBins(data):
    dataSorted = sorted(data)
    # create bins for each term
    n = len(dataSorted)

    bins = []
    probBins = []

    # put the data into appropraite bins
    for i in range(len(dataSorted)):
        if(i > 0 and dataSorted[i - 1] == dataSorted[i]):
            bins[len(bins) - 1].append(dataSorted[i])
        else:
            bins.append([])
            bins[len(bins) - 1].append(dataSorted[i])
    
    # determine probability of each bin
    for i in range(len(bins)):
        bin_i = bins[i]
        probBins.append(len(bins[i]) / n)

    return (bins, probBins)

def getBinIndiciesToMerge(bins):
    if(len(bins) <= 1):
        return (-1, -1)
    minDiff = float('inf')
    first_i = -1
    second_i = -1

    for i in range(len(bins) - 1):
        first = statistics.mean(bins[i])
        second = statistics.mean(bins[i+1])
        diff = abs(first - second)
        if(diff < minDiff):
            first_i = i
            second_i = i + 1
            minDiff = diff
    
    return (first_i, second_i)

def mergeBins(bins, i, j, n):
    newBins = []
    probBins = []

    bins_i = list(bins[i])
    bins_j = list(bins[j])
    combined = bins_i + bins_j

    k = 0
    while(k < len(bins)):
        if(k == i):
            newBins.append(list(combined))
            k+=2
        else:
            newBins.append(list(bins[k]))
            k+=1

    for m in range(len(newBins)):
        probBins.append(len(newBins[m]) / n)
    
    return (newBins, probBins)

def printFinalBins(bins, probabilites):
    print("________________________________")
    print("After decretization: ")
    print("Final Bins State: ",  bins)
    print("Final Bin Probability Distribution: ", probabilites)

# Get the average information loss if we discretize all elements to a single bin
# This will be used as a reference point for determining when to stop "bin combination"
def getAverageInformationLoss(bins, probBins, n):
    print("___ Getting Average Information Loss to Determine Stopping Point ___")
    print("We first determine the average information loss by combing all bins into a single bin and tracking the information loss at each iteration")
    print("This will determine our stoping point for discretization if the info loss at a certain point exceeds the average\n")

    infoLoss = []
    binsCopy = copy.deepcopy(bins)
    probBinsCopy = copy.deepcopy(probBins)

    while(len(binsCopy) > 1):
        # determing (i, j) bins to merge
        (i, j) = getBinIndiciesToMerge(binsCopy)
        probFirstBinForMerge = probBinsCopy[i]
        probSecondBinForMerge = probBinsCopy[j]
        sumOfProbBins = probFirstBinForMerge + probSecondBinForMerge
        
        # get info loss and new bins after merging
        print("Information Loss:" )
        print("Sum of Bins " + str(i+1) + " and " + str(j+1) +" probabilites = " + str(probFirstBinForMerge) + " + " + str(probSecondBinForMerge) + " = " + str(sumOfProbBins))
        print("H_2 for bins "+ str(i+1) + " and " + str(j+1) + ": ")
        infoLoss_i = (probFirstBinForMerge + probSecondBinForMerge) * getEntropy([ probFirstBinForMerge / sumOfProbBins, probSecondBinForMerge/ sumOfProbBins ])
        print("Information loss = Sum of bins probabilities *  H_2 = " + str(infoLoss_i))
        print("\n")
        newBins, newProb = mergeBins(binsCopy, i, j, n)

        # add info loss to list, and set new bins/prob dist
        infoLoss.append(infoLoss_i)
        binsCopy = newBins
        probBinsCopy = newProb
    
    avg = sum(infoLoss) / len(infoLoss)
    print("Average info loss is: " + str(avg))
    print("_______________________________\n")
    return avg

def main(data):
    print("Data: ", data)

    n = len(data)
    bins, probBins = createBins(data)

    print("Initial Bins State: ", bins)
    print("Probability distribution of bins: ", probBins)
    print("\n")
    
    c = 1
    # H_n
    entropy = None
    averageInfoLoss = getAverageInformationLoss(bins, probBins, n)

    print("___ Discretization Process ___")
    while(len(bins) > 1):
        # calculate entropy
        print("Entropy with current bins:" )
        entropy = getEntropy(probBins)
        print("\n")

        c += 1
        
        # keep track of probabilty of bins to be merged to calulcate info loss
        (i, j) = getBinIndiciesToMerge(bins)
        probFirstBinForMerge = probBins[i]
        probSecondBinForMerge = probBins[j]
        sumOfProbBins = probFirstBinForMerge + probSecondBinForMerge

        # merge two consecutive bins with the lesat difference in mean value
        newBins, newProb = mergeBins(bins, i, j, n)
        previousBinsState = list(bins)
        previousProbBins = list(probBins)
        bins = newBins
        probBins = newProb


        # calculate info loss
        print("Information Loss:" )
        print("Sum of Bins " + str(i+1) + " and " + str(j+1) +" probabilites = " + str(probFirstBinForMerge) + " + " + str(probSecondBinForMerge) + " = " + str(sumOfProbBins))
        print("H_2 for bins "+ str(i+1) + " and " + str(j+1) + ": ")
        infoLoss_i = (probFirstBinForMerge + probSecondBinForMerge) * getEntropy([ probFirstBinForMerge / sumOfProbBins, probSecondBinForMerge/ sumOfProbBins ])
        print("Information loss = Sum of bins probabilities *  H_2 = " + str(infoLoss_i))
        print("\n")

        
        # Check if info loss rate of change starts to accelerate
        # in subsequent iterations. If so, we stop
        if(infoLoss_i > averageInfoLoss):
            print("Info loss exceeds average info loss, so we still stop the iteration process\n")
            printFinalBins(bins,probBins)
            return

        print("----------------------------------")
        print("After " + str(c-1) + " iterations")
        print("Bins " + str(i + 1) + " and " +  str(j + 1) + " will be combined since they have smallest difference in means")
        print("Previous Bins State: ", previousBinsState)
        print("Next Bins State: ", bins)
        print("Previous Probability Distribution of bins: ", previousProbBins)
        print("Next Probability distribution of bins: ", probBins)
        print("----------------------------------\n")
    
    printFinalBins(bins,probBins)
    return

data=[
0.5593630989,
0.9247589927,
0.9286284722,
0.7990522484,
0.8204542483,
0.6647561547,
0.9765970904,
0.9745135246,
0.8049680872,
0.7426005405,
0.6156845434,


0.9758901663,
0.645977965,
#N/A
#N/A
0.9194581,
]


main(data)