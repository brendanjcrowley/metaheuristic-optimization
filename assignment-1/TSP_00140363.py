

"""
Author: Brendan J Crowley
file: TSP_R00140363.py
Rename this file to TSP_x.py where x is your student number 
"""

import random
from Individual import *
import sys

myStudentNum = 140363 # Replace 12345 with your student number
random.seed(myStudentNum)

class BasicTSP:
    def __init__(self, _fName, _popSize, _mutationRate, _maxIterations):
        """
        Parameters and general variables
        """

        self.population     = []
        self.matingPool     = []
        self.best           = None
        self.popSize        = _popSize
        self.genSize        = None
        self.mutationRate   = _mutationRate
        self.maxIterations  = _maxIterations
        self.iteration      = 0
        self.fName          = _fName
        self.data           = {}

        self.readInstance()
        self.initPopulation()


    def readInstance(self):
        """
        Reading an instance from fName
        """
        file = open(self.fName, 'r')
        self.genSize = int(file.readline())
        self.data = {}
        for line in file:
            (id, x, y) = line.split()
            self.data[int(id)] = (int(x), int(y))
        file.close()

    def initPopulation(self):
        """
        Creating random individuals in the population
        """
        for i in range(0, self.popSize):
            individual = Individual(self.genSize, self.data)
            individual.computeFitness()
            self.population.append(individual)

        self.best = self.population[0].copy()
        for ind_i in self.population:
            if self.best.getFitness() > ind_i.getFitness():
                self.best = ind_i.copy()
        print ("Best initial sol: ",self.best.getFitness())

    def updateBest(self, candidate):
        if self.best == None or candidate.getFitness() < self.best.getFitness():
            self.best = candidate.copy()
            print ("iteration: ",self.iteration, "best: ",self.best.getFitness())

    def randomSelection(self):
        """
        Random (uniform) selection of two individuals
        """
        indA = self.matingPool[ random.randint(0, self.popSize-1) ]
        indB = self.matingPool[ random.randint(0, self.popSize-1) ]
        return [indA, indB]

    def stochasticUniversalSampling(self):
        """
        Your stochastic universal sampling Selection Implementation
        """
        pass

    def uniformCrossover(self, indA, indB):
        pass

    def pmxCrossover(self, indA, indB):
        """
        Your PMX Crossover Implementation
        """
        pass
    
    def reciprocalExchangeMutation(self, ind):
        # Generate random integer (index of first gene to swap)
        gene1 = random.randint(0, len(ind.genSize))

        # Generate random integer (index of second gene to swap)
        gene2 = random.randint(0, len(ind.genSize))

        # If both point to the same index, change the second index
        while gene1 == gene2:
            gene2 = random.randint(0, len(ind.genSize))

        # Store the value of the gene at the first gene index
        gene_holder = ind.genes[gene1]

        # Set the value of the gene at the first gene index to the value of the gene at the second gene index
        ind.genes[gene1] = ind.genes[gene2]

        # Set the value of the gene at the first gene index to the value in the gene_holder object

        return ind

    def inversionMutation(self, ind):
        # Randomly generate an integer (index of one of the swapping points)
        gene1 = None
        gene2 = None
        # If both index points are None or the same, re-generate new index
        # points
        while gene1 == gene2:
            # Randomly generate an integer (index of one of the swapping points)
            gene1 = random.randint(0, len(ind.genSize))

            # Randomly generate an integer (index of one of the swapping points)
            gene2 = random.randint(0, len(ind.genSize))
        i = 0
        j = 0
        # If the first number generated is greater than the 2nd, set i = 2nd
        # index generated
        if gene1 > gene2:
            i = gene2
            j = gene1
        else:
            i = gene1
            j = gene2
            
        # Double-ended search
        # While the lower index is less than the higher index:
        while i < j:
            # Store the value of the child's ith gene
            gene_holder = ind.genes[i]

            # Replace the value of the child's ith gene with the value in the
            # child's jth gene
            ind.genes[i] = ind.genes[j]

            # Replace the value of the child's jth gene with the value in the
            # gene holder (the original value in the child's ith gene)
            ind.genes[j] = gene_holder

            #Increment i, decrement j
            i += 1
            j -= 1
        
        return ind

    def crossover(self, indA, indB):
        """
        Executes a 1 order crossover and returns a new individual
        """
        child = []
        tmp = {}

        indexA = random.randint(0, self.genSize-1)
        indexB = random.randint(0, self.genSize-1)

        for i in range(0, self.genSize):
            if i >= min(indexA, indexB) and i <= max(indexA, indexB):
                tmp[indA.genes[i]] = False
            else:
                tmp[indA.genes[i]] = True
        aux = []
        for i in range(0, self.genSize):
            if not tmp[indB.genes[i]]:
                child.append(indB.genes[i])
            else:
                aux.append(indB.genes[i])
        child += aux
        return child

    def mutation(self, ind):
        """
        Mutate an individual by swaping two cities with certain probability (i.e., mutation rate)
        """
        if random.random() > self.mutationRate:
            return
        indexA = random.randint(0, self.genSize-1)
        indexB = random.randint(0, self.genSize-1)

        tmp = ind.genes[indexA]
        ind.genes[indexA] = ind.genes[indexB]
        ind.genes[indexB] = tmp

        ind.computeFitness()
        self.updateBest(ind)

    def updateMatingPool(self):
        """
        Updating the mating pool before creating a new generation
        """
        self.matingPool = []
        for ind_i in self.population:
            self.matingPool.append( ind_i.copy() )

    def newGeneration(self):
        """
        Creating a new generation
        1. Selection
        2. Crossover
        3. Mutation
        """
        for i in range(0, len(self.population)):
            """
            Depending of your experiment you need to use the most suitable algorithms for:
            1. Select two candidates
            2. Apply Crossover
            3. Apply Mutation
            """

    def GAStep(self):
        """
        One step in the GA main algorithm
        1. Updating mating pool with current population
        2. Creating a new Generation
        """

        self.updateMatingPool()
        self.newGeneration()

    def search(self):
        """
        General search template.
        Iterates for a given number of steps
        """
        self.iteration = 0
        while self.iteration < self.maxIterations:
            self.GAStep()
            self.iteration += 1

        print ("Total iterations: ",self.iteration)
        print ("Best Solution: ", self.best.getFitness())

if len(sys.argv) < 2:
    print ("Error - Incorrect input")
    print ("Expecting python BasicTSP.py [instance] ")
    sys.exit(0)


problem_file = sys.argv[1]

ga = BasicTSP(sys.argv[1], 300, 0.1, 500)
ga.search()
