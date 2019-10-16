

"""
Author: Brendan J Crowley
file: TSP_R00140363.py
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
        """Executes a uniform crossover and returns a new individual
        :param ind1: The first parent (or individual)
        :param ind2: The second parent (or individual)
        :returns: A new individual"""

        child = Individual(self.genSize, self.data)

        # Select random parent to inherit initial genes from
        parent = random.randint(0,1)
        parent1 = indA
        parent2 = indB
        if parent == 0:
            parent1 = indB
            parent2 = indA

        # Keep track of values already checked
        values_crossed = []
        
        # For each gene in the selected parent
        for i in range(len(parent1.genes)):
            
            # Randomly decide if the current gene is being kept
            cross = random.randint(0,1)

            # If the gene is selected:
            if cross == 1:

                # Set the child's gene as the same as the parent, and add
                # the value to the tracker
                child.genes[i] = parent1.genes[i]
                values_crossed += [parent1.genes[i]]
                
            else:
                # Otherwise, set the current gene as none
                child.genes[i] = None
        
        j = 0
        i = 0
        # For each gene in the other parent
        while i < len(parent2.genes):

            # If the current gene has a value, continue
            if child.genes[i] != None:
                i += 1

            # Otherwise
            elif child.genes[i] == None:

                # If the current value hasn't been selected previously
                over = parent2.genes[j]
                if over not in values_crossed:

                # Add the value to the checked values, and update the child
                # with the new gene
                    values_crossed += [over]
                    child.genes[i] = over
                    i += 1
                # Continue
                j += 1
        return child

    def pmxCrossover(self, indA, indB):

        child = Individual(self.genSize, self.data)
        
        parent1 = indA
        parent2 = indB

        # Keep track of values already checked, and their map
        indices_checked = []
        values_checked = []
        value_maps = dict()

        index1 = random.randint(0, len(parent1.genes)-1)
        index2 = random.randint(0, len(parent1.genes)-1)

        if index2 < index1:
            tmp = index1
            index1 = index2
            index2 = tmp

        for i in range(len(parent1.genes)):
            if i < index1:
                child.genes[i] = None
            elif i > index2:
                child.genes[i] = None
            else: #index is between 2 selected
                child.genes[i] = parent1.genes[i]
                indices_checked += [i]
                values_checked += [parent1.genes[i]]
                value_maps[parent1.genes[i]] = [parent2.genes[i]]
        print(child.genes)
                
        while len(indices_checked) > 0:
            
            current = indices_checked[0]
            value1 = parent1.genes[current]
            value2 = parent2.genes[current]
            
            if value1 == value2:
                indices_checked.remove(current)
                
            else:
                
                cur_val2 = parent2.genes.index(value1)
                mapped = value_maps[value1]
                
                while True:
                    
                    indices_checked.remove(current)
                    
                    if mapped not in list(value_maps.keys()):
                        
                        child.genes[cur_val2] = value2
                        values_checked += [value2]
                        break
                    
                    mapped = value_maps[mapped]
                    
        j = 0
        i = 0
        # For each gene in the other parent
        while i < len(parent2.genes):

            # If the current gene has a value, continue
            if child.genes[i] != None:
                i += 1

            # Otherwise
            elif child.genes[i] == None:

                # If the current value hasn't been selected previously
                over = parent2.genes[j]
                if over not in values_checked:

                # Add the value to the checked values, and update the child
                # with the new gene
                    values_checked += [over]
                    child.genes[i] = over
                    i += 1
                # Continue
                j += 1
        return child
    
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
            gene1 = random.randint(0, ind.genSize)

            # Randomly generate an integer (index of one of the swapping points)
            gene2 = random.randint(0, ind.genSize)
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

"""if len(sys.argv) < 2:
    print ("Error - Incorrect input")
    print ("Expecting python BasicTSP.py [instance] ")
    sys.exit(0)


problem_file = sys.argv[1]"""

ga = BasicTSP("inst-13.tsp", 100, 0.1, 500)
parent1 = ga.population[2]
parent2 = ga.population[1]
r = ga.pmxCrossover(parent1, parent2)
