"""
Written by Merve Kalpar
For constructing a genetic algorithm 
Date: June 28, 2021
"""



#!/usr/bin/env python3
import numpy as np
import sys

"""
Construct a random distance matrix between cities / genes
"""
def distance_matrix(num_of_cities):
    b = np.random.randint(low=0, high=100, size=(num_of_cities,num_of_cities))
    b_symm = (b + b.T)/2
    np.fill_diagonal(b_symm, 0)
    return b_symm
    
"""
Construct a random initial population (e.g., 1000 individuals) 
(Random solutions including list of cities together with 
their corresponding total distances which is the entity we want 
to minimize)
"""
def init_population(population_size,num_of_cities):
    #each row of population matrix represents a sequence of gene of an individual, which has 4 genes, or 4 cities to travel to
    population_matrix = np.zeros((population_size,num_of_cities),dtype=int)
    city_matrix = [city for city in range(num_of_cities)]
    for i in range(population_size):
        np.random.shuffle(city_matrix)
        population_matrix[i,:]=city_matrix
    return population_matrix    
    
"""
Create unique pairs until no pair repeat itself
so that no nonsensical matches where the pair is one individual occurs
"""
def pair_population(population_size):  
    pool = np.arange(population_size)
    new_pool = pool.copy()
    while True:
        np.random.shuffle(new_pool)
        for a,b in zip(pool,new_pool):
            if a==b:
                break
        else:
            return new_pool    

"""
Create offspring with single crossover operator
"""
def crossover(population_matrix, genes, parent1, parent2, crossover_point):
    genes = genes[:crossover_point]
    parent1 = population_matrix[parent1,:]
    parent2 = population_matrix[parent2,:]
    tmp = parent1[np.where(np.isin(parent1,genes,invert=True))]
    parent1[np.where(np.isin(parent1,genes,invert=True))] = parent2[np.where(np.isin(parent2,genes,invert=True))]
    parent2[np.where(np.isin(parent2,genes,invert=True))] = tmp
    return parent1,parent2

"""
Apply mutation
"""
def swapgene_mutation(chromosome,GENE_SIZE):
    # Mutation changes a single gene in each offspring randomly.
    swap_gene = np.random.randint(0,GENE_SIZE,2)
    tmp = chromosome[swap_gene[0]]
    chromosome[swap_gene[0]] = chromosome[swap_gene[1]]
    chromosome[swap_gene[1]] = tmp
    return chromosome

"""
Fitness calculation of a chromosome
"""
def fitness_calculation(chromosome,distance_matrix,GENE_SIZE):
    fitness_score = 0
    for i, j in zip(range(GENE_SIZE - 1), range(1, GENE_SIZE)):
        fitness_score = fitness_score + distance_matrix[chromosome[i], chromosome[j]]
    
    # last element to the first element bc salesman has to return to initial city
    fitness_score = fitness_score + distance_matrix[chromosome[GENE_SIZE - 1],chromosome[0]]  
    return fitness_score
 
 
"""
GENETIC ALGORITHM ITERATIONS START HERE. 
"""  
def genetic_algo_iterations(POPULATION_SIZE, GENE_SIZE, CROSSOVER_POINT, NUM_OF_GENERATIONS):
    distance_mat = distance_matrix(GENE_SIZE)
    population_matrix = init_population(POPULATION_SIZE,GENE_SIZE) 
    genes = np.arange(0,GENE_SIZE)
    
    print("DISTANCE MATRIX:")
    for i in range(len(distance_mat)):
        print('{} {}'.format(i, distance_mat[i]))
        
    for generation in range(NUM_OF_GENERATIONS):
        print("GENERATION NO:",generation) 
        
        #get the current size of population matrix and create pairs 
        corresponding_parent = pair_population(population_matrix.shape[0])  
        for i,j in zip(np.arange(0,POPULATION_SIZE),corresponding_parent):
            #print(i,j)
            offspring1, offspring2 = crossover(population_matrix, genes, i, j, CROSSOVER_POINT)
            offspring1 = swapgene_mutation(offspring1,GENE_SIZE)
            offspring2 = swapgene_mutation(offspring2,GENE_SIZE)
            population_matrix = np.vstack([population_matrix,offspring1])
            population_matrix = np.vstack([population_matrix,offspring2])
        
        #fit score
        scores = []
        for i in population_matrix:
            score = fitness_calculation(i,distance_mat,GENE_SIZE)
            scores.append(score)
        scores = np.array(scores)
        #fitness_calculation(offspring1,distance_mat,GENE_SIZE)
        #fitness_calculation(offspring2,distance_mat,GENE_SIZE)
        
        print("New population size:",population_matrix.shape[0])
        print("{} th generation: shortest path of distance {} on individual no {}, with chromosome ".format(generation,np.sort(scores)[0]
        ,np.where(scores.argsort()[::1]==0)[0][0]))
    

                   
if len(sys.argv)<5 or len(sys.argv)>5:
    print("Arguments are wrong.\n" \
          "Correct usage: genetic.py population-size gene-size crossover-gene-point, num-of-generations(num-of-iterations) \n" \
          "Sample usage: genetic.py 100 10 4 5")
    exit(0)
else:
    pop_size  = int(sys.argv[1])
    gene_size  = int(sys.argv[2])  
    crossover_gene_point  = int(sys.argv[3]) 
    num_of_generations  = int(sys.argv[4]) 
    genetic_algo_iterations(pop_size,gene_size,crossover_gene_point,num_of_generations)
