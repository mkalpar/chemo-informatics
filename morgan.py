#!/usr/bin/env python3
import numpy as np
import sys

def print_degrees_of_node(connTable, atomLabel):
    for i in range(len(connTable)):
        print('{} {}'.format(atomLabel[i], connTable[i]))
        
def degree_calculator(fileName):
    text = open(fileName)
    lines = text.readlines()[3:]
    structure = lines[0].split()
    atomCnt = int(structure[0])
    bondCnt = int(structure[1])
    bondBlc = lines[atomCnt+1:atomCnt+1+bondCnt]            #bond block
    adjacencyTab = np.zeros((atomCnt,atomCnt),dtype = int)  #connection table
    bondTab = np.zeros((atomCnt,atomCnt),dtype = int)       #bond table
    
    for line in bondBlc:
        i = int(line.split()[0])-1
        j = int(line.split()[1])-1
        k = int(line.split()[2])
        bondTab[i][j] = k
        bondTab[j][i] = k
        if(k!=0):
            adjacencyTab[i][j] = 1 
            adjacencyTab[j][i] = 1
    
    atomLabels = []
    for line in lines[1:atomCnt+1]:
        atomLabels.append(line.split()[3])
    atomLabels = np.asarray(atomLabels)

    return adjacencyTab, bondTab, atomLabels

def morgan_algo(fileName,iteration_cnt):
    adjacencytab, bondtab, atomlabels = degree_calculator(fileName)
    degrees = np.sum(adjacencytab!=0,axis=0) # each index of the array stands for atom in that index of the atomlabel
    
    neighboring_indexes = []
    for i in range(len(atomlabels)):
        neighboring_indexes.append(list(np.where(adjacencytab[i]>0)[0]))
    
    print("neighboring atoms for each atom")
    for i in range(len(atomlabels)):
        print('{} {}'.format(atomlabels[i], neighboring_indexes[i]))
          
    for iteration in range(iteration_cnt):
        neighboring_degrees=[[degrees[j] for j in i] for i in neighboring_indexes]
        degrees = [sum(i) for i in neighboring_degrees]
        print(len(np.unique(degrees)),"different values of degrees after iteration",iteration+1,":",np.unique(degrees))
    
    # turn into numpy array to apply argsort function to it
    degrees = np.array(degrees)
    
    # morgan list will store indexes of new node numbers for each index at the atom label list
    # therefore, the node at index 0 in morgan, is the new node number of the 0th index in the atomlabel list
    morgan = []
    
    # start with finding the first node; the node with the highest degree
    morgan.append(np.where(degrees.argsort()[::-1]==0)[0].tolist())
    
    # fill morgan list with new items
    while len(morgan)<len(atomlabels):
        for i in morgan:
            neighboring_degrees = [] # will store the degrees of the neighboring nodes
            cp_index = [] # will store neighboring indexes for that one atom we are focusing on 
            for j in neighboring_indexes[i[0]]:
                if j not in morgan:
                    neighboring_degrees.append(degrees[j])
                    cp_index.append(j)
                    
            # indexes of atoms sorted by their degrees, 
            # so that the numbering can follow that order
            sorted_cp_index =[]
            for index in np.array(neighboring_degrees).argsort()[::-1]:
                sorted_cp_index.append(cp_index[index])
            
            # go through each neighbor of the atom we are focusing on
            for index in sorted_cp_index:
                # since it's already sorted, all it checks is whether that atom is not already in morgan list
                if index not in morgan:
                    morgan.append([index])
                    #neighboring_degrees.remove(sorted(neighboring_degrees,reverse=True)[0]
        print("FINAL MORGAN CANONICALIZATION:")

        original_index =0
        for i in range(len(atomlabels)):
            print('{}{} :{}'.format(atomlabels[i], original_index,morgan[i]))
            original_index +=1
          

                    
if len(sys.argv)<3 or len(sys.argv)>3:
    print("Arguments are wrong.\n" \
          "Correct usage: morgan.py file-path-of-mol-file iteration_count")
    exit(0)
else:
    molFile  = sys.argv[1]  
    iterCount  = int(sys.argv[2])  
    morgan_algo(molFile,iterCount)
