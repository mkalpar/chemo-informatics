# chemo-informatics

Chem(o)informatics is 
concerned with the application of computational methods to tackle chemical 
problems. It is clear that chemistry is a scientific discipline that is largely 
built on experimental observations and data. The amount of data and 
information accumulated is, however, enormous, and the size of this mountain is increasing with increasing speed. The problem is, then, to extract knowledge from these data and this information, and use this knowledge to make predictions. This is where chem(o)informatics is useful.

This repo includes some solutions to problems faced in chemoinformatics.
- Morgan algorithm : The representation of a chemical structure by a CTAB is neither unambiguous nor unique. A structure with n atoms can be numbered in n! different manners. Deriving a unique way of numbering is called canonicalization. The most popular canonizer is the Morgan algorithm.
- Ring detection (WIP)
- Isomorphism detection (WIP)
- SDF structure file parser (WIP)
- Degree of structure calculator (WIP)
- Symmetry of structure detection (WIP)
- Genetic algorithm: One of the main problems in chemoinfo is conformational analysis which studies the difference a different conformation of a molecule causes to occur. Each conformer corresponds to a single point on the potential surface, the lowest energy would be a global minimum. Evidently, it's computationally very complex to figure all confirmations of a molecule. Genetic algos are stochastic methods used to find global minimum based on the mechanics of natural genetics to evolve solutions. In GA, solutions are expressed as sequences of values. Each sequence is called a “chromosome” or a “string”, and each parameter within the chromosome is defined as a “gene”. 
    This code only does one-point crossover to evolve genes and finds as much as the input that is entered.
