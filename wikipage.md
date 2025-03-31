# Random-Key Optimizer (RKO)

The Random-Key Optimizer (RKO) is a metaheuristic optimization algorithm that employs a continuous random-key representation to encode solutions. It is designed to solve combinatorial and continuous optimization problems by leveraging evolutionary and trajectory algorithm principles.

## Overview
RKO is based on the concept of random-key encoding, where solutions are represented as vectors of real numbers within a predefined range (typically [0,1]). These vectors are then decoded into feasible solutions using a problem-specific decoding function. The method is inspired by the Biased Random-Key Genetic Algorithm (BRKGA) but incorporates modifications to enhance efficiency and adaptability.

## History
The RKO method was introduced as a general-purpose optimization framework for solving complex decision-making problems. It builds upon the foundation of random-key genetic algorithms (RKGA), first introduced in the late 1990s, and extends their applicability by incorporating adaptive strategies and heuristic refinements.

## Algorithm
The general structure of the RKO involves the following steps:

1. **Initialization**: Generate a population of random-key encoded solutions.
2. **Decoding**: Convert random-key vectors into feasible solutions using a problem-specific decoding function.
3. **Evaluation**: Assess the quality of each solution based on a predefined objective function.
4. **Selection and Variation**:
   - Apply selection mechanisms to retain high-quality solutions.
   - Implement crossover and mutation operations to explore the search space.
5. **Termination**: Repeat the process until a stopping criterion is met (e.g., maximum iterations or convergence threshold).

## Applications
RKO has been applied to a variety of combinatorial optimization problems, including:
- Scheduling and resource allocation
- Network design and logistics
- Multi-plant capacitated lot sizing problems (MPCLSP)
- Machine learning hyperparameter tuning

## Advantages
- **Scalability**: Efficiently handles large-scale optimization problems.
- **Robustness**: Adapts to different problem structures with minimal customization.
- **Flexibility**: Can be integrated with other heuristic and exact methods for hybrid approaches.

## References
- Londe, L. R., et al. (2024). "A Systematic Review on Biased Random-Key Genetic Algorithms." *Journal of Optimization Research*.
- Gonçalves, J. F., & Resende, M. G. C. (2011). "Biased Random-Key Genetic Algorithms for Combinatorial Optimization." *Mathematical Programming*.

## External Links
- [Official Repository or Implementation (if available)]
- [Related Research Papers]
