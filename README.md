# Poisson-Channel-Capacity-Bounds
MATLAB codes for the ISIT 2023 paper "Capacity Bounds for the Poisson-Repeat Channel"

The codes are commented, with explanations about all employed input, output and intermediate variables.

The codes are:
- Bound_Comp (the main file to be executed): outputs the proposed upper bound on the Poisson-Repeat Channel versus the provided range of Poisson parameter \lambda

- Transition_Matrix_RepCh_allRv3_par: generates the transition probability matrix for a given Poisson-Repeat Channel

- TM_RepCh_allR_diffLambda: updates the transition probability matrix for a new \lambda value without calculating it from beginning (to reduce the computational cost)

- BAA (and its parfor version BAA_par): performs the Blahut-Arimoto Algorithm (BAA) for numerical computation of DMC capacity
