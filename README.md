# Improved Particle Swarm Optimization (IPSO) algorithm

## Overview

The IPSO algorithm is a variant of the Particle Swarm Optimization (PSO) algorithm, which is a population-based stochastic optimization technique inspired by the social behavior of bird flocking or fish schooling. IPSO improves upon traditional PSO by incorporating additional mechanisms for better exploration and exploitation of the search space.

## Features

- Implementation of the IPSO algorithm in C++.
- Easy-to-use interface for optimizing user-defined objective functions.
- Customizable parameters for fine-tuning the optimization process.
- Example code demonstrating the usage of the IPSO algorithm for solving optimization problems.

## Getting Started

Follow these steps to get started with the IPSO algorithm:

1. **Clone the Repository**: Clone this repository to your local machine using the following command:

    ```bash
    git clone https://github.com/sippathamm/IPSO.git
    ```

2. **Build the Project**: Navigate to the root directory and create a build directory. Then, run CMake to configure the project and generate build files:

    ```bash
    cd IPSO
    mkdir build
    cmake -S . -B build
    ```

3. **Compile the Code**: Once the build files are generated, build the IPSO executable by running:

    ```bash
    cmake --build build --target IPSO
    ```

4. **Run the Executable**: Execute the IPSO algorithm by navigating to the build directory and running the generated executable:

    ```bash
    cd build
    ./IPSO
    ```
   
## Usage

Replace the `ObjectiveFunction` implementation with your own objective function logic. This function should take a `std::vector<double>` representing the position of the particle in the search space and return the value of the objective function at that position.

Do not forget to call `IPSO.SetObjectiveFunction(ObjectiveFunction)` before calling `IPSO.Run()`

## Example

This example demonstrates the usage of the IPSO algorithm to find the global minimum with the following parameters:
- Maximum iteration: 1000
- Number of population: 50
- Number of dimensions: 30
- Lowerbound: -100
- Upperbound: 100
- Social coefficient (c1): 2.0
- Cognitive coefficient (c2): 1.3

The objective function used in this example is a Sphere function.

```cpp
#include <iostream>

#include "IPSO.h"

// Define your objective function here
double ObjectiveFunction (const std::vector<double> &Position) 
{
    // This function should return the value of the objective function at the given position
    
    double Sum = 0.0;

    for (const double &i : Position)
    {
        Sum += i * i;
    }

    return Sum;
}

int main () 
{
    // Initialize parameters
    int MaximumIteration = 1000;
    int NPopulation = 50;
    int NVariable = 30;
    std::vector<double> LowerBound = std::vector<double> (NVariable, -100);
    std::vector<double> UpperBound = std::vector<double> (NVariable, 100);
    double SocialCoefficient = 2.0; // c1
    double CognitiveCoefficient = 1.3; // c2
    double VelocityFactor = 0.5; // Factor for limiting velocity update
    int VelocityConfinement = MTH::IPSO::VELOCITY_CONFINEMENT::HYPERBOLIC;

    // Initialize IPSO algorithm
    MTH::IPSO::AIPSO<double> IPSO(LowerBound, UpperBound,
                                  MaximumIteration, NPopulation, NVariable,
                                  SocialCoefficient, CognitiveCoefficient,
                                  VelocityFactor,
                                  VelocityConfinement,
                                  true);
    
    // Set objective function for the algorithm
    IPSO.SetObjectiveFunction(ObjectiveFunction); 
    
    if (IPSO.Run()) // If the algorithm runs successfully
    {
        auto GlobalBestPosition = IPSO.GetGlobalBestPosition();

        std::cout << "Global Best Position:\t";
        std::for_each(GlobalBestPosition.begin(), GlobalBestPosition.end(), [](const auto &i) { std::cout << i << "\t"; });
        std::cout << std::endl;

        double GlobalBestCost = IPSO.GetGlobalBestCost();
        std::cout << "Global Best Cost:\t" << GlobalBestCost << std::endl;
    }
    else // If the algorithm fails to run
    {
        break; // Do nothing
    }
        
    return 0;
}
```

## Feedback and Bugs Report

If you have any feedback, suggestions, or encounter bugs while using the IPSO algorithm, please feel free to open an issue on the [GitHub repository](https://github.com/sippathamm/IPSO/issues).

## Author

This repository is maintained by Sippawit Thammawiset. You can contact the author at sippawit.t@kkumail.com
