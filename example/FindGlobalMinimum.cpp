//
// Created by Sippawit Thammawiset on 25/1/2024 AD.
//

/* TODO:    - Add a convenience way to test benchmark function
 */

#include <iostream>

#include "PSO.h"
#include "BenchmarkFunction.h"

//#define SPHERE
//#define SCHWEFEL_S_2_22
//#define SCHWEFEL_S_1_20
//#define ROSENBROCK
#define STEP
//#define QUARTIC_NOISE

//#define SCHWEFEL_S_2_26
//#define RASTRIGIN
//#define ACKLEY
//#define GRIEWANK

double GetMean (const std::vector<double> &Sample);
double GetVariance (const std::vector<double> &Sample);

double FitnessFunction (const std::vector<double> &Position)
{
    // Define your fitness function here

#ifdef SPHERE

    return Benchmark::Function::Sphere(Position);

#endif // SPHERE

#ifdef SCHWEFEL_S_2_22

    return Benchmark::Function::Schwefel_s_2_22(Position);

#endif // SCHWEFEL_S_2_22

#ifdef SCHWEFEL_S_1_20

    return Benchmark::Function::Schwefel_s_1_20(Position);

#endif // SCHWEFEL_S_1_20

#ifdef ROSENBROCK

    return Benchmark::Function::Rosenbrock(Position);

#endif // ROSENBROCK

#ifdef STEP

    return Benchmark::Function::Step(Position);

#endif // STEP

#ifdef QUARTIC_NOISE
    
    return Benchmark::Function::QuarticNoise(Position);

#endif // QUARTIC_NOISE

#ifdef SCHWEFEL_S_2_26

    return Benchmark::Function::Schwefel_s_2_26(Position);

#endif // SCHWEFEL_S_2_26

#ifdef RASTRIGIN

    return Benchmark::Function::Rastrigin(Position);

#endif // RASTRIGIN

#ifdef ACKLEY

    return Benchmark::Function::Ackley(Position);

#endif // ACKLEY

#ifdef GRIEWANK

    return Benchmark::Function::Griewank(Position);

#endif // GRIEWANK
}

int main() {
    std::vector<double> LowerBound, UpperBound;

    int MaximumIteration, NPopulation, NVariable;
    double SocialCoefficient, CognitiveCoefficient;
    double VelocityFactor;

#ifdef SPHERE

    Benchmark::Condition::Sphere(LowerBound, UpperBound,
                                 MaximumIteration, NPopulation, NVariable,
                                 SocialCoefficient, CognitiveCoefficient,
                                 VelocityFactor);

#endif // SPHERE

#ifdef SCHWEFEL_S_2_22

    Benchmark::Condition::Schwefel_s_2_22(LowerBound, UpperBound,
                                          MaximumIteration, NPopulation, NVariable,
                                          SocialCoefficient, CognitiveCoefficient,
                                          VelocityFactor);

#endif // SCHWEFEL_S_2_22

#ifdef SCHWEFEL_S_1_20

    Benchmark::Condition::Schwefel_s_1_20(LowerBound, UpperBound,
                                          MaximumIteration, NPopulation, NVariable,
                                          SocialCoefficient, CognitiveCoefficient,
                                          VelocityFactor);

#endif // SCHWEFEL_S_1_20

#ifdef ROSENBROCK

    Benchmark::Condition::Rosenbrock(LowerBound, UpperBound,
                                     MaximumIteration, NPopulation, NVariable,
                                     SocialCoefficient, CognitiveCoefficient,
                                     VelocityFactor);

#endif // ROSENBROCK

#ifdef STEP

    Benchmark::Condition::Step(LowerBound, UpperBound,
                               MaximumIteration, NPopulation, NVariable,
                               SocialCoefficient, CognitiveCoefficient,
                               VelocityFactor);

#endif // STEP

#ifdef QUARTIC_NOISE

    Benchmark::Condition::QuarticNoise(LowerBound, UpperBound,
                                       MaximumIteration, NPopulation, NVariable,
                                       SocialCoefficient, CognitiveCoefficient,
                                       VelocityFactor);

#endif // QUARTIC_NOISE

#ifdef SCHWEFEL_S_2_26

    Benchmark::Condition::Schwefel_s_2_26(LowerBound, UpperBound,
                                          MaximumIteration, NPopulation, NVariable,
                                          SocialCoefficient, CognitiveCoefficient,
                                          VelocityFactor);

#endif // SCHWEFEL_S_2_26

#ifdef RASTRIGIN

    Benchmark::Condition::Rastrigin(LowerBound, UpperBound,
                                    MaximumIteration, NPopulation, NVariable,
                                    SocialCoefficient, CognitiveCoefficient,
                                    VelocityFactor);

#endif // RASTRIGIN

#ifdef ACKLEY

    Benchmark::Condition::Ackley(LowerBound, UpperBound,
                                 MaximumIteration, NPopulation, NVariable,
                                 SocialCoefficient, CognitiveCoefficient,
                                 VelocityFactor);

#endif // ACKLEY

#ifdef GRIEWANK

    Benchmark::Condition::Griewank(LowerBound, UpperBound,
                                   MaximumIteration, NPopulation, NVariable,
                                   SocialCoefficient, CognitiveCoefficient,
                                   VelocityFactor);

#endif // GRIEWANK

    int NRun = 30;

    // Store the result
    double Maximum = -INFINITY;
    double Minimum = INFINITY;
    std::vector<double> Sample;

    for (int Run = 0; Run < NRun; Run++)
    {
        std::cout << "Run:\t" << Run << std::endl;

        Optimizer::APSO PSO(LowerBound, UpperBound,
                            MaximumIteration, NPopulation, NVariable,
                            SocialCoefficient, CognitiveCoefficient,
                            VelocityFactor,
                            false);

        PSO.SetFitnessFunction(FitnessFunction);

        if (PSO.Run())
        {
            auto GlobalBestPosition = PSO.GetGlobalBestPosition();

            std::cout << "Global Best Position:\t";
            for (const auto &i : GlobalBestPosition)
            {
                std::cout << i << "\t";
            }
            std::cout << std::endl;

            double GlobalBestFitnessValue = PSO.GetGlobalBestFitnessValue();

            std::cout << "Global Best Fitness Value:\t" << GlobalBestFitnessValue << std::endl;

            Maximum = std::max(Maximum, GlobalBestFitnessValue);
            Minimum = std::min(Minimum, GlobalBestFitnessValue);

            Sample.push_back(GlobalBestFitnessValue);
        }
        else
        {
            break;
        }

        std::cout << "------------------------" << std::endl;
    }

    std::cout << "Benchmark Result" << std::endl;

    double Mean = GetMean(Sample);
    double Variance = GetVariance(Sample);
    double SD = sqrtf(Variance);

    std::cout << "Maximum:\t" << Maximum << std::endl;
    std::cout << "Minimum:\t" << Minimum << std::endl;
    std::cout << "Mean:\t" << Mean << std::endl;
    std::cout << "SD:\t" << SD << std::endl;

    return 0;
}

double GetMean (const std::vector<double> &Sample)
{
    double Sum = 0.0f;

    for (const auto &i : Sample)
    {
        Sum += i;
    }

    return Sum / (double)Sample.size();
}

double GetVariance (const std::vector<double> &Sample)
{
    double Mean = GetMean(Sample);
    double Variance = 0.0f;

    for (const auto &i : Sample)
    {
        Variance += powf(i - Mean, 2);
    }

    return Variance / ((double)Sample.size() - 1);
}