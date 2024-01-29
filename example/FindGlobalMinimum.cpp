//
// Created by Sippawit Thammawiset on 25/1/2024 AD.
//

/* TODO:    - Add a convenience way to test benchmark function
 */

#include <iostream>

#include "PSO.h"
#include "BenchmarkFunction.h"

double GetMean (const std::vector<double> &Sample);
double GetVariance (const std::vector<double> &Sample);

double FitnessFunction (const std::vector<double> &Position)
{
    // Define your fitness function here

    return Benchmark::BenchmarkFunction(SPHERE, Position);
}

int main() {
    std::vector<double> LowerBound, UpperBound;

    int MaximumIteration, NPopulation, NVariable;
    double SocialCoefficient, CognitiveCoefficient;
    double VelocityFactor;

    Benchmark::BenchmarkCondition(SPHERE,
                                  LowerBound, UpperBound,
                                  MaximumIteration, NPopulation, NVariable,
                                  SocialCoefficient, CognitiveCoefficient,
                                  VelocityFactor);

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

    for (const double &i : Sample)
    {
        std::cout << i << std::endl;
    }

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