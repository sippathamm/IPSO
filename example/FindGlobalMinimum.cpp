//
// Created by Sippawit Thammawiset on 25/1/2024 AD.
//

#include <iostream>

#include "BenchmarkFunction.h"
#include "PSO.h"

double GetMean (const std::vector<double> &Sample);
double GetVariance (const std::vector<double> &Sample);

double ObjectiveFunction (const std::vector<double> &Position)
{
    // Define your objective function here

    return Benchmark::BenchmarkFunction(Benchmark::SPHERE, Position);
}

int main() {
    std::vector<double> LowerBound, UpperBound;

    int MaximumIteration, NPopulation, NVariable;
    double SocialCoefficient = 1.5f, CognitiveCoefficient = 1.5f;
    double VelocityFactor = 0.5f;

    Benchmark::BenchmarkCondition(Benchmark::SPHERE,
                                  LowerBound, UpperBound,
                                  MaximumIteration, NPopulation, NVariable);

    int NRun = 30;

    // Save the results
    double Maximum = -INFINITY;
    double Minimum = INFINITY;
    std::vector<double> Sample;

    for (int Run = 1; Run <= NRun; Run++)
    {
        std::cout << "Run:\t" << Run << std::endl;

        Optimizer::APSO PSO(LowerBound, UpperBound,
                            MaximumIteration, NPopulation, NVariable,
                            SocialCoefficient, CognitiveCoefficient,
                            VelocityFactor,
                            false);

        PSO.SetObjectiveFunction(ObjectiveFunction);

        if (PSO.Run())
        {
            auto GlobalBestPosition = PSO.GetGlobalBestPosition();

            std::cout << "Global Best Position:\t";
            for (const auto &i : GlobalBestPosition)
            {
                std::cout << i << "\t";
            }
            std::cout << std::endl;

            double GlobalBestCost = PSO.GetGlobalBestCost();

            std::cout << "Global Best Cost:\t" << GlobalBestCost << std::endl;

            Maximum = std::max(Maximum, GlobalBestCost);
            Minimum = std::min(Minimum, GlobalBestCost);

            Sample.push_back(GlobalBestCost);
        }
        else
        {
            break;
        }

        std::cout << "------------------------" << std::endl;
    }

    std::cout << "Benchmark Result" << std::endl;

    for (const auto &i : Sample)
    {
        std::cout << i << "\t";
    }
    std::cout << std::endl;

    double Mean = GetMean(Sample);
    double Variance = GetVariance(Sample);
    double SD = sqrt(Variance);

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
        Variance += pow(i - Mean, 2);
    }

    return Variance / ((double)Sample.size() - 1);
}