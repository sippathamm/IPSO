#include <iostream>

#include "PSO.h"

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

int main() {
    std::vector<double> LowerBound = {-15, -3};
    std::vector<double> UpperBound = {-5, 3};

    int NRun = 30;

    int MaxIteration = 5000;
    int NPopulation = 50;
    int NVariable = 2;

    double Phi = 2.05f;
    double VelocityFactor = 0.5f;

    double Maximum = -INFINITY;
    double Minimum = INFINITY;
    std::vector<double> Sample;

    for (int Run = 0; Run < NRun; Run++)
    {
        std::cout << "Run: " << Run << std::endl;

        Optimizer::APSO PSO(LowerBound, UpperBound,
                            MaxIteration, NPopulation, NVariable,
                            Phi, VelocityFactor);

        if (PSO.Run())
        {
            auto GlobalBestPosition = PSO.GetGlobalBestPosition();

            std::cout << "Global Best Position: ";
            for (const auto &i : GlobalBestPosition)
            {
                std::cout << i << " ";
            }
            std::cout << std::endl;

            double GlobalBestFitnessValue = PSO.GetGlobalBestFitnessValue();

            std::cout << "Global Best Fitness Value: " << GlobalBestFitnessValue << std::endl;

            Maximum = std::max(Maximum, GlobalBestFitnessValue);
            Minimum = std::min(Minimum, GlobalBestFitnessValue);

            Sample.push_back(GlobalBestFitnessValue);
        }

        std::cout << "------------------------" << std::endl;
    }

    std::cout << "Result" << std::endl;

    std::cout << "Maximum: " << Maximum << std::endl;
    std::cout << "Minimum: " << Minimum << std::endl;
    std::cout << "Mean: " << GetMean(Sample) << std::endl;
    std::cout << "Variance: " << GetVariance(Sample) << std::endl;

    return 0;
}
