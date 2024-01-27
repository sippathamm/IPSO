#include <iostream>

#include "PSO.h"

double GetMean (const std::vector<double> &Sample);
double GetVariance (const std::vector<double> &Sample);

double FitnessFunction (const std::vector<double> &Position)
{
    // Define your fitness function here

    double X1 = Position[0];
    double X2 = Position[1];

//    return X1 * X2;
//    return X1 * X1 + X2 * X2;
//    return pow(X1 - 3.14f, 2) + pow(X2 - 2.72f, 2) + sin(3 * X1 + 1.41f) + sin(4 * X2 - 1.73f);

//    return 2.0f * X1 * X1 - 1.05f * pow(X1, 4) + pow(X1, 6) / 6 + X1 * X2 + X2 * X2;

    double Term1 = 100.0f * sqrtf(fabs(X2 - (0.01f * X1 * X1)));
    double Term2 = 0.01f * fabs(X1 + 10.0f);

    return Term1 + Term2;
}

int main() {
    int NRun = 30;

    int MaxIteration = 1000;
    int NPopulation = 50;
    int NVariable = 2;

    std::vector<double> LowerBound = {-15, -3};
    std::vector<double> UpperBound = {-5, 3};

    double InertialWeight = 0.9f;
    double SocialCoefficient = 0.2f;
    double CognitiveCoefficient = 0.3f;
    double VelocityFactor = 0.5f;

    double Maximum = -INFINITY;
    double Minimum = INFINITY;
    std::vector<double> Sample;

    for (int Run = 0; Run < NRun; Run++)
    {
        std::cout << "Run:\t" << Run << std::endl;

        Optimizer::APSO PSO(LowerBound, UpperBound,
                            MaxIteration, NPopulation, NVariable,
                            InertialWeight, SocialCoefficient, CognitiveCoefficient,
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

    std::cout << "Result" << std::endl;

    std::cout << "Maximum:\t" << Maximum << std::endl;
    std::cout << "Minimum:\t" << Minimum << std::endl;
    std::cout << "Mean:\t" << GetMean(Sample) << std::endl;
    std::cout << "Variance:\t" << GetVariance(Sample) << std::endl;

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