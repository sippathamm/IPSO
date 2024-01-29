//
// Created by Sippawit Thammawiset on 29/1/2024 AD.
//

#ifndef BENCHMARK_FUNCTION_H
#define BENCHMARK_FUNCTION_H

enum
{
    SPHERE = 0,
    SCHWEFEL_S_2_22 = 1,
    SCHWEFEL_S_1_20 = 2,
    ROSENBROCK = 3,
    STEP = 4,
    QUARTIC_NOISE = 5,
    SCHWEFEL_S_2_26 = 6,
    RASTRIGIN = 7,
    ACKLEY = 8,
    GRIEWANK = 9
};

namespace Benchmark
{
    namespace Function
    {
        double Sphere (const std::vector<double> &Position)
        {
            double Sum = 0.0f;

            for (const double &i : Position)
            {
                Sum += i * i;
            }

            return Sum;
        }

        double Schwefel_s_2_22 (const std::vector<double> &Position)
        {
            double Term1 = 0.0f;
            double Term2 = 1.0f;

            for (const double &i : Position)
            {
                Term1 += std::abs(i);
                Term2 *= std::abs(i);
            }

            return Term1 + Term2;
        }

        double Schwefel_s_1_20 (const std::vector<double> &Position)
        {
            double Sum = 0.0;

            for (int i = 0; i < Position.size(); i++)
            {
                double InnerSum = 0.0f;

                for (int j = 0; j <= i; j++)
                {
                    InnerSum += Position[j];
                }

                Sum += InnerSum * InnerSum;
            }

            return Sum;
        }

        double Rosenbrock (const std::vector<double> &Position)
        {
            double Sum = 0.0f;

            for (int i = 0; i < Position.size() - 1; i++)
            {
                double Term1 = 100 * (Position[i + 1] - Position[i] * Position[i]) * (Position[i + 1] - Position[i] * Position[i]);
                double Term2 = (Position[i] - 1) * (Position[i] - 1);

                Sum += Term1 + Term2;
            }

            return Sum;
        }

        double Step (const std::vector<double> &Position)
        {
            double Sum = 0.0;

            for (const double &i : Position)
            {
                Sum += (i + 0.5f) * (i + 0.5f);
            }

            return Sum;
        }

        double QuarticNoise (const std::vector<double> &Position)
        {
            double Sum = 0.0f;

            for (int i = 0; i < Position.size(); i++)
            {
                Sum += i * Position[i] * Position[i] * Position[i] * Position[i];
            }

            return Sum + Optimizer::GenerateRandom(0.0f, 1.0f);
        }

        double Schwefel_s_2_26 (const std::vector<double> &Position)
        {
            double Sum = 0.0;

            for (const double &i : Position)
            {
                Sum += i * sinf(sqrtf(abs(i)));
            }

            return -Sum;
        }

        double Rastrigin (const std::vector<double> &Position)
        {
            double Sum = 0.0;

            for (const double &i : Position) {
                Sum += i * i - 10 * cosf(2.0f * M_PI * i) + 10.0f;
            }

            return Sum;
        }

        double Ackley (const std::vector<double> &Position)
        {
            int Dimension = Position.size();

            double SumSquare = 0.0;
            double SumCosine = 0.0;

            for (const double &i : Position)
            {
                SumSquare += i * i;
                SumCosine += cosf(2.0f * M_PI * i);
            }

            double Term1 = -20.0f * expf(-0.2f * sqrtf(SumSquare / Dimension));
            double Term2 = -expf(SumCosine / Dimension);

            return Term1 + Term2 + 20.0f + expf(1.0f);
        }

        double Griewank (const std::vector<double> &Position)
        {
            double SumSquare = 0.0;
            double ProductCosine = 1.0;

            for (int i = 0; i < Position.size(); ++i) {
                SumSquare += Position[i] * Position[i];
                ProductCosine *= cosf(Position[i] / sqrtf(i + 1));
            }

            return SumSquare / 4000 - ProductCosine + 1;
        }
    } // Benchmark::Function

    double BenchmarkFunction (int FUNCTION, const std::vector<double> &Position)
    {
        switch (FUNCTION)
        {
            case SPHERE:
                return Function::Sphere(Position);
            case SCHWEFEL_S_2_22:
                return Function::Schwefel_s_2_22(Position);
            case SCHWEFEL_S_1_20:
                return Function::Schwefel_s_1_20(Position);
            case ROSENBROCK:
                return Function::Rosenbrock(Position);
            case STEP:
                return Function::Step(Position);
            case QUARTIC_NOISE:
                return Function::QuarticNoise(Position);
            case SCHWEFEL_S_2_26:
                return Function::Schwefel_s_2_26(Position);
            case RASTRIGIN:
                return Function::Rastrigin(Position);
            case ACKLEY:
                return Function::Ackley(Position);
            case GRIEWANK:
                return Function::Griewank(Position);

            default:
                return Function::Sphere(Position);
        }

        return -1;
    }
}

namespace Benchmark
{
    namespace Condition
    {
        void Sphere(std::vector<double> &LowerBound,
                    std::vector<double> &UpperBound,
                    int &MaximumIteration, int &NPopulation, int &NVariable,
                    double &SocialCoefficient, double &CognitiveCoefficient,
                    double &VelocityFactor) {
            MaximumIteration = 1000;
            NPopulation = 50;

            NVariable = 30;

            LowerBound = std::vector<double>(NVariable, -100);
            UpperBound = std::vector<double>(NVariable, 100);

            SocialCoefficient = 1.5f;
            CognitiveCoefficient = 1.5f;
            VelocityFactor = 0.5f;
        }

        void Schwefel_s_2_22(std::vector<double> &LowerBound,
                             std::vector<double> &UpperBound,
                             int &MaximumIteration, int &NPopulation, int &NVariable,
                             double &SocialCoefficient, double &CognitiveCoefficient,
                             double &VelocityFactor) {
            MaximumIteration = 1000;
            NPopulation = 50;

            NVariable = 30;

            LowerBound = std::vector<double>(NVariable, -10);
            UpperBound = std::vector<double>(NVariable, 10);

            SocialCoefficient = 1.5f;
            CognitiveCoefficient = 1.5f;
            VelocityFactor = 0.5f;
        }

        void Schwefel_s_1_20(std::vector<double> &LowerBound,
                             std::vector<double> &UpperBound,
                             int &MaximumIteration, int &NPopulation, int &NVariable,
                             double &SocialCoefficient, double &CognitiveCoefficient,
                             double &VelocityFactor) {
            MaximumIteration = 1000;
            NPopulation = 50;

            NVariable = 30;

            LowerBound = std::vector<double>(NVariable, -100);
            UpperBound = std::vector<double>(NVariable, 100);

            SocialCoefficient = 1.5f;
            CognitiveCoefficient = 1.5f;
            VelocityFactor = 0.5f;
        }

        void Rosenbrock(std::vector<double> &LowerBound,
                        std::vector<double> &UpperBound,
                        int &MaximumIteration, int &NPopulation, int &NVariable,
                        double &SocialCoefficient, double &CognitiveCoefficient,
                        double &VelocityFactor) {
            MaximumIteration = 1000;
            NPopulation = 50;

            NVariable = 30;

            LowerBound = std::vector<double>(NVariable, -30);
            UpperBound = std::vector<double>(NVariable, 30);

            SocialCoefficient = 1.5f;
            CognitiveCoefficient = 1.5f;
            VelocityFactor = 0.5f;
        }

        void Step(std::vector<double> &LowerBound,
                  std::vector<double> &UpperBound,
                  int &MaximumIteration, int &NPopulation, int &NVariable,
                  double &SocialCoefficient, double &CognitiveCoefficient,
                  double &VelocityFactor) {
            MaximumIteration = 1000;
            NPopulation = 50;

            NVariable = 30;

            LowerBound = std::vector<double>(NVariable, -100);
            UpperBound = std::vector<double>(NVariable, 100);

            SocialCoefficient = 1.5f;
            CognitiveCoefficient = 1.5f;
            VelocityFactor = 0.5f;
        }

        void QuarticNoise(std::vector<double> &LowerBound,
                          std::vector<double> &UpperBound,
                          int &MaximumIteration, int &NPopulation, int &NVariable,
                          double &SocialCoefficient, double &CognitiveCoefficient,
                          double &VelocityFactor) {
            MaximumIteration = 1000;
            NPopulation = 50;

            NVariable = 30;

            LowerBound = std::vector<double>(NVariable, -1.28f);
            UpperBound = std::vector<double>(NVariable, 1.28f);

            SocialCoefficient = 1.5f;
            CognitiveCoefficient = 1.5f;
            VelocityFactor = 0.5f;
        }

        void Schwefel_s_2_26(std::vector<double> &LowerBound,
                             std::vector<double> &UpperBound,
                             int &MaximumIteration, int &NPopulation, int &NVariable,
                             double &SocialCoefficient, double &CognitiveCoefficient,
                             double &VelocityFactor) {
            MaximumIteration = 1000;
            NPopulation = 50;

            NVariable = 30;

            LowerBound = std::vector<double>(NVariable, -500);
            UpperBound = std::vector<double>(NVariable, 500);

            SocialCoefficient = 1.5f;
            CognitiveCoefficient = 1.5f;
            VelocityFactor = 0.5f;
        }

        void Rastrigin(std::vector<double> &LowerBound,
                       std::vector<double> &UpperBound,
                       int &MaximumIteration, int &NPopulation, int &NVariable,
                       double &SocialCoefficient, double &CognitiveCoefficient,
                       double &VelocityFactor) {
            MaximumIteration = 1000;
            NPopulation = 50;

            NVariable = 30;

            LowerBound = std::vector<double>(NVariable, -5.12f);
            UpperBound = std::vector<double>(NVariable, 5.12f);

            SocialCoefficient = 1.5f;
            CognitiveCoefficient = 1.5f;
            VelocityFactor = 0.5f;
        }

        void Ackley(std::vector<double> &LowerBound,
                    std::vector<double> &UpperBound,
                    int &MaximumIteration, int &NPopulation, int &NVariable,
                    double &SocialCoefficient, double &CognitiveCoefficient,
                    double &VelocityFactor) {
            MaximumIteration = 1000;
            NPopulation = 50;

            NVariable = 30;

            LowerBound = std::vector<double>(NVariable, -32);
            UpperBound = std::vector<double>(NVariable, 32);

            SocialCoefficient = 1.5f;
            CognitiveCoefficient = 1.5f;
            VelocityFactor = 0.5f;
        }

        void Griewank(std::vector<double> &LowerBound,
                      std::vector<double> &UpperBound,
                      int &MaximumIteration, int &NPopulation, int &NVariable,
                      double &SocialCoefficient, double &CognitiveCoefficient,
                      double &VelocityFactor) {
            MaximumIteration = 1000;
            NPopulation = 50;

            NVariable = 30;

            LowerBound = std::vector<double>(NVariable, -600);
            UpperBound = std::vector<double>(NVariable, 600);

            SocialCoefficient = 1.5f;
            CognitiveCoefficient = 1.5f;
            VelocityFactor = 0.5f;
        }
    }

    void BenchmarkCondition(int FUNCTION,
                            std::vector<double> &LowerBound,
                            std::vector<double> &UpperBound,
                            int &MaximumIteration, int &NPopulation, int &NVariable,
                            double &SocialCoefficient, double &CognitiveCoefficient,
                            double &VelocityFactor)
    {
        switch (FUNCTION)
        {
            case SPHERE:
                return Condition::Sphere(LowerBound,
                                         UpperBound,
                                         MaximumIteration, NPopulation, NVariable,
                                         SocialCoefficient, CognitiveCoefficient,
                                         VelocityFactor);
            case SCHWEFEL_S_2_22:
                return Condition::Schwefel_s_2_22(LowerBound,
                                                  UpperBound,
                                                  MaximumIteration, NPopulation, NVariable,
                                                  SocialCoefficient, CognitiveCoefficient,
                                                  VelocityFactor);
            case SCHWEFEL_S_1_20:
                return Condition::Schwefel_s_1_20(LowerBound,
                                                  UpperBound,
                                                  MaximumIteration, NPopulation, NVariable,
                                                  SocialCoefficient, CognitiveCoefficient,
                                                  VelocityFactor);
            case ROSENBROCK:
                return Condition::Rosenbrock(LowerBound,
                                             UpperBound,
                                             MaximumIteration, NPopulation, NVariable,
                                             SocialCoefficient, CognitiveCoefficient,
                                             VelocityFactor);
            case STEP:
                return Condition::Step(LowerBound,
                                       UpperBound,
                                       MaximumIteration, NPopulation, NVariable,
                                       SocialCoefficient, CognitiveCoefficient,
                                       VelocityFactor);
            case QUARTIC_NOISE:
                return Condition::QuarticNoise(LowerBound,
                                               UpperBound,
                                               MaximumIteration, NPopulation, NVariable,
                                               SocialCoefficient, CognitiveCoefficient,
                                               VelocityFactor);
            case SCHWEFEL_S_2_26:
                return Condition::Schwefel_s_2_26(LowerBound,
                                                  UpperBound,
                                                  MaximumIteration, NPopulation, NVariable,
                                                  SocialCoefficient, CognitiveCoefficient,
                                                  VelocityFactor);
            case RASTRIGIN:
                return Condition::Rastrigin(LowerBound,
                                            UpperBound,
                                            MaximumIteration, NPopulation, NVariable,
                                            SocialCoefficient, CognitiveCoefficient,
                                            VelocityFactor);
            case ACKLEY:
                return Condition::Ackley(LowerBound,
                                         UpperBound,
                                         MaximumIteration, NPopulation, NVariable,
                                         SocialCoefficient, CognitiveCoefficient,
                                         VelocityFactor);
            case GRIEWANK:
                return Condition::Griewank(LowerBound,
                                           UpperBound,
                                           MaximumIteration, NPopulation, NVariable,
                                           SocialCoefficient, CognitiveCoefficient,
                                           VelocityFactor);
            default:
                return Condition::Sphere(LowerBound,
                                         UpperBound,
                                         MaximumIteration, NPopulation, NVariable,
                                         SocialCoefficient, CognitiveCoefficient,
                                         VelocityFactor);
        }
    }
} // Benchmark::Condition

#endif //BENCHMARK_FUNCTION_H