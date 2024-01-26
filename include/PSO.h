//
// Created by Sippawit Thammawiset on 25/1/2024 AD.
//

#ifndef PSO_H
#define PSO_H

#include <iostream>
#include <vector>
#include <random>
#include <cmath>

#define CLAMP(X, MIN, MAX)        std::max(MIN, std::min(MAX, X))

namespace Optimizer
{
    double GenerateRandom (double LowerBound = 0.0f, double UpperBound = 1.0f)
    {
        std::random_device Engine;
        std::uniform_real_distribution<double> RandomDistribution(0.0f, 1.0f);
        return LowerBound + RandomDistribution(Engine) * (UpperBound - LowerBound);
    }

    typedef struct AParticle
    {
        AParticle () : BestFitnessValue(INFINITY) {}

        std::vector<double> Position;
        std::vector<double> Velocity;

        std::vector<double> BestPosition;
        double BestFitnessValue;
    } AParticle;

    class APSO
    {
    public:
        APSO (const std::vector<double>& LowerBound, const std::vector<double>& UpperBound,
              int MaxIteration, int NPopulation, int NVariable,
              double Phi = 2.05f, double VelocityFactor = 0.5f) :
              LowerBound_(LowerBound),
              UpperBound_(UpperBound),
              MaxIteration_(MaxIteration),
              NPopulation_(NPopulation),
              NVariable_(NVariable),
              Phi_(Phi),
              VelocityFactor_(VelocityFactor)
        {

        }

        ~APSO () = default;

        static double FitnessFunction (const std::vector<double> &Position)
        {
            double X1 = Position[0];
            double X2 = Position[1];

//            return X1 * X2;
//            return X1 * X1 + X2 * X2;
//            return pow(X1 - 3.14f, 2) + pow(X2 - 2.72f, 2) + sin(3 * X1 + 1.41f) + sin(4 * X2 - 1.73f);

            double Term1 = 100.0f * sqrtf(fabs(X2 - (0.01f * X1 * X1)));
            double Term2 = 0.01f * fabs(X1 + 10.0f);

            return Term1 + Term2;
        }

        bool Run ()
        {
            // Initialize
            this->Population_ = std::vector<AParticle> (NPopulation_);

            this->MaximumVelocity_ = std::vector<double> (this->NVariable_);
            this->MinimumVelocity_ = std::vector<double> (this->NVariable_);

            this->W_ = 1.0f / (this->Phi_ - 1.0f + sqrt(this->Phi_ * this->Phi_ - 2.0f * this->Phi_));
            this->CMax_ = this->W_ * this->Phi_;

            for (int VariableIndex = 0; VariableIndex < this->NVariable_; VariableIndex++)
            {
                this->MaximumVelocity_[VariableIndex] = this->VelocityFactor_ * (this->UpperBound_[VariableIndex] - this->LowerBound_[VariableIndex]);
                this->MinimumVelocity_[VariableIndex] = -this->MaximumVelocity_[VariableIndex];
            }

            for (int PopulationIndex = 0; PopulationIndex < this->NPopulation_; PopulationIndex++)
            {
                auto *CurrentPopulation = &this->Population_[PopulationIndex];

                std::vector<double> Position(this->NVariable_);
                std::vector<double> Velocity(this->NVariable_);

                for (int VariableIndex = 0; VariableIndex < this->NVariable_; VariableIndex++)
                {
                    double RandomPosition = GenerateRandom(this->LowerBound_[VariableIndex], this->UpperBound_[VariableIndex]);
                    double RandomVelocity = (this->LowerBound_[VariableIndex] - RandomPosition) +
                                            GenerateRandom(0.0f, 1.0f) * (this->UpperBound_[VariableIndex] - this->LowerBound_[VariableIndex]);

                    RandomVelocity = CLAMP(RandomVelocity, this->MinimumVelocity_[VariableIndex], this->MaximumVelocity_[VariableIndex]);

                    Position[VariableIndex] = RandomPosition;
                    Velocity[VariableIndex] = RandomVelocity;
                }

                CurrentPopulation->Position = Position;
                CurrentPopulation->Velocity = Velocity;

                double FitnessValue = FitnessFunction(CurrentPopulation->Position);

                CurrentPopulation->BestPosition = CurrentPopulation->Position;
                CurrentPopulation->BestFitnessValue = FitnessValue;

                if (FitnessValue < this->GlobalBestFitnessValue_)
                {
                    this->GlobalBestPosition_ = CurrentPopulation->Position;
                    this->GlobalBestFitnessValue_ = FitnessValue;
                    this->GlobalBestIndex_ = PopulationIndex;
                }
            }
            
            // Optimize
            for (int Iteration = 0; Iteration < this->MaxIteration_; Iteration++)
            {
                for (int PopulationIndex = 0; PopulationIndex < this->NPopulation_; PopulationIndex++)
                {
                    auto *CurrentPopulation = &this->Population_[PopulationIndex];

                    Optimize(CurrentPopulation, PopulationIndex);
                }

                std::cout << "[INFO] Iteration: " << Iteration << " >>> " << "Best fitness value: " << this->GlobalBestFitnessValue_ << std::endl;
            }

            std::cout << "[INFO] Completed." << std::endl;

            return true;
        }

        void Optimize (AParticle *CurrentPopulation, int PopulationIndex)
        {
            std::vector<double> UpdatedPosition(this->NVariable_);
            std::vector<double> UpdatedVelocity(this->NVariable_);

            double Radius = 0.0f;
            double Norm = 0.0f;
            std::vector<double> SphereCenter(this->NVariable_);
            std::vector<double> Random(this->NVariable_);

            // Hypersphere
            for (int VariableIndex = 0; VariableIndex < this->NVariable_; VariableIndex++)
            {
                double Gr = CurrentPopulation->Position[VariableIndex] + (1.0f / 3.0f) * this->CMax_ *
                            (CurrentPopulation->BestPosition[VariableIndex] +
                            this->GlobalBestPosition_[VariableIndex] - 2.0f * CurrentPopulation->Position[VariableIndex]) *
                            (PopulationIndex == this->GlobalBestIndex_ ? 0.75f : 1.0f);

                SphereCenter[VariableIndex] = Gr;

                Radius += powf(Gr - CurrentPopulation->Position[VariableIndex], 2);

                Random[VariableIndex] = GenerateRandom(0.0f, 1.0f);
                Norm += Random[VariableIndex] * Random[VariableIndex];
            }

            double MaximumRadius = sqrt(Radius);
            Norm = sqrt(Norm);

            for (int VariableIndex = 0; VariableIndex < this->NVariable_; VariableIndex++)
            {
                double XSphere = Random[VariableIndex] * (GenerateRandom(0.0f, MaximumRadius) / Norm);

                // Update Velocity
                double NewVelocity = this->W_ * CurrentPopulation->Velocity[VariableIndex] +
                                     SphereCenter[VariableIndex] + XSphere - CurrentPopulation->Position[VariableIndex];

                NewVelocity = CLAMP(NewVelocity, this->MinimumVelocity_[VariableIndex], this->MaximumVelocity_[VariableIndex]);

                // Temporary Update Position
                double NewPosition = CurrentPopulation->Position[VariableIndex] + NewVelocity;

                if (NewPosition < this->LowerBound_[VariableIndex] || NewPosition > this->UpperBound_[VariableIndex])
                {
                    NewVelocity = -GenerateRandom(0.0f, 1.0f) * NewVelocity;
                }

                // Update Position
                NewPosition = CurrentPopulation->Position[VariableIndex] + NewVelocity;

                NewPosition = CLAMP(NewPosition, this->LowerBound_[VariableIndex], this->UpperBound_[VariableIndex]);

                UpdatedPosition[VariableIndex] = NewPosition;
                UpdatedVelocity[VariableIndex] = NewVelocity;
            }

            // Evaluate Fitness Value
            double FitnessValue = FitnessFunction(CurrentPopulation->Position);

            // Update PBest
            if (FitnessValue < CurrentPopulation->BestFitnessValue)
            {
                CurrentPopulation->BestPosition = CurrentPopulation->Position;
                CurrentPopulation->BestFitnessValue = FitnessValue;
            }

            // Update GBest
            if (FitnessValue < this->GlobalBestFitnessValue_)
            {
                this->GlobalBestPosition_ = CurrentPopulation->Position;
                this->GlobalBestFitnessValue_ = FitnessValue;
                this->GlobalBestIndex_ = PopulationIndex;
            }

            CurrentPopulation->Position = UpdatedPosition;
            CurrentPopulation->Velocity = UpdatedVelocity;
        }

        std::vector<double> GetGlobalBestPosition ()
        {
            return this->GlobalBestPosition_;
        }

        double GetGlobalBestFitnessValue () const
        {
            return this->GlobalBestFitnessValue_;
        }

    private:
        std::vector<double> LowerBound_, UpperBound_;
        int MaxIteration_, NPopulation_, NVariable_;
        double Phi_, W_{}, CMax_{};
        double VelocityFactor_;

        std::vector<AParticle> Population_;
        std::vector<double> MinimumVelocity_, MaximumVelocity_;

        std::vector<double> GlobalBestPosition_;
        double GlobalBestFitnessValue_ = INFINITY;
        int GlobalBestIndex_{};
    };
}

#endif // PSO_H
