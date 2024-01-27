//
// Created by Sippawit Thammawiset on 25/1/2024 AD.
//

// TODO: Add documentation

#ifndef PSO_H
#define PSO_H

#include <iostream>
#include <vector>
#include <random>
#include <cmath>

#define CLAMP(X, MIN, MAX)        std::max(MIN, std::min(MAX, X))

enum
{
    FAILED = 0,
    OK = 1
};

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
              double InertialWeight, double SocialCoefficient, double CognitiveCoefficient,
              double VelocityFactor = 0.5f,
              bool Log = true) :
              LowerBound_(LowerBound),
              UpperBound_(UpperBound),
              MaxIteration_(MaxIteration),
              NPopulation_(NPopulation),
              NVariable_(NVariable),
              InertialWeight_(InertialWeight),
              SocialCoefficient_(SocialCoefficient),
              CognitiveCoefficient_(CognitiveCoefficient),
              VelocityFactor_(VelocityFactor),
              Log_(Log)
        {

        }

        ~APSO () = default;

        void SetFitnessFunction (double (*UserFitnessFunction)(const std::vector<double> &Position))
        {
            this->FitnessFunction = UserFitnessFunction;
        }

        double UpdateVelocity (const AParticle *CurrentPopulation, int VariableIndex)
        {
            double NewVelocity =
                    this->InertialWeight_ * CurrentPopulation->Velocity[VariableIndex] +
                    this->SocialCoefficient_ * GenerateRandom(0.0f, 1.0f) * (CurrentPopulation->BestPosition[VariableIndex] - CurrentPopulation->Position[VariableIndex])+
                    this->CognitiveCoefficient_ * GenerateRandom(0.0f, 1.0f) * (this->GlobalBestPosition_[VariableIndex] - CurrentPopulation->Position[VariableIndex]);

            NewVelocity = CLAMP(NewVelocity, this->MinimumVelocity_[VariableIndex], this->MaximumVelocity_[VariableIndex]);

            return NewVelocity;
        }

        double UpdatePosition (const AParticle *CurrentPopulation, int VariableIndex)
        {
            double NewPosition = CurrentPopulation->Position[VariableIndex] + CurrentPopulation->Velocity[VariableIndex];

            NewPosition = CLAMP(NewPosition, this->LowerBound_[VariableIndex], this->UpperBound_[VariableIndex]);

            return NewPosition;
        }

        bool Run ()
        {
            if (FitnessFunction == nullptr)
            {
                std::cerr << "Fitness function is not defined. Please use SetFitnessFunction(UserFitnessFunction) before calling Run()." << std::endl;

                return FAILED;
            }

            if (this->LowerBound_.size() != NVariable_ || this->UpperBound_.size() != NVariable_)
            {
                std::cerr << "Size of lowerbound or upperbound does not match with number of variables." << std::endl;

                return FAILED;
            }

            // Initialize
            this->Population_ = std::vector<AParticle> (this->NPopulation_);

            this->MaximumVelocity_ = std::vector<double> (this->NVariable_);
            this->MinimumVelocity_ = std::vector<double> (this->NVariable_);

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
                    double RandomVelocity = GenerateRandom(0.0f, 1.0f);

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
                }
            }
            
            // Optimize
            for (int Iteration = 0; Iteration < this->MaxIteration_; Iteration++)
            {
                for (int PopulationIndex = 0; PopulationIndex < this->NPopulation_; PopulationIndex++)
                {
                    auto *CurrentPopulation = &this->Population_[PopulationIndex];

                    Optimize(CurrentPopulation);
                }

                if (this->Log_)
                {
                    std::cout << "[INFO] Iteration: " << Iteration << " >>> " << "Best fitness value: " << this->GlobalBestFitnessValue_ << std::endl;
                }
            }

            std::cout << "[INFO] Completed." << std::endl;

            return OK;
        }

        void Optimize (AParticle *CurrentPopulation)
        {
            for (int VariableIndex = 0; VariableIndex < this->NVariable_; VariableIndex++)
            {
                // Update Velocity
                double NewVelocity = UpdateVelocity(CurrentPopulation, VariableIndex);
                CurrentPopulation->Velocity[VariableIndex] = NewVelocity;
            }

            for (int VariableIndex = 0; VariableIndex < this->NVariable_; VariableIndex++)
            {
                // Update Position
                double NewPosition = UpdatePosition(CurrentPopulation, VariableIndex);
                CurrentPopulation->Position[VariableIndex] = NewPosition;
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
            }
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
        double InertialWeight_, SocialCoefficient_, CognitiveCoefficient_;
        double VelocityFactor_;

        double (*FitnessFunction)(const std::vector<double> &Position) = nullptr;

        std::vector<AParticle> Population_;
        std::vector<double> MaximumVelocity_, MinimumVelocity_;

        std::vector<double> GlobalBestPosition_;
        double GlobalBestFitnessValue_ = INFINITY;

        bool Log_;
    };
} // Optimizer

#endif // PSO_H
