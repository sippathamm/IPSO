//
// Created by Sippawit Thammawiset on 25/1/2024 AD.
//

/* TODO:    - Add documentation
 *          - Add more velocity confinements
 */

#ifndef PSO_H
#define PSO_H

#include <iostream>
#include <vector>
#include <random>
#include <cmath>

#define CLAMP(X, MIN, MAX)                      std::max(MIN, std::min(MAX, X))
#define IS_OUT_OF_BOUND(X, MIN, MAX)            X < MIN || X > MAX

namespace Optimizer
{
    enum
    {
        FAILED = 0,
        SUCCESS = 1
    };

    double GenerateRandom (double LowerBound = 0.0f, double UpperBound = 1.0f)
    {
        std::random_device Engine;
        std::uniform_real_distribution<double> RandomDistribution(0.0f, 1.0f);
        return LowerBound + RandomDistribution(Engine) * (UpperBound - LowerBound);
    }

    typedef struct AParticle
    {
        AParticle () : BestFitnessValue(INFINITY), FitnessValue(0.0f) {}

        std::vector<double> Position;
        std::vector<double> Velocity;
        double FitnessValue;

        std::vector<double> BestPosition;
        double BestFitnessValue;

        std::vector<double> Feedback;
    } AParticle;

    class APSO
    {
    public:
        APSO (const std::vector<double> &LowerBound, const std::vector<double> &UpperBound,
              int MaximumIteration, int NPopulation, int NVariable,
              double SocialCoefficient = 1.5f, double CognitiveCoefficient = 1.5f,
              double VelocityFactor = 0.5f,
              bool Log = true) :
              LowerBound_(LowerBound),
              UpperBound_(UpperBound),
              MaximumIteration_(MaximumIteration),
              NPopulation_(NPopulation),
              NVariable_(NVariable),
              SocialCoefficient_(SocialCoefficient),
              CognitiveCoefficient_(CognitiveCoefficient),
              VelocityFactor_(VelocityFactor),
              Log_(Log)
        {

        }

        ~APSO () = default;

        void SetFitnessFunction (double (*UserFitnessFunction)(const std::vector<double> &Position))
        {
            this->FitnessFunction_ = UserFitnessFunction;
        }

        double UpdateVelocity (const AParticle *CurrentPopulation, int VariableIndex)
        {
            double NewVelocity = this->InertialWeight_ * CurrentPopulation->Velocity[VariableIndex] +
                                 this->SocialCoefficient_ * GenerateRandom(0.0f, 1.0f) * (CurrentPopulation->BestPosition[VariableIndex] - CurrentPopulation->Position[VariableIndex])+
                                 this->CognitiveCoefficient_ * GenerateRandom(0.0f, 1.0f) * (this->GlobalBestPosition_[VariableIndex] - CurrentPopulation->Position[VariableIndex]) +
                                 CurrentPopulation->Feedback[VariableIndex];

            NewVelocity = CLAMP(NewVelocity, this->MinimumVelocity_[VariableIndex], this->MaximumVelocity_[VariableIndex]);

            return NewVelocity;
        }

        double UpdatePosition (AParticle *CurrentPopulation, int VariableIndex)
        {
            double TemporaryNewPosition = CurrentPopulation->Position[VariableIndex] + CurrentPopulation->Velocity[VariableIndex];

            if (IS_OUT_OF_BOUND(TemporaryNewPosition, this->LowerBound_[VariableIndex], this->UpperBound_[VariableIndex]))
            {
                double VelocityConfinement = -GenerateRandom(0.0f, 1.0f) * CurrentPopulation->Velocity[VariableIndex];

                CurrentPopulation->Velocity[VariableIndex] = VelocityConfinement;
            }

            double NewPosition = CurrentPopulation->Position[VariableIndex] + CurrentPopulation->Velocity[VariableIndex];

            NewPosition = CLAMP(NewPosition, this->LowerBound_[VariableIndex], this->UpperBound_[VariableIndex]);

            return NewPosition;
        }

        bool Run ()
        {
            if (FitnessFunction_ == nullptr)
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
                    double RandomVelocity = (this->LowerBound_[VariableIndex] - RandomPosition) +
                                             GenerateRandom(0.0f, 1.0f) * (this->UpperBound_[VariableIndex] - this->LowerBound_[VariableIndex]);
                    
                    Position[VariableIndex] = RandomPosition;
                    Velocity[VariableIndex] = RandomVelocity;
                }

                CurrentPopulation->Position = Position;
                CurrentPopulation->Velocity = Velocity;

                double FitnessValue = FitnessFunction_(CurrentPopulation->Position);
                CurrentPopulation->FitnessValue = FitnessValue;
                
                this->AverageFitnessValue_ += FitnessValue;

                CurrentPopulation->BestPosition = CurrentPopulation->Position;
                CurrentPopulation->BestFitnessValue = FitnessValue;
                CurrentPopulation->Feedback = std::vector<double> (this->NVariable_);

                if (FitnessValue < this->GlobalBestFitnessValue_)
                {
                    this->GlobalBestPosition_ = CurrentPopulation->Position;
                    this->GlobalBestFitnessValue_ = FitnessValue;
                }
            }

            this->AverageFitnessValue_ /= this->NPopulation_;

            // Optimize
            for (int Iteration = 0; Iteration < this->MaximumIteration_; Iteration++)
            {
                for (int PopulationIndex = 0; PopulationIndex < this->NPopulation_; PopulationIndex++)
                {
                    auto *CurrentPopulation = &this->Population_[PopulationIndex];

                    Optimize(Iteration, CurrentPopulation);
                }

                this->NextAverageFitnessValue_ /= this->NPopulation_;
                this->AverageFitnessValue_ = this->NextAverageFitnessValue_;
                this->NextAverageFitnessValue_ = 0.0f;

                if (this->Log_)
                {
                    std::cout << "[INFO] Iteration: " << Iteration << " >>> " << "Best fitness value: " << this->GlobalBestFitnessValue_ << std::endl;
                }
            }

            std::cout << "[INFO] Completed." << std::endl;

            return SUCCESS;
        }

        void CalculateAdaptiveInertialWeight (AParticle *CurrentPopulation)
        {
            if (CurrentPopulation->FitnessValue <= this->AverageFitnessValue_)
            {
                this->InertialWeight_ = this->MinimumInertialWeight_ + (this->MaximumInertialWeight_ - this->MinimumInertialWeight_) *
                                        ((CurrentPopulation->FitnessValue - this->GlobalBestFitnessValue_) / (this->AverageFitnessValue_ - this->GlobalBestFitnessValue_));
            }
            else
            {
                this->InertialWeight_ = this->MinimumInertialWeight_ + (this->MaximumInertialWeight_ - this->MinimumInertialWeight_) *
                                        ((this->AverageFitnessValue_ - this->GlobalBestFitnessValue_) / (CurrentPopulation->FitnessValue - this->GlobalBestFitnessValue_));
            }
        }

        void Optimize (int Iteration, AParticle *CurrentPopulation)
        {
            CalculateAdaptiveInertialWeight(CurrentPopulation);

            // Update Velocity
            for (int VariableIndex = 0; VariableIndex < this->NVariable_; VariableIndex++)
            {
                double NewVelocity = UpdateVelocity(CurrentPopulation, VariableIndex);

                CurrentPopulation->Velocity[VariableIndex] = NewVelocity;
            }

            // Update Position
            for (int VariableIndex = 0; VariableIndex < this->NVariable_; VariableIndex++)
            {
                double NewPosition = UpdatePosition(CurrentPopulation, VariableIndex);

                CurrentPopulation->Position[VariableIndex] = NewPosition;
            }

            // Evaluate Fitness Value
            double FitnessValue = FitnessFunction_(CurrentPopulation->Position);
            CurrentPopulation->FitnessValue = FitnessValue;

            this->NextAverageFitnessValue_ += FitnessValue;

            // Update Previous Best
            if (FitnessValue < CurrentPopulation->BestFitnessValue)
            {
                CurrentPopulation->BestPosition = CurrentPopulation->Position;
                CurrentPopulation->BestFitnessValue = FitnessValue;

                // This kind of awful!
                for (int VariableIndex = 0; VariableIndex < this->NVariable_; VariableIndex++)
                {
                    CurrentPopulation->Feedback[VariableIndex] = (1.0f / static_cast<double>(Iteration + 1)) * CurrentPopulation->Feedback[VariableIndex] +
                                                                 (CurrentPopulation->Velocity[VariableIndex]) *
                                                                 GenerateRandom(0.0f, 1.0f);
                }
            }
            else
            {
                for (int VariableIndex = 0; VariableIndex < this->NVariable_; VariableIndex++)
                {
                    CurrentPopulation->Feedback[VariableIndex] = (1.0f / static_cast<double>(Iteration + 1)) * CurrentPopulation->Feedback[VariableIndex] -
                                                                 (CurrentPopulation->Velocity[VariableIndex]) *
                                                                 GenerateRandom(0.0f, 1.0f);
                }
            }

            // Update Global Best
            if (FitnessValue < this->GlobalBestFitnessValue_)
            {
                this->GlobalBestPosition_ = CurrentPopulation->Position;
                this->GlobalBestFitnessValue_ = FitnessValue;
            }
        }

        std::vector<double> GetGlobalBestPosition () const
        {
            return this->GlobalBestPosition_;
        }

        double GetGlobalBestFitnessValue () const
        {
            return this->GlobalBestFitnessValue_;
        }

    private:
        std::vector<double> LowerBound_, UpperBound_;
        int MaximumIteration_, NPopulation_, NVariable_;
        double InertialWeight_{}, SocialCoefficient_, CognitiveCoefficient_;
        double MaximumInertialWeight_ = 0.9f, MinimumInertialWeight_ = 0.4;
        double VelocityFactor_;

        double (*FitnessFunction_)(const std::vector<double> &Position) = nullptr;

        std::vector<AParticle> Population_;
        std::vector<double> MaximumVelocity_, MinimumVelocity_;

        std::vector<double> GlobalBestPosition_;
        double GlobalBestFitnessValue_ = INFINITY;
        double AverageFitnessValue_ = 0.0f;
        double NextAverageFitnessValue_ = 0.0f;

        bool Log_;
    };
} // Optimizer

#endif // PSO_H
