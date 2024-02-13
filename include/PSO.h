//
// Created by Sippawit Thammawiset on 25/1/2024 AD.
//

/* TODO:    - Add comments and documentation
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
        AParticle () : Cost(0.0f), BestCost(INFINITY) {}

        std::vector<double> Position;
        std::vector<double> Velocity;
        double Cost;

        std::vector<double> BestPosition;
        double BestCost;

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

        void SetObjectiveFunction (double (*UserObjectiveFunction)(const std::vector<double> &Position))
        {
            this->ObjectiveFunction_ = UserObjectiveFunction;
        }

        bool Run ()
        {
            if (ObjectiveFunction_ == nullptr)
            {
                std::cerr << "Objective function is not defined. Please use SetObjectiveFunction(UserObjectiveFunction) before calling Run()." << std::endl;

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

                double Cost = ObjectiveFunction_(CurrentPopulation->Position);
                CurrentPopulation->Cost = Cost;
                
                this->AverageCost_ += Cost;

                CurrentPopulation->BestPosition = CurrentPopulation->Position;
                CurrentPopulation->BestCost = Cost;
                CurrentPopulation->Feedback = std::vector<double> (this->NVariable_);

                if (Cost < this->GlobalBestCost_)
                {
                    this->GlobalBestPosition_ = CurrentPopulation->Position;
                    this->GlobalBestCost_ = Cost;
                }
            }

            this->AverageCost_ /= this->NPopulation_;

            // Optimize
            for (int Iteration = 1; Iteration <= this->MaximumIteration_; Iteration++)
            {
                for (int PopulationIndex = 0; PopulationIndex < this->NPopulation_; PopulationIndex++)
                {
                    auto *CurrentPopulation = &this->Population_[PopulationIndex];

                    Optimize(Iteration, CurrentPopulation);
                }

                this->NextAverageCost_ /= this->NPopulation_;
                this->AverageCost_ = this->NextAverageCost_;
                this->NextAverageCost_ = 0.0f;

                if (this->Log_)
                {
                    std::cout << "[INFO] Iteration: " << Iteration << " >>> "
                              << "Best Cost: " << this->GlobalBestCost_ << std::endl;
                }
            }

            std::cout << "[INFO] Completed." << std::endl;

            return SUCCESS;
        }

        std::vector<double> GetGlobalBestPosition () const
        {
            return this->GlobalBestPosition_;
        }

        double GetGlobalBestCost () const
        {
            return this->GlobalBestCost_;
        }

    private:
        std::vector<double> LowerBound_, UpperBound_;
        int MaximumIteration_, NPopulation_, NVariable_;
        double InertialWeight_{}, SocialCoefficient_, CognitiveCoefficient_;
        double MaximumInertialWeight_ = 0.9f, MinimumInertialWeight_ = 0.4;
        double VelocityFactor_;

        double (*ObjectiveFunction_)(const std::vector<double> &Position) = nullptr;

        std::vector<AParticle> Population_;
        std::vector<double> MaximumVelocity_, MinimumVelocity_;

        std::vector<double> GlobalBestPosition_;
        double GlobalBestCost_ = INFINITY;
        double AverageCost_ = 0.0f;
        double NextAverageCost_ = 0.0f;

        bool Log_;

        void CalculateAdaptiveInertialWeight (AParticle *CurrentPopulation)
        {
            if (CurrentPopulation->Cost <= this->AverageCost_)
            {
                this->InertialWeight_ = this->MinimumInertialWeight_ + (this->MaximumInertialWeight_ - this->MinimumInertialWeight_) *
                                        ((CurrentPopulation->Cost - this->GlobalBestCost_) / (this->AverageCost_ - this->GlobalBestCost_));
            }
            else
            {
                this->InertialWeight_ = this->MinimumInertialWeight_ + (this->MaximumInertialWeight_ - this->MinimumInertialWeight_) *
                                        ((this->AverageCost_ - this->GlobalBestCost_) / (CurrentPopulation->Cost - this->GlobalBestCost_));
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

            // Evaluate Cost
            double Cost = ObjectiveFunction_(CurrentPopulation->Position);
            CurrentPopulation->Cost = Cost;

            this->NextAverageCost_ += Cost;

            // Update Previous Best
            if (Cost < CurrentPopulation->BestCost)
            {
                CurrentPopulation->BestPosition = CurrentPopulation->Position;
                CurrentPopulation->BestCost = Cost;

                for (int VariableIndex = 0; VariableIndex < this->NVariable_; VariableIndex++)
                {
                    CurrentPopulation->Feedback[VariableIndex] = (1.0f / static_cast<double>(Iteration)) * CurrentPopulation->Feedback[VariableIndex] +
                                                                 (CurrentPopulation->Velocity[VariableIndex]) *
                                                                 GenerateRandom(0.0f, 1.0f);
                }
            }
            else
            {
                for (int VariableIndex = 0; VariableIndex < this->NVariable_; VariableIndex++)
                {
                    CurrentPopulation->Feedback[VariableIndex] = (1.0f / static_cast<double>(Iteration)) * CurrentPopulation->Feedback[VariableIndex] -
                                                                 (CurrentPopulation->Velocity[VariableIndex]) *
                                                                 GenerateRandom(0.0f, 1.0f);
                }
            }

            // Update Global Best
            if (Cost < this->GlobalBestCost_)
            {
                this->GlobalBestPosition_ = CurrentPopulation->Position;
                this->GlobalBestCost_ = Cost;
            }
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
    };
} // Optimizer

#endif // PSO_H
