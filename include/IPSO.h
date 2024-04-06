//
// Created by Sippawit Thammawiset on 25/1/2024 AD.
//

#ifndef IPSO_H
#define IPSO_H

#include <iostream>
#include <vector>
#include <random>
#include <cmath>

/**
 * @brief A macro to clamp a value between a minimum and a maximum.
 */
#define CLAMP(X, MIN, MAX)                      std::max(MIN, std::min(MAX, X))

/**
 * @brief A macro to check if a value is out of bounds.
 */
#define IS_OUT_OF_BOUND(X, MIN, MAX)            X < MIN || X > MAX

namespace MTH::IPSO
{
    namespace STATE
    {
        /**
         * @brief An enum defining the states of the algorithm.
         */
        enum
        {
            FAILED = 0,
            SUCCESS = 1
        };
    }

    namespace VELOCITY_CONFINEMENT
    {
        /**
         * @brief An enum defining different types of velocity confinement methods.
         */
        enum VELOCITY_CONFINEMENT
        {
            RANDOM_BACK = 0,
            HYPERBOLIC = 1,
            MIXED = 2
        };
    }

    /**
     * @brief Generate a random number within a specified range.
     *
     * @param LowerBound Lower bound of the range.
     * @param UpperBound Upper bound of the range.
     * @return Random number within the specified range.
     */
    double GenerateRandom (double LowerBound = 0.0f, double UpperBound = 1.0f)
    {
        std::random_device Engine;
        std::uniform_real_distribution<double> RandomDistribution(0.0f, 1.0f);
        return LowerBound + RandomDistribution(Engine) * (UpperBound - LowerBound);
    }

    /**
     * @brief A Struct representing a particle in the IPSO algorithm.
     */
    typedef struct AParticle
    {
        AParticle () : BestCost(2e10), Cost(0.0f) {}

        std::vector<double> Position;
        std::vector<double> Velocity;
        double Cost;

        std::vector<double> BestPosition;
        double BestCost;

        std::vector<double> Feedback;
    } AParticle;

    /**
     * @brief A class representing the Improved Particle Swarm Optimization (IPSO) algorithm.
     */
    class AIPSO
    {
    public:
        /**
         * @brief Constructor
         *
         * @param LowerBound Lower bound of the search space.
         * @param UpperBound Upper bound of the search space.
         * @param MaximumIteration Maximum number of iterations.
         * @param NPopulation Population size.
         * @param NVariable Number of variables.
         * @param SocialCoefficient Social coefficient for velocity update.
         * @param CognitiveCoefficient Cognitive coefficient for velocity update.
         * @param VelocityFactor Factor for limiting velocity update.
         * @param VelocityConfinement Type of velocity confinement method.
         * @param Log Flag indicating whether to log information during optimization.
         */
        inline AIPSO (const std::vector<double> &LowerBound, const std::vector<double> &UpperBound,
               int MaximumIteration, int NPopulation, int NVariable,
               double SocialCoefficient = 1.5f, double CognitiveCoefficient = 1.5f,
               double VelocityFactor = 0.5f,
               int VelocityConfinement = VELOCITY_CONFINEMENT::RANDOM_BACK,
               bool Log = true) :
               LowerBound_(LowerBound),
               UpperBound_(UpperBound),
               MaximumIteration_(MaximumIteration),
               NPopulation_(NPopulation),
               NVariable_(NVariable),
               SocialCoefficient_(SocialCoefficient),
               CognitiveCoefficient_(CognitiveCoefficient),
               VelocityFactor_(VelocityFactor),
               VelocityConfinement_(VelocityConfinement),
               Log_(Log)
        {

        }

        inline ~AIPSO () = default;

        /**
         * @brief Set the objective function for optimization.
         *
         * @param UserObjectiveFunction Pointer to the objective function.
         */
        void SetObjectiveFunction (double (*UserObjectiveFunction)(const std::vector<double> &Position))
        {
            this->ObjectiveFunction_ = UserObjectiveFunction;
        }

        /**
         * @brief Run the IPSO algorithm.
         *
         * @return Success or failure flag.
         */
        bool Run ()
        {
            // Check if the objective function is defined
            if (ObjectiveFunction_ == nullptr)
            {
                std::cerr << "Objective function is not defined. Please use SetObjectiveFunction(UserObjectiveFunction) before calling Run()." << std::endl;

                return STATE::FAILED;
            }

            // Check if the size of lower bound or upper bound matches the number of variables
            if (this->LowerBound_.size() != NVariable_ || this->UpperBound_.size() != NVariable_)
            {
                std::cerr << "Size of lowerbound or upperbound does not match with number of variables." << std::endl;

                return STATE::FAILED;
            }

            // Initialize population and velocity bounds
            this->Population_ = std::vector<AParticle> (this->NPopulation_);
            this->MaximumVelocity_ = std::vector<double> (this->NVariable_);
            this->MinimumVelocity_ = std::vector<double> (this->NVariable_);

            // Compute maximum and minimum velocity bounds for each variable
            for (int VariableIndex = 0; VariableIndex < this->NVariable_; ++VariableIndex)
            {
                this->MaximumVelocity_[VariableIndex] = this->VelocityFactor_ * (this->UpperBound_[VariableIndex] - this->LowerBound_[VariableIndex]);
                this->MinimumVelocity_[VariableIndex] = -this->MaximumVelocity_[VariableIndex];
            }

            // Initialize particles with random positions and velocities
            for (int PopulationIndex = 0; PopulationIndex < this->NPopulation_; ++PopulationIndex)
            {
                auto *CurrentPopulation = &this->Population_[PopulationIndex];

                std::vector<double> Position(this->NVariable_);
                std::vector<double> Velocity(this->NVariable_);

                // Generate random positions and velocities within bounds for each variable
                for (int VariableIndex = 0; VariableIndex < this->NVariable_; ++VariableIndex)
                {
                    double RandomPosition = GenerateRandom(this->LowerBound_[VariableIndex], this->UpperBound_[VariableIndex]);

                    // The initialized velocity is derived from Equation (3.4) of "Standard Particle Swarm Optimisation" by Maurice Clerc. Link: https://hal.science/hal-00764996/document
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

                // Update personal best position and cost
                CurrentPopulation->BestPosition = CurrentPopulation->Position;
                CurrentPopulation->BestCost = Cost;
                CurrentPopulation->Feedback = std::vector<double> (this->NVariable_);

                // Update global best position and cost
                if (Cost < this->GlobalBestCost_)
                {
                    this->GlobalBestPosition_ = CurrentPopulation->Position;
                    this->GlobalBestCost_ = Cost;
                }
            }

            // Compute average cost of the initial population
            this->AverageCost_ /= this->NPopulation_;

            // Run optimization process for a maximum number of iterations
            for (int Iteration = 1; Iteration <= this->MaximumIteration_; ++Iteration)
            {
                // Optimize particle positions and velocities
                Optimize(Iteration);

                // Update average cost for the next iteration
                this->NextAverageCost_ /= this->NPopulation_;
                this->AverageCost_ = this->NextAverageCost_;
                this->NextAverageCost_ = 0.0f;

                if (this->Log_)
                {
                    std::cout << "[INFO] Iteration: " << Iteration << " >>> \n\t"
                              << "Best Cost: " << this->GlobalBestCost_ <<
                              std::endl;;
                }
            }

            std::cout << "[INFO] Completed." << std::endl;

            return STATE::SUCCESS;
        }

        /**
         * @brief Get the global best position found by the algorithm.
         *
         * @return Global best position.
         */
        std::vector<double> GetGlobalBestPosition () const
        {
            return this->GlobalBestPosition_;
        }

        /**
         * @brief Get the global best cost found by the algorithm.
         *
         * @return Global best cost.
         */
        double GetGlobalBestCost () const
        {
            return this->GlobalBestCost_;
        }

    private:
        std::vector<double> LowerBound_; /**< Lower bound of the search space. */
        std::vector<double> UpperBound_; /**< Upper bound of the search space. */
        int MaximumIteration_; /**< Maximum number of iterations. */
        int NPopulation_; /**< Population size. */
        int NVariable_; /**< Number of variables. */
        double InertialWeight_{}; /**< Current inertial weight. */
        double SocialCoefficient_; /**< Social coefficient for velocity update. */
        double CognitiveCoefficient_; /**< Cognitive coefficient for velocity update. */
        double MaximumInertialWeight_ = 0.9f; /**< Maximum value of inertial weight. */
        double MinimumInertialWeight_ = 0.4f; /**< Minimum value of inertial weight. */
        double VelocityFactor_; /**< Factor for limiting velocity update. */

        double (*ObjectiveFunction_)(const std::vector<double> &Position) = nullptr; /**< Pointer to the objective function. */

        std::vector<AParticle> Population_; /**< Vector containing particles. */
        std::vector<double> MaximumVelocity_; /**< Vector containing maximum velocity for each variable. */
        std::vector<double> MinimumVelocity_; /**< Vector containing minimum velocity for each variable. */

        std::vector<double> GlobalBestPosition_; /**< Global best position found by the algorithm. */
        double GlobalBestCost_ = INFINITY; /**< Global best cost found by the algorithm. */
        double AverageCost_ = 0.0f; /**< Average cost of the population. */
        double NextAverageCost_ = 0.0f; /**< Next average cost of the population. */

        int VelocityConfinement_; /**< Type of velocity confinement method. */
        bool Log_; /**< Flag indicating whether to log information during optimization. */

        /**
         * @brief Optimize the IPSO algorithm for the current iteration.
         *
         * This function optimizes the IPSO algorithm for the current iteration.
         * It updates the velocity and position of each particle, evaluates the cost of the updated position,
         * updates the best position of each particle, and updates the global best position.
         *
         * @note A random positive feedback factor is added into the velocity update formula of the particles, which enhances the ability to find a local optimum.
         *       This information is extracted from the paper "A cubic spline method combing improved particle swarm optimization for robot path planning in dynamic uncertain environment"
         *       by Wen Li, Mao Tan, Ling Wang, and Qiuzhen Wang. The feedback update is derived from Equation (7) and (8) of their paper.
         *       Link to the paper: https://journals.sagepub.com/doi/10.1177/1729881419891661
         */
        void Optimize (int Iteration)
        {
            // Iterate over each particle in the population
            for (int PopulationIndex = 0; PopulationIndex < this->NPopulation_; ++PopulationIndex)
            {
                auto CurrentPopulation = &this->Population_[PopulationIndex];

                // Calculate adaptive inertial weight for the current particle
                CalculateAdaptiveInertialWeight(CurrentPopulation);

                // Update velocity for each dimension of the particle
                for (int VariableIndex = 0; VariableIndex < this->NVariable_; ++VariableIndex)
                {
                    double NewVelocity = UpdateVelocity(CurrentPopulation, VariableIndex);

                    CurrentPopulation->Velocity[VariableIndex] = NewVelocity;
                }

                // Update position for each dimension of the particle
                for (int VariableIndex = 0; VariableIndex < this->NVariable_; ++VariableIndex)
                {
                    double NewPosition = UpdatePosition(CurrentPopulation, VariableIndex);

                    CurrentPopulation->Position[VariableIndex] = NewPosition;
                }

                // Evaluate cost of the updated position for the particle
                double Cost = ObjectiveFunction_(CurrentPopulation->Position);
                CurrentPopulation->Cost = Cost;

                // Accumulate cost for calculating average cost
                this->NextAverageCost_ += Cost;

                // Update personal best position, cost, and feedback of the particle if applicable
                if (Cost < CurrentPopulation->BestCost)
                {
                    CurrentPopulation->BestPosition = CurrentPopulation->Position;
                    CurrentPopulation->BestCost = Cost;

                    for (int VariableIndex = 0; VariableIndex < this->NVariable_; ++VariableIndex)
                    {
                        // Equation (7)
                        CurrentPopulation->Feedback[VariableIndex] = (1.0f / static_cast<double>(Iteration)) * CurrentPopulation->Feedback[VariableIndex] +
                                                                     (CurrentPopulation->Velocity[VariableIndex]) *
                                                                     GenerateRandom(0.0f, 1.0f);
                    }
                }
                else
                {
                    for (int VariableIndex = 0; VariableIndex < this->NVariable_; ++VariableIndex)
                    {
                        // Equation (8)
                        CurrentPopulation->Feedback[VariableIndex] = (1.0f / static_cast<double>(Iteration)) * CurrentPopulation->Feedback[VariableIndex] -
                                                                     (CurrentPopulation->Velocity[VariableIndex]) *
                                                                     GenerateRandom(0.0f, 1.0f);
                    }
                }

                // Update global best position and cost if applicable
                if (Cost < this->GlobalBestCost_)
                {
                    this->GlobalBestPosition_ = CurrentPopulation->Position;
                    this->GlobalBestCost_ = Cost;
                }
            }
        }

        /**
         * @brief Calculate the adaptive inertial weight based on the current population's cost.
         *
         * This function calculates the adaptive inertial weight based on the current population's cost,
         * the average cost of the population, and the global best cost.
         * The inertial weight is adjusted based on the difference between the current cost and the global best cost and
         * the difference between the average cost of the population and the global best cost.
         *
         * @param CurrentPopulation Pointer to the current population.
         *
         * @note This adaptive inertia weight adjustment mechanism is inspired by the work of Zhenjian Yang, Ning Li, Yunjie Zhang, and Jin Li,
         *       who introduced adaptive inertia weight to balance the exploration ability.
         *       Specifically, they proposed an adaptive inertia weight approach in their paper "Mobile Robot Path Planning Based on Improved Particle Swarm Optimization and Improved Dynamic Window Approach".
         *       The inertia weight adjustment method described here is derived from Equation (8) of their paper.
         *       Link to the paper: https://www.hindawi.com/journals/jr/2023/6619841/
         */
        void CalculateAdaptiveInertialWeight(AParticle *CurrentPopulation)
        {
            if (CurrentPopulation->Cost <= this->AverageCost_)
            {
                this->InertialWeight_ = this->MinimumInertialWeight_ +
                                        (this->MaximumInertialWeight_ - this->MinimumInertialWeight_) *
                                        ((CurrentPopulation->Cost - this->GlobalBestCost_) /
                                         (this->AverageCost_ - this->GlobalBestCost_));
            }
            else
            {
                this->InertialWeight_ = this->MinimumInertialWeight_ +
                                        (this->MaximumInertialWeight_ - this->MinimumInertialWeight_) *
                                        ((this->AverageCost_ - this->GlobalBestCost_) /
                                         (CurrentPopulation->Cost - this->GlobalBestCost_));
            }
        }

        /**
         * @brief Update the velocity of a particle for a specific dimension.
         *
         * This function updates the velocity of a particle for a specific dimension.
         * It considers the inertial weight, social and cognitive coefficients, as well as a random positive feedback factor
         * to adjust the velocity.
         *
         * @param CurrentPopulation Pointer to the current particle.
         * @param VariableIndex The index of the dimension for which the velocity is updated.
         *
         * @return The updated velocity for the specified dimension of the particle.
         */
        double UpdateVelocity(const AParticle *CurrentPopulation, int VariableIndex)
        {
            // Calculate the new velocity
            double NewVelocity = this->InertialWeight_ * CurrentPopulation->Velocity[VariableIndex] +
                                 this->SocialCoefficient_ * GenerateRandom(0.0f, 1.0f) * (CurrentPopulation->BestPosition[VariableIndex] - CurrentPopulation->Position[VariableIndex]) +
                                 this->CognitiveCoefficient_ * GenerateRandom(0.0f, 1.0f) * (this->GlobalBestPosition_[VariableIndex] - CurrentPopulation->Position[VariableIndex]) +
                                 CurrentPopulation->Feedback[VariableIndex];

            // Clamp the new velocity to stay within specified bounds
            NewVelocity = CLAMP(NewVelocity, this->MinimumVelocity_[VariableIndex], this->MaximumVelocity_[VariableIndex]);

            return NewVelocity;
        }

        /**
         * @brief Update the position of a particle for a specific dimension.
         *
         * This function updates the position of a particle for a specific dimension.
         * It considers the particle's current position and velocity, applies velocity confinement if necessary,
         * and clamps the new position within specified bounds.
         *
         * @param CurrentPopulation Pointer to the current particle.
         * @param VariableIndex The index of the dimension for which the position is updated.
         *
         * @return The updated position for the specified dimension of the particle.
         */
        double UpdatePosition(AParticle *CurrentPopulation, int VariableIndex)
        {
            // Calculate the temporary new position by adding the updated velocity to the current position
            double TemporaryNewPosition = CurrentPopulation->Position[VariableIndex] + CurrentPopulation->Velocity[VariableIndex];

            // Apply velocity confinement if the temporary new position is out of bounds
            if (IS_OUT_OF_BOUND(TemporaryNewPosition, this->LowerBound_[VariableIndex], this->UpperBound_[VariableIndex]))
            {
                double Confinement;

                switch (this->VelocityConfinement_)
                {
                    case VELOCITY_CONFINEMENT::RANDOM_BACK:
                        Confinement = RandomBackConfinement(CurrentPopulation->Velocity[VariableIndex]);
                        break;

                    case VELOCITY_CONFINEMENT::HYPERBOLIC:
                        Confinement = HyperbolicConfinement(this->LowerBound_[VariableIndex],
                                                            this->UpperBound_[VariableIndex],
                                                            CurrentPopulation->Position[VariableIndex],
                                                            CurrentPopulation->Velocity[VariableIndex]);
                        break;

                    case VELOCITY_CONFINEMENT::MIXED:
                        Confinement = MixedConfinement(this->LowerBound_[VariableIndex],
                                                       this->UpperBound_[VariableIndex],
                                                       CurrentPopulation->Position[VariableIndex],
                                                       CurrentPopulation->Velocity[VariableIndex]);
                        break;

                    default:
                        Confinement = RandomBackConfinement(CurrentPopulation->Velocity[VariableIndex]);
                }

                // Apply the confinement to the velocity
                CurrentPopulation->Velocity[VariableIndex] = Confinement;
            }

            // Calculate the new position by adding the updated velocity to the current position
            double NewPosition = CurrentPopulation->Position[VariableIndex] + CurrentPopulation->Velocity[VariableIndex];

            // Clamp the new position to stay within specified bounds
            NewPosition = CLAMP(NewPosition, this->LowerBound_[VariableIndex], this->UpperBound_[VariableIndex]);

            return NewPosition;
        }

        /**
         * @brief Apply random back velocity confinement.
         *
         * This function applies random back velocity confinement to a given velocity.
         * It generates a random factor within the range [0, 1] and multiplies it with the given velocity.
         *
         * @param Velocity The velocity to which random back confinement is applied.
         *
         * @return The velocity after applying random back confinement.
         */
        static double RandomBackConfinement (const double Velocity)
        {
            double VelocityConfinement = -GenerateRandom(0.0f, 1.0f) * Velocity;

            return VelocityConfinement;
        }

        /**
         * @brief Apply hyperbolic velocity confinement.
         *
         * This function applies hyperbolic velocity confinement to a given velocity.
         * It calculates the velocity confinement based on the position and velocity of the particle,
         * as well as the specified lower and upper bounds for the position.
         *
         * @param LowerBound The lower bound for the position (search space).
         * @param UpperBound The upper bound for the position (search space).
         * @param Position The current position of the particle.
         * @param Velocity The velocity to which hyperbolic confinement is applied.
         *
         * @return The velocity after applying hyperbolic confinement.
         */
        static double HyperbolicConfinement(const double LowerBound,
                                            const double UpperBound,
                                            const double Position,
                                            const double Velocity)
        {
            double VelocityConfinement;

            if (Velocity > 0)
            {
                VelocityConfinement = Velocity / (1.0f + abs(Velocity / (UpperBound - Position)));
            }
            else
            {
                VelocityConfinement = Velocity / (1.0f + abs(Velocity / (Position - LowerBound)));
            }

            return VelocityConfinement;
        }

        /**
          * @brief Apply mixed velocity confinement.
          *
          * This function applies mixed velocity confinement to a given velocity.
          * It randomly selects between hyperbolic and random back velocity confinement based on a generated random factor.
          *
          * @param LowerBound The lower bound for the position (search space).
          * @param UpperBound The upper bound for the position (search space).
          * @param Position The current position of the particle.
          * @param Velocity The velocity to which mixed confinement is applied.
          *
          * @return The velocity after applying mixed confinement.
          */
        static double MixedConfinement(const double LowerBound,
                                       const double UpperBound,
                                       const double Position,
                                       const double Velocity)
        {
            double VelocityConfinement;

            if (GenerateRandom(0.0f, 1.0) >= 0.5f)
            {
                VelocityConfinement = HyperbolicConfinement(LowerBound, UpperBound, Position, Velocity);
            }
            else
            {
                VelocityConfinement = RandomBackConfinement(Velocity);
            }

            return VelocityConfinement;
        }
    };
} // MTH

#endif // IPSO_H
