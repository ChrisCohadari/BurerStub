#ifndef HEURISTICS_BASE_SOLUTION_H_
#define HEURISTICS_BASE_SOLUTION_H_

#include <vector>

namespace Burer
{
    class BaseSolution {
     public:
      // Constructor takes the number of nodes in the graph / vars in the problem
      // and the initial value for the assignments_ vector in the init_assignment
      // parameter.
      BaseSolution(int N, int init_assignment);
      
      // Constructor takes the assignments, weight, and problem.
      BaseSolution(const std::vector<int>& assignments, double weight);

      // Return the Hamming distance to another solution.
      int SymmetricDifference(const BaseSolution& other) const;
      
      // Return the Hamming distance to another solution, and populate a vector of
      // the indices that differ.
      int SymmetricDifference(const BaseSolution& other,
			      std::vector<int>* diff) const;

      // Return the Hamming distance to another solution, and populate a vector of
      // the indices that differ and indices in common
      int SymmetricDifference(const BaseSolution& other,
			      std::vector<int>* diff,
			      std::vector<int>* common) const;

      // Identify if one solution improves over another solution. These functions take
      // into account the fact that our problems are stated with floating point
      // numbers instead of integers as weights.
      static bool ImprovesOver(double weight1, double weight2);
      bool ImprovesOver(double other_weight) const {
	return ImprovesOver(weight_, other_weight);
      }
      bool ImprovesOver(const BaseSolution& other) const {
	return ImprovesOver(weight_, other.weight_);
      }

      // Print out the decision variables for a solutio
      void PrintSolution() const;

      // Getters
      double get_weight() const {  return weight_; }
      const std::vector<int>& get_assignments() const {  return assignments_; }

      // Assignment operator
      BaseSolution& operator=(const BaseSolution &rhs);

      // Copy constructor
      BaseSolution(const BaseSolution& x);

      // Equals operator (checks if the solutions have identical assignments_)
      bool operator==(const BaseSolution& other) const;

      // Not equals operator
      bool operator!=(const BaseSolution& other) const {  return !(*this == other); }

      // >, <, >=, <=: compare based on weight_
      bool operator>(const BaseSolution& other) const {
	return weight_ > other.weight_;
      }
      bool operator<(const BaseSolution& other) const {
	return weight_ < other.weight_;
      }
      bool operator>=(const BaseSolution& other) const {
	return weight_ >= other.weight_;
      }
      bool operator<=(const BaseSolution& other) const {
	return weight_ <= other.weight_;
      }

     protected:
      // Assignment of each node (-1/1 in MAXCUT, 0/1 in QUBO)
      std::vector<int> assignments_;
      // Objective value
      double weight_;
      // The number of nodes in the graph / variables in the problem
      int N_;

     private:
      // Default constructor disabled
      BaseSolution();
    };


    class ExtendedSolution : public BaseSolution {
     public:
      // Constructor takes the number of nodes in the graph / vars in the problem
      // and the initial value for the assignments_ vector in the init_assignment
      // parameter.
      ExtendedSolution(int N, int init_assignment);

      // 1-swap functions: startpos means you won't consider 1-swaps for any index
      // less than this parameter.

      // Perform all available 1-moves, at each step selecting the most valuable.
      void AllBest1Swap(int startpos = 0);

      // Perform all 1-moves, at each step selecting the first improving move
      void AllFirst1Swap(int startpos = 0);

      // Perform all 1-moves, iteratively taking the first improving move on a
      // shuffled list of indices.
      void AllShuffle1Swap(int startpos = 0);

      // A popular simulated annealing update mechanism uses the standard deviation
      // of the diff_weights_ values, so we'll provide an accessor for this value.
      double DiffWeightStandardDeviation() const;

      // Identify if a move sufficiently improves the quality of this solution to be
      // taken.
      bool ImprovingMove(int index) const {
	return ImprovesOver(weight_ + diff_weights_[index], weight_);
      }
      bool ImprovingMove(double diff_weight) const {
	return ImprovesOver(weight_ + diff_weight, weight_);
      }
      bool ImprovesOverAfterMove(double weight2, int index) const {
	return ImprovesOver(weight_ + diff_weights_[index], weight2);
      }
      bool ImprovesOverAfterMove(const BaseSolution& other, int index) const {
	return ImprovesOver(weight_ + diff_weights_[index], other.get_weight());
      }
      bool NonDetrimentalMove(int index) const {
	return !ImprovesOver(weight_, weight_ + diff_weights_[index]);
      }

      // Getter
      const std::vector<double>& get_diff_weights() const {  return diff_weights_; }

      // Assignment operator
      ExtendedSolution& operator=(const ExtendedSolution &rhs);

      // Copy constructor
      ExtendedSolution(const ExtendedSolution& x);

     protected:
      // Switch the set of update_index, updating the assignments (x), the
      // diff_weights, and the objective. Extending classes must implement.
      virtual void UpdateCutValues(int update_index, std::vector<int>* x,
				   std::vector<double>* diff_weights,
				   double *objective) const = 0;

      // Switch the set of update_index, updating class variables assignments_,
      // diff_weights_, and weight_
      void UpdateCutValues(int update_index) {
	UpdateCutValues(update_index, &assignments_, &diff_weights_, &weight_);
      }

      // Amount you would gain from switching sets / flipping variable
      std::vector<double> diff_weights_;

     private:
      // Default constructor disabled
      ExtendedSolution();
    };


    class MaxCutSimpleSolution : public BaseSolution {
     public:
      // Initialize a completely random solution (p=0.5 for each vertex)
      static MaxCutSimpleSolution RandomSolution() {
	return MaxCutSimpleSolution(0, 0);  // Private constructor
      }

      // Initialize a random solution (p for each vertex provided)
      static MaxCutSimpleSolution RandomSolution(
						 const std::vector<double>& p
						 ) {
	return MaxCutSimpleSolution(p);
      }

      // Initialize mi and heuristic members but nothing else. This is a light-
      // weight constructor that labels the solution as having weight 0.
      MaxCutSimpleSolution(int initVal = -1);

      // Initialize from mi, heuristic, assignments, and weight. This is a light-
      // weight constructor that does not re-compute the weight.
      MaxCutSimpleSolution(
			   const std::vector<int>& assignments, double weight);

      // Assignment operator
      MaxCutSimpleSolution& operator=(const MaxCutSimpleSolution &rhs);

      // Copy constructor
      MaxCutSimpleSolution(const MaxCutSimpleSolution& x);

      // This function initializes weight_, given that assignments_ is populated.
      void PopulateFromAssignments();

      // Initialize a completely random solution
      MaxCutSimpleSolution(
			   int ignored1, int ignored2);

      // Random solution given vertex probabilities
      MaxCutSimpleSolution(const std::vector<double>& p);
    };

    // Some constructive procedures for MAXCUT assign edges one at a time
    // and maintain the gains associated with assigning a node either to S
    // (assignment value 1) or NS (assignment value -1). The MaxCutPartialSolution
    // supports such procedures, handling variable updates. The MaxCutPartialSolution
    // can be converted (cheaply) to a MaxCutSolution through a constructor in
    // MaxCutSolution.
    class MaxCutPartialSolution {
     public:
      // Getters
      double get_weight() const {  return weight_; }
      int get_num_unassigned() const {  return num_unassigned_; }
      const std::vector<int>& get_assignments() const {  return assignments_; }
      const std::vector<double>& get_gainS() const {  return gainS_; }
      const std::vector<double>& get_gainNS() const {  return gainNS_; }

     protected:
      // Initialize a solution to all unassigned
      MaxCutPartialSolution();

      // Update gainS_, gainNS_, num_unassigned_, and weight_ from assignments_
      void PopulateFromAssignments();

      // Set the value of update_index to the specified new_value, updating
      // gainS_, gainNS_, assignments_, num_unassigned_, and weight_.
      void UpdateCutValues(int update_index, int new_value);

      // Convenience variable for number of variables in the problem
      int N_;
      // Assignment of each node (1 means S, -1 means NS, 0 means unassigned).
      std::vector<int> assignments_;
      // Change in objective from assigning each vertex to set S (value 1)
      std::vector<double> gainS_;
      // Change in objective from assigning each vertex to set NS (value -1)
      std::vector<double> gainNS_;
      // Number of assignments that are unassigned (value 0)
      int num_unassigned_;
      // Objective value
      double weight_;
    };

    class MaxCutSolution : public ExtendedSolution {
     public:
      // Copy over a MaxCutPartialSolution with no unassigned values
      MaxCutSolution(const MaxCutPartialSolution& x);

      // Initialize a completely random solution (p=0.5 for each vertex)
      static MaxCutSolution RandomSolution() {
	return MaxCutSolution(0, 0);  // Private constructor
      }

      // Initialize a random solution (p for each vertex provided)
      static MaxCutSolution RandomSolution(
					   const std::vector<double>& p
    ) {
	return MaxCutSolution(p);
      }

      // Initialize from vector of assignments
      MaxCutSolution(const std::vector<int>& assignments);

      // Perform all available 2-moves, at each step selecting the most valuable.
      void AllBest2Swap()  {  AllBest2Swap(0);  }

      // Perform all 2-moves, at each step selecting the first improving move
      void AllFirst2Swap()  {  AllFirst2Swap(0);  }

      // Assignment operator
      MaxCutSolution& operator=(const MaxCutSolution &rhs);

      // Copy constructor
      MaxCutSolution(const MaxCutSolution& x);

      // This function initializes diff_weights_ and weight_, given that
      // assignments_ is populated.
      void PopulateFromAssignments();

     protected:
      // Initialize mi and heuristic members but nothing else
      // ********************** IMPORTANT NOTE ***************************
      // diff_weights_ is not properly set after using this constructor, so you
      // must either set that vector manually or by using PopulateFromAssignments
      // *****************************************************************
      MaxCutSolution(
		     int initVal = -1);

      // Switch the set of update_index, updating, the assignments (x), the
      // diff_weights, and the objective.
      void UpdateCutValues(int update_index, std::vector<int>* x,
			   std::vector<double>* diff_weights,
			   double *objective) const;

      using ExtendedSolution::UpdateCutValues;  // Unhide single-argument version

      // Perform all available 2-moves, at each step selecting the most valuable.
      void AllBest2Swap(int startpos);

      // Perform all 2-moves, at each step selecting the first improving move
      void AllFirst2Swap(int startpos);

     private:
      // Initialize a completely random solution
      MaxCutSolution(
		     int ignored1, int ignored2);

      // Random solution given vertex probabilities
      MaxCutSolution(const std::vector<double>& p);
    };

    class FirstFixedMaxCutSolution : public MaxCutSolution {
     public:
      // Initialize a completely random solution (p=0.5 for each vertex except 1st)
      static FirstFixedMaxCutSolution RandomSolution(int fixedVal) {
	return FirstFixedMaxCutSolution(fixedVal, 0);
      }

      // Initialize a random solution (p for each vertex provided)
      static FirstFixedMaxCutSolution RandomSolution(
						     const std::vector<double>& p,

						     int fixedVal) {
	return FirstFixedMaxCutSolution(p, fixedVal);
      }

      // Perform all available 1-moves, at each step selecting the most valuable.
      // We set local search start pos to 1 so we don't flip the first index.
      void AllBest1Swap() {  ExtendedSolution::AllBest1Swap(1);  }

      // Perform all 1-moves, at each step selecting the first improving move.
      // We set local search start pos to 1 so we don't flip the first index.
      void AllFirst1Swap() {  ExtendedSolution::AllFirst1Swap(1);  }

      // Perform all 1-moves, iteratively taking the first improving move on a
      // shuffled list of indices. We set local search start pos to 1 so we don't
      // flip the first index.
      void AllShuffle1Swap() {  ExtendedSolution::AllShuffle1Swap(1);  }

      // Perform all available 2-moves, at each step selecting the most valuable.
      // We set local search start pos to 1 so we don't flip the first index.
      void AllBest2Swap() {  AllBest2Swap(1);  }

      // Perform all 2-moves, at each step selecting the first improving move.
      // We set local search start pos to 1 so we don't flip the first index.
      void AllFirst2Swap() {  AllFirst2Swap(1);  }

      // Assignment operator
      FirstFixedMaxCutSolution& operator=(const FirstFixedMaxCutSolution &rhs);

      // Copy constructor
      FirstFixedMaxCutSolution(const FirstFixedMaxCutSolution& x);

     protected:
      // Initialize mi, heuristic, and fixedVal_ members but nothing else
      FirstFixedMaxCutSolution(
			       int fixedVal);

      // This function initializes diff_weights_ and weight_, given that
      // assignments_ is populated.
      void PopulateFromAssignments();

      // Switch the set of update_index, updating, the assignments (x), the
      // diff_weights, and the objective.
      void UpdateCutValues(int update_index, std::vector<int>* x,
			   std::vector<double>* diff_weights,
			   double *objective) const;

      using ExtendedSolution::UpdateCutValues;  // Unhide single-argument version
      using MaxCutSolution::AllBest2Swap;  // Unhide 1-argument version
      using MaxCutSolution::AllFirst2Swap;  // Unhide 1-argument version

      // fixedVal_ is the fixed value for the first node (1 or -1)
      int fixedVal_;

     private:
      // Initialize a completely random solution
      FirstFixedMaxCutSolution(
			       int fixedVal, int ignored);

      // Random solution given vertex probabilities
      FirstFixedMaxCutSolution(
			       const std::vector<double>& p,
			       int fixedVal);
    };
}
#endif
