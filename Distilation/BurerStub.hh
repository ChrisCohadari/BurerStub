#ifndef BURERSTUB_GUARD_PROCEDURECUT_HH
#define BURERSTUB_GUARD_PROCEDURECUT_HH 1

#include <vector>

namespace Burer
{

  class Solution {
  public:
    Solution(int N, int init_assignment); //USED
    
    // Constructor takes the assignments, weight, and problem.
    Solution(const std::vector<int>& assignments, double weight); //USED

      static Solution Rank2Cut(double w1norm, std::vector<double>* theta, int G_n) {
	return Solution(w1norm, theta,G_n);
      }

	  Solution(double w1norm, std::vector<double>* theta, int G_n); //USED for in Algortihm1

    void UpdateCutValues(int update_index, std::vector<int>* x,
	  			     std::vector<double>* diff_weights,
	  			     double *objective); 

    void UpdateCutValues(int update_index) {
      UpdateCutValues(update_index, &assignments_, &diff_weights_, &weight_);
    }

    void PopulateFromAssignments();

    static bool ImprovesOver(double weight1, double weight2);

    bool ImprovingMove(int index) const {
	    return ImprovesOver(weight_ + diff_weights_[index], weight_);
    }

    // Print out the decision variables for a solutio
    void PrintSolution() const;

    // Getters
    double get_weight() const {  return weight_; }
    const std::vector<int>& get_assignments() const {  return assignments_; }
    const std::vector<double>& get_diff_weights() const {  return diff_weights_; }

    // Assignment operator
    Solution& operator=(const Solution &rhs);

    // Copy constructor
    Solution(const Solution& x);

    // Equals operator (checks if the solutions have identical assignments_)
    bool operator==(const Solution& other) const;

    // Not equals operator
    bool operator!=(const Solution& other) const {  return !(*this == other); }

    // >, <, >=, <=: compare based on weight_
    bool operator>(const Solution& other) const {
	return weight_ > other.weight_;
    }
    bool operator<(const Solution& other) const {
	return weight_ < other.weight_;
    }
    bool operator>=(const Solution& other) const {
	return weight_ >= other.weight_;
    }
    bool operator<=(const Solution& other) const {
	return weight_ <= other.weight_;
    }
    //Black box local minimisation algo from MQLIB

    void All1Swap(double tolerance);
    void All2Swap(double tolerance);

    double LoadNewTheta(const std::vector<double>& theta,
				       std::vector<double>* cos_theta,
				       std::vector<double>* sin_theta,
				       std::vector<double>* dH,
              int G_n,
              int G_m,
              std::vector<int> G_first,
              std::vector<int> G_second,
              std::vector<double> G_weight);

    // std::vector<double> MQLIB_m_inf(double w1norm, std::vector<double>* theta);

    //Performs ProcudureCUT as in paper
    // std::vector<int> ProcedureCUT(std::vector<double>* theta);

  protected:
    // Assignment of each node (-1/1 in MAXCUT, 0/1 in QUBO)
    std::vector<int> assignments_;
    // Objective value
    double weight_;
    // The number of nodes in the graph / variables in the problem
    int N_;

    std::vector<double> diff_weights_;
  };

class Burer2002{
 public:
  // Solves Max-Cut on edge-weighted graph mi, reporting each new best
  // solution to the reporter.
  //Perform algorithm 1 from paper
  Burer2002(int n, int m, int * f, int * s, double * w);
};

}
#endif
