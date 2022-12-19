#ifndef BURERSTUB_GUARD_PROCEDURECUT_HH
#define BURERSTUB_GUARD_PROCEDURECUT_HH 1

#include <vector>

namespace Burer
{

  class Solution1 {
  public:
    void UpdateCutValues(int update_index, std::vector<int>* x,
	  			     std::vector<double>* diff_weights,
	  			     double *objective); 

    void PopulateFromAssignments();

    bool ImprovesOver(double weight1, double weight2);

    //Black box local minimisation algo from MQLIB

    double LoadNewTheta(const std::vector<double>& theta,
				       std::vector<double>* cos_theta,
				       std::vector<double>* sin_theta,
				       std::vector<double>* dH);

    std::vector<double> MQLIB_m_inf(double w1norm, std::vector<double>* theta);

    //Performs ProcudureCUT as in paper
    std::vector<int> ProcedureCUT(std::vector<double>* theta);

    //Perform algorithm 1 from paper
    // std::vector<int> Algorithm1(int N, std::vector<double> theta);
  protected:
    // Assignment of each node (-1/1 in MAXCUT, 0/1 in QUBO)
    std::vector<int> assignments_;
    // Objective value
    double weight_;
    // The number of nodes in the graph / variables in the problem
    int N_;

    std::vector<double> diff_weights_;
  };
}
#endif
