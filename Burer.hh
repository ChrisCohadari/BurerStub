#ifndef HEURISTICS_MAXCUT_BURER_2002_H_
#define HEURISTICS_MAXCUT_BURER_2002_H_

#include <vector>
#include "BurerSolution.hh"

namespace Burer
{
    class Burer2002Solution : public MaxCutSolution {
     public:
      static Burer2002Solution Rank2Cut(double w1norm, std::vector<double>* theta) {
	return Burer2002Solution(w1norm, theta);
      }


      void All1Swap(double tolerance);
      void All2Swap(double tolerance);

     private:
      // Obtain an initial cut using gradient descent to minimize f(theta) from a
      // random initial theta, followed by Procedure-CUT. w1norm is the 1-norm of
      // the matrix of weights, and theta is passed as the initial theta vector
      // for the search and returned as the final one.
      Burer2002Solution(double w1norm, std::vector<double>* theta);

      // Given a new theta vector, compute the element-wise cosine and sine, as well
      // as weighted sum of the cosines of the angle differences between paired
      // nodes (returned through command-line arguments). Return objective.
      double LoadNewTheta(const std::vector<double>& theta,
			  std::vector<double>* cos_theta,
			  std::vector<double>* sin_theta, std::vector<double>* dH);
    };

    class Burer2002{
     public:
      // Solves Max-Cut on edge-weighted graph mi, reporting each new best
      // solution to the reporter.
      Burer2002(int n, int m, int * f, int * s, double * w);
    };
}
#endif
