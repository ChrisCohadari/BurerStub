#include <math.h>
#include <string.h>
#include <algorithm>
#include <iostream>
#include <limits>
#include <vector>

#include "BurerStub.hh"


//Definiton of the Global Vars which were in BurerGlobal.hh are now in the block below

int G_n;
int G_m;
std::vector<int> G_first;
std::vector<int> G_second;
std::vector<double> G_weight;
std::vector<int> adjList;
std::vector<int> offset;
std::vector<int> edgeCorrForAdjList;

//Beginning of functions concerning BurerStub.hh
Solution::Solution(int N, int init_assignment) :
  assignments_(N, init_assignment),
  weight_(0.0),
  N_(N) {}

Solution::Solution(const std::vector<int>& assignments, double weight) :
  assignments_(assignments),
  weight_(weight),
  N_(assignments.size()) {}

Solution::Solution(double w1norm, std::vector<double>* theta, int G_n, std::vector<int> &adjList, std::vector<int> &offset, std::vector<int> &edgeCorrForAdjList){
  assignments_.resize(G_n,1);
  N_ = G_n;
  // ***Parameters for construction method
  // Stopping value for relative change in f(theta) when optimizing f
  const double f_tolerance = 1e-4;
  // Stopping value for normalized value of g(theta) when optimizing f
  const double g_tolerance = 1e-4;
  // Maximum number of backtracks before termination
  const int maxback = 10;
  // Backtrack multiplier
  const double tau = 0.5;
  // The "beta" value in the Armijo line-search (required fraction of 1st order
  // reduction) -- see http://people.maths.ox.ac.uk/hauser/hauser_lecture2.pdf
  // slide 9
  const double gamma = 0.01;
  // The initial step size in the line search
  const double alpha_init = 1.0;
  // The maximum step size in the line search
  const double max_alpha = 4.0;
  // Maximum number of gradient descent steps when optimizing f(theta)
  const int max_opt_iter = 200;

  // Minimize f(theta) using gradient descent with backtracking Armijo
  // line-search.

  // *** Setup variables
  // The backtracking alpha value for each iteration's line search
  double bt_alpha = alpha_init;
  // negative of cosine differences incident to each node
  std::vector<double> dH(N_);
  // element-wise cosine of the theta vector
  std::vector<double> cos_theta(N_);
  // element-wise sine of the theta vector
  std::vector<double> sin_theta(N_);
  // f: the current objective value of nonlinear optimization problem
  double f = LoadNewTheta(*theta, &cos_theta, &sin_theta, &dH, G_n, G_m, G_first, G_second, G_weight);


  for (int opt_iter = 0; opt_iter < max_opt_iter; ++opt_iter) {
    // Compute the current gradient
    std::vector<double> g(N_, 0.0);  // Gradient of f(theta)

      for (int e(0) ; e < G_m ; ++e)
      {
         int i = G_first[e];
         int j = G_second[e];
         double w_ij = G_weight[e];
      double scaled_sin_diff =
	w_ij * (sin_theta[j]*cos_theta[i] - cos_theta[j]*sin_theta[i]);
      g[i] += scaled_sin_diff;
      // sin(theta[i] - theta[j]) = -sin(theta[j] - theta[i])
      g[j] -= scaled_sin_diff;
    }
    
    // Stop the gradient descent if the gradient is too close to 0
    double norm_gradient = 0.0;
    for (int ct=0; ct < N_; ++ct) {
      norm_gradient += g[ct] * g[ct];
    }
    if (norm_gradient / w1norm < g_tolerance) {
      break;
    }
    
    // Determine the descent direction via scaling; we determined this behavior
    // by actually looking at the circut code as this is not mentioned in the
    // Burer2002 paper
    std::vector<double> desc(N_);  // Descent direction
    double g_times_desc = 0.0;
    double divisor = 1.0;
    for (int ct=0; ct < N_; ++ct) {
      divisor = std::max(divisor, dH[ct]);
    }
    for (int ct=0; ct < N_; ++ct) {
      desc[ct] = -g[ct] / divisor;
      g_times_desc += desc[ct] * g[ct];
    }
    
    // Use backtracking Armijo line-search to determine a good step size
    std::vector<double> new_theta(N_);
    int numback;
    double recent_f = -1.0;
    for (numback=1; numback <= maxback; ++numback) {
      // Compute the new theta with step size alpha
      for (int ct=0; ct < N_; ++ct) {
	new_theta[ct] = (*theta)[ct] + bt_alpha * desc[ct];
      }
      
      // Update the cosine and sine vectors, as well as dH and the
      // objective. Exit on Armijo condition.
      recent_f = LoadNewTheta(new_theta, &cos_theta, &sin_theta, &dH, G_n, G_m, G_first, G_second, G_weight);
      if (recent_f <= f + gamma * bt_alpha * g_times_desc) {
	break;
      }
      
      // Exponential scaling on step size
      bt_alpha *= tau;
    }
    double f_prev = f;
    f = recent_f;
    *theta = new_theta;  // Copy over the current theta
    
    // Stop the gradient descent if there was not enough change in the
    // objective value
    double rel_change = fabs(f - f_prev) / (1.0 + fabs(f_prev));
    if (rel_change < f_tolerance) {
      break;
    }
    
    // Step size update lifted from circut code; not mentioned in paper...
    if (numback <= 2) {
      if (max_alpha < 2.0 * bt_alpha) {
	bt_alpha = max_alpha;
      } else {
	bt_alpha *= 2.0;
      }
    }
  }
  
  // Modulo the angles to be between 0 and 2*PI, and add to a vector of
  // index/angle pairs. Sort wrt angle.
  std::vector<std::pair<double, int> > angles;
  for (int ct=0; ct < N_; ++ct) {
    (*theta)[ct] -= 2 * M_PI * floor((*theta)[ct] / (2*M_PI));
    angles.push_back(std::pair<double, int>((*theta)[ct], ct)); //Q: Do we want to pushback every time? Vector has fixed length N_
    //OPt. would be to init size N angles vector 
    //then no new memory will need to be allocated
  }
  angles.push_back(std::pair<double, int>(2*M_PI, N_)); //Theta_n+1 from paper
  std::sort(angles.begin(), angles.end()); //sorts by first entry (of pair) by default

  // Determine initial set inclusion, and setup variables to be updated
  // later.
  for (int ct=0; ct < N_; ++ct) {
    assignments_[ct] = -1; 
  }
  std::vector<std::pair<double, int> >::iterator first_it = angles.begin();
  std::vector<std::pair<double, int> >::iterator second_it =  angles.begin();
  while (second_it->first <= M_PI) {
    assignments_[second_it->second] = 1;  // The first set has angles in [0, pi]
    ++second_it;
  }

  // Now, first_it points at the beginning, and second_it points at the first
  // node with theta > pi.
  // Correspondence to paper:
  // i is first_it and j is second_it

  // Fill in diff_weights_ and weight_ from assignments_, and copy them over to
  // the temporary solution we'll be updating in the search for a better
  // solution.
  PopulateFromAssignments();
  double curr_weight = weight_; //This is \gamma of the cut given by \alpha = 0
  //curr_weight can be thought of as \Gamma in the paper
  std::vector<int> curr_assignments(assignments_);
  std::vector<double> curr_diff_weights(diff_weights_);

  //Until this point we effectively did the computation for \alpha = 0
  //Next we need to decide on the new \alpha
  //We are effictively in the first iteration of the while loop with \alpha = 0 about to start step 3.
  
  //Q: How does this work? update_index == N_ must somehow correspond to alpha > pi
  //But why 

  // Compute the optimal cut by exhaustively searching through the possible cuts
  while (1) {
    // Update the cut angle
    int update_index;
    //note: first_it->first corresponds to \theta_i in paper
    //second_it->first corresponds to \theta_j in paper
    if (first_it->first <= second_it->first - M_PI) { 
      update_index = first_it->second;  // We will exclude this guy.
      ++first_it;
    } else {
      update_index = second_it->second;  // We will include this guy
      ++second_it;
    }
    if (update_index == N_) {
      break;  // We have exhausted all the possible cuts
    }
    
    // See if this is a new best cut; if so, copy as the current solution
    UpdateCutValues(update_index, &curr_assignments, &curr_diff_weights,
		    &curr_weight, adjList, offset, edgeCorrForAdjList);
    if (ImprovesOver(curr_weight, weight_)) {
      weight_ = curr_weight;
      assignments_ = curr_assignments;
      diff_weights_ = curr_diff_weights;
    }
  }
}

void Solution::UpdateCutValues(int update_index, std::vector<int>* x,
				     std::vector<double>* diff_weights,
				     double *objective,
             const std::vector<int> &adjList,
             const std::vector<int> &offset,
             const std::vector<int> &edgeCorrForAdjList) 
{
  *objective += (*diff_weights)[update_index];
  (*x)[update_index] = -(*x)[update_index]; //node ass. to update_index changes ass. set
  (*diff_weights)[update_index] = -(*diff_weights)[update_index]; //This is what would happen in PopulateFromAssignments

  // Iterate the set of all neighbors for node update_index
#ifdef false
  for (auto iter = mi_.get_edges_begin(update_index);
       iter != mi_.get_edges_end(update_index); ++iter) {
    int j = iter->first;  // There is an edge (update_index, j) in the graph
    double w_ij = iter->second;  // The edge has weight w_ij
    (*diff_weights)[j] += 2.0 * (*x)[update_index] * (*x)[j] * w_ij;
  }
#endif

//  printf("update index is %d\n", update_index);

// This is the only loop that gets ``inefficient'' because we do not have
// adjacency lists in this implementation
  // for (int e(0) ; e < G_m ; ++e)
  // {
  //    int i = G_first[e];
  //    int j = G_second[e];
  //
  //    //gamma changes only around the neighbors of the node update_index
  //    //This decreases the number of computations
  //    if (i == update_index)
  //    {
  //         //The factor 2 flips the sign of the summand
  //        (*diff_weights)[j] += 2.0 * (*x)[update_index] * (*x)[j] * G_weight[e]; 
  //    }
  //    else if (j == update_index)
  //    {
  //        (*diff_weights)[i] += 2.0 * (*x)[update_index] * (*x)[i] * G_weight[e];
  //    }
  // }
  for (int i = offset[update_index]  ; i < offset[update_index +1] ; i++)
  {
    //j corresponds to the "i-th" neighbour of update_index 
    int j = adjList[i];

    //gamma changes only around the neighbors of the node update_index
    //This decreases the number of computations
    //The factor 2 flips the sign of the summand
    (*diff_weights)[j] += 2.0 * (*x)[update_index] * (*x)[j] * G_weight[edgeCorrForAdjList[i]]; 
  }
}


//Given an assignment this function fills in the diff_weights_ vector and computes the weight_ (aka the value of \gamma(assignments_))
void Solution::PopulateFromAssignments() 
{
  weight_ = 0.0;
  diff_weights_.assign(N_, 0.0); //Note: assigns all entries 0.0

  for (int e(0) ; e < G_m ; ++e)
  {
     int i = G_first[e];
     int j = G_second[e];

    if (assignments_[i] == assignments_[j]) {
      // These two nodes are in the same set
      diff_weights_[i] += G_weight[e];
      diff_weights_[j] += G_weight[e];
    } else {
      // These two nodes are in different sets
      weight_ += G_weight[e];
      diff_weights_[i] -= G_weight[e];
      diff_weights_[j] -= G_weight[e];
    }
  }
}


bool Solution::ImprovesOver(double weight1, double weight2) { 
  return weight1 - weight2 > 1e-6 &&
    (weight2 <= 0.0 || (weight1 / weight2 - 1.0) >= 1e-10);
}

// Perform all first available 1-swaps. This is nearly identical to
// ExtendedSolution::AllFirst1Swap, except it has a tolerance.
void Solution::All1Swap(double tolerance) {
  bool move_made = true;
  while (move_made) {
    move_made = false;
    for (int i=0; i < N_; ++i) {
      if (diff_weights_[i] > tolerance && ImprovingMove(i)) {
	UpdateCutValues(i, adjList, offset, edgeCorrForAdjList);
	move_made = true;
	break;
      }
    }
  }
}

// Perform all first available 2-swaps. This is nearly identical to
// MaxCutSolution::AllFirst2Swap, except it has a tolerance.
void Solution::All2Swap(double tolerance) {
  bool move_made = true;
  while (move_made) {
    move_made = false;

      for (int e(0) ; e < G_m ; ++e)
      {
         int i = G_first[e];
         int j = G_second[e];
         double w_ij = G_weight[e];
      double benefit = diff_weights_[i] + diff_weights_[j] -
	2.0 * assignments_[i] * assignments_[j] * w_ij;
      if (benefit > tolerance && ImprovingMove(benefit)) {
	UpdateCutValues(i, adjList, offset, edgeCorrForAdjList);
	UpdateCutValues(j, adjList, offset, edgeCorrForAdjList);
	move_made = true;
	break;
      }
    }
  }
}

//MQLIB minimize f starts here
double Solution::LoadNewTheta(const std::vector<double>& theta,
				       std::vector<double>* cos_theta,
				       std::vector<double>* sin_theta,
				       std::vector<double>* dH,
              int G_n,
              int G_m,
              std::vector<int>& G_first,
              std::vector<int>& G_second,
              std::vector<double>& G_weight 
               ) {
  for (int i=0; i < N_; ++i) {
    (*cos_theta)[i] = cos(theta[i]);
    (*sin_theta)[i] = sin(theta[i]);
    (*dH)[i] = 0.0;
  }
  double objective = 0.0;

      for (int e(0) ; e < G_m ; ++e)
      {
         int i = G_first[e];
         int j = G_second[e];
         double w_ij = G_weight[e];

    double scaled_cos_diff = w_ij *
      ((*cos_theta)[i] * (*cos_theta)[j] + (*sin_theta)[i] * (*sin_theta)[j]);
    objective += scaled_cos_diff;
    (*dH)[i] -= scaled_cos_diff;
    (*dH)[j] -= scaled_cos_diff;
  }
  return objective;
}

// std::vector<double> Solution::MQLIB_m_inf(double w1norm, std::vector<double>* theta) 
// {
//   // ***Parameters for construction method
//   // Stopping value for relative change in f(theta) when optimizing f
//   const double f_tolerance = 1e-4;
//   // Stopping value for normalized value of g(theta) when optimizing f
//   const double g_tolerance = 1e-4;
//   // Maximum number of backtracks before termination
//   const int maxback = 10;
//   // Backtrack multiplier
//   const double tau = 0.5;
//   // The "beta" value in the Armijo line-search (required fraction of 1st order
//   // reduction) -- see http://people.maths.ox.ac.uk/hauser/hauser_lecture2.pdf
//   // slide 9
//   const double gamma = 0.01;
//   // The initial step size in the line search
//   const double alpha_init = 1.0;
//   // The maximum step size in the line search
//   const double max_alpha = 4.0;
//   // Maximum number of gradient descent steps when optimizing f(theta)
//   const int max_opt_iter = 200;
//
//   // Minimize f(theta) using gradient descent with backtracking Armijo
//   // line-search.
//
//   // *** Setup variables
//   // The backtracking alpha value for each iteration's line search
//   double bt_alpha = alpha_init;
//   // negative of cosine differences incident to each node
//   std::vector<double> dH(N_);
//   // element-wise cosine of the theta vector
//   std::vector<double> cos_theta(N_);
//   // element-wise sine of the theta vector
//   std::vector<double> sin_theta(N_);
//   // f: the current objective value of nonlinear optimization problem
//   double f = LoadNewTheta(*theta, &cos_theta, &sin_theta, &dH, G_n, G_m, G_first, G_second, G_weight);
//   for (int opt_iter = 0; opt_iter < max_opt_iter; ++opt_iter) {
//     // Compute the current gradient
//     std::vector<double> g(N_, 0.0);  // Gradient of f(theta)
//
//       for (int e(0) ; e < G_m ; ++e)
//       {
//          int i = G_first[e];
//          int j = G_second[e];
//          double w_ij = G_weight[e];
//       double scaled_sin_diff =
// 	w_ij * (sin_theta[j]*cos_theta[i] - cos_theta[j]*sin_theta[i]);
//       g[i] += scaled_sin_diff;
//       // sin(theta[i] - theta[j]) = -sin(theta[j] - theta[i])
//       g[j] -= scaled_sin_diff;
//     }
//     
//     // Stop the gradient descent if the gradient is too close to 0
//     double norm_gradient = 0.0;
//     for (int ct=0; ct < N_; ++ct) {
//       norm_gradient += g[ct] * g[ct];
//     }
//     if (norm_gradient / w1norm < g_tolerance) {
//       break;
//     }
//     
//     // Determine the descent direction via scaling; we determined this behavior
//     // by actually looking at the circut code as this is not mentioned in the
//     // Burer2002 paper
//     std::vector<double> desc(N_);  // Descent direction
//     double g_times_desc = 0.0;
//     double divisor = 1.0;
//     for (int ct=0; ct < N_; ++ct) {
//       divisor = std::max(divisor, dH[ct]);
//     }
//     for (int ct=0; ct < N_; ++ct) {
//       desc[ct] = -g[ct] / divisor;
//       g_times_desc += desc[ct] * g[ct];
//     }
//     
//     // Use backtracking Armijo line-search to determine a good step size
//     std::vector<double> new_theta(N_);
//     int numback;
//     double recent_f = -1.0;
//     for (numback=1; numback <= maxback; ++numback) {
//       // Compute the new theta with step size alpha
//       for (int ct=0; ct < N_; ++ct) {
// 	new_theta[ct] = (*theta)[ct] + bt_alpha * desc[ct];
//       }
//       
//       // Update the cosine and sine vectors, as well as dH and the
//       // objective. Exit on Armijo condition.
//       recent_f = LoadNewTheta(new_theta, &cos_theta, &sin_theta, &dH, G_n, G_m, G_first, G_second, G_weight);
//       if (recent_f <= f + gamma * bt_alpha * g_times_desc) {
// 	break;
//       }
//       
//       // Exponential scaling on step size
//       bt_alpha *= tau;
//     }
//     double f_prev = f;
//     f = recent_f;
//     *theta = new_theta;  // Copy over the current theta
//     
//     // Stop the gradient descent if there was not enough change in the
//     // objective value
//     double rel_change = fabs(f - f_prev) / (1.0 + fabs(f_prev));
//     if (rel_change < f_tolerance) {
//       break;
//     }
//     
//     // Step size update lifted from circut code; not mentioned in paper...
//     if (numback <= 2) {
//       if (max_alpha < 2.0 * bt_alpha) {
// 	bt_alpha = max_alpha;
//       } else {
// 	bt_alpha *= 2.0;
//       }
//     }
//   }
//   
//   return *theta;
// }

//MQLIB minimize f ends here

// std::vector<int> Solution::ProcedureCUT(std::vector<double>* theta)
// {
//   // Modulo the angles to be between 0 and 2*PI, and add to a vector of
//   // index/angle pairs. Sort wrt angle.
//   std::vector<std::pair<double, int> > angles;
//   for (int ct=0; ct < N_; ++ct) {
//     (*theta)[ct] -= 2 * M_PI * floor((*theta)[ct] / (2*M_PI));
//     angles.push_back(std::pair<double, int>((*theta)[ct], ct)); //Q: Do we want to pushback every time? Vector has fixed length N_
//     //OPt. would be to init size N angles vector 
//     //then no new memory will need to be allocated
//   }
//   angles.push_back(std::pair<double, int>(2*M_PI, N_)); //Theta_n+1 from paper
//   std::sort(angles.begin(), angles.end()); //sorts by first entry (of pair) by default
//
//   // Determine initial set inclusion, and setup variables to be updated
//   // later.
//   for (int ct=0; ct < N_; ++ct) {
//     assignments_[ct] = -1;
//   }
//   std::vector<std::pair<double, int> >::iterator first_it = angles.begin();
//   std::vector<std::pair<double, int> >::iterator second_it =  angles.begin();
//   while (second_it->first <= M_PI) {
//     assignments_[second_it->second] = 1;  // The first set has angles in [0, pi]
//     ++second_it;
//   }
//
//   // Now, first_it points at the beginning, and second_it points at the first
//   // node with theta > pi.
//   // Correspondence to paper:
//   // i is first_it and j is second_it
//
//   // Fill in diff_weights_ and weight_ from assignments_, and copy them over to
//   // the temporary solution we'll be updating in the search for a better
//   // solution.
//   PopulateFromAssignments();
//   double curr_weight = weight_; //This is \gamma of the cut given by \alpha = 0
//   //curr_weight can be thought of as \Gamma in the paper
//   std::vector<int> curr_assignments(assignments_);
//   std::vector<double> curr_diff_weights(diff_weights_);
//
//   //Until this point we effectively did the computation for \alpha = 0
//   //Next we need to decide on the new \alpha
//   //We are effictively in the first iteration of the while loop with \alpha = 0 about to start step 3.
//   
//   //Q: How does this work? update_index == N_ must somehow correspond to alpha > pi
//   //But why 
//
//   // Compute the optimal cut by exhaustively searching through the possible cuts
//   while (1) {
//     // Update the cut angle
//     int update_index;
//     //note: first_it->first corresponds to \theta_i in paper
//     //second_it->first corresponds to \theta_j in paper
//     if (first_it->first <= second_it->first - M_PI) { 
//       update_index = first_it->second;  // We will exclude this guy.
//       ++first_it;
//     } else {
//       update_index = second_it->second;  // We will include this guy
//       ++second_it;
//     }
//     if (update_index == N_) {
//       break;  // We have exhausted all the possible cuts
//     }
//     
//     // See if this is a new best cut; if so, copy as the current solution
//     UpdateCutValues(update_index, &curr_assignments, &curr_diff_weights,
// 		    &curr_weight);
//     if (ImprovesOver(curr_weight, weight_)) {
//       weight_ = curr_weight;
//       assignments_ = curr_assignments;
//       diff_weights_ = curr_diff_weights;
//     }
//   }
//   return assignments_;
// }

Burer2002::Burer2002(int n, int m, int * f, int * s, double * w, std::vector<int> &adjList, std::vector<int> &offset, std::vector<int> &edgeCorrForAdjList)
{
    // GLobal Vars in BurerGlobal (!)

    G_n = n;
    G_m = m;

    for (int i = 0 ; i < m ; ++i)
    {
      G_first.push_back(f[i]);
      G_second.push_back(s[i]);
      G_weight.push_back(w[i]);
    }

    //debugging
    adjList = adjList;
    offset = offset;
    edgeCorrForAdjList = edgeCorrForAdjList;
  // Parameters
  // New: Number of outer iterations
  const int M = 3; //used to be 1
  // Number of permitted non-improving perturbations to optimal theta before
  // search is stopped. This was set to a few different values in the
  // computational results of burer2002 (0, 4, and 8 for torus set; 0 and 10
  // for G-set; and 10 and 50 for spin-glass dataset) so we'll use 50 since
  // it was the most common choice (and an intermediate value).
  const int N = 10;
  // Whether we perform greedy 1- and 2-moves after 
  const int local_search = 0;
  // Amount of improvement required to do a 1-move
  const double one_move_tolerance = 0.01;
  // Amount of improvement required to do a 2-move
  const double two_move_tolerance = 0.1;
  // Perturbation in the range [-perturbation*PI, perturbation*PI]
  const double perturbation = 0.2;

  // Compute 1-norm of vec(W), the vectorization of the weight matrix.
  double w1norm = 0.0;

  for (int e(0) ; e < G_m ; ++e)
  {
      //Note: Factor 2 comes form Weight matrix being symmetric
      w1norm += 2.0 * fabs(G_weight[e]); // Count both directions of edge
  }

  srand(time(NULL)); //setting seed to investigate different results


  //Keep track of best solution over all outer iterations (of M)
  Solution heur_sol(G_n, 1);

  std::vector<double> theta(G_n);

  for (int iter=0; iter < M ; ++iter) {  // Random restart until termination criterion
    // Generate random starting set of angles
    // std::vector<double> theta(G_n);
    for (int ct=0; ct < G_n; ++ct) {
      theta[ct] = (((double)rand()) / (((long)RAND_MAX)+1)) * 2 * M_PI; //RHS is double in [0,2PI)
    }


    // Algorithm 1 from Section 4
    double best_weight = -std::numeric_limits<double>::max();
    int k = 0;  // # of runs without an improvement in the best solution
    while (k <= N) {
      // Rank2Cut minimizes f(theta) and then uses Procedure-CUT to get the
      // best cut associated with this theta.
      
      Solution x = Solution::Rank2Cut(w1norm, &theta, G_n, adjList, offset, edgeCorrForAdjList);

      // Perform local searches (all 1- and 2-moves better than the tolerance)
      if (local_search) {
        x.All1Swap(one_move_tolerance);
        x.All2Swap(two_move_tolerance);
      }

      // The following is taken out of function: Termination only when iteration
      // counts exhausted
      // Check termination criterion (runtime on non-validation runs; iteration
      // count on validation runs).
      //if (!Report(x, iter)) {
      //	return;
      //      }
      
      // Update counter keeping track of iterations without improvement
      if (Solution::ImprovesOver(x.get_weight(),best_weight)) {
        best_weight = x.get_weight();
        k = 0;
        //If solution higher weight than all previous solutions then assigns as new candidate
        if(Solution::ImprovesOver(x.get_weight(),heur_sol.get_weight()) )
          heur_sol = x;
      } else {
        ++k;
      }

      printf("%d:%d: Got solution of weight %lf\n", iter, k, x.get_weight()); 

      // FILE *fp;
      // fp = fopen("output_fixed.log", "a");
      // fprintf(fp,"%d:%d: Got solution of weight %lf\n", iter, k, x.get_weight()); 
      // fclose(fp);

      // Perturb the angles associated with the current solution
      for (int ct=0; ct < G_n; ++ct) {
        theta[ct] = M_PI / 2.0 * (1.0 - x.get_assignments()[ct]) +
          perturbation * (2 * M_PI * (((double)rand()) / (((long)RAND_MAX)+1)) - M_PI);
      }
    }
  }

  printf("Heuristic solution of weight %lf\n", heur_sol.get_weight()); 

  // FILE *fp;
  // fp = fopen("output_fixed.log", "a");
  // fprintf(fp,"Heuristic solution of weight %lf\n", heur_sol.get_weight()); 
  // fclose(fp);

  // delete[] G_first;
  // delete[] G_second;
  // delete[] G_weight;

}

