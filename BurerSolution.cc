#include "BurerGlobal.hh"
#include "BurerSolution.hh"
#include "BurerRandom.hh"

#include <stdlib.h>
#include <string.h>
#include <algorithm>
#include <math.h>
#include <iostream>
#include <limits>
#include <vector>

using namespace Burer;

BaseSolution::BaseSolution(int N, int init_assignment) :
  assignments_(N, init_assignment),
  weight_(0.0),
  N_(N) {}

BaseSolution::BaseSolution(const std::vector<int>& assignments, double weight) :
  assignments_(assignments),
  weight_(weight),
  N_(assignments.size()) {}

int BaseSolution::SymmetricDifference(const BaseSolution& other) const {
  int num_different = 0;
  for (int i=0; i < N_; ++i) {
    num_different += (assignments_[i] != other.assignments_[i]);
  }
  return num_different;
}

int BaseSolution::SymmetricDifference(const BaseSolution& other,
				      std::vector<int>* diff) const {
  diff->clear();
  for (int i=0; i < N_; ++i) {
    if (assignments_[i] != other.assignments_[i]) {
      diff->push_back(i);
    }
  }
  return diff->size();
}

int BaseSolution::SymmetricDifference(const BaseSolution& other,
				      std::vector<int>* diff,
				      std::vector<int>* common) const {
  diff->clear();
  common->clear();
  for (int i=0; i < N_; ++i) {
    if (assignments_[i] != other.assignments_[i]) {
      diff->push_back(i);
    } else {
      common->push_back(i);
    }
  }
  return diff->size();
}

bool BaseSolution::operator==(const BaseSolution& other) const {
  // If the weights are different, no need to check the assignments
  if (ImprovesOver(other) || other.ImprovesOver(*this)) {
    return false;
  }

  // Check the assignments
  for (int i=0; i < N_; ++i) {
    if (assignments_[i] != other.assignments_[i]) {
      return false;
    }
  }
  return true;
}

BaseSolution& BaseSolution::operator=(const BaseSolution &rhs) {
  assignments_ = rhs.assignments_;
  weight_ = rhs.weight_;
  return *this;
}

BaseSolution::BaseSolution(const BaseSolution& x) :
  assignments_(x.assignments_),
  weight_(x.weight_),
  N_(x.N_) {}

bool BaseSolution::ImprovesOver(double weight1, double weight2) {
  return weight1 - weight2 > 1e-6 &&
    (weight2 <= 0.0 || (weight1 / weight2 - 1.0) >= 1e-10);
}

void BaseSolution::PrintSolution() const {
  for (int i=0; i < N_; ++i) {
    if (i != 0) {
      std::cout << " ";
    }
    std::cout << assignments_[i];
  }
  std::cout << std::endl;
}

ExtendedSolution::ExtendedSolution(int N, int init_assignment) :
  BaseSolution(N, init_assignment),
  diff_weights_(N, 0.0) {}

void ExtendedSolution::AllBest1Swap(int startpos) {
  // Take all profitable 1-moves, taking the most profitable first.
  while (true) {
    double best_move = 0.0;
    int best_pos = -1;
    for (int i=startpos; i < N_; ++i) {
      if (diff_weights_[i] > best_move) {
	best_move = diff_weights_[i];
	best_pos = i;
      }
    }
    if (best_pos < 0 || !ImprovingMove(best_pos)) {
      // No more profitable moves
      break;
    }
    
    // Update the diff_weights_ variable and objective
    UpdateCutValues(best_pos);
  }
}

void ExtendedSolution::AllFirst1Swap(int startpos) {
  // Take all profitable 1-moves, taking the first one we find when scanning
  // the nodes/variables sequentially.
  bool move_made = true;
  while (move_made) {
    move_made = false;
    for (int i=startpos; i < N_; ++i) {
      if (ImprovingMove(i)) {
	UpdateCutValues(i);
	move_made = true;
	break;
      }
    }
  }
}

void ExtendedSolution::AllShuffle1Swap(int startpos) {
  // First, build a list of indices to consider and shuffle them
  std::vector<int> indices;
  for (int idx=startpos; idx < N_; ++idx) {
    indices.push_back(idx);
  }
  std::random_shuffle(indices.begin(), indices.end());

  // Take all profitable 1-moves, taking the first one we find when scanning the
  // nodes/variables sequentially.
  bool move_made = true;
  while (move_made) {
    move_made = false;
    for (auto iter=indices.begin(); iter != indices.end(); ++iter) {
      if (ImprovingMove(*iter)) {
	UpdateCutValues(*iter);
	move_made = true;
	break;
      }
    }
  }
}

double ExtendedSolution::DiffWeightStandardDeviation() const {
  // Compute the standard deviation in one pass:
  // http://www.strchr.com/standard_deviation_in_one_pass
  double sum = 0.0;
  double sq_sum = 0.0;
  for (int i=0; i < N_; ++i) {
    sum += diff_weights_[i];
    sq_sum += diff_weights_[i] * diff_weights_[i];
  }
  double mean = sum / N_;
  return sqrt(sq_sum / N_ - mean * mean);
}

ExtendedSolution& ExtendedSolution::operator=(const ExtendedSolution &rhs) {
  BaseSolution::operator=(rhs);
  diff_weights_ = rhs.diff_weights_;
  return *this;
}

ExtendedSolution::ExtendedSolution(const ExtendedSolution& x) :
  BaseSolution(x),
  diff_weights_(x.diff_weights_) {}


// Empty solution
MaxCutSimpleSolution::MaxCutSimpleSolution(
                                           int initVal) :
  BaseSolution(G_n, initVal)
{}

// Solution with provided assignments and weight
MaxCutSimpleSolution::MaxCutSimpleSolution(const std::vector<int>& assignments,
                                           double weight) :
  BaseSolution(assignments, weight)
{}

// Random solution (p=0.5 for each vertex)
MaxCutSimpleSolution::MaxCutSimpleSolution(
                                           int ignored1, int ignored2):
  BaseSolution(G_n, -1){
  // Random assignments in {-1, 1}
  for (int i=0; i < G_n; ++i) {
    assignments_[i] = 2 * Random::RandInt(0, 1) - 1;
  }

  // Obtain weight_
  PopulateFromAssignments();
}

// Random solution (p provided for each vertex)
MaxCutSimpleSolution::MaxCutSimpleSolution(const std::vector<double>& p) :
  BaseSolution(G_n, -1){
  // Random assignments in {-1, 1} using probabilities p of being set 1
  for (int i=0; i < G_n; ++i) {
    assignments_[i] = (Random::RandDouble() <= p[i]) ? 1 : -1;
  }

  // Obtain weights_
  PopulateFromAssignments();
}

MaxCutSimpleSolution& MaxCutSimpleSolution::operator=(const MaxCutSimpleSolution &rhs) {
  BaseSolution::operator=(rhs);
  // Can't copy over mi_ because it's a reference
  return *this;
}

MaxCutSimpleSolution::MaxCutSimpleSolution(const MaxCutSimpleSolution &x) :
  BaseSolution(x)
{}

void MaxCutSimpleSolution::PopulateFromAssignments() {
  weight_ = 0.0;

  for (int e(0) ; e < G_m ; ++e)
  {
      int t = G_first[e];
      int h = G_second[e];

    if (assignments_[t] != assignments_[h]) {
      // These two nodes are in different sets
      weight_ += G_weight[e];
    }
  }
}


// Empty solution
MaxCutSolution::MaxCutSolution(
			       int initVal) :
  ExtendedSolution(G_n, initVal)
{}

// Random solution (p=0.5 for each vertex)
MaxCutSolution::MaxCutSolution(int ignored1,
			       int ignored2):
  ExtendedSolution(G_n, -1)
{
  // Random assignments in {-1, 1}
  for (int i=0; i < G_n; ++i) {
    assignments_[i] = 2 * Random::RandInt(0, 1) - 1;
  }

  // Obtain weight_ and diff_weights_
  PopulateFromAssignments();
}

// Random solution (p provided for each vertex)
MaxCutSolution::MaxCutSolution(
			       const std::vector<double>& p
			      ) :
  ExtendedSolution(G_n, -1) {
  // Random assignments in {-1, 1} using probabilities p of being set 1
  for (int i=0; i < G_n; ++i) {
    assignments_[i] = (Random::RandDouble() <= p[i]) ? 1 : -1;
  }

  // Obtain weights_ and diff_weights_
  PopulateFromAssignments();
}

// Initialize from vector of assignments
MaxCutSolution::MaxCutSolution(const std::vector<int>& assignments
                               ) :
  ExtendedSolution(G_n, -1) {
  // Copy over assignments
  assignments_ = assignments;

  // Obtain weights_ and diff_weights_
  PopulateFromAssignments();
}


// Does all the updates required when the node at update_index has its cut set
// inclusion flipped. x is updated to show this flipped inclusion, and the
// difference in weights to own set and other set are updated.
void MaxCutSolution::UpdateCutValues(int update_index, std::vector<int>* x,
				     std::vector<double>* diff_weights,
				     double *objective) const {
  *objective += (*diff_weights)[update_index];
  (*x)[update_index] = -(*x)[update_index];
  (*diff_weights)[update_index] = -(*diff_weights)[update_index];

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
  for (int e(0) ; e < G_m ; ++e)
  {
     int i = G_first[e];
     int j = G_second[e];

     if (i == update_index)
     {
         (*diff_weights)[j] += 2.0 * (*x)[update_index] * (*x)[j] * G_weight[e];
     }
     else if (j == update_index)
     {
         (*diff_weights)[i] += 2.0 * (*x)[update_index] * (*x)[i] * G_weight[e];
     }
   }
}


void MaxCutSolution::AllBest2Swap(int startpos) {
  // Take all profitable 2-moves, taking the most profitable first.
  while (true) {
    double best_move = 0.0;
    int best_i = -1;
    int best_j = -1;

      for (int e(0) ; e < G_m ; ++e)
      {
         int i = G_first[e];
         int j = G_second[e];
         double w_ij = G_weight[e];

      double benefit = diff_weights_[i] + diff_weights_[j] -
	2.0 * assignments_[i] * assignments_[j] * w_ij;
      // Only note move if it improves on best found so far and i and j
      // are both large enough. Note that there's limited efficiency penalty
      // compared to writing two separate functions (one without startpos and
      // one with it) because the last two conditions are only checked
      // the few times when we've found an improving 2-move. The same argument
      // also applies for AllFirst2Swap() below.
      if (benefit > best_move && i >= startpos && j >= startpos) {
	best_move = benefit;
	best_i = i;
	best_j = j;
      }
    }
    if (best_i < 0 || !ImprovingMove(best_move)) {
      // No more profitable moves
      break;
    }

    // Update the diff_weights_ variables and objective
    UpdateCutValues(best_i);
    UpdateCutValues(best_j);
  }
}

void MaxCutSolution::AllFirst2Swap(int startpos) {
  // Take all profitable 2-moves, taking the first one we find when scanning
  // the edges sequentially.
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
      if (ImprovingMove(benefit) && i >= startpos && j >= startpos) {
	UpdateCutValues(i);
	UpdateCutValues(j);
	move_made = true;
	break;
      }
    }
  }
}

MaxCutSolution& MaxCutSolution::operator=(const MaxCutSolution &rhs) {
  ExtendedSolution::operator=(rhs);
  // Can't copy over mi_ because it's a reference
  return *this;
}

MaxCutSolution::MaxCutSolution(const MaxCutSolution &x) :
  ExtendedSolution(x)
{}

void MaxCutSolution::PopulateFromAssignments() {
  weight_ = 0.0;
  diff_weights_.assign(N_, 0.0);

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

// Empty solution
FirstFixedMaxCutSolution::FirstFixedMaxCutSolution(
						   int fixedVal) :
  MaxCutSolution(fixedVal),
  fixedVal_(fixedVal) {}

// Random solution
FirstFixedMaxCutSolution::FirstFixedMaxCutSolution(
						   int fixedVal, int ignored2):
  MaxCutSolution(fixedVal),
  fixedVal_(fixedVal) {
  // Random assignments in {-1, 1} (don't assign the first one; already handled
  // by the MaxCutSolution initializer).
  for (int i=1; i < G_n; ++i) {
    assignments_[i] = 2 * Random::RandInt(0, 1) - 1;
  }

  // Obtain weight_ and diff_weights_
  PopulateFromAssignments();
}

// Random solution from vertex probabilities
FirstFixedMaxCutSolution::FirstFixedMaxCutSolution(
						   const std::vector<double>& p,
						   int fixedVal) :
  MaxCutSolution(fixedVal),
  fixedVal_(fixedVal) {
  // Random assignments in {-1, 1} using probabilities p of being set 1 (don't
  // assign the first one; already handled by MaxCutSolution initializer)
  for (int i=1; i < G_n; ++i) {
    assignments_[i] = (Random::RandDouble() <= p[i]) ? 1 : -1;
  }

  // Obtain weights_ and diff_weights_
  PopulateFromAssignments();
}

// Does all the updates required when the node at update_index has its cut set
// inclusion flipped. x is updated to show this flipped inclusion, and the
// difference in weights to own set and other set are updated.
void FirstFixedMaxCutSolution::UpdateCutValues(int update_index,
					       std::vector<int>* x,
					       std::vector<double>* diff_wts,
					       double *objective) const {
  if (update_index == 0) {
    std::cout << "Error: flipping first index of a FirstFixedMaxCutSolution" <<
      std::endl;
    exit(0);
  }
  MaxCutSolution::UpdateCutValues(update_index, x, diff_wts, objective);
}

FirstFixedMaxCutSolution&
FirstFixedMaxCutSolution::operator=(const FirstFixedMaxCutSolution &rhs) {
  MaxCutSolution::operator=(rhs);
  fixedVal_ = rhs.fixedVal_;
  return *this;
}

FirstFixedMaxCutSolution::FirstFixedMaxCutSolution(const
						   FirstFixedMaxCutSolution &x)
  : MaxCutSolution(x),
    fixedVal_(x.fixedVal_) {}

void FirstFixedMaxCutSolution::PopulateFromAssignments() {
  if (assignments_[0] != fixedVal_) {
    std::cout << "Error: wrong start val in PopulateFromAssignments" <<
      std::endl;
    exit(0);
  }
  MaxCutSolution::PopulateFromAssignments();
}

#ifdef false
// Copy over a MaxCutPartialSolution to be a full MaxCutSolution
MaxCutSolution::MaxCutSolution(const MaxCutPartialSolution& x) :
  ExtendedSolution(x.get_mi().get_size(), 0),
  mi_(x.get_mi()),
  heuristic_(x.get_heuristic()) {
  if (x.get_num_unassigned() > 0) {
    std::cout << "Cannot copy over MaxCutPartialSolution with unassigned nodes"
	      << std::endl;
    exit(0);
  }

  // Copy over assignments and objective value from partial solution
  assignments_.assign(x.get_assignments().begin(), x.get_assignments().end());
  weight_ = x.get_weight();

  // Assign diff_weights_ from gainS and gainNS of the partial solution
  for (int i=0; i < N_; ++i) {
    if (assignments_[i] == 1) {
      diff_weights_[i] = x.get_gainNS()[i];
    } else {
      diff_weights_[i] = x.get_gainS()[i];
    }
  }
}
#endif
