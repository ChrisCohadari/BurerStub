#include <vector>

class Graph {
  int G_n;
  int G_m;
  int * G_first;
  int * G_second;
  double * G_weight;
  int * adjList;
  int * offset;

  //returns a tuple where the first vector is the adjacency list and the second is the list with the offset vector 
  //The neighbours of vertex i are the entries from offset[i] to (excluding) offset[i+1]
  std::vector<int> create_adjList(int G_n, int G_m, std::vector<int> G_first, std::vector<int> G_second);
};
