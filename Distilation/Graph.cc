#include <vector>
#include <iostream>

#include "Graph.hh"

//TODO: This can probably be implemented without the cumCounts vector
std::pair<std::vector<int>, std::vector<int>> create_adjList(int G_n, int G_m, std::vector<int> G_first, std::vector<int> G_second){

  std::vector<int> adjList(2*G_m, -1);
  std::vector<int> offset(G_n, 0);

  //acc to .mc doc every edge is listed once, thus i think that this counts correctly
  //count occurences
  for(int i = 0; i < G_m; i++){
    offset[G_first[i]]++;
    offset[G_second[i]]++;
  }

  std::cout << "offset\n";
  for(int i =0; i < offset.size(); i++)
    std::cout << offset[i] << " ";
  std::cout << "\n";


  //an implementation without cumCounts would be nice
  // for(int i = 1; i < G_n; i++)
  //   offset[i] += offset[i-1];

  //cummulative counts
  std::vector<int> cumCounts(G_n, offset[0]);
  // cumCounts[0] = offset[0];
  for(int i = 1 ; i< G_n; i++)
    cumCounts[i] = offset[i] + cumCounts[i-1];

  std::cout << "cumCounts\n";
  for(int i =0; i < cumCounts.size(); i++)
    std::cout << cumCounts[i] << " ";
  std::cout << "\n";

  for(int i = 0; i<G_n; i++)
    offset[i] = cumCounts[i] - offset[i];

  int pos;
  for(int i = 0; i < G_m; i++){
    //insert neigh of first index
    pos = offset[G_first[i]];
    adjList[pos] = G_second[i];
    offset[G_first[i]]++;

    //insert neigh of second index
    pos = offset[G_second[i]];
    adjList[pos] = G_first[i];
    offset[G_second[i]]++;
  }

  //Note that now the i-th entry of indexvector contains the startindex of the neighbours of vertex i+1
  //the number of neighbours of vertex 
  for(int i = G_n -1; i > 0; i--)
  {
    std::cout << "i = " << i << "\n";
    offset[i] = offset[i-1];
  }
  offset[0] = 0;

  return std::pair<std::vector<int>,std::vector<int>>(adjList,offset);
}


// int main(int argc, char *argv[])
// {
//         int n = 4;
//         int m = 4;
//         std::vector<int> frst = {0,0,0,1};
//         std::vector<int> scnd = {1,2,3,2};
//
//         std::pair<std::vector<int>,std::vector<int>> adjList = create_adjList(n,m,frst,scnd);
//         // if (ok)
//         // {
//         //     Burer::Burer2002 * b = new Burer::Burer2002(n, m, frst, scnd, wght);
//         //
//         //     delete b;
//         //
//         //     for (int i(0) ; i < m ; ++i)
//         //         delete[] wght_str[i];
//         // }
//         for( int i = 0 ; i < 2*m; i++)
//           std::cout << adjList.first[i] << " ";
//         std::cout << "\n";
//
//         for( int i = 0 ; i < n; i++)
//           std::cout << adjList.second[i] << " ";
//         std::cout << "\n";
//
// }
