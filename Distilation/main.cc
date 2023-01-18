/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/* Sven Mallach (2022) */

#include <vector>
#include <common/parsedata.hh>
#include <common/lesspain.hh>
#include "BurerStub.hh"
// #include "BurerGlobal.hh"

#include <iostream>
void print_usage()
{
   fprintf(stderr, "Usage: ./burer <file in mc or bq format>\n");
}

//Implemntation of adjacency list avoiding a graph class
void create_adjList(int G_n, int G_m, int * G_first, int * G_second, std::vector<int> &adjList, std::vector<int> &offset, std::vector<int> &edgeCorrForAdjList){

  // std::vector<int> adjList(2*G_m, -1);
  // std::vector<int> offset(G_n, 0);

  //acc to .mc doc every edge is listed once, thus i think that this counts correctly
  //count occurences
  for(int i = 0; i < G_m; i++){
    offset[G_first[i]]++;
    offset[G_second[i]]++;
  }

  // std::cout << "offset\n";
  // for(int i =0; i < offset.size(); i++)
  //   std::cout << offset[i] << " ";
  // std::cout << "\n";


  //an implementation without cumCounts would be nice
  // for(int i = 1; i < G_n; i++)
  //   offset[i] += offset[i-1];

  //cummulative counts
  std::vector<int> cumCounts(G_n, offset[0]);
  // cumCounts[0] = offset[0];
  for(int i = 1 ; i< G_n; i++)
    cumCounts[i] = offset[i] + cumCounts[i-1];

  // std::cout << "cumCounts\n";
  // for(int i =0; i < cumCounts.size(); i++)
  //   std::cout << cumCounts[i] << " ";
  // std::cout << "\n";

  for(int i = 0; i<G_n; i++)
    offset[i] = cumCounts[i] - offset[i];

  int pos;
  for(int i = 0; i < G_m; i++){
    //insert neigh of first index
    pos = offset[G_first[i]];
    adjList[pos] = G_second[i];
    offset[G_first[i]]++;
    edgeCorrForAdjList[pos]=i;

    //insert neigh of second index
    pos = offset[G_second[i]];
    adjList[pos] = G_first[i];
    offset[G_second[i]]++;
    edgeCorrForAdjList[pos]=i;
  }

  //Note that now the i-th entry of indexvector contains the startindex of the neighbours of vertex i+1
  //the number of neighbours of vertex 
  for(int i = G_n -1; i > 0; i--)
  {
    // std::cout << "i = " << i << "\n";
    offset[i] = offset[i-1];
  }
  offset[0] = 0;
  // offset[G_n] = G_n;//debugging
}


int main(int argc, char *argv[])
{

    if (argc != 2)
    {
        print_usage();
        return -1;
    }

    char * const filename = argv[1];

    const int fnlen(strlen(filename));

    if (fnlen < 3)
    {
        print_usage();
        return -1;
    }

    const bool MC(0 == strncmp(filename + fnlen-3, ".mc", 3));
    const bool BQ(0 == strncmp(filename + fnlen-3, ".bq", 3));

    if (! (MC || BQ))
    {
        print_usage();
        return -1;
    }

    if (! check_readability(filename))
    {
        fprintf(stderr, "Could not open file %s for reading. Exiting.\n", filename);
        return -1;
    }

    if (MC)
    {
        int n, m;
        int * frst(nullptr);
        int * scnd(nullptr);
        double * wght(nullptr);
        char ** wght_str(nullptr);
        char * cmmnts(nullptr);

        bool ok = parseMC_extern(filename, n, m, frst, scnd, wght, wght_str, cmmnts);

        std::vector<int> adjList(2*m,-1);
        std::vector<int> offset(n+1,0);
        std::vector<int> edgeCorrForAdjList(2*m,-1);
        
        // create_adjList(n,m,frst,scnd, adjList, offset, edgeCorrForAdjList);


        if (ok)
        {
            Burer2002 * b = new Burer2002(n, m, frst, scnd, wght, adjList, offset, edgeCorrForAdjList);

            delete b;

            for (int i(0) ; i < m ; ++i)
                delete[] wght_str[i];
        }

        delete[] frst;
        delete[] scnd;
        delete[] wght;
        delete[] wght_str;
        delete[] cmmnts;
    }
    // else if (BQ)
    // {
    //     int n, m;
    //     int * frst(nullptr);
    //     int * scnd(nullptr);
    //     double * nnz(nullptr);
    //     char ** nnz_str(nullptr);
    //     char * cmmnts(nullptr);
    //
    //     bool ok = parseBQ_extern(filename, n, m, frst, scnd, nnz, nnz_str, cmmnts);
    //
    //     if (ok)
    //     {
    //         for (int e(0) ; e < m ; ++e)
    //         {
    //             for (int f(e+1) ; f < m ; ++f)
    //             {
    //                 if ((frst[e] == scnd[f]) && (scnd[e] == frst[f]))
    //                 {
    //                     nnz[e] += nnz[f];
    //                     frst[f] = frst[m - 1];
    //                     scnd[f] = scnd[m - 1];
    //                     nnz[f] = nnz[m - 1];
    //                     --m;
    //                 }
    //             }
    //         }
    //
    //         int mMC = m;
    //
    //         double * sune = new double[n];
    //
    //         for (int i(0) ; i < n ; ++i)
    //             sune[i] = 0.0;
    //
    //         for (int e(0) ; e < m ; ++e)
    //         {
    //             if (frst[e] != scnd[e])
    //             {
    //                 double half = nnz[e] / 2.0;
    //                sune[frst[e]] -= half;
    //                sune[scnd[e]] -= half;
    //             }
    //             else
    //             {
    //                 sune[frst[e]] -= nnz[e];
    //                 --mMC;
    //             }
    //         }
    //
    //         // If there is NO sun-edge (can happen e.g. when converting back and forth),
    //         // there is no need to add an isolated sun node
    //
    //         int sune1(0);
    //
    //         for (int i(0) ; i < n ; ++i)
    //         {
    //             if (sune[i] != 0.0)
    //             {
    //                 sune1 = 1;
    //                 ++mMC;
    //             }
    //         }
    //
    //         int nMC = n + sune1;
    //
    //         int * frstMC = new int[mMC];
    //         int * scndMC = new int[mMC];
    //         double * wghtMC = new double[mMC];
    //
    //         int pos(0);
    //         for (int i(0) ; i < n ; ++i)
    //         {
    //             if (sune[i] != 0.0)
    //             {
    //                 frstMC[pos] = 0;
    //                 scndMC[pos] = i + sune1;
    //                 wghtMC[pos] = sune[i];
    //                 ++pos;
    //             }
    //         }
    //
    //         delete[] sune;
    //
    //         for (int e(0) ; e < m ; ++e)
    //         {
    //             if (frst[e] != scnd[e])
    //             {
    //                 frstMC[pos] = frst[e] + sune1;
    //                 scndMC[pos] = scnd[e] + sune1;
    //                 wghtMC[pos] = nnz[e] / 2.0;
    //                 ++pos;
    //             }
    //         }
    //
    //         Burer2002 * b = new Burer2002(nMC, mMC, frstMC, scndMC, wghtMC);
    //
    //         delete b;
    //
    //         delete[] frstMC;
    //         delete[] scndMC;
    //         delete[] wghtMC;
    //
    //         for (int i(0) ; i < m ; ++i)
    //             delete[] nnz_str[i];
    //     }
    //
    //     delete[] frst;
    //     delete[] scnd;
    //     delete[] nnz;
    //     delete[] nnz_str;
    //     delete[] cmmnts;
    // }
}
