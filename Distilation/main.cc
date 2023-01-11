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

        if (ok)
        {
            Burer2002 * b = new Burer2002(n, m, frst, scnd, wght);

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
    else if (BQ)
    {
        int n, m;
        int * frst(nullptr);
        int * scnd(nullptr);
        double * nnz(nullptr);
        char ** nnz_str(nullptr);
        char * cmmnts(nullptr);

        bool ok = parseBQ_extern(filename, n, m, frst, scnd, nnz, nnz_str, cmmnts);

        if (ok)
        {
            for (int e(0) ; e < m ; ++e)
            {
                for (int f(e+1) ; f < m ; ++f)
                {
                    if ((frst[e] == scnd[f]) && (scnd[e] == frst[f]))
                    {
                        nnz[e] += nnz[f];
                        frst[f] = frst[m - 1];
                        scnd[f] = scnd[m - 1];
                        nnz[f] = nnz[m - 1];
                        --m;
                    }
                }
            }

            int mMC = m;

            double * sune = new double[n];

            for (int i(0) ; i < n ; ++i)
                sune[i] = 0.0;

            for (int e(0) ; e < m ; ++e)
            {
                if (frst[e] != scnd[e])
                {
                    double half = nnz[e] / 2.0;
                   sune[frst[e]] -= half;
                   sune[scnd[e]] -= half;
                }
                else
                {
                    sune[frst[e]] -= nnz[e];
                    --mMC;
                }
            }

            // If there is NO sun-edge (can happen e.g. when converting back and forth),
            // there is no need to add an isolated sun node

            int sune1(0);

            for (int i(0) ; i < n ; ++i)
            {
                if (sune[i] != 0.0)
                {
                    sune1 = 1;
                    ++mMC;
                }
            }

            int nMC = n + sune1;

            int * frstMC = new int[mMC];
            int * scndMC = new int[mMC];
            double * wghtMC = new double[mMC];

            int pos(0);
            for (int i(0) ; i < n ; ++i)
            {
                if (sune[i] != 0.0)
                {
                    frstMC[pos] = 0;
                    scndMC[pos] = i + sune1;
                    wghtMC[pos] = sune[i];
                    ++pos;
                }
            }

            delete[] sune;

            for (int e(0) ; e < m ; ++e)
            {
                if (frst[e] != scnd[e])
                {
                    frstMC[pos] = frst[e] + sune1;
                    scndMC[pos] = scnd[e] + sune1;
                    wghtMC[pos] = nnz[e] / 2.0;
                    ++pos;
                }
            }

            Burer2002 * b = new Burer2002(nMC, mMC, frstMC, scndMC, wghtMC);

            delete b;

            delete[] frstMC;
            delete[] scndMC;
            delete[] wghtMC;

            for (int i(0) ; i < m ; ++i)
                delete[] nnz_str[i];
        }

        delete[] frst;
        delete[] scnd;
        delete[] nnz;
        delete[] nnz_str;
        delete[] cmmnts;
    }
}
