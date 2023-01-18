/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/* Sven Mallach (2021) */

#ifndef MCBQ_PARSEDAT_HH
#define MCBQ_PARSEDAT_HH 1

#define MAX_WEIGHT_LENGTH 32

// #include <common/lesspain.hh>
#include "lesspain.hh" //changed this

/* Call this only for VALIDATED files! */

// Returns true if everything worked fine, false otherwise.
// Does allocation with new - if true is returned, calling delete[] is your task
bool parseBQ_extern(char * const filename, int & n, int & m, int *& frst, int *& scnd, double *& nnz, char **& nnz_str, char *& cmmnts)
{
    FILE * fp = fopen(filename, "r");

    buffer buf;
    buf.cap = 32;
    buf.line = (char *) malloc(buf.cap);

    if (buf.line == NULL)
    {
        fprintf(stderr, "Fatal error while allocating small buffer\n");
        return false;
    }

    buffer dummy;
    dummy.cap = 32;
    dummy.line = (char *) malloc(dummy.cap);

    if (dummy.line == NULL)
    {
        fprintf(stderr, "Fatal error while allocating small buffer\n");
        free(buf.line);
        return false;
    }

    buffer comments;
    comments.cap = 32;
    comments.line = (char *) malloc(comments.cap);

    long long int commentchars(0);

    if (comments.line == NULL)
    {
        fprintf(stderr, "Fatal error while allocating small buffer\n");
        free(buf.line);
        free(dummy.line);
        return false;
    }

    bool header(false);

    while (! header && buffer_line(fp, &buf))
    {
        if (buf.line[0] != '#')
        {
            if (strlen(buf.line) > dummy.cap)
            {
                dummy.cap = 1 + strlen(buf.line);

                void * rptr = realloc(dummy.line, dummy.cap);

                if (NULL == rptr)
                {
                    fprintf(stderr, "Fatal error while enlarging buffer\n");
                    free(buf.line);
                    free(dummy.line);
                    free(comments.line);
                    return false;
                }

                dummy.line = (char *) rptr;

            }

            int r = sscanf(buf.line, "%d %d%s", &n, &m, dummy.line);

            if (r >= 2 && n >= 1 && m >= 0)
                header = true;

            break;
        }
        else
        {
            if (commentchars + strlen(buf.line) > comments.cap)
            {
                comments.cap = commentchars + 1 + strlen(buf.line);

                void * rptr = realloc(comments.line, comments.cap);

                if (NULL == rptr)
                {
                    fprintf(stderr, "Fatal error while enlarging buffer\n");
                    free(buf.line);
                    free(dummy.line);
                    free(comments.line);
                    return false;
                }

                comments.line = (char *) rptr;
            }

            // omit terminating character
            memcpy(comments.line + commentchars, buf.line, strlen(buf.line));
            commentchars += strlen(buf.line);
        }

    }

    if (! header)
    {
        fclose(fp);
        free(buf.line);
        free(dummy.line);
        free(comments.line);

        fprintf(stdout, "This should not happen: First non-comment line has invalid form. Did you use a non-validated file?\n");

        return false;
    }

    frst = new int[1+m];
    scnd = new int[1+m];
    nnz = new double[1+m];
    nnz_str = new char*[1+m];

    int eread(0);

    // Detect malformed lines
    while (buffer_line(fp, &buf))
    {
        if (buf.line[0] != '#')
        {
            if (strlen(buf.line) > dummy.cap)
            {
                dummy.cap = 1 + strlen(buf.line);

                void * rptr = realloc(dummy.line, dummy.cap);

                if (NULL == rptr)
                {
                    fprintf(stderr, "Fatal error while enlarging buffer\n");
                    free(buf.line);
                    free(dummy.line);
                    delete[] frst;
                    delete[] scnd;
                    delete[] nnz;
                    delete[] nnz_str;
                    return false;
                }

                dummy.line = (char *) rptr;
            }

            char * str = new char[MAX_WEIGHT_LENGTH];

            int r = sscanf(buf.line, "%d %d %s%s", frst + eread, scnd + eread, str, dummy.line);

            if (r < 3)
            {
                delete[] str;
                break;
            }

            nnz_str[eread] = str;
            r = sscanf(str, "%lf", nnz + eread);

            if (r != 1)
            {
                delete[] str;
                break;
            }

            --frst[eread];
            --scnd[eread];

            ++eread;
        }
        else
        {
            if (commentchars + strlen(buf.line) > comments.cap)
            {
                comments.cap = commentchars + 1 + strlen(buf.line);

                void * rptr = realloc(comments.line, comments.cap);

                if (NULL == rptr)
                {
                    fprintf(stderr, "Fatal error while enlarging buffer\n");
                    free(buf.line);
                    free(dummy.line);
                    free(comments.line);
                    delete[] frst;
                    delete[] scnd;
                    delete[] nnz;
                    delete[] nnz_str;
                    return false;
                }

                comments.line = (char *) rptr;
            }

            // omit terminating character
            memcpy(comments.line + commentchars, buf.line, strlen(buf.line));
            commentchars += strlen(buf.line);
        }
    }

    fclose(fp);

    // Add termination to comments, should always
    // work as we always reserve one character of space
    // in addition

    comments.line[commentchars] = '\0';
    cmmnts = new char[commentchars + 1];
    strcpy(cmmnts, comments.line);

    free(buf.line);
    free(dummy.line);
    free(comments.line);

    if (eread != m)
    {
        delete[] frst;
        delete[] scnd;
        delete[] nnz;
        delete[] nnz_str;
        delete[] cmmnts;

        fprintf(stdout, "This should not happen: Number of non-zero entries does not match the specified number. This can also be due to malformed lines. Did you use a non-validated file?\n");

        return false;
    }

    return true;
}

// Returns true if everything worked fine, false otherwise.
// Does allocation with new - if true is returned, calling delete[] is your task
bool parseMC_extern(char * const filename, int & n, int & m, int *& frst, int *& scnd, double *& wght, char **& wght_str, char *& cmmnts)
{
    FILE * fp = fopen(filename, "r");

    buffer buf;
    buf.cap = 32;
    buf.line = (char *) malloc(buf.cap);

    if (buf.line == NULL)
    {
        fprintf(stderr, "Fatal error while allocating small buffer\n");
        return false;
    }

    buffer dummy;
    dummy.cap = 32;
    dummy.line = (char *) malloc(dummy.cap);

    if (dummy.line == NULL)
    {
        fprintf(stderr, "Fatal error while allocating small buffer\n");
        free(buf.line);
        return false;
    }

    buffer comments;
    comments.cap = 32;
    comments.line = (char *) malloc(comments.cap);

    long long int commentchars(0);

    if (comments.line == NULL)
    {
        fprintf(stderr, "Fatal error while allocating small buffer\n");
        free(buf.line);
        free(dummy.line);
        return false;
    }

    bool header(false);

    while (! header && buffer_line(fp, &buf))
    {
        if (buf.line[0] != '#')
        {
            if (strlen(buf.line) > dummy.cap)
            {
                dummy.cap = 1 + strlen(buf.line);

                void * rptr = realloc(dummy.line, dummy.cap);

                if (NULL == rptr)
                {
                    fprintf(stderr, "Fatal error while enlarging buffer\n");
                    free(buf.line);
                    free(dummy.line);
                    free(comments.line);

                    return false;
                }

                dummy.line = (char *) rptr;

            }

            int r = sscanf(buf.line, "%d %d%s", &n, &m, dummy.line);

            if (r >= 2 && n >= 1 && m >= 0)
                header = true;

            break;
        }
        else
        {
            if (commentchars + strlen(buf.line) > comments.cap)
            {
                comments.cap = commentchars + 1 + strlen(buf.line);

                void * rptr = realloc(comments.line, comments.cap);

                if (NULL == rptr)
                {
                    fprintf(stderr, "Fatal error while enlarging buffer\n");
                    free(buf.line);
                    free(dummy.line);
                    free(comments.line);
                    return false;
                }

                comments.line = (char *) rptr;
            }

            // omit terminating character
            memcpy(comments.line + commentchars, buf.line, strlen(buf.line));
            commentchars += strlen(buf.line);
        }

    }

    if (! header)
    {
        fclose(fp);
        free(buf.line);
        free(dummy.line);
        free(comments.line);

        fprintf(stdout, "This should not happen: First non-comment line has invalid form. Did you use a non-validated file?\n");

        return false;
    }

    frst = new int[1+m];
    scnd = new int[1+m];
    wght = new double[1+m];
    wght_str = new char*[1+m];

    int mread(0);

    // Detect malformed lines and self-loops
    while (buffer_line(fp, &buf))
    {
        if (buf.line[0] != '#')
        {
            if (strlen(buf.line) > dummy.cap)
            {
                dummy.cap = 1 + strlen(buf.line);

                void * rptr = realloc(dummy.line, dummy.cap);

                if (NULL == rptr)
                {
                    fprintf(stderr, "Fatal error while enlarging buffer\n");
                    free(buf.line);
                    free(dummy.line);
                    delete[] frst;
                    delete[] scnd;
                    delete[] wght;
                    delete[] wght_str;
                    return false;
                }

                dummy.line = (char *) rptr;
            }

            char * str = new char[MAX_WEIGHT_LENGTH];

            int r = sscanf(buf.line, "%d %d %s%s", frst + mread, scnd + mread, str, dummy.line);

            if (r < 3 || frst[mread] == scnd[mread])
            {
                delete[] str;
                break;
            }

            wght_str[mread] = str;
            r = sscanf(str, "%lf", wght + mread);

            if (r != 1)
            {
                delete[] str;
                break;
            }

            --frst[mread];
            --scnd[mread];

            ++mread;
        }
        else
        {
            if (commentchars + strlen(buf.line) > comments.cap)
            {
                comments.cap = commentchars + 1 + strlen(buf.line);

                void * rptr = realloc(comments.line, comments.cap);

                if (NULL == rptr)
                {
                    fprintf(stderr, "Fatal error while enlarging buffer\n");
                    free(buf.line);
                    free(dummy.line);
                    free(comments.line);
                    delete[] frst;
                    delete[] scnd;
                    delete[] wght;
                    delete[] wght_str;
                    return false;
                }

                comments.line = (char *) rptr;
            }

            // omit terminating character
            memcpy(comments.line + commentchars, buf.line, strlen(buf.line));
            commentchars += strlen(buf.line);
        }
    }

    fclose(fp);

    // Add termination to comments, should always
    // work as we always reserve one character of space
    // in addition

    comments.line[commentchars] = '\0';
    cmmnts = new char[commentchars + 1];
    strcpy(cmmnts, comments.line);

    free(buf.line);
    free(dummy.line);
    free(comments.line);

    if (mread != m)
    {
        delete[] frst;
        delete[] scnd;
        delete[] wght;
        delete[] wght_str;
        delete[] cmmnts;

        fprintf(stdout, "This should not happen: Number of valid edge lines does not match the specified number of edges. This may also be due to present self-loops. Did you use a non-validated file?\n");

        return false;
    }

    return true;
}
#endif
