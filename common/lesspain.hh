/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/* Sven Mallach (2021) */

#ifndef MCBQ_LESSPAIN_HH
#define MCBQ_LESSPAIN_HH 1

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>

struct buffer
{
    char * line;
    size_t cap;
};

bool find_option(int argc, char * argv[], const char * option)
{
    for (int i(1) ; i < argc ; ++i)
    {
        if (0 == strcmp(argv[i], option))
        {
            return true;
        }
    }

    return false;
}

// Attention, this does not give a NEW string, just a position!
char * remove_path_from_filename(char * filename)
{
    // Strip path off from filename
    char * fns = strtok(filename, "/");
    char * fn;

    do
    {
        fn = fns;
        fns = strtok(NULL, "/");
    }
    while (fns != NULL);

    return fn;
}

bool check_readability(char * filename)
{
    FILE * fp = fopen(filename, "r");

    if( fp == NULL )
    {
        return false;
    }
    else
    {
        fclose(fp);
        return true;
    }
}

char * get_tempfile()
{
    char * tmpfn = new char[100];
    sprintf(tmpfn, "/tmp/tool.XXXXXX");

    if (mkstemp(tmpfn) == -1)
    {
        fprintf(stderr, "Could not create temporary file\n");
        return nullptr;
    }
    else
    {
        return tmpfn;
    }
}

void remove_tempfile(char * tmpfn)
{
    unlink(tmpfn);
    delete[] tmpfn;
}

bool buffer_line(FILE * fp, buffer * b)
{
    if (NULL == fgets(b->line, b->cap, fp))
        return false;

    // We can assume here that at least one
    // character was read, so strlen(b->line) is at least 1.

    unsigned int offset = 0;

    while (! (feof(fp) || b->line[strlen(b->line) - 1] == '\n') )
    {
        offset = b->cap;
        b->cap <<= 1;

        void * rptr = realloc(b->line, b->cap);

        if (NULL == rptr)
        {
            fprintf(stderr, "Fatal error while enlarging buffer\n");
            return false;
        }

        b->line = (char *) rptr;

        // Now append into "newly allocated half" thereby overwriting
        // the null termination of the former fgets call
        if (NULL == fgets(b->line + offset - 1, offset + 1, fp))
            return false;
    }

    return true;
}
#endif
