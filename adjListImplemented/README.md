# Overview

The project is an implementation of the [BMZ heuristic]{https://doi.org/10.1137/S1052623400382467} as described in the report.
The code in this project is a distillation from the code found in the [MQLIB library]{https://github.com/MQLib/MQLib/blob/master/src/heuristics/maxcut/burer2002.cpp}.
Besides releasing the code from the class hierarchy an adjacency list has been implemented increasing code performance.

# Using the code

## Compiling and running the program

The project can be compiled using the "make" command.
The program can be run using:
./OUTPUTNAME INSTANCE
where INSTANCE is the graph on which the heuristic is to be run. 
It is assumed that INSTANCE is in the "mc" or "bq" format as specified in [format specifications]{http://bqp.cs.uni-bonn.de/library/html/formatspecifications.html}
The parameters are fixed after compilation and need to be changed as described in the next subsection.

## Changing parameters

The parameter M of outer iterations and N of inner iterations need to be changed in the sourcecode.
(The formulations are inner and outer iterations are chosen in accordance with the report "report.pdf".)
They are be changed in the file BurerStub.cc in the constructor Burer2002, i.e. lines 407 and 413.
