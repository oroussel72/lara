#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
using namespace std;

#include "Global.h"
#include "Vector.h"
#include "Matrix.h"
#include "Timer.h"
#include "Parameters.h"
#include "Cell.h"
#include "RegularMesh.h"
#include "AdaptiveMesh.h"

int main(int argc, char** argv)
{
  // General init
  
  initParameters();

//  Verbose = Restart = false;
  
  switch(argc)
  {
    case 1:
      break;
      
    case 2:
      if (argv[1][1] == 'v')
	Verbose = true;
      if (argv[1][1] == 'r')
	Restart = true;
      break;
      
    case 3:
      Verbose = Restart = true;
      break;
  };
    
  // Choose multiresolution or fine grid computation
  
  if (Multiresolution)
  {
    AdaptiveMesh *mesh;

    // Init

    mesh = new AdaptiveMesh();
    mesh->init();
    
    // Iterate

    while (RemainingTime > 0.)
    {
      IterationNo++;
      mesh->iterate();
    }

    // End computation

    mesh->end();

    // Deallocate mesh

    delete mesh;
  }
  else
  {
    RegularMesh *mesh;

    // Init

    mesh = new RegularMesh();
    mesh->init();
    
    // Iterate

    while (RemainingTime > 0.)
    {
      IterationNo++;
      mesh->iterate();
    }

    // End computation

    mesh->end();

    // Deallocate mesh

    delete mesh;
  }
  
  exit(1);
}

