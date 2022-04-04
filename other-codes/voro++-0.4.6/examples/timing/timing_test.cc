// Timing test example code
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011
#include <vector>
#include <ctime>
using namespace std;

#include "voro++.hh"
//using namespace voro;

// Set up constants for the container geometry
const double xmin=-1,xmax=1;
const double ymin=-1,ymax=1;
const double zmin=-1,zmax=1;

// Set up the number of blocks that the container is divided into. If the
// preprocessor variable NNN hasn't been passed to the code, then initialize it
// to a good value. Otherwise, use the value that has been passed.
#ifndef NNN
#define NNN 26
#endif
const int n_x=NNN,n_y=NNN,n_z=NNN;

// Set the number of particles that are going to be randomly introduced
const int particles=100;

// This function returns a random double between 0 and 1
double rnd() {return double(rand())/RAND_MAX;}

int main() {
	clock_t start,end;
	int i;

	// Create a container with the geometry given above, and make it
	// periodic in each of the three coordinates. Allocate space for eight
	// particles within each computational block.
	const int numNeighbors = 100;
	//double *x, *y, *z;
	std::vector<double> x,y,z;
	for( int i = 0; i < numNeighbors; i++) {
	  x.push_back((xmin+rnd()*(xmax-xmin))*.8);
	  y.push_back((ymin+rnd()*(ymax-ymin))*.8);
	  z.push_back((zmin+rnd()*(zmax-zmin))*.8);
	  //    con.put(i,x,y,z);
	}
	
	// Store the initial clock time
	start=clock();

	voro::container *con;
	
	for( int iter = 0; iter < 10; iter++) {
	  con = new voro::container(xmin,xmax,ymin,ymax,zmin,zmax,n_x,n_y,n_z,
			      false,false,false,8);
	  
	  //Randomly add particles into the container
	  for(i=0;i<numNeighbors;i++) {
	    con->put(i,x[i],y[i],z[i]);
	  }
	  // Carry out a dummy computation of all cells in the entire container
	  con->compute_all_cells();
	  delete(con);
	}


	// Calculate the elapsed time and print it
	end=clock();
	double runtime=double(end-start)/CLOCKS_PER_SEC;
	printf("%g\n",runtime);
}
