#include <iostream>
#include <iomanip>
#include <cmath>
#include <chrono>
#include <thread>
#include <mpi.h>

using namespace std;

double g(int i, int j, int N)
{
  if (i == N/4 && j == N*3/4) return pow(N,2) * 2*M_PI;
  return 0;
}

double boundary(int i, int j, int N)
{
  if (!(i==0 || i==N-1 || j==0 || j==N-1)) {
    throw invalid_argument("i or j should be at the boundary!");
  }
  int i0 = N/4;
  int j0 = N*3/4;

  return log(sqrt(pow(i-i0,2)+pow(j-j0,2)));
}

double f( double f, double l, double r, double d, double u, double g, int N )
{ 
  double gam = 0.1;
  return (1-gam)*f + gam*0.25*(l+r+u+d - g/pow(N,2) );
}

// Update the left/right border of the subgrid
void update_lr(double** M, int N, int id, int x_min, int x_max, int id_next)
{
  if (id%2 == 0) {
    MPI_Send(M[x_max-1], N, MPI_DOUBLE, id_next, 1, MPI_COMM_WORLD);
    MPI_Recv(M[x_max], N, MPI_DOUBLE, id_next, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  } else {
    MPI_Recv(M[x_min-1], N, MPI_DOUBLE, id_next, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Send(M[x_min], N, MPI_DOUBLE, id_next, 1, MPI_COMM_WORLD);
  } 
}

// Update the top/bottom border of the subgrid
void update_tb(double** M, int N, int id, int x_min, int x_max, int y_min, int y_max, int id_next)
{
  if (id < 2) {
    for (int i=x_min; i<x_max; i++) {
      MPI_Send(&(M[i][y_max-1]), 1, MPI_DOUBLE, id_next, 1, MPI_COMM_WORLD);
      MPI_Recv(&(M[i][y_max]), 1, MPI_DOUBLE, id_next, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
  } else {
    for (int i=x_min; i<x_max; i++) {
      MPI_Recv(&(M[i][y_min-1]), 1, MPI_DOUBLE, id_next, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Send(&(M[i][y_min]), 1, MPI_DOUBLE, id_next, 1, MPI_COMM_WORLD);
    } 
  } 
} 

void red_black( double** M, int N, int iters, int id, int ntasks )
{
  if ( !(ntasks==1 || ntasks==2 || ntasks==4 ) ) {
    throw invalid_argument("Number of processes should be 1, 2 or 4");
  }
  int x_min = 1;
  int y_min = 1;
  int x_max = N-1;
  int y_max = N-1;
  int id_lr, id_tb, id_top, id_bottom;

  // Setting boundaries for each
  // parallel process
  if (ntasks != 1) {
    if ( id%2 == 0 ) {
      x_max = x_max / 2;
      id_lr = id+1;
    } else {
      x_min = x_max / 2;
      id_lr = id-1;
    }
    if (ntasks == 4) {
      if ( id<2 ) {
        y_max = y_max / 2;
        id_bottom = id+2;
        id_tb = id+2;
      } else {
        y_min = y_max / 2;
        id_top = id-2;
        id_tb = id-2;
      }
    }
  }
  
  for ( int k=0; k<iters; k++ ) {
    // RED sublattice
    for ( int i=x_min; i<x_max ; i++ ) {
      for ( int j=i%2+y_min ; j<y_max; j+=2 ) {
        M[i][j] = f(M[i][j], M[i-1][j], M[i+1][j], M[i][j-1], M[i][j+1], g(i,j,N), N);
      }
    }
    // Sync
    if (ntasks != 1) {
      update_lr(M,N,id,x_min,x_max,id_lr);
      if (ntasks == 4) update_tb(M,N,id,x_min,x_max,y_min,y_max,id_tb);
    } 

    // BLACK sublattice
    for ( int i=x_min; i<x_max; i++ ) {
      for ( int j=(i+1)%2+y_min; j<y_max; j+=2 ) {
        M[i][j] = f(M[i][j], M[i-1][j], M[i+1][j], M[i][j-1], M[i][j+1], g(i,j,N), N);
      }
    }
    // Sync
    if (ntasks != 1) {
      update_lr(M,N,id,x_min,x_max,id_lr);
      if (ntasks == 4) update_tb(M,N,id,x_min,x_max,y_min,y_max,id_tb);
    } 
  }

  // MERGE THE RESULTS
  
  // For four subprocesses, first the top subgrid 
  // is merged with the corresponding bottom subgrid
  if (ntasks == 4) {
    if (id < 2) {
      int from_id = id+2;
      for (int j=y_max; j<N-1; j++) {
        for (int i=x_min; i<x_max; i++) {
          MPI_Recv(&M[i][j], 1, MPI_DOUBLE, from_id, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
      }
    } else {
      int to_id = id-2;
      for (int j=y_min; j<y_max; j++) {
        for (int i=x_min; i<x_max; i++) {
          MPI_Send(&M[i][j], 1, MPI_DOUBLE, to_id, 1, MPI_COMM_WORLD);
        }
      }
    }
  }

  // Left and right subgrids merged
  if (ntasks != 1) {
    if (id == 0) {
      int from_id = 1;
      for (int i=x_max; i<N-1; i++) {
        MPI_Recv(M[i], N, MPI_DOUBLE, from_id, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }
    } else if (id == 1) {
      for (int i=x_min; i<x_max; i++) {
        MPI_Send(M[i], N, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
      }
    }
  }
}

int main(int argc, char *argv[])
{
  auto t0 = chrono::system_clock::now();
  int rc = MPI_Init(&argc, &argv);
  if (rc!=MPI_SUCCESS) {
    cout << "MPI_Init failed!" << endl;
    return 1;
  }

  int ntasks;
  MPI_Comm_size(MPI_COMM_WORLD, &ntasks);

  int id;
  MPI_Comm_rank(MPI_COMM_WORLD, &id);

  int N = 300;
  double** M = new double*[N];
  for (int i = 0; i < N; i++) {
    // Declare a memory block of size n
    M[i] = new double[N];
  }
  for (int j = 0; j < N; j++) {
    for (int i = 0; i < N; i++) {
      M[j][i] = 0;
    }
  }
  for (int i = 0; i < N; i++) {
    M[1][i] = 0;
  }

  for (int i=0; i<N; i++) {
    M[i][0] = boundary(i,0,N);
    M[i][N-1] = boundary(i,N-1,N);
    M[0][i] = boundary(0,i,N);
    M[N-1][i] = boundary(N-1,i,N);
  }

//  cout << fixed << setprecision(2);
//  for (int i=0; i<N; i++) {
//    for (int j=0; j<N; j++) {
//      cout << M[i][j] << " ";
//    }
//    cout << endl;
//  }

  red_black( M, N, 10000, id, ntasks );
  auto t1 = chrono::system_clock::now();
  auto wct = chrono::duration_cast<chrono::milliseconds>(t1-t0);
  if (id == 0) cerr << wct.count()/1000.0 << endl;
  rc = MPI_Finalize();
  return 0;

  this_thread::sleep_for(chrono::seconds(id));
  cout << endl;
  cout <<id << endl;
  for (int i=0; i<N; i++) {
    for (int j=0; j<N; j++) {
      cout << M[i][j] << " ";
    }
    cout << endl;
  }

  cout << endl;
  cout << endl;

  rc = MPI_Finalize();
  return 0;
}

