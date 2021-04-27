#include <omp.h>
#include <cstdio>
#include <vector>
#include <math.h>

using Vector = std::vector<double>;
using Matrix = std::vector<std::vector<double>>;

Vector solve_system(Matrix A, Vector b, Vector x_0, double tol);
void print_vector(Vector x);
void print_matrix(Matrix A);

int main()
{
  Matrix A = {{1,0,1},
              {0,1,0},
              {1,0,1}};
  Vector b = {10,1,10};
  Vector x_0 = {0,0,0};
  double tol = 0.001;
  printf("Starting CG algorithm.\n");
  printf("For A = \n");
  print_matrix(A);
  printf("b = \n");
  print_vector(b);
  printf("x_0 = \n");
  print_vector(x_0);
  Vector x = solve_system(A, b, x_0, tol);
  printf("\nSolution\n");
  print_vector(x);
}

void print_vector(Vector x)
{
  for(auto const &num: x)
    printf("%f ", num);
  printf("\n");
}

void print_matrix(Matrix A)
{
  for(auto const &row: A)
  {
    for(auto const &element: row)
    {
      printf("%f ", element);
    }
    printf("\n");
  }
}



Vector solve_system(Matrix A, Vector b, Vector x_0, double tol)
{
  int n = A.size();
  int i;
  Vector h(n);
  Vector x, r, p;
  x.resize(n);
  r.resize(n);
  p.resize(n);
  double alpha, beta, gamma, delta, zeta, sigma;
  double error = 0.0;
  //printf("Init\n");
  #pragma omp parallel for private(i) default(shared) reduction(+:error)
  for(i = 0; i < n; i++)
  {
    for(int j = 0; j < n; j++)
    {
      h[i] += A[i][j] * x_0[j];
    }
    r[i] = b[i] - h[i];
    p[i] = r[i];
    x[i] = x_0[i];
    error += r[i]*r[i];
  }
  error = std::sqrt(error);
  //print_vector(r);
  //print_vector(p);
  //print_vector(x);

  //printf("Setup\n");
  while(error > tol)
  {
    alpha = 0.0;
    beta = 0.0;
    gamma = 0.0;
    delta = 0.0;
    zeta = 0.0;
    sigma = 0.0;

    #pragma omp parallel for private(i) default(shared) reduction(+:gamma,beta,delta,zeta)
    for(i = 0; i < n; i++)
    {
      h[i] = 0;
      for(int j = 0; j < n; j++)
      {
        h[i] += A[i][j] * p[j];
      }
      gamma += p[i]*h[i];
      beta += r[i]*r[i];
      delta += r[i]*h[i];
      zeta += h[i]*h[i];
    }
    alpha = beta / gamma;
    sigma = beta - 2*alpha*delta + alpha*alpha*zeta;
    error = 0.0;
    #pragma omp parallel for private(i) default(shared) reduction(+:error)
    for(i = 0; i < n; i++){
      x[i] = x[i] + alpha*p[i];
      r[i] = r[i] - alpha*h[i];
      p[i] = r[i] + sigma/beta*p[i];
      error += r[i]*r[i];
    }
    error = std::sqrt(error);
    //printf("x: ");
    //print_vector(x);
    //printf("h: ");
    //print_vector(h);
    //printf("r: ");
    //print_vector(r);
    //printf("p: ");
    //print_vector(p);
  }
  return x;
}

