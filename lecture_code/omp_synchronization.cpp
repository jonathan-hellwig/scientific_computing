#include <omp.h>
#include <cstdio>

int main()
{
  int x;
  x = 0;
  #pragma omp parallel shared(x)
  {
    #pragma omp critical
    x = x + 1;
    printf("x: %d\n", x);
  }
}
