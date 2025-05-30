#include "test.h"

int main () {
  double x;
  int i, j;
  long long_i;

  for (j = 0; j < 4; j++) {
    printf("j = %d\n", j);
  }

  x = M_PI;
  printf("x = %24.16e\n", x);

  /* test of 2^32 = 4294967296: */
  i      =  (int)(pow(2.0,32));
  long_i = (long)(pow(2.0,32));
  printf("i      = 2^32: %%21d = %21d, %%21ld = %21ld\n",      i,      i);
  printf("long_i = 2^32: %%21d = %21d, %%21ld = %21ld\n", long_i, long_i);
  printf("sizeof(int)           = %lu\n", sizeof(int));
  printf("sizeof(long)          = %lu\n", sizeof(long));
  printf("sizeof(long int)      = %lu\n", sizeof(long int));
  printf("sizeof(long long)     = %lu\n", sizeof(long long));
  printf("sizeof(long long int) = %lu\n", sizeof(long long int));

  i = (int)(pow(2.0,31) - 4.0);
  printf("2^31 - 4 = \n");
  for (j = 0; j < 6; j++) {
    i++;
    printf("i++ = %11d = %8x\n", i, i);
  }

  return 0;
}

