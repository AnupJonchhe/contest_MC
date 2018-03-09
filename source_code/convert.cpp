const int NUM_SIG = 476251;

#include <cstdio>
#include <cassert>
#include <cstdlib>

int main(int argc, char **argv) {
  if (argc != 4) {
    fprintf(stderr, "Usage %s num_queries binaryfile csvfile\n", argv[0]);
    exit(1);
  }
  int q = atoi(argv[1]);
  FILE *fin = fopen(argv[2], "r");
  FILE *fout = fopen(argv[3], "w");
  float truth[NUM_SIG];
  for (int i = 0; i < NUM_SIG; ++i) {
    int x = fread(truth, sizeof(float), q, fin);
    //    fprintf(stderr, "x=%d\n", x);
    assert(x == q);
    for (int j = 0; j < q; ++j) {
      if (j) fprintf(fout, ",");
      fprintf(fout, "%.6lf", truth[j]);
    }
    fprintf(fout, "\n");
  }
  return 0;
}
