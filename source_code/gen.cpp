#include <cstdio>
#include <cstdlib>
#include <vector>

#define REP(i,n) for (int i=0; i<(n); ++i)
#define SIZE(x) (int)x.size()

using namespace std;

vector<int> v;

void gen(const char *name, int n) {
  FILE *f = fopen(name, "w");
  REP(i,n) {
    fprintf(f, ",");
    REP(j,100) fprintf(f,",%d", v[rand() % SIZE(v)]);
    fprintf(f, "\n");
  }
  fclose(f);
}

int main() {
  FILE *f = fopen("data/gene_ids.txt", "r");
  int x;
  while (fscanf(f, "%d", &x) == 1) v.push_back(x);
  int n = 250;
  gen("up_random3.csv", n);
  gen("down_random3.csv", n);
}
