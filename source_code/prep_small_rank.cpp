const int NUM_GENES = 10174;
const int NUM_SIG = 476251;

#define ALIGNED __attribute__ ((aligned(16)))

#define SSELOADF(a)     _mm256_load_ps((__m256*)&a)

#include <bits/stdtr1c++.h>
#include <unistd.h>
#include <sys/time.h>

//#include <avxintrin.h>
//#include <immintrin.h>
#include <x86intrin.h>

using namespace std;
string TEST_DESCRIPTION;
double MYSCORE;

#ifdef LOCAL
const double TIME_SCALE = 0.8;
#else
const double TIME_SCALE = 1.0;
#endif
const double TIME_LIMIT1 = TIME_SCALE * 16.0;
const double TIME_LIMIT2 = TIME_SCALE * 17.0;
const double TIME_LIMIT_LAST = TIME_SCALE * 19;

typedef vector <int> VI;
typedef vector <VI> VVI;
typedef long long LL;
typedef vector <LL> VLL;
typedef vector <double> VD;
typedef vector <VD> VVD;
typedef vector <string> VS;
typedef vector <VS> VVS;
typedef pair<int,int> PII;
typedef vector <PII> VPII;
typedef istringstream ISS;

#define ALL(x) x.begin(),x.end()
#define REP(i,n) for (int i=0; i<(n); ++i)
#define FOR(var,pocz,koniec) for (int var=(pocz); var<=(koniec); ++var)
#define FORD(var,pocz,koniec) for (int var=(pocz); var>=(koniec); --var)
#define FOREACH(it, X) for(__typeof((X).begin()) it = (X).begin(); it != (X).end(); ++it)
#define PB push_back
#define PF push_front
#define MP(a,b) make_pair(a,b)
#define ST first
#define ND second
#define SIZE(x) (int)x.size()

//#ifndef DEBUG
//#define assert(X) ;
//#endif

#ifndef LOCAL
#define assert(X) ;
#endif

template<class T> string i2s(T x) {ostringstream o; o << x; return o.str();}
template<class T1,class T2> ostream& operator<<(ostream &os, pair<T1,T2> &p) {os << "(" << p.first << "," << p.second << ")"; return os;}
template<class T> ostream& operator<<(ostream &os, vector<T> &v) {os << "{"; REP(i, (int)v.size()) {if (i) os << ", "; os << v[i];} os << "}"; return os;}
#define DB(a) {cerr << #a << ": " << (a) << endl; fflush(stderr); }
//#define DB2(a,b) {cerr << (a) << ": " << (b) << endl; fflush(stderr); }

namespace Time{
  double start_time;
  static double last_call = 0;
   
  double get_time() {
    timeval tv;
    gettimeofday(&tv, 0);
    return tv.tv_sec+tv.tv_usec*1e-6;
  }

  void print_time(string s) {
    double x = get_time();
    fprintf(stderr,"%s cur=%.6lf lap=%.6lf\n",s.c_str(),x,x-last_call);
    last_call = x;
  }

  void init_time() {
    start_time = get_time();
    last_call = start_time;
  }
}

#define STOPER(name) name(#name)

struct Stoper;
vector<Stoper*> stoper_pointers;

struct Stoper {
  double used_time;
  double last_call;
  string name;
  void start() {
#ifdef LOCAL
    last_call = Time::get_time();
#endif
  }
  void stop() {
#ifdef LOCAL
    used_time += Time::get_time() - last_call;
#endif
  }

  Stoper(string s="") {
    used_time = 0.0;
    name=s;
    stoper_pointers.PB(this);
  }
}; 
Stoper STOPER(st_whole);
Stoper STOPER(st_pow);

/************************************************************************/
/************************ Code starts here ******************************/
/************************************************************************/


const int MAX_LINE = 1000000;
char line[MAX_LINE];

string name1, name2;

void load_rank(int idx, int *output) {
  FILE *f = fopen(name1.c_str(), "r");
  int err = fseek(f, (long)idx * NUM_GENES * sizeof(int), SEEK_SET);
  assert(err == 0);
  int elts_read = fread(output, sizeof(int), NUM_GENES, f);
  assert(elts_read == NUM_GENES);
  fclose(f);
}

void load_signature(int idx, double *output) {
  FILE *f = fopen(name2.c_str(), "r");
  int err = fseek(f, (long)idx * NUM_GENES * sizeof(double), SEEK_SET);
  assert(err == 0);
  int elts_read = fread(output, sizeof(double), NUM_GENES, f);
  assert(elts_read == NUM_GENES);
  fclose(f);
}

int main(int argc, char **argv){
  if (argc != 3) {
    fprintf(stderr, "Usage %s orig_rank_trans new_rank_trans orig_score_trans new_score_trans\n", argv[0]);
    exit(1);
  }
  name1 = argv[1];
//  name2 = argv[3];

  {
    FILE *f2 = fopen(argv[2], "w");

    int rank[NUM_GENES];
    short int short_order[NUM_GENES];
    REP(i, NUM_SIG) {
      load_rank(i, rank);
      REP(j, NUM_GENES) short_order[rank[j]-1] = j;
      int x = fwrite(short_order, sizeof(short int), NUM_GENES, f2);
      assert(x == NUM_GENES);
      if (i % 1000 == 0) DB(i);
    }
    fclose(f2);
  }
/*
  {
    FILE *f2 = fopen(argv[4], "w");

    double score[NUM_GENES];
    float float_score[NUM_GENES];
    REP(i, NUM_SIG) {
      load_signature(i, score);
      REP(j, NUM_GENES) float_score[j] = fabs(score[j]);
      int x = fwrite(float_score, sizeof(float), NUM_GENES, f2);
      assert(x == NUM_GENES);
      if (i % 1000 == 0) DB(i);
    }
    fclose(f2);
  }
*/
  return 0;
}
