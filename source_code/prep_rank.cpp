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

/*
void multiply_and_add(const float* a, const float* b, const float* c, float* d) {  
  for(int i=0; i<8; i++) {
    d[i] = a[i] * b[i];
    d[i] = d[i] + c[i];
    printf("d[%d]=%.6lf\n", i, d[i]);
  }
}

inline __m256 multiply_and_add(__m256 a, __m256 b, __m256 c) {
  return _mm256_add_ps(_mm256_mul_ps(a, b), c);
}
*/

/* loads scores for a given signature */
#ifdef LOCAL
void load_signature(int idx, double *output) {
  FILE *f = fopen("data/score_trans.bin", "r");
  int err = fseek(f, (long)idx * (long)NUM_GENES * sizeof(double), SEEK_SET);
  assert(err == 0);
  if (err) {
    DB(err);
    DB(strerror(err));
  }
  int elts_read = fread(output, sizeof(double), NUM_GENES, f);
  assert(elts_read == NUM_GENES);
  fclose(f);
}
#else
void load_signature(int idx, double *output) {
  //loadFromDoubleFile("scoresBySig", idx * 10174, 10174);
}
#endif

map<int,int> gene_map; //maping ids to range [0,NUM_GENES-1]

VVI parse(VS vs) {
  VVI res;
  REP(i, SIZE(vs)) {
    for (auto &c: vs[i]) if (c == ',') c = ' ';
    istringstream in(vs[i]);
    int x;
    VI v;
    while (in >> x) {
      assert(gene_map.count(x));
      v.PB(gene_map[x]);
    }
    res.PB(v);
    //if (i < 5) DB(res.back());
  }
  return res;
}

double cur_sig[NUM_GENES+10];

void solve_sig(int sig, double *output, const VVI &vall) {
  int q = SIZE(vall) / 2;
  load_signature(sig, cur_sig);
  REP(i, SIZE(vall)) {
  }
}

class CMAP2Updated {
  public:

    int init(VI gene_list) {
      REP(i, SIZE(gene_list)) gene_map[gene_list[i]] = i;
      return 0;
    }

    vector <double> getWTKScomb(vector <string> up, vector <string> down) {
      REP(i,5) printf("up[%d]=%s\n", i, up[i].substr(0, 50).c_str());
      REP(i,5) printf("down[%d]=%s\n", i, down[i].substr(0, 50).c_str());

      int q = SIZE(up);
      assert(SIZE(up) == SIZE(down));
      VVI vup = parse(up);
      VVI vdown = parse(down);
      VVI vall = vup;
      vall.insert(vall.end(), ALL(vdown));
      assert(SIZE(vall) == 2 * q);

      VD res(NUM_SIG * q, 0.0);
      double tmp[q];
      REP(sig, NUM_SIG) {
        solve_sig(sig, tmp, vall);
        REP(i, q) {
          if ((tmp[i] > 0 && tmp[i+q] > 0) || (tmp[i] < 0 && tmp[q+i] < 0)) continue;
          res[q * sig + i] = (tmp[i] - tmp[q+i]) / 2.0;
          if (sig < 10 && i < 10) fprintf(stderr, "res for sig=%d query=%d is %.6lf\n", sig, i, res[q * sig + i]);
        }
      }
      /*
      float a[8], b[8], c[8], d[8] ALIGNED;
      REP(i,8) {
        a[i] = drand48();
        b[i] = drand48();
        c[i] = drand48();
      }
      multiply_and_add(a, b, c, d);
      __m256 aa = _mm256_load_ps(a);
      __m256 bb = _mm256_load_ps(b);
      __m256 cc = _mm256_load_ps(c);
      __m256 x = multiply_and_add(aa, bb, cc);
      float *f = (float*)&x;
      puts("");
      puts("");
      REP(i,8) printf("d[%d]=%.6lf\n", i, f[i]);
      */

      return VD();
    }
};


#ifdef LOCAL

const int MAX_LINE = 1000000;
char line[MAX_LINE];

VS load_queries(string filename) {
  FILE *f = fopen(filename.c_str(), "r");
  VS res;
  while (fgets(line, MAX_LINE, f)) {
    int pos = 0;
    while (line[pos] != ',') pos++;
    pos++;
    assert(line[pos] == ',');
    pos++;
    assert(isdigit(line[pos]));
    res.PB(line+pos);
  }
  return res;
}

int main(int argc, char **argv){

  double tmp[1000000];
  vector<pair<double,int>> v;
  FILE *f = fopen("data/rank_trans.bin", "w");
  int rank[NUM_GENES];
  REP(i, NUM_SIG) {
    load_signature(i, tmp);
    v.clear();
    REP(j, NUM_GENES) v.PB(MP(tmp[j], -j));
    sort(ALL(v));
    reverse(ALL(v));
    REP(j, NUM_GENES) rank[-v[j].ND] = j+1;
    int x = fwrite(rank, sizeof(int), NUM_GENES, f);
    //DB(x);
    //DB(NUM_GENES);
    assert(x == NUM_GENES);
    DB(i);
  }
  fclose(f);
  return 0;
}

#endif

