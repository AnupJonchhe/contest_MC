const int NUM_GENES = 10174;
const int NUM_SIG = 476251;
const int NUM_Q = 2013;
const int MAX_THREADS = 20;
const int BUFFER = 60000;

int NUM_THREADS = 1;
int BATCH_SIZE = 10;
int TEST = 0;
int CALC_SIG = NUM_SIG; 

#include <string>
std::string PATH = "data";
std::string PATH2 = "data";
std::string DOWN_FILE = "offline_query_down_n250.csv";
std::string UP_FILE = "offline_query_up_n250.csv";
std::string OUT_FILE = "fp32_output.bin";

#define ALIGNED __attribute__ ((aligned(16)))

#include <bits/stdtr1c++.h>
#include <unistd.h>
#include <sys/time.h>
#include <thread>

//#include <avxintrin.h>
//#include <immintrin.h>
#include <x86intrin.h>

using namespace std;
string TEST_DESCRIPTION;
double MYSCORE;

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
  string name;
  mutex m;

  void start() {
#ifdef LOCAL
    double x = Time::get_time();
    m.lock();
    used_time -= x;
    m.unlock();
#endif
  }
  void stop() {
#ifdef LOCAL
    double x = Time::get_time();
    m.lock();
    used_time += x;
    m.unlock();
#endif
  }

  Stoper(string s="") {
    used_time = 0.0;
    name=s;
    stoper_pointers.PB(this);
  }
}; 
Stoper STOPER(st_disk);
Stoper STOPER(st_prep);
Stoper STOPER(st_zero);
Stoper STOPER(st_rs);
Stoper STOPER(st_sort);
Stoper STOPER(st_load_rank);
Stoper STOPER(st_load_sig);
Stoper STOPER(st_write);
Stoper STOPER(st_full);

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

map<int,int> gene_map; //maping ids to range [0,NUM_GENES-1]
const int MAX_LINE = 1000000;
char line[MAX_LINE];
VI vall[NUM_Q];

int parse(const VS &vs, int offset) {
  int q = SIZE(vs);
  REP(i, q) {
    vall[i+offset].reserve(100);
    int pos = 0;
    for (auto &c: vs[i]) {
      if (c == ',') {
        line[pos] = 0;
        vall[offset+i].PB(gene_map[atoi(line)]);
        pos = 0;
      } else line[pos++] = c;
    }
    if (pos) {
      line[pos] = 0;
      vall[offset+i].PB(gene_map[atoi(line)]);
    }
  }
  return q;
}

//For each gene store its list of queries
int qtab[NUM_Q * (NUM_GENES + 53)];
int qbeg[NUM_GENES+10];
int qend[NUM_GENES+10];

int query_size[NUM_Q];

//For each query store its list of scores and ranks

// i - index queries
void prepare_queries(int q) {

  int sum = 0;
  REP(i, q) {
    query_size[i] = SIZE(vall[i]);
    sum += SIZE(vall[i]);
  }
  sum += NUM_GENES * 3;
  REP(i, sum) qtab[i] = q;

  REP(i, NUM_GENES+1) qbeg[i] = qend[i] = 0;
  REP(i, q) for (auto x : vall[i]) qend[x]++;

  //REP(i, NUM_GENES) fprintf(stderr, "occ[%d]=%d\n", i, qend[i]);
  REP(i, NUM_GENES) {
    while (qend[i] & 3) qend[i]++;
    qend[i+1] += qend[i];
    qbeg[i+1] += qend[i];
  }
  REP(i, q) for (auto x : vall[i]) qtab[qbeg[x]++] = i;
  qbeg[0] = 0;
  REP(i, NUM_GENES) qbeg[i+1] = qend[i];

  REP(i, 5) DB(SIZE(vall[i]));

}

float tab_sig[BUFFER+5][NUM_GENES+500];
short int tab_order[BUFFER+5][NUM_GENES+500];

class Solver{
  //short int order[NUM_GENES+10];

  //float cur_sig[NUM_GENES+10];
  //int cur_rank[NUM_GENES+10];

  //modifiable
  float ssum_tab[NUM_Q];
//  short int rtab[NUM_Q * NUM_GENES];
  struct RS{
    float rs;
    float res;
  } rs_tab[NUM_Q];

  public:
  void solve_sig2(int bid, int sig, double *output, int q) {
    if (sig < 10 || sig % 1000 == 0) DB(sig);
    //DB(sig);
    short int *order = tab_order[bid];
    float *cur_sig = tab_sig[bid];

    //TODO: optimize with memcpy/memset
    st_zero.start();
    REP(i, q) {
      ssum_tab[i] = 0.0;
      rs_tab[i].rs = rs_tab[i].res = 0.0;
    }
    st_zero.stop();

    // make a pass in normal order
    st_load_rank.start();
    REP(x, NUM_GENES) {
      float score = cur_sig[x];
      for (int a = qbeg[x]; a < qend[x]; ++a) {
        int i = qtab[a];
        ssum_tab[i] += score;
      }
    }

    //REP(i, q) ssum_tab[i] /= NUM_GENES - query_size[i];
    st_load_sig.start();
    REP(i, q) ssum_tab[i] = (NUM_GENES - query_size[i]) / ssum_tab[i];
    st_load_sig.stop();
    st_load_rank.stop();

    typedef float D;
    st_rs.start();
    D ones[] ALIGNED = {1., 1., 1., 1.};
    __m128 m_ones = _mm_load_ps(ones);
    __m128 m_zeros = _mm_set_ps1(0.0);

    __m128i minus1 = _mm_set1_epi32(-1);
    __m128 abs_mask = _mm_castsi128_ps(_mm_srli_epi32(minus1, 1));

    REP(xx, NUM_GENES) {
      __builtin_prefetch(qtab+qbeg[order[xx+2]]);
      int x = order[xx];
      float score = cur_sig[x];
      D rank = xx;
      int a = qbeg[x];

      while (a < qend[x]) {
        D scores[] ALIGNED = {score, score, score, score};
        __m128 m_scores = _mm_load_ps(scores);


        D ranks[] ALIGNED = {rank, rank, rank, rank};
        __m128 m_ranks = _mm_load_ps(ranks);

        int i1 = qtab[a++]; // query id
        int i2 = qtab[a++];
        int i3 = qtab[a++];
        int i4 = qtab[a++];

        D ssumi1 = ssum_tab[i1];
        D ssumi2 = ssum_tab[i2];
        D ssumi3 = ssum_tab[i3];
        D ssumi4 = ssum_tab[i4];
        D sums[] ALIGNED = {ssumi1, ssumi2, ssumi3, ssumi4};
        __m128 m_sums = _mm_load_ps(sums);

        D rss[] ALIGNED = {rs_tab[i1].rs, rs_tab[i2].rs, rs_tab[i3].rs, rs_tab[i4].rs};
        __m128 m_rss = _mm_load_ps(rss);

        D rs_maxs[] ALIGNED = {rs_tab[i1].res, rs_tab[i2].res, rs_tab[i3].res, rs_tab[i4].res}; //TODO: change the two into one load by playing with doubles and converting back?
        __m128 m_rs_maxs = _mm_load_ps(rs_maxs);
        __m128 m_abs = _mm_and_ps(abs_mask, m_rs_maxs); // abs(res)

        __m128 tmp = _mm_sub_ps(m_rss, m_ranks);
        __m128 tmp_abs = _mm_and_ps(abs_mask, tmp); // abs(tmp)

        __m128 cmpres = _mm_cmpgt_ps(tmp_abs, m_abs);
        //m_rs_maxs = _mm_or_ps(_mm_and_ps(cmpres, tmp), _mm_andnot_ps(cmpres, m_rs_maxs));
        m_rs_maxs = _mm_blendv_ps(m_rs_maxs, tmp, cmpres);

        /*
           rs_tab[i1].rs += score * ssumi1;
         */
        __m128 tmp2 = _mm_mul_ps(m_scores, m_sums);
        m_rss = _mm_add_ps(m_rss, tmp2);
        tmp = _mm_add_ps(tmp, tmp2);

        m_abs = _mm_and_ps(abs_mask, m_rs_maxs); // abs(res)
        tmp_abs = _mm_and_ps(abs_mask, tmp); // abs(tmp)

        cmpres = _mm_cmpgt_ps(tmp_abs, m_abs);
        //m_rs_maxs = _mm_or_ps(_mm_and_ps(cmpres, tmp), _mm_andnot_ps(cmpres, m_rs_maxs));
        m_rs_maxs = _mm_blendv_ps(m_rs_maxs, tmp, cmpres);

        _mm_store_ps(rs_maxs, m_rs_maxs);
        rs_tab[i1].res = rs_maxs[0];
        rs_tab[i2].res = rs_maxs[1];
        rs_tab[i3].res = rs_maxs[2];
        rs_tab[i4].res = rs_maxs[3];

        _mm_store_ps(rss, _mm_add_ps(m_rss, m_ones));
        rs_tab[i1].rs = rss[0];
        rs_tab[i2].rs = rss[1];
        rs_tab[i3].rs = rss[2];
        rs_tab[i4].rs = rss[3];
      }
    }

    //REP(i, q) output[i] = (rs_tab[i].max_rs >= -rs_tab[i].min_rs ? rs_tab[i].max_rs : rs_tab[i].min_rs) / (NUM_GENES - query_size[i]);
    REP(i, q) output[i] = rs_tab[i].res / (NUM_GENES - query_size[i]);
    st_rs.stop();
  }
} solvers[MAX_THREADS];

void init_genes(VI gene_list) {
  REP(i, SIZE(gene_list)) gene_map[gene_list[i]] = i;
  gene_list.resize(10);
  DB(gene_list);
}

float *response;

class Buffer {
  private:
  PII q_buffer[BUFFER]; 
  int q_beg, q_end, q_cnt;

  mutex lock;
  condition_variable not_full;
  condition_variable not_empty;

  public:
  Buffer() {
    q_beg = q_end = q_cnt = 0;
  }

  void add(PII p){
    std::unique_lock<std::mutex> l(lock);

    not_full.wait(l, [this](){return q_cnt != BUFFER; });

    q_buffer[q_end] = p;
    q_end = (q_end + 1) % BUFFER;
    ++q_cnt;

    not_empty.notify_one();
  }

  PII get(){
    std::unique_lock<std::mutex> l(lock);

    not_empty.wait(l, [this](){return q_cnt != 0; });

    PII res = q_buffer[q_beg];
    q_beg = (q_beg + 1) % BUFFER;
    --q_cnt;

    not_full.notify_one();
    return res;
  }
} buffer_sig, buffer_ids1, buffer_ids2;

int reading_done;
mutex sig_mutex;


struct Reader {
  FILE *fsig;
  FILE *frank;
  int start, finish;

  void load_signature_new(int idx, float *output) {
    int elts_read = fread(output, sizeof(float), NUM_GENES, fsig);
    cerr << "Number of elements read = " << (elts_read) << endl;
    assert(elts_read == NUM_GENES);
  }

  void load_rank(int idx, short int *output) {
    cerr << "output:" << (output) << ". idx: " << (idx) << endl;
    int elts_read = fread(output, sizeof(short int), NUM_GENES, frank);
    cerr << "Number of rank elements read = " << (elts_read) << endl;
    cerr << "feof(frank) returns: " << (feof(frank)) << " for index: " << (idx) << endl;
    cerr << "ferror(frank) returns: " << (ferror(frank)) << " for index: " << (idx) << endl;
    assert(elts_read == NUM_GENES);
  }

  Reader(int _start, int _finish) {
    start = _start;
    finish = _finish;
    fsig = fopen(((start == 0 ? PATH : PATH2)+"/score_trans_small.bin").c_str(), "r");
    frank = fopen(((start == 0 ? PATH : PATH2)+"/rank_trans_small.bin").c_str(), "r");
    cerr << ((start == 0 ? PATH : PATH2)+"/rank_trans_small.bin").c_str() << endl;
    if (start != 0) {
      int err = fseek(fsig, (long)start * NUM_GENES * sizeof(float), SEEK_SET);
      assert(err == 0);
      err = fseek(frank, (long)start * NUM_GENES * sizeof(short int), SEEK_SET);
      assert(err == 0);
    }
  }

  void go() {
    PII p;

    while (start < finish) {
      if (finish == CALC_SIG) p = buffer_ids2.get();
      else p = buffer_ids1.get();

      p.ND = start++;

      st_disk.start();
      load_signature_new(p.ND, tab_sig[p.ST]);
      cerr << "p.ND value:" << (p) << endl;
      cerr << "tab_order[p.ST] :" << (*tab_order[p.ST]) << endl;
      load_rank(p.ND, tab_order[p.ST]);
      st_disk.stop();

      buffer_sig.add(p);
    }
    sig_mutex.lock();
    reading_done += 1;
    int last = (reading_done == 2);
    sig_mutex.unlock();
    if (last) {
      REP(foo, NUM_THREADS) buffer_sig.add(MP(-1,-1));
    }
  }
};

void thread_wrapper(int id) {
  assert(id >= 0 && id < 2);
  int a = 0;
  int b = (CALC_SIG / 2) * 0.95;
  if (id == 1) {
    a = b;
    b = CALC_SIG;
  } 
  Reader reader(a, b);
  reader.go();
}

int written_flag[NUM_SIG+1000];
mutex lock_writer;
condition_variable cond_writer;
int writer_pos;

void thread_writer(int q) {
  FILE *f = fopen(OUT_FILE.c_str(), "w");
  while (writer_pos < CALC_SIG) {
    {
      std::unique_lock<std::mutex> l(lock_writer);
      cond_writer.wait(l, [](){return written_flag[writer_pos];});
    }
    if (writer_pos % 1000 == 0) DB(writer_pos);
    st_write.start();
    int x = fwrite(response + writer_pos * q, sizeof(float), q, f);
    assert(x == q);
    st_write.stop();
    writer_pos++;
  }
  fclose(f);
}

void thread_solve(int id, int q) {
  double tmp[2 * q];

  while (true) {
    PII p = buffer_sig.get();
    if (p.ST < 0) break;

    int bid = p.ST;
    int sig = p.ND;

    solvers[id].solve_sig2(bid, sig, tmp, 2 * q);
    REP(i, q) {
      if ((tmp[i] > 0 && tmp[i+q] > 0) || (tmp[i] < 0 && tmp[q+i] < 0)) response[q * sig + i] = 0.0;
      else response[q * sig + i] = (tmp[i] - tmp[q+i]) / 2.0;
    }
    if (bid < BUFFER/2) buffer_ids1.add(MP(bid, bid));
    else buffer_ids2.add(MP(bid, bid));
    {
      std::unique_lock<std::mutex> l(lock_writer);
      written_flag[sig] = 1;
      cond_writer.notify_one();
    }
  }
  double finish_time = Time::get_time();
  fprintf(stderr, "thread %d finished %.6lf\n", id, finish_time);
}

void getWTKScomb(vector <string> up, vector <string> down) {
  int q = SIZE(up);
  assert(SIZE(up) == SIZE(down));
  int q1 = parse(up, 0);
  assert(q1 == q);
  int q2 = parse(down, q);
  assert(q2 == q);

  response = new float[CALC_SIG * q]; //TODO is the memory set to 0?

  prepare_queries(2 * q);
  st_prep.stop();

  FOR(bid, 1, BUFFER/2) buffer_ids1.add(MP(bid,bid));
  FOR(bid, BUFFER/2+1, BUFFER) buffer_ids2.add(MP(bid,bid));
  thread reader1 = thread(thread_wrapper, 0);
  thread reader2 = thread(thread_wrapper, 1);

  if (NUM_THREADS == 1) {
    cerr << "We only have one thread." << endl;
    thread_solve(0, q);
    reader1.join();
    reader2.join();
  } else {
    cerr << "We have more than 1 thread." << endl;
    thread threads[MAX_THREADS];
    REP(i, NUM_THREADS) threads[i] = thread(thread_solve, i, q);
    reader1.join();
    thread writer = thread(thread_writer, q);
    reader2.join();
    writer.join();
    REP(i, NUM_THREADS) threads[i].join();
  }
}


VS load_queries(string filename) {
  DB(filename);
  FILE *f = fopen(filename.c_str(), "r");
  cerr << "Made it past the fopen call" << endl;
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

VD load_ground_truth(string filename) {
  FILE *f = fopen(filename.c_str(), "r");
  VD res;
  fgets(line, MAX_LINE, f);
  while (fgets(line, MAX_LINE, f)) {
    int pos = 0;
    while (line[pos] != ',') pos++;
    pos++;
    int i = pos;
    while (line[i]) {
      if (line[i] == ',') line[i] = ' ';
      i++;
    }
    istringstream in(line+pos);
    double x;
    while (in >> x) res.PB(x);
  }
  assert(SIZE(res) == 250 * NUM_SIG);
  return res;
}

float truth[250 * NUM_SIG];
float mytmp[1000000];

int main(int argc, char **argv){
  st_full.start();
  st_prep.start();

  for (int i = 1; i < argc; ++i) {
    if (string(argv[i]) == "--calc-sig") {
      ++i;
      CALC_SIG = atoi(argv[i]);
    } else if (string(argv[i]) == "--test") {
      TEST = 1;
    } else if (string(argv[i]) == "--out") {
      ++i;
      OUT_FILE = argv[i];
    } else if (string(argv[i]) == "--down") {
      ++i;
      DOWN_FILE = argv[i];
    } else if (string(argv[i]) == "--up") {
      ++i;
      UP_FILE = argv[i];
    } else if (string(argv[i]) == "--path2") {
      ++i;
      PATH2 = argv[i];
    } else if (string(argv[i]) == "--path") {
      ++i;
      PATH = argv[i];
    } else if (string(argv[i]) == "--batch-size") {
      ++i;
      BATCH_SIZE = atoi(argv[i]);
    } else if (string(argv[i]) == "--num-threads") {
      ++i;
      NUM_THREADS = atoi(argv[i]);
    } else {
      fprintf(stderr, "Unknown option %s\n", argv[i]);
      exit(0);
    }
  }
  DB(NUM_THREADS);
  DB(BATCH_SIZE);
  DB(PATH);
  DB(UP_FILE);
  DB(DOWN_FILE);
  DB(OUT_FILE);
  DB(CALC_SIG);
  assert(NUM_THREADS <= MAX_THREADS);
  cerr << "Made it past the assert call" << endl;
  FILE *f = fopen((PATH+"/gene_ids.txt").c_str(), "r");
  cerr << "Made it past the file open gene_ids.txt call" << endl;
  VI v;
  int x;
  cerr << "Made it to the fscan call" << endl;
  while (fscanf(f, "%d", &x) == 1) v.PB(x);
  assert(SIZE(v) == NUM_GENES);
  
  VS up = load_queries(UP_FILE);
  VS down = load_queries(DOWN_FILE);
  cerr << "Made it past the load_queries call" << endl;
  
if (!TEST) {
    cerr << "Apparently we are testing?" << endl;
    f = fopen((PATH+"/fp32_truth250_data.bin").c_str(), "r");
    x = fread(truth, sizeof(float), 250 * NUM_SIG, f);
    assert(x == 250 * NUM_SIG);
    REP(i,5) REP(j,5) fprintf(stderr, "i=%d j=%d truth=%.6lf\n", i, j, truth[i*250+j]);
  }
  
  cerr << "Starting 'init_genes' function" << endl;
  init_genes(v);
  cerr << "Finished 'init_genes' function. Starting 'getWTKScomb'..." << endl;
  getWTKScomb(up, down);
  cerr << "getWTKScomb completed." << endl;
  st_full.stop();
  if (!TEST) {
    int good = 0;
    int bad = 0;
    int res_size = CALC_SIG * SIZE(up);
    REP(i, res_size) {
      int qq = i / SIZE(up);
      int qi = (i % 250);
      float tr = truth[qq * 250 + qi];
      int ok = fabs(response[i] - tr) < 0.001;
//      if (!ok) fprintf(stderr, "discrepancy for res=%.6lf truth=%.6lf\n", response[i], truth[i]);
      if (!ok && response[i] != 0.0 && tr != 0.0) bad++;
      good += ok;
    }
    DB(bad);
    DB(good);
    DB(100.0 * good / res_size);
  }

  DB(SIZE(stoper_pointers));
  for (auto s : stoper_pointers) fprintf(stderr, "%s time %.6lf\n", s->name.c_str(), s->used_time);
  fflush(stderr);
  return 0;
}
