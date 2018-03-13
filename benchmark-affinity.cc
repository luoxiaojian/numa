#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include <float.h>
#include <limits.h>
#include <malloc.h>
#include <numa.h>
#include <pthread.h>
#include <sched.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <unistd.h>

#include <thread>
#include <vector>

#ifndef NUM_THREADS
#define NUM_THREADS 8
#endif
#ifndef MEM_OFF
#define MEM_OFF 0
#endif
#ifndef STRIDE
#define STRIDE (64 / sizeof(double))
#endif
#ifndef N
#define N 4800000
#endif
#ifndef NTIMES
#define NTIMES 10
#endif
#define REPEAT (5 * STRIDE)

#define HLINE "-------------------------------------------------------------\n"

#ifndef MIN
#define MIN(x, y) ((x) < (y) ? (x) : (y))
#endif
#ifndef MAX
#define MAX(x, y) ((x) > (y) ? (x) : (y))
#endif

static double **a, **b, **c;
static double avgtime[5] = {0}, maxtime[5] = {0},
              mintime[5] = {FLT_MAX, FLT_MAX, FLT_MAX, FLT_MAX, FLT_MAX};
static char *label[5] = {"Read:      ", "Copy:      ", "Scale:     ",
                         "Add:       ", "Triad:     "};
static double bytes[5] = {
    1 * sizeof(double) * N * NUM_THREADS * REPEAT / STRIDE,
    2 * sizeof(double) * N *NUM_THREADS *REPEAT / STRIDE,
    2 * sizeof(double) * N *NUM_THREADS *REPEAT / STRIDE,
    3 * sizeof(double) * N *NUM_THREADS *REPEAT / STRIDE,
    3 * sizeof(double) * N *NUM_THREADS *REPEAT / STRIDE};

double results[NUM_THREADS];

typedef struct Arg_T {
  int proc;
  int allocNode;
  int tid;
  int cpuid;
} arg_t;

void readProc(int tid, int me) {
  double *a2 = a[tid];
  double total = 0.0;

  for (int i = 0; i < REPEAT; i++) {
    for (int j = 0; j < N; j += STRIDE) {
      total += a2[j];
    }
  }

  results[tid] = total;
}

void copyProc(int tid, int me) {
  double *a2 = a[tid];
  double *c2 = c[tid];

  for (int i = 0; i < REPEAT; i++) {
    for (int j = 0; j < N; j += STRIDE) {
      c2[j] = a2[j];
    }
  }
}

void scaleProc(int tid, int me) {
  double *b2 = b[tid];
  double *c2 = c[tid];

  for (int i = 0; i < REPEAT; i++) {
    for (int j = 0; j < N; j += STRIDE) {
      b2[j] = 3.0 * c2[j];
    }
  }
}

void addProc(int tid, int me) {
  double *a2 = a[tid];
  double *b2 = b[tid];
  double *c2 = c[tid];

  for (int i = 0; i < REPEAT; i++) {
    for (int j = 0; j < N; j += STRIDE) {
      c2[j] = a2[j] + b2[j];
    }
  }
}

void triadProc(int tid, int me) {
  double *a2 = a[tid];
  double *b2 = b[tid];
  double *c2 = c[tid];

  for (int i = 0; i < REPEAT; i++) {
    for (int j = 0; j < N; j += STRIDE) {
      a2[j] = b2[j] + 3.0 * c2[j];
    }
  }
}

double mysecond() {
  struct timeval tp;
  gettimeofday(&tp, NULL);
  return ((double)tp.tv_sec + (double)tp.tv_usec * 1.e-6);
}

#define M 20

int checktick() {
  int i, minDelta, Delta;
  double t1, t2, timesfound[M];

  for (i = 0; i < M; i++) {
    t1 = mysecond();
    while (((t2 = mysecond()) - t1) < 1.0E-6)
      ;
    timesfound[i] = t2 - t1;
  }

  minDelta = 1000000;
  for (i = 1; i < M; i++) {
    Delta = (int)(1.0E6 * (timesfound[i] - timesfound[i - 1]));
    minDelta = MIN(minDelta, MAX(Delta, 0));
  }

  return minDelta;
}

int main(int argc, char **argv) {
  int BytesPerWord, quantum;
  register int j, k;
  double t, times[5][NTIMES];

  printf(HLINE);
  BytesPerWord = sizeof(double);
  printf("This system uses %d bytes per DOUBLE PRECISION word.\n",
         BytesPerWord);

  printf(HLINE);
  printf("Array size = %llu.\n", (unsigned long long)N);
  printf("Total memory required = %.1f MB.\n",
         (3.0 * BytesPerWord) * ((double)N / 1048576.0) * NUM_THREADS);
  printf(
      "Each test is run %d times, but only the best time for each is used.\n",
      NTIMES);
  printf("Number of threads requested = %d\n", NUM_THREADS);

  printf(HLINE);
  if (numa_available() < 0) {
    printf("System does not support NUMA.\n");
    exit(1);
  }

  int num_nodes = numa_max_node() + 1;
  printf("Number of available nodes = %d\n", num_nodes);

  std::vector<std::thread> threads(NUM_THREADS);
  std::vector<arg_t> p(NUM_THREADS);

  int cpuid1 = 0;
  int cpuid2 = 12;
  for (int i = 0; i < NUM_THREADS; i++) {
    p[i].proc = i % num_nodes;
    p[i].allocNode = ((i + MEM_OFF) % num_nodes);
    p[i].tid = i;
    p[i].cpuid = (p[i].proc == 0) ? cpuid1++ : cpuid2++;
  }

  a = (double **)malloc(sizeof(double *) * NUM_THREADS);
  b = (double **)malloc(sizeof(double *) * NUM_THREADS);
  c = (double **)malloc(sizeof(double *) * NUM_THREADS);

  for (int i = 0; i < NUM_THREADS; i++) {
    a[i] = (double *)numa_alloc_onnode(sizeof(double) * N, p[i].allocNode);
    b[i] = (double *)numa_alloc_onnode(sizeof(double) * N, p[i].allocNode);
    c[i] = (double *)numa_alloc_onnode(sizeof(double) * N, p[i].allocNode);
    if ((a[i] == NULL) || (b[i] == NULL) || (c[i] == NULL)) {
      printf("Failed to allocate %d bytes on node %d\n",
             (int)sizeof(double) * N, p[i].allocNode);
      exit(-1);
    }
    for (int j = 0; j < N; j++) {
      a[i][j] = 1.0;
      b[i][j] = 2.0;
      c[i][j] = 3.0;
    }
  }

  if ((quantum = checktick()) >= 1) {
    printf("Your clock granularity/precision appears to be %d microseconds.\n",
           quantum);
  } else {
    printf("Your clock granularity appears to be less than one microsecond.\n");
    quantum = 1;
  }
  printf(HLINE);
  cpu_set_t cpuset;
  for (k = 0; k < NTIMES; k++) {
    times[0][k] = mysecond();
    for (int i = 0; i < NUM_THREADS; i++) {
      threads[i] = std::thread([&](int i) { readProc(i, p[i].proc); }, i);
      CPU_ZERO(&cpuset);
      CPU_SET(p[i].cpuid, &cpuset);
      int rc = pthread_setaffinity_np(threads[i].native_handle(),
                                      sizeof(cpu_set_t), &cpuset);
    }
    for (auto &thread : threads) {
      thread.join();
    }
    times[0][k] = mysecond() - times[0][k];

    times[1][k] = mysecond();
    for (int i = 0; i < NUM_THREADS; i++) {
      threads[i] = std::thread([&](int i) { copyProc(i, p[i].proc); }, i);
      CPU_ZERO(&cpuset);
      CPU_SET(p[i].cpuid, &cpuset);
      int rc = pthread_setaffinity_np(threads[i].native_handle(),
                                      sizeof(cpu_set_t), &cpuset);
    }
    for (auto &thread : threads) {
      thread.join();
    }
    times[1][k] = mysecond() - times[1][k];

    times[2][k] = mysecond();
    for (int i = 0; i < NUM_THREADS; i++) {
      threads[i] = std::thread([&](int i) { scaleProc(i, p[i].proc); }, i);
      CPU_ZERO(&cpuset);
      CPU_SET(p[i].cpuid, &cpuset);
      int rc = pthread_setaffinity_np(threads[i].native_handle(),
                                      sizeof(cpu_set_t), &cpuset);
    }
    for (auto &thread : threads) {
      thread.join();
    }
    times[2][k] = mysecond() - times[2][k];

    times[3][k] = mysecond();
    for (int i = 0; i < NUM_THREADS; i++) {
      threads[i] = std::thread([&](int i) { addProc(i, p[i].proc); }, i);
      CPU_ZERO(&cpuset);
      CPU_SET(p[i].cpuid, &cpuset);
      int rc = pthread_setaffinity_np(threads[i].native_handle(),
                                      sizeof(cpu_set_t), &cpuset);
    }
    for (auto &thread : threads) {
      thread.join();
    }
    times[3][k] = mysecond() - times[3][k];

    times[4][k] = mysecond();
    for (int i = 0; i < NUM_THREADS; i++) {
      threads[i] = std::thread([&](int i) { triadProc(i, p[i].proc); }, i);
      CPU_ZERO(&cpuset);
      CPU_SET(p[i].cpuid, &cpuset);
      int rc = pthread_setaffinity_np(threads[i].native_handle(),
                                      sizeof(cpu_set_t), &cpuset);
    }
    for (auto &thread : threads) {
      thread.join();
    }
    times[4][k] = mysecond() - times[4][k];
  }

  for (int i = 0; i < NUM_THREADS; i++) {
    numa_free(a[i], sizeof(double) * N);
    numa_free(b[i], sizeof(double) * N);
    numa_free(c[i], sizeof(double) * N);
  }
  free(a);
  free(b);
  free(c);

  for (k = 1; k < NTIMES; k++) {
    for (j = 0; j < 5; j++) {
      avgtime[j] = avgtime[j] + times[j][k];
      mintime[j] = MIN(mintime[j], times[j][k]);
      maxtime[j] = MAX(maxtime[j], times[j][k]);
    }
  }
  printf(
      "Function      Rate (MB/s)   Latency(ns)   Avg time     Min time     Max "
      "time\n");
  for (j = 0; j < 5; j++) {
    avgtime[j] = avgtime[j] / (double)(NTIMES - 1);
    printf("%s%11.4f  %11.4f  %11.4f  %11.4f  %11.4f\n", label[j],
           1.0E-06 * bytes[j] / avgtime[j],
           (avgtime[j] / (bytes[j] / (sizeof(double) * NUM_THREADS))) * 1.0E9,
           avgtime[j], mintime[j], maxtime[j]);
  }
  printf(HLINE);
  return 0;
}
