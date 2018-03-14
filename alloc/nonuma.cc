#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include <float.h>
#include <limits.h>
#include <malloc.h>
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

void readProc(std::vector<double> &a, double &result) {
  double total = 0.0;

  for (int i = 0; i < REPEAT; i++) {
    for (int j = 0; j < N; j += STRIDE) {
      total += a[j];
    }
  }

  result = total;
}

void copyProc(std::vector<double> &a, std::vector<double> &c) {
  for (int i = 0; i < REPEAT; i++) {
    for (int j = 0; j < N; j += STRIDE) {
      c[j] = a[j];
    }
  }
}

void scaleProc(std::vector<double> &b, std::vector<double> &c) {
  for (int i = 0; i < REPEAT; i++) {
    for (int j = 0; j < N; j += STRIDE) {
      b[j] = 3.0 * c[j];
    }
  }
}

void addProc(std::vector<double> &a, std::vector<double> &b,
             std::vector<double> &c) {
  for (int i = 0; i < REPEAT; i++) {
    for (int j = 0; j < N; j += STRIDE) {
      c[j] = a[j] + b[j];
    }
  }
}

void triadProc(std::vector<double> &a, std::vector<double> &b,
               std::vector<double> &c) {
  for (int i = 0; i < REPEAT; i++) {
    for (int j = 0; j < N; j += STRIDE) {
      a[j] = b[j] + 3.0 * c[j];
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

  std::vector<double> results(NUM_THREADS);

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

  std::vector<std::thread> threads(NUM_THREADS);
  std::vector<std::vector<double>> a, b, c;

  for (int i = 0; i < NUM_THREADS; i++) {
    a.emplace_back(N, 1.0);
    b.emplace_back(N, 2.0);
    c.emplace_back(N, 3.0);
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
      threads[i] = std::thread([&](int i) { readProc(a[i], results[i]); }, i);
    }
    for (auto &thread : threads) {
      thread.join();
    }
    times[0][k] = mysecond() - times[0][k];

    times[1][k] = mysecond();
    for (int i = 0; i < NUM_THREADS; i++) {
      threads[i] = std::thread([&](int i) { copyProc(a[i], c[i]); }, i);
    }
    for (auto &thread : threads) {
      thread.join();
    }
    times[1][k] = mysecond() - times[1][k];

    times[2][k] = mysecond();
    for (int i = 0; i < NUM_THREADS; i++) {
      threads[i] = std::thread([&](int i) { scaleProc(b[i], c[i]); }, i);
    }
    for (auto &thread : threads) {
      thread.join();
    }
    times[2][k] = mysecond() - times[2][k];

    times[3][k] = mysecond();
    for (int i = 0; i < NUM_THREADS; i++) {
      threads[i] = std::thread([&](int i) { addProc(a[i], b[i], c[i]); }, i);
    }
    for (auto &thread : threads) {
      thread.join();
    }
    times[3][k] = mysecond() - times[3][k];

    times[4][k] = mysecond();
    for (int i = 0; i < NUM_THREADS; i++) {
      threads[i] = std::thread([&](int i) { triadProc(a[i], b[i], c[i]); }, i);
    }
    for (auto &thread : threads) {
      thread.join();
    }
    times[4][k] = mysecond() - times[4][k];
  }

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
