/**
 * Name: Giang Duong, ID 014533857
 * a. Solve the Zodiac 408 by HMM with 1000 random restart of HMM and 200 iterations
 * b. Repeat a but use 10,000 restarts
 * c. Repeat b but use 100,000 restarts
 * d. Repeat c but use 1,000,000 restarts
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// pi[N], A[N][N] and B[N][M]
#define N 26
#define M 54
//#define M 26

// character set
// Note: Size of character set must be M
#define ALPHABET {"abcdefghijklmnopqrstuvwxyz"} // Removing space as it is not needed
//#define ALPHABET {"abcdefghijklmnopqrstuvwxyz"}

// maximum characters per line
#define MAX_CHARS 500

// other
#define EPSILON 0.00001
#define DABS(x) ((x) < (0.0) ? -(x) : (x))

// debugging and/or printing
//#define PRINT_OBS
//#define PRINT_GET_T
//#define CHECK_GAMMAS
#define PRINT_REESTIMATES

struct stepStruct
{
    int obs;
    double c;
    double alpha[N];
    double beta[N];
    double gamma[N];
    double diGamma[N][N];
};

void alphaPass(struct stepStruct *step,
               double pi[],
               double A[][N],
               double B[][M],
               int T);

void betaPass(struct stepStruct *step,
              double pi[],
              double A[][N],
              double B[][M],
              int T);

void computeGammas(struct stepStruct *step,
                   double pi[],
                   double A[][N],
                   double B[][M],
                   int T);

void reestimatePi(struct stepStruct *step,
                  double piBar[]);

void reestimateA(struct stepStruct *step,
                 double Abar[][N],
                 int T);

void reestimateB(struct stepStruct *step,
                 double Bbar[][M],
                 int T);

void initMatrices(double pi[],
                  double B[][M],
                  int seed);

int GetT(char fname[],
         int startPos,
         int startChar,
         int maxChars);

int GetObservations(char fname[],
                    struct stepStruct *step,
                    int T,
                    int startPos,
                    int startChar,
                    int maxChars,
                    int flag);

void printPi(double pi[]);

void printA(double A[][N]);

void printBT(double B[][M]);
