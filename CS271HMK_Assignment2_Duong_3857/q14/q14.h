//
// Created by Giang Duong on 9/11/20.
//
/*
 * Name: Giang Duong, student ID 014533857
 * Homework 14: CS 271 Fall 2020
 * 14.	Since HMM training is a hill climb, we are only assured of reaching a local maximum. And, as with any hill climb,
 * the speciﬁc local maximum that we ﬁnd will depend on our choice of initial values. Therefore, by training a hidden
 * Markov model multiple times with diﬀerent initial values, we would expect to obtain better results than when training
 * only once.
    In the paper [16], the authors use an expectation maximization (EM) approach with multiple random restarts as a
    means of attacking ho-mophonic substitution ciphers. An analogous HMM-based technique is analyzed in the report
    [158], where the eﬀectiveness of multiple ran-dom restarts on simple substitution cryptanalysis is explored in
    detail. Multiple random restarts are especially helpful in the most challenging cases, that is, when little data
    (i.e., ciphertext) is available. However, the tradeoﬀ is that the work factor can be high, since the number of
    restarts required may be very large (millions of random restarts are required in some cases).
    a)	Obtain an English plaintext message consisting of 1000 plaintext characters, consisting only of lower case a
    through z (i.e., remove all punctuation, special characters, and spaces, and convert all upper case letters to lower
    case). Encrypt this plaintext using a randomly selected shift of the alphabet. Remember the key. Also generate a
    digraph frequency matrix 󂗀, as discussed in part c) of Problem 11.
b)	Train 󂗀 HMMs, for each of 󂗀 = 1, 󂗀 = 10, 󂗀 = 100, and 󂗀 = 1000, following the same process as in Problem 11, part
 d), but using the 󂗀 = 1000 observations generated in part a) of this problem. For a given 󂗀 select the best result
 based on the model scores and give the fraction of the putative key that is correct, calculated as in Problem 11, part
 d).
c)	Repeat part b), but only use the ﬁrst 󂗀 = 400 observations.
d)	Repeat part c), but only use the ﬁrst 󂗀 = 300 observations.

 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// pi[N], A[N][N] and B[N][M]
#define N 26
#define M 26
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
