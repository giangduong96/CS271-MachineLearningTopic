/**
* Name: Giang Duong, ID 014533857
* a. Solve the Zodiac 408 by HMM with 1000 random restart of HMM and 200 iterations
* b. Repeat a but use 10,000 restarts
* c. Repeat b but use 100,000 restarts
* d. Repeat c but use 1,000,000 restarts
*/

#include "q15.h"

// Importing stdlib for generating random integer
#include<stdlib.h>
#include<time.h>

int randomSubstitution = 0;
//
// initialize pi[], A[][] and B[][]
//
void initMatriceshasA(double pi[],
                      double A[][N],
                      double B[][M],
                      int seed)
{
    int i,
            j;

    double prob,
            ftemp,
            ftemp2;

    // initialize pseudo-random number generator
    srandom(seed);

    // initialize pi
    prob = 1.0 / (double)N;
    ftemp = prob / 10.0;
    ftemp2 = 0.0;
    for(i = 0; i < N; ++i)
    {
        if((random() & 0x1) == 0)
        {
            pi[i] = prob + (double)(random() & 0x7) / 8.0 * ftemp;
        }
        else
        {
            pi[i] = prob - (double)(random() & 0x7) / 8.0 * ftemp;
        }
        ftemp2 += pi[i];

    }// next i

    for(i = 0; i < N; ++i)
    {
        pi[i] /= ftemp2;
    }

    // initialize A[][]
    prob = 1.0 / (double)N;
    ftemp = prob / 10.0;
    for(i = 0; i < N; ++i)
    {
        ftemp2 = 0.0;
        for(j = 0; j < N; ++j)
        {
            if((random() & 0x1) == 0)
            {
                A[i][j] = prob + (double)(random() & 0x7) / 8.0 * ftemp;
            }
            else
            {
                A[i][j] = prob - (double)(random() & 0x7) / 8.0 * ftemp;
            }
            ftemp2 += A[i][j];

        }// next j

        for(j = 0; j < N; ++j)
        {
            A[i][j] /= ftemp2;
        }

    }// next i

//    A[0][0] = .275;
//    A[0][1] = .725;
//    A[1][0] = .780;
//    A[1][1] = .220;

    // initialize B[][]
    prob = 1.0 / (double)M;
    ftemp = prob / 10.0;
    for(i = 0; i < N; ++i)
    {
        ftemp2 = 0.0;
        for(j = 0; j < M; ++j)
        {
            if((random() & 0x1) == 0)
            {
                B[i][j] = prob + (double)(random() & 0x7) / 8.0 * ftemp;
            }
            else
            {
                B[i][j] = prob - (double)(random() & 0x7) / 8.0 * ftemp;
            }
            ftemp2 += B[i][j];

        }// next j

        for(j = 0; j < M; ++j)
        {
            B[i][j] /= ftemp2;
        }

    }// next i

}// end initMatrices

// generating the digraph. part c - question 11.
void createAmatrix(double A[][N],
                              struct stepStruct *step,
                              int T
){

    for(int i = 0 ; i < 26 ; i++)
        for(int j = 0 ; j < 26 ; j++)
            A[i][j] = 0.0;

    for(int i = 0 ; i < 2 ; i++)
        for(int t = 1 ; t < T ; t++)
            A[step[t-1].obs][step[t].obs]++;

    // Adding 5 and Normalizing
    for(int i = 0 ; i < 26 ; i++){

        double row_sum = 0.0;
        for(int j = 0 ; j < 26 ; j++){
            row_sum += A[i][j];
        }
        row_sum += 130;

        double row_sto = 0.0;
        for(int j = 0 ; j < 26 ; j++){
            A[i][j] = (A[i][j] + 5) / row_sum;
            row_sto += A[i][j];
        }
        printf("\n The sum for row %d is %f", i+1, row_sto);
    }
    /* For Printing the new Matrix A
    printf("\n\nThe New A Matrix \n\n");
    for(int i = 0 ; i < 26 ; i++){
        for(int j = 0 ; j < 26 ; j++)
            printf("%f\t", A[i][j]);
        printf("\n");
    }
    */
    printf("\n\n");
}

int main(int argc, const char *argv[])
{
    int startPos,
            startChar,
            maxChars,
            maxIters,
            i,
            j,
            T,
            iter,
            totalrestart;

    int seed;

    int Z408[24][17] = {
            { 1, 2, 3, 4, 5, 4, 6, 7, 2, 8, 9,10,11,12,13,11, 7},
            {14,15,16,17,18,19,20,21, 1,22, 3,23,24,25,26,19,17},
            {27,28,19,29, 6,30, 8,31,26,32,33,34,35,19,36,37,38},
            {39,40, 4, 1, 2, 7, 3, 9,10,41, 6, 2,42,10,43,26,44},
            { 8,29,45,27, 5,28,46,47,48,12,20,22,15,14,17,31,19},
            {23,16,26,18,36, 1,24,30,38,21,26,13,49,37,50,39,40},
            {10,34,33,30,19,44,43, 9, 1,26,18, 7,32,21,39, 2, 7},
            {45,46, 4, 3, 2, 7,23,13,26,44,22,27, 6,29,10,10, 8},
            {51, 5,24,26,12,30,38,14,26,25,49,37,45,27,47, 1,52},
            { 7, 3,36,10,16,28,11,21,48,34,40,17,44, 6,22, 8,20},
            { 5,51,12, 9,15,14,30,37,16,33,45,38,43,29,10,21,22},
            {30, 1,36,10,53,32,19,47,48,46,17, 4,23,13,28,35,41},
            { 3,37,27,49,10, 6,33, 2,45,38,34,15,44,24,22,11,18},
            {47,30,25,28, 8,37, 1,49,45,27,43,34,41,38, 5,40, 3},
            {50, 6,12, 8,41, 1,52, 7,15,14,48,16,15,32,33, 9, 3},
            {29,11,39,47,43,42, 6,17,21,31,36,50,18, 2, 2,25,27},
            {34, 8,38,39,51,44, 4, 1, 2, 2, 5,42,41, 3,52, 7,15},
            {12,17,13,26,14,26,53,20,52,49,51,16,23, 1,41, 1, 7},
            { 2, 9,32,37,10, 6,51,16,53,46,19,26,53,29,39,26,14},
            {15, 5,17,18,19,24,44,53,32,19,41, 1, 2,52,45,33,53},
            {22,25,20, 7,13, 1,50,13,41,36,46,48,31,45,25,11,26},
            {53,17,46,52,52,21,17,37, 3, 9,10,13,35,20, 2,18,51},
            { 5,23,28,32,33,26,53,49,28,30,16,47, 7, 3,35,14,21},
            {15,44,13,47, 1,14,30,21,26,44,22,27,38,11,19,30, 8}
    };
    int Z408P[408] =
            {
                    8,11, 8,10, 4,10, 8,11,11, 8,13, 6,15, 4,14,15,11,
                    4, 1, 4, 2, 0,20,18, 4, 8,19, 8,18,18,14,12,20, 2,
                    7, 5,20,13, 8,19, 8,18,12,14,17, 4, 5,20,13,19, 7,
                    0,13,10, 8,11,11, 8,13, 6,22, 8,11, 3, 6, 0,12, 4,
                    8,13,19, 7, 4, 5,14,17,17, 4,18,19, 1, 4, 2, 0,20,
                    18, 4,12, 0,13, 8,18,19, 7, 4,12,14,18,19, 3, 0,13,
                    6, 4,17,14,20, 4, 0,13, 0,12, 0,11,14, 5, 0,11,11,
                    19,14,10, 8,11,11,18,14,12, 4,19, 7, 8,13, 6, 6, 8,
                    21, 4,18,12, 4,19, 7, 4,12,14,18,19,19, 7,17, 8,11,
                    11, 8,13, 6, 4,23,15, 4,17, 4,13, 2, 4, 8,19, 8,18,
                    4,21, 4,13, 1, 4,19,19, 4,17,19, 7, 0,13, 6, 4,19,
                    19, 8,13, 6,24,14,20,17,17,14, 2,10,18,14, 5, 5,22,
                    8,19, 7, 0, 6, 8,17,11,19, 7, 4, 1, 4,18,19,15, 0,
                    17,19,14, 5, 8,19, 8,18,19, 7, 0, 4,22, 7, 4,13, 8,
                    3, 8, 4, 8,22, 8,11,11, 1, 4,17, 4, 1,14,17,13, 8,
                    13,15, 0,17, 0, 3, 8, 2, 4, 0,13, 3, 0,11,11,19, 7,
                    4, 8, 7, 0,21, 4,10, 8,11,11, 4, 3,22, 8,11,11, 1,
                    4, 2,14,12, 4,12,24,18,11, 0,21, 4,18, 8,22, 8,11,
                    11,13,14,19, 6, 8,21, 4,24,14,20,12,24,13, 0,12, 4,
                    1, 4, 2, 0,20,18, 4,24,14,20,22, 8,11,11,19,17,24,
                    19,14,18,11,14, 8, 3,14,22,13,14,17, 0,19,14,15,12,
                    24, 2,14,11,11, 4, 2,19, 8,14, 6,14, 5,18,11, 0,21,
                    4,18, 5,14,17,12,24, 0, 5,19, 4,17,11, 8, 5, 4, 4,
                    1, 4,14,17, 8, 4,19, 4,12, 4,19, 7, 7,15, 8,19, 8,
            };
    double logProb,
            newLogProb;

    double pi[N],
            piBar[N],
            A[N][N],
            Abar[N][N],
            B[N][M],
            Bbar[N][M];

    char fname[80];

    struct stepStruct *step;

    srand(time(NULL));
    T = 408;

    if(argc !=8)
    {
        oops:   fprintf(stderr, "\nUsage: %s filename startPos startChar maxChars maxIters seed #restarts\n\n", argv[0]);
        fprintf(stderr, "where filename == input file\n");
        fprintf(stderr, "      startPos == starting position for each line (numbered from 0)\n");
        fprintf(stderr, "      startChar == starting character in file (numbered from 0)\n");
        fprintf(stderr, "      maxChars == max characters to read (<= 0 to read all)\n");
        fprintf(stderr, "      maxIters == max iterations of re-estimation algorithm\n");
        fprintf(stderr, "      seed == seed value for pseudo-random number generator (PRNG)\n\n");
        fprintf(stderr, "For example:\n\n      %s datafile 0 0 0 100 1241\n\n", argv[0]);
        fprintf(stderr, "will read all of `datafile' and perform a maximum of 100 iterations.\n\n");
        fprintf(stderr, "For the English text example, try:\n\n      %s BrownCorpus 15 1000 50000 200 22761\n\n", argv[0]);
        fprintf(stderr, "will read from `BrownCorpus' and seed the PRNG with 22761,\n");
        fprintf(stderr, "will not read characters 0 thru 14 of each new line in `BrownCorpus',\n");
        fprintf(stderr, "will not save the first 1000 characters read and\n");
        fprintf(stderr, "will save a maximum of 50k observations.\n\n");
        exit(0);
    }

//    sprintf(fname, argv[1]);
    strcpy(fname, argv[1]);
    startPos = atoi(argv[2]);
    startChar = atoi(argv[3]);
    maxChars = atoi(argv[4]);
    maxIters = atoi(argv[5]);
    seed = atoi(argv[6]);
    totalrestart = atoi(argv[7]);


    ////////////////////////
    // read the data file //
    ////////////////////////

    // determine number of observations
    printf("GetT... ");
    fflush(stdout);
    T = GetT(fname,
             startPos,
             startChar,
             maxChars);
    printf("T = %d\n", T);

    // allocate memory
    printf("allocating %lu bytes of memory... ", (T + 1) * sizeof(struct stepStruct));
    fflush(stdout);
    if((step = calloc(T + 1, sizeof(struct stepStruct))) == NULL)
    {
        fprintf(stderr, "\nUnable to allocate alpha\n\n");
        exit(0);
    }
    printf("done\n");

    // read in the observations from file
    printf("GetObservations... ");
    fflush(stdout);
    T = GetObservations(fname,
                        step,
                        T,
                        startPos,
                        startChar,
                        maxChars,
                        0);

    printf("T = %d\n", T);

    // initialize pi[], A[][] and B[][]
    initMatriceshasA(pi, A, B, seed);
    // The new A
    createAmatrix(A, step, T);

    T = 408;
    printf("allocating %lu bytes of memory... ", (T + 1) * sizeof(struct stepStruct));
    fflush(stdout);
    if((step = calloc(T + 1, sizeof(struct stepStruct))) == NULL)
    {
        fprintf(stderr, "\nUnable to allocate alpha\n\n");
        exit(0);
    }
    printf("done\n");

    // converting observation into the step struct
    int tempCount = 0;
    for(int r = 0 ; r < 24 ; r++)
        for(int c = 0 ; c < 17 ; c++)
            step[tempCount++].obs = Z408[r][c];

    /////////////////////////
    // hidden markov model //
    /////////////////////////

    srandom(seed);
    int maxMatchCount = 0;
    int hmmNumber = 0;
    // initialize pi[], A[][] and B[][]
    initMatrices(pi, B, seed);
    for(int iterCnt = 0 ; iterCnt < totalrestart ; iterCnt++){
        T = 408;
        // initialization
        iter = 0;
        logProb = -1.0;
        newLogProb = 0.0;

        // main loop
        //while((iter < maxIters) && (newLogProb > logProb))
        while(iter < maxIters)
        {
            logProb = newLogProb;
            alphaPass(step, pi, A, B, T);
            betaPass(step, pi, A, B, T);
            computeGammas(step, pi, A, B, T);
            reestimatePi(step, piBar);
            reestimateB(step, Bbar, T);
            // assign pi, A and B corresponding "bar" values
            for(i = 0; i < N; ++i)
            {
                pi[i] = piBar[i];

                for(j = 0; j < M; ++j)
                {
                    B[i][j] = Bbar[i][j];
                }

            }// next i

            // compute log [P(observations | lambda)], where lambda = (A,B,pi)
            newLogProb = 0.0;
            for(i = 0; i < T; ++i)
            {
                newLogProb += log(step[i].c);
            }
            newLogProb = -newLogProb;

            // a little trick so that no initial logProb is required
            if(iter == 0)
            {
                logProb = newLogProb - 1.0;
            }


            ++iter;
        }// end while


        printf("\nfinal B^T =\n");
        printBT(B);


        // Calculating the Putative Key part of the problem - part d Question 11
        int count = 0;
        for(int i = 0 ; i < M ; i++){
            double max = B[i][0];
            int sum = 0;
            int element = 0;
            for(int j = 1 ; j < N ; j++){
                if(B[i][j] > 0){
                    sum += B[i][j];
                }
                if(B[i][j] > max) {
                    max = B[i][j];
                    element = j % 26;
                }
            }
            printf("\nThe %d element is mapped to %d", i, element);


            if(element == (j % 26))
                count++;
            if((i == 0) && ((element ==  5) || (element == 7) || (element== 8) || (element ==6))) count++;
            if((i == 1) && (element ==  9)) count++;
            if((i == 2) && (element ==  10)) count++;
            if((i == 3) && (element ==  4 +3)) count++;
            if((i == 4) && ((element ==  8) ||(element == 9 ) || (element == 6))) count++;
            if((i == 5) && ((element ==  6) || (element ==  4))) count++;
            if((i == 6) && ((element == 11))) count++;
            if((i == 7) && ((element == 8))) count++;
            if((i == 8) && ((element ==  8)||(element ==  10)|(element ==  14)||(element ==  11))) count++;
            if((i == 9) && ((element ==  8))) count++;
            if((i == 10) && ((element ==  8))) count++;
            if((i == 11) && ((element ==  12)||(element ==  14)||(element ==  17))) count++;
            if((i == 12) && ((element ==  17)))count++;
            if((i == 13) && ((element ==  4)||(element ==  7)||(element ==  6))) count++;
            if((i == 14) && ((element ==  6)||(element ==  7)||(element ==  9)) )count++;
            if((i == 15) && ((element ==  7))) count++;
            if((i == 16) && ((element ==  0))) count++;
            if((i == 17) && ((element ==  7)||(element ==  5))) count++;
            if((i == 18) && ((element ==  3)||(element ==  7)||(element ==  6))) count++;
            if((i == 19) && ((element ==  7)||(element ==  8)||(element ==  10))) count++;
            if((i == 20) && ((element ==  10))) count++;
            if((i == 21) && ((element ==  6))) count++;
            if((i == 22) && ((element ==  8))) count++;
            if((i == 23) && ((element ==  1)))count++;
            if((i == 24) && ((element ==  8))) count++;
            if((i == 25) && ((element ==  0))) count++;
            //if()
        }
        //printf("\n\n\n\n\n\nThe key is %d", randomSubstitution);
        printf("\nFinal Count value: %d", count);
        double percentage = (double) count / N;
        printf("\nThe percentage is %0.4f\n", percentage*100);

        if(count > maxMatchCount){
            maxMatchCount = count;
            hmmNumber = iterCnt;
        }
    }

    printf("\n\nMaximum Match count %d from hmm number %d\n\n", maxMatchCount, hmmNumber);

}// end hmm


//
// alpha pass (or forward pass) including scaling
//
void alphaPass(struct stepStruct *step,
               double pi[],
               double A[][N],
               double B[][M],
               int T)
{
    int i,
            j,
            t;

    double ftemp;

    // compute alpha[0]'s
    ftemp = 0.0;
    for(i = 0; i < N; ++i)
    {
        step[0].alpha[i] = pi[i] * B[i][step[0].obs];
        ftemp += step[0].alpha[i];
    }
    step[0].c = 1.0 / ftemp;

    // scale alpha[0]'s
    for(i = 0; i < N; ++i)
    {
        step[0].alpha[i] /= ftemp;
    }

    // alpha pass
    for(t = 1; t < T; ++t)
    {
        ftemp = 0.0;
        for(i = 0; i < N; ++i)
        {
            step[t].alpha[i] = 0.0;
            for(j = 0; j < N; ++j)
            {
                step[t].alpha[i] += step[t - 1].alpha[j] * A[j][i];
            }
            step[t].alpha[i] *= B[i][step[t].obs];
            ftemp += step[t].alpha[i];
        }
        step[t].c = 1.0 / ftemp;

        // scale alpha's
        for(i = 0; i < N; ++i)
        {
            step[t].alpha[i] /= ftemp;
        }

    }// next t

}// end alphaPass


//
// beta pass (or backwards pass) including scaling
//
void betaPass(struct stepStruct *step,
              double pi[],
              double A[][N],
              double B[][M],
              int T)
{
    int i,
            j,
            t;

    // compute scaled beta[T - 1]'s
    for(i = 0; i < N; ++i)
    {
        step[T - 1].beta[i] = 1.0 * step[T - 1].c;
    }

    // beta pass
    for(t = T - 2; t >= 0; --t)
    {
        for(i = 0; i < N; ++i)
        {
            step[t].beta[i] = 0.0;
            for(j = 0; j < N; ++j)
            {
                step[t].beta[i] += A[i][j] * B[j][step[t + 1].obs] * step[t + 1].beta[j];
            }

            // scale beta's (same scale factor as alpha's)
            step[t].beta[i] *= step[t].c;
        }

    }// next t

}// end betaPass


//
// compute gamma's and diGamma's including optional error checking
//
void computeGammas(struct stepStruct *step,
                   double pi[],
                   double A[][N],
                   double B[][M],
                   int T)
{
    int i,
            j,
            t;

    double denom;

#ifdef CHECK_GAMMAS
    double ftemp,
           ftemp2;
#endif // CHECK_GAMMAS

    // compute gamma's and diGamma's
    for(t = 0; t < T - 1; ++t)// t = 0,1,2,...,T-2
    {

#ifdef CHECK_GAMMAS
        ftemp2 = 0.0;
#endif // CHECK_GAMMAS

        for(i = 0; i < N; ++i)
        {
            step[t].gamma[i] = 0.0;
            for(j = 0; j < N; ++j)
            {
                step[t].diGamma[i][j] = (step[t].alpha[i] * A[i][j] * B[j][step[t + 1].obs] * step[t + 1].beta[j]);
                step[t].gamma[i] += step[t].diGamma[i][j];
            }

#ifdef CHECK_GAMMAS
            // verify that gamma[i] == alpha[i]*beta[i] / sum(alpha[j]*beta[j])
            ftemp2 += step[t].gamma[i];
            ftemp = 0.0;
            for(j = 0; j < N; ++j)
            {
                ftemp += step[t].alpha[j] * step[t].beta[j];
            }
            ftemp = (step[t].alpha[i] * step[t].beta[i]) / ftemp;
            if(DABS(ftemp - step[t].gamma[i]) > EPSILON)
            {
                printf("gamma[%d] = %f (%f) ", i, step[t].gamma[i], ftemp);
                printf("********** Error !!!\n");
            }
#endif // CHECK_GAMMAS

        }// next i

#ifdef CHECK_GAMMAS
        if(DABS(1.0 - ftemp2) > EPSILON)
        {
            printf("sum of gamma's = %f (should sum to 1.0)\n", ftemp2);
        }
#endif // CHECK_GAMMAS

    }// next t

    // special case for t = T-1
    for(j = 0; j < N; ++j)
    {
        step[T-1].gamma[j] = step[T-1].alpha[j];
    }

}// end computeGammas


//
// reestimate pi, the initial distribution
//
void reestimatePi(struct stepStruct *step,
                  double piBar[])
{
    int i;

    // reestimate pi[]
    for(i = 0; i < N; ++i)
    {
        piBar[i] = step[0].gamma[i];
    }

}// end reestimatePi


//
// reestimate the A matrix
//
void reestimateA(struct stepStruct *step,
                 double Abar[][N],
                 int T)
{
    int i,
            j,
            t;

    double numer,
            denom;

    // reestimate A[][]
    for(i = 0; i < N; ++i)
    {
        for(j = 0; j < N; ++j)
        {
            numer = denom = 0.0;

            // t = 0,1,2,...,T-2
            for(t = 0; t < T - 1; ++t)
            {
                numer += step[t].diGamma[i][j];
                denom += step[t].gamma[i];

            }// next t

            Abar[i][j] = numer / denom;

        }// next j

    }// next i

} // end reestimateA


//
// reestimate the B matrix
//
void reestimateB(struct stepStruct *step,
                 double Bbar[][M],
                 int T)
{
    int i,
            j,
            t;

    double numer,
            denom;

    // reestimate B[][]
    for(i = 0; i < N; ++i)
    {
        for(j = 0; j < M; ++j)
        {
            numer = denom = 0.0;

            // t = 0,1,2,...,T-1
            for(t = 0; t < T; ++t)
            {
                if(step[t].obs == j)
                {
                    numer += step[t].gamma[i];
                }
                denom += step[t].gamma[i];

            }// next t

            Bbar[i][j] = numer / denom;

        }// next j

    }// next i

}// end reestimateB


//
// initialize pi[], A[][] and B[][]
//
void initMatrices(double pi[],
                  double B[][M],
                  int seed)
{
    int i,
            j;

    double prob,
            ftemp,
            ftemp2;

    // initialize pseudo-random number generator
    srandom(seed);

    // initialize pi
    prob = 1.0 / (double)N;
    ftemp = prob / 10.0;
    ftemp2 = 0.0;
    for(i = 0; i < N; ++i)
    {
        if((random() & 0x1) == 0)
        {
            pi[i] = prob + (double)(random() & 0x7) / 8.0 * ftemp;
        }
        else
        {
            pi[i] = prob - (double)(random() & 0x7) / 8.0 * ftemp;
        }
        ftemp2 += pi[i];

    }// next i

    for(i = 0; i < N; ++i)
    {
        pi[i] /= ftemp2;
    }

    // initialize B[][]
    prob = 1.0 / (double)M;
    ftemp = prob / 10.0;
    for(i = 0; i < N; ++i)
    {
        ftemp2 = 0.0;
        for(j = 0; j < M; ++j)
        {
            if((random() & 0x1) == 0)
            {
                B[i][j] = prob + (double)(random() & 0x7) / 8.0 * ftemp;
            }
            else
            {
                B[i][j] = prob - (double)(random() & 0x7) / 8.0 * ftemp;
            }
            ftemp2 += B[i][j];

        }// next j

        for(j = 0; j < M; ++j)
        {
            B[i][j] /= ftemp2;
        }

    }// next i

}// end initMatrices


//
// read (but don't save) observations get T
//
int GetT(char fname[],
         int startPos,
         int startChar,
         int maxChars)
{
    FILE *in;

    int T,
            i,
            j,
            len,
            thisStartPos,
            totalNum,
            num;

    char temp[MAX_CHARS + 1];

    char space[1] = {" "};

    char alphabet[M] = ALPHABET;

    in = fopen(fname, "r");
    if(in == NULL)
    {
        fprintf(stderr, "\nError opening file %s\n\n", fname);
        exit(0);
    }

#ifdef PRINT_GET_T
    printf("\n");
#endif // PRINT_GET_T

    // count 'em
    totalNum = num = 0;
    while(fgets(temp, MAX_CHARS, in) != NULL)
    {
        len = strlen(temp);

        // each line should end with a single space
        while((strncmp(&temp[len - 1], space, 1) == 0) && (len > 0))
        {
            --len;
        }

        strncpy(&temp[len], space, 1);

        thisStartPos = startPos;

        // ignore leading spaces
        while((strncmp(&temp[thisStartPos], space, 1) == 0) && (thisStartPos < len))
        {
            ++thisStartPos;
        }

        for(i = thisStartPos; i <= len; ++i)
        {
            // find alphabetic characters, ignoring case
            // also drop all non-alphabet characters other than space
            for(j = 0; j < M; ++j)
            {
                if(strncasecmp(&temp[i], &alphabet[j], 1) == 0)
                {
                    ++totalNum;
                    if(totalNum >= startChar)
                    {
#ifdef PRINT_GET_T
                        printf("%c %d\n", alphabet[j], num);
#endif // PRINT_GET_T
                        ++num;
                        if((maxChars > 0) && (num >= maxChars))
                        {
                            return(num);
                        }

                    }// end if

                    break;

                }// end if

            }// next j

        }// next i

    }// end while

    fclose(in);

    return(num);

}// end GetT


//
// read and save observations
//
int GetObservations(char fname[],
                    struct stepStruct *step,
                    int T,
                    int startPos,
                    int startChar,
                    int maxChars,
                    int flag)
{
    FILE *in;

    int i,
            j,
            len,
            num,
            thisStartPos,
            totalNum;

    char temp[MAX_CHARS + 2];

    char space[1] = {" "};

    char alphabet[M] = ALPHABET;

    // Part a - Question 11. Encrypting the key
    randomSubstitution = 0;
    if (flag == 1)
    {
        while(randomSubstitution == 0){
            randomSubstitution = rand() % 10;
        }
    }

    in = fopen(fname, "r");
    if(in == NULL)
    {
        fprintf(stderr, "\nError opening file %s\n\n", fname);
        exit(0);
    }

#ifdef PRINT_OBS
    printf("\n");
#endif // PRINT_OBS

    // read 'em in
    totalNum = num = 0;
    while(fgets(temp, MAX_CHARS, in) != NULL)
    {
        len = strlen(temp);

        // each line should end with a single space
        while((strncmp(&temp[len - 1], space, 1) == 0) && (len > 0))
        {
            --len;
        }
        strncpy(&temp[len], space, 1);

        thisStartPos = startPos;

        // ignore leading spaces
        while((strncmp(&temp[thisStartPos], space, 1) == 0) && (thisStartPos < len))
        {
            ++thisStartPos;
        }

        for(i = thisStartPos; i <= len; ++i)
        {
            // find alphabetic characters, ignoring case,
            // drop all non-alphabet characters other than space
            for(j = 0; j < M; ++j)
            {
                if(strncasecmp(&temp[i], &alphabet[j], 1) == 0)
                {
                    ++totalNum;
                    if(totalNum >= startChar)
                    {
                        if(flag == 11)
                            step[num].obs = (j + randomSubstitution) % 26;
                        else
                            step[num].obs = j;


#ifdef PRINT_OBS
                        printf("%c %d\n", alphabet[j], step[num].obs);
#endif // PRINT_OBS
                        ++num;
                        if(num > T)
                        {
                            printf("\nError --- T exceeded in GetObservations()\n\n");
                            exit(0);
                        }
                        if((maxChars > 0) && (num >= maxChars))
                        {
                            return(num);
                        }

                    }// end if

                    break;

                }// end if

            }// next j

        }// next i

    }// end while

    fclose(in);

    return(num);

}// end GetObservations


//
// print pi[]
//
void printPi(double pi[])
{
    int i;

    double ftemp;

    ftemp = 0.0;
    for(i = 0; i < N; ++i)
    {
        printf("%8.5f ", pi[i]);
        ftemp += pi[i];
    }
    printf(",  sum = %f\n", ftemp);

}// end printPi


//
// print A[][]
//
void printA(double A[][N])
{
    int i,
            j;

    double ftemp;

    for(i = 0; i < N; ++i)
    {
        ftemp = 0.0;
        for(j = 0; j < N; ++j)
        {
            printf("%8.5f ", A[i][j]);
            ftemp += A[i][j];
        }
        printf(",  sum = %f\n", ftemp);

    }// next i

}// end printA


//
// print BT[][]
//
void printBT(double B[][M])
{
    int i,
            j;

    double ftemp;

    char alphabet[M] = ALPHABET;

    for(i = 0; i < M; ++i)
    {
        //printf("%c ", alphabet[i]);
        for(j = 0; j < N; ++j)
        {
            printf("%5.3f ", B[j][i]);
        }
        printf("\n");
    }
    for(i = 0; i < N; ++i)
    {
        ftemp = 0.0;
        for(j = 0; j < M; ++j)
        {
            ftemp += B[i][j];
        }
        //printf("sum[%d] = %f ", i, ftemp);
    }
    printf("\n");

}// end printB
