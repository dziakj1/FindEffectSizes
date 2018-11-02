%MACRO FindEffectSizes(gammas=,
                       rhos=,
                       nstarts=200);
   /***************************************************************
    | FindEffectSizes macro                                                 
    | By John Dziak
    | Performs effect size calculations for the bootstrap likelihood 
    | ratio test in the case of latent class analysis with dichotomous
    | indicators and no covariates, as described in
    | Dziak, Lanza and Tan (submitted).
    *----------------------------------------------------------------
    | Copyright: 
    | (c) 2012 The Pennsylvania State University               
    *----------------------------------------------------------------
    | License:
    | This program is free software; you can redistribute it and/or
    | modify it under the terms of the GNU General Public License as
    | published by the Free Software Foundation; either version 2 of
    | the License, or (at your option) any later version.
    |
    | This program is distributed in the hope that it will be useful,
    | but WITHOUT ANY WARRANTY; without even the implied warranty of
    | MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
    | General Public License for more details.
    *----------------------------------------------------------------
    | Example of syntax:
    |           %INCLUDE "C:\Documents\FindEffectSizes.sas";
    |            DATA gammas;
    |                INPUT gamma;
    |                DATALINES;
    |                .6
    |                .3
    |                .1
    |                ;
    |            RUN;
    |            DATA rhos;
    |                INPUT rho01 rho02 rho03;
    |                DATALINES;
    |                .2 .8 .8
    |                .2 .8 .8
    |                .2 .8 .8
    |                .2 .2 .8
    |                .2 .2 .8    
    |            ;
    |            %FindEffectSizes(gammas=gammas,
    |                             rhos=rhos); 
    *----------------------------------------------------------------*/
     PROC IML;
          USE &gammas;
              READ ALL INTO gammas;
          CLOSE &gammas; 
          USE &rhos;
              READ ALL INTO rhos;
          CLOSE &rhos;  * Assumes dichotomous items only!;
          * Determine alternative hypothesis cell probabilities;
          nClasses = NROW(gammas);
          nItems = NROW(rhos);
          nCells = 2**nItems; 
          CALL SYMPUT("FindEffectSizesH0Classes",CHAR(nClasses-1));
          CALL SYMPUT("FindEffectSizesH1Classes",CHAR(nClasses));
          CALL SYMPUT("FindEffectSizesNItems",CHAR(nItems));
          %LET FindEffectSizesCategories = ;
          %LET FindEffectSizesItemNames = ;
          %DO FindEffectSizesTemp = 1 %TO &FindEffectSizesNItems;
             %IF %EVAL(&FindEffectSizesTemp>9) %THEN %DO;
                 %LET FindEffectSizesNumber = &FindEffectSizesTemp;
             %END;
             %ELSE %DO; 
                 %LET FindEffectSizesNumber = 0&FindEffectSizesTemp;
             %END;
             %LET FindEffectSizesCategories = &FindEffectSizesCategories 2;
             %LET FindEffectSizesItemNames = &FindEffectSizesItemNames I&FindEffectSizesNumber;
         %END;
          allCells = J(nCells, nItems, 0);
          DO j = 1 TO nItems; 
               allCells[,j] = SHAPE(REPEAT({1, 2}, 2**(j - 1), 2**(nItems - j)),nCells) ;
          END;
          cellsProbForH1GivenClass = J(nCells,nClasses,0);
          cellsProbForH1 = J(nCells,1,0);
          DO i = 1 TO nCells;
               DO j = 1 TO nClasses; 
                    m1 = rhos[,j] || (1-rhos[,j]);
                    m2 = (allCells[i,]=1)` || (allCells[i,]=2)`; 
                    cellsProbForH1GivenClass[i,j] = EXP(SUM(LOG(1E-30+((m1#m2)[,+]))));
                     /*  The above three lines are a confusing but concise way to say this:
                         The probability under the true H1, of having response profile i,
                         for an individual in class j, is the product of the rhos associated
                         with each of the responses.  For example, if response profile i is
                         (2,1,2,2,1), and rho[k] is the probability of answering 1 to
                         dichotomous item k, then we want 
                          (1-rho[1])*(rho[2])*(1-rho[3])*(1-rho[4])*(rho[5]). */
               END;
               cellsProbForH1[i] = cellsProbForH1GivenClass[i,] * gammas; * dot product, i.e., sum of cross-products;
          END;
          mockPopulation = allCells || cellsProbForH1; 
          CREATE FindEffectSizesMockPopulation VAR {&FindEffectSizesItemNames cellsProbForH1};
              APPEND FROM mockPopulation;
          CLOSE FindEffectSizesMockPopulation; 
          CREATE FindEffectSizesCellsProbForH1 FROM cellsProbForH1[COLNAME="ExpectedH1"];
               APPEND FROM cellsProbForH1;
          CLOSE FindEffectSizesCellsProbForH1;
          CREATE FindEffectSizesCellsProbH1Class FROM cellsProbForH1GivenClass;
               APPEND FROM cellsProbForH1GivenClass;
          CLOSE FindEffectSizesCellsProbH1Class;
     QUIT;
     DATA FindEffectSizesMockPopulation;
         SET FindEffectSizesMockPopulation;
         WHERE cellsProbForH1>0;
         mockCount = 1000*cellsProbForH1;
         cell = _N_;
     RUN;
     * Determine null hypothesis cell probabilities;
     PROC LCA DATA=FindEffectSizesMockPopulation 
              OUTPARAM=FindEffectSizesH0Params 
              OUTPOST=FindEffectSizesPostH0 NOPRINT;
        NCLASS &FindEffectSizesH0Classes;
        ITEMS &FindEffectSizesItemNames;
        CATEGORIES &FindEffectSizesCategories; 
        ESTIMATION EM;
        SEED 100000;  
        MAXITER 5000;          
        FREQ mockCount;
        CRITERION 0.0000000001;
        ID cell;
        NSTARTS &nstarts;
     RUN; 
     PROC IML;
          USE FindEffectSizesH0Params;
               READ ALL INTO H0Params;
          CLOSE FindEffectSizesH0Params;
          gammas = H0Params[1,3:NCOL(H0Params)]`;
          nItems = (NROW(H0Params)-1)/2; * only want the rho1s, not the gammas or rho2s;
          rhos = H0Params[2:(1+nItems),3:NCOL(H0Params)];
          nClasses = NROW(gammas);
          nItems = NROW(rhos);
          nCells = 2**nItems;
          allCells = J(nCells, nItems, 0);
          DO j = 1 TO nItems; 
               allCells[,j] = SHAPE(REPEAT({1, 2}, 2**(j - 1), 2**(nItems - j)),nCells) ;
          END;
          cellsProbForH0GivenClass = J(nCells,nClasses,0);
          cellsProbForH0 = J(nCells,1,0); 
          DO i = 1 TO nCells;
               DO j = 1 TO nClasses;
                    m1 = rhos[,j] || (1-rhos[,j]);
                    m2 = (allCells[i,]=1)` || (allCells[i,]=2)`;
                    cellsProbForH0GivenClass[i,j] = EXP(SUM(LOG(1e-30+((m1#m2)[,+]))));
                     /*  The above three lines are a confusing but concise way to say this:
                         The probability under the true H0, of having response profile i,
                         for an individual in class j, is the product of the rhos associated
                         with each of the responses.  For example, if response profile i is
                         (2,1,2,2,1), and rho[k] is the probability of answering 1 to
                         dichotomous item k, then we want 
                          (1-rho[1])*(rho[2])*(1-rho[3])*(1-rho[4])*(rho[5]). */
               END;
               cellsProbForH0[i] = cellsProbForH0GivenClass[i,] * gammas; * dot product, i.e., sum of cross-products;
          END; 
          CREATE FindEffectSizesAllCells FROM allCells;
               APPEND FROM allCells;
          CLOSE FindEffectSizesAllCells;
          CREATE FindEffectSizesCellsProbForH0 FROM cellsProbForH0[COLNAME="ExpectedH0"];
               APPEND FROM cellsProbForH0;
          CLOSE FindEffectSizesCellsProbForH0;
          CREATE FindEffectSizesCellsProbH0Class FROM cellsProbForH0GivenClass;
               APPEND FROM cellsProbForH0GivenClass;
          CLOSE FindEffectSizesCellsProbH0Class;
     QUIT;
     * Determine effect size;
     DATA FindEffectSizesExpectedProbs;
          MERGE FindEffectSizesCellsProbForH0 FindEffectSizesCellsProbForH1; 
          CellDifference = ((ExpectedH1 - ExpectedH0)**2)/(ExpectedH0+1E-30);
          CellKLDivergence = (LOG(1e-30+ExpectedH1)-LOG(1e-30+ExpectedH0))*ExpectedH1;
     RUN;
     PROC IML;
          USE FindEffectSizesExpectedProbs;
               READ ALL VAR {CellDifference CellKLDivergence };
          CLOSE FindEffectSizesExpectedProbs; 
          EffectSize = SQRT(SUM(CellDifference));
          EffectSizeKL = SUM(CellKLDivergence); 
          J = &FindEffectSizesNItems;
          K = &FindEffectSizesH1Classes;
          EffectSizes = J || K || EffectSize || EffectSizeKL  ;
          cols = 'J' || 'K' || 'W' || 'KL' ;
          CREATE FindEffectSizesEffectSize FROM EffectSizes [COLNAME = cols];
               APPEND  FROM EffectSizes;
          CLOSE FindEffectSizesEffectSize;
     QUIT; 
          * Determine population entropy measure for H0;
          PROC IML;
               USE FindEffectSizesH0Params;
                    READ ALL INTO H0Params;
               CLOSE FindEffectSizesH0Params; 
               gammas = H0Params[1,3:2+&FindEffectSizesH0classes];  
               USE FindEffectSizesCellsProbForH0;
                    READ ALL INTO ExpectedH0;
               CLOSE FindEffectSizesCellsProbForH0; 
               USE FindEffectSizesCellsProbH0Class;
                    READ ALL INTO ExpectedH0GivenClass;
               CLOSE FindEffectSizesCellsProbH0Class;
               PostProbs = J(NROW(ExpectedH0GivenClass),NCOL(ExpectedH0GivenClass),0);
               DO j = 1 TO NCOL(ExpectedH0GivenClass);
                    PostProbs[,j] = gammas[j]*ExpectedH0GivenClass[,j]/ExpectedH0;
               END;
               EntropyContribByCell = (1/LOG(&FindEffectSizesH0classes+1E-30))*( PostProbs # LOG(PostProbs+1E-30))[,+]; 
               RescaledEntropy = 1 + SUM( ExpectedH0 # EntropyContribByCell );  
               CREATE FindEffectSizesH0RescaledEntropy FROM RescaledEntropy [COLNAME = 'H0RescaledEntropy'];
                    APPEND  FROM RescaledEntropy;
               CLOSE FindEffectSizesH0RescaledEntropy;
          QUIT;
          * Determine population entropy measure for H1;
          PROC IML;
               USE &gammas;
                    READ ALL INTO gammas;
               CLOSE &gammas;  
               USE FindEffectSizesCellsProbForH1;
                    READ ALL INTO ExpectedH1;
               CLOSE FindEffectSizesCellsProbForH1; 
               USE FindEffectSizesCellsProbH1Class;
                    READ ALL INTO ExpectedH1GivenClass;
               CLOSE FindEffectSizesCellsProbH1Class;
               PostProbs = J(NROW(ExpectedH1GivenClass),NCOL(ExpectedH1GivenClass),0);
               DO j = 1 TO NCOL(ExpectedH1GivenClass);
                    PostProbs[,j] = gammas[j]*ExpectedH1GivenClass[,j]/ExpectedH1;
               END;
               EntropyContribByCell = (1/LOG(&FindEffectSizesH1classes+1E-30))*( PostProbs # LOG(PostProbs+1E-30))[,+]; 
               RescaledEntropy = 1 + SUM( ExpectedH1 # EntropyContribByCell );
               CREATE FindEffectSizesPostH1 FROM PostProbs;
                    APPEND FROM PostProbs;
               CLOSE FindEffectSizesPostH1;
               CREATE FindEffectSizesH1RescaledEntropy FROM RescaledEntropy [COLNAME = 'H1RescaledEntropy'];
                    APPEND  FROM RescaledEntropy;
               CLOSE FindEffectSizesH1RescaledEntropy;
          QUIT;
     DATA FindEffectSizesAnswers;
          MERGE FindEffectSizesEffectSize FindEffectSizesH0RescaledEntropy FindEffectSizesH1RescaledEntropy;
     RUN;
     PROC DATASETS;
        DELETE FindEffectSizesAllCells 
               FindEffectSizesCellsProbForH0 
               FindEffectSizesCellsProbH0Class
               FindEffectSizesCellsProbForH1
               FindEffectSizesCellsProbH1Class 
               FindEffectSizesEffectSize 
               FindEffectSizesExpectedProbs
               FindEffectSizesH0Params
               FindEffectSizesH0RescaledEntropy
               FindEffectSizesH1RescaledEntropy
               FindEffectSizesMockPopulation 
               FindEffectSizesPostH0
               FindEffectSizesPostH1 ;
     RUN;
     PROC PRINT DATA=FindEffectSizesAnswers; RUN;
%MEND;
