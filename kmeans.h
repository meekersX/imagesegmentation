#include <vector>
using namespace std;

void OPTRA(	unsigned int *IC1,
			int *IC2,
			int *NC,
			int *NCP,
			int *ITRAN,
			int *LIVE,
			vector<vector<double> > &A,
			double *D,
			vector<vector<double> > &C,
			double *AN1,
			double *AN2,
			int M,
			int N,
			int K,
			int &INDX );

void QTRAN(	unsigned int *IC1,
			int *IC2,
			int *NC,
			int *NCP,
			int *ITRAN,
			vector<vector<double> > &A,
			double *D,
			vector<vector<double> > &C,
			double *AN1,
			double *AN2,
			int M,
			int N,
			int K,
			int &INDX );

pair<vector<vector<double> >, vector<unsigned int> > kmeans(const vector<vector<double> > &fvs, vector<vector<double> > means, unsigned int k, double &d)
{
	//      SUBROUTINE KMNS(A, M, N, C, K, IC1, IC2, NC, AN1, AN2, NCP, D,
	//     *    ITRAN, LIVE, ITER, WSS, IFAULT)
	//C
	//C     ALGORITHM AS 136  APPL. STATIST. (1979) VOL.28, NO.1
	//C
	//C     Divide M points in N-dimensional space into K clusters so that
	//C     the within cluster sum of squares is minimized.
	//C
	//      INTEGER IC1(M), IC2(M), NC(K), NCP(K), ITRAN(K), LIVE(K)
	//      REAL    A(M,N), D(M), C(K,N), AN1(K), AN2(K), WSS(K), DT(2)
	//      REAL    ZERO, ONE
	//C
	//C     Define BIG to be a very large positive number
	//C
	//      DATA BIG /1.E30/, ZERO /0.0/, ONE /1.0/
	//C
	//      IFAULT = 3
	//      IF (K .LE. 1 .OR. K .GE. M) RETURN
	//C
	//C     For each point I, find its two closest centres, IC1(I) and
	//C     IC2(I).     Assign it to IC1(I).
	//C

	int M = fvs.size();
	int K = k;
	int N = 6;
	int ITER = 100;
	unsigned int *IC1 = new unsigned int[M];
	int *IC2 = new int[M];
	int *NC = new int[K];
	int *NCP = new int[K];
	int *ITRAN = new int[K];
	int *LIVE = new int[K];
	vector<vector<double> > A(fvs);
	double *D = new double[M];
	vector<vector<double> > C(means);
	double *AN1 = new double[K];
	double *AN2 = new double[K];
	double *WSS = new double[K];
	double *DT = new double[2];
	int I;
	int IL;
	int J;
	int IJ;
	int L;
	int INDX;
	int IFAULT;
	int II;
	double AA;
	double DA;
	double DB;
	double DC;
	double TEMP;

	for (I = 0; I < M; I++)											//	DO 50 I = 1, M
	{
		IC1[I] = 0;													//	  IC1(I) = 1
		IC2[I] = 1;													//	  IC2(I) = 2

		for (IL = 0; IL < 2; IL++)									//	DO 10 IL = 1, 2
		{
			DT[IL] = 0;												//	  DT(IL) = ZERO
			for (J = 0; J < N; J++)									//	  DO 10 J = 1, N
			{
				DA = A[I][J] - C[IL][J];							//	    DA = A(I,J) - C(IL,J)
				DT[IL] = DT[IL] + DA*DA;							//	    DT(IL) = DT(IL) + DA*DA
			}
		}															//   10   CONTINUE

		if (DT[0] > DT[1])											//	IF (DT(1) .GT. DT(2)) THEN
		{
			IC1[I] = 1;												//	  IC1(I) = 2
			IC2[I] = 0;												//	  IC2(I) = 1
			TEMP = DT[0];											//	  TEMP = DT(1)
			DT[0] = DT[1];											//	  DT(1) = DT(2)
			DT[1] = TEMP;											//	  DT(2) = TEMP
		}															//	END IF

		for (L = 2; L < K; L++)										//	DO 50 L = 3, K
		{
			DB = 0.0;												//	  DB = ZERO

			for (J = 0; J < N; J++)									//	  DO 30 J = 1, N
			{
				DC = A[I][J] - C[L][J];								//	    DC = A(I,J) - C(L,J)
				DB = DB + DC*DC;									//	    DB = DB + DC*DC

				if (DB >= DT[1])									//	    IF (DB .GE. DT(2)) GO TO 50
				{
					goto fifty;
				}
			}														//   30     CONTINUE

			if (DB >= DT[0])										//	  IF (DB .LT. DT(1)) GO TO 40
			{
				DT[1] = DB;											//	  DT(2) = DB
				IC2[I] = L;											//	  IC2(I) = L
			}														//	  GO TO 50
			else
			{
				DT[1] = DT[0];										//   40     DT(2) = DT(1)
				IC2[I] = IC1[I];									//	  IC2(I) = IC1(I)
				DT[0] = DB;											//	  DT(1) = DB
				IC1[I] = L;											//	  IC1(I) = L
			}
fifty:
			;
		}															//   50 CONTINUE
	}

	//C
	//C     Update cluster centres to be the average of points contained
	//C     within them.
	//C

	for (L = 0; L < K; L++)											//	DO 70 L = 1, K
	{
		NC[L] = 0;													//	  NC(L) = 0

		for (J = 0; J < N; J++)										//	  DO 60 J = 1, N
		{
			C[L][J] = 0.0;											//	60   C(L,J) = ZERO
		}
	}																//	70 CONTINUE

	for (I = 0; I < M; I++)											//      DO 90 I = 1, M
	{
		L = IC1[I];													//	L = IC1(I)
		NC[L] = NC[L] + 1;											//	NC(L) = NC(L) + 1

		for (J = 0; J < N; J++)										//	DO 80 J = 1, N
		{
			C[L][J] = C[L][J] + A[I][J];							//   80   C(L,J) = C(L,J) + A(I,J)
		}
	}																//   90 CONTINUE

	//C
	//C     Check to see if there is any empty cluster at this stage
	//C

	for (L = 0; L < K; L++)											//      DO 120 L = 1, K
	{
		if (NC[L] == 0)												//	IF (NC(L) .EQ. 0) THEN
		{
			IFAULT = 1;												//	  IFAULT = 1
			goto Error;												//	  RETURN
		}															//	END IF

		AA = NC[L];													//	AA = NC(L)

		for (J = 0; J < N; J++)										//	DO 110 J = 1, N
		{
			C[L][J] = C[L][J] / AA;									//  110   C(L,J) = C(L,J) / AA
		}

		//C
		//C     Initialize AN1, AN2, ITRAN & NCP
		//C     AN1(L) = NC(L) / (NC(L) - 1)
		//C     AN2(L) = NC(L) / (NC(L) + 1)
		//C     ITRAN(L) = 1 if cluster L is updated in the quick-transfer stage,
		//C              = 0 otherwise
		//C     In the optimal-transfer stage, NCP(L) stores the step at which
		//C     cluster L is last updated.
		//C     In the quick-transfer stage, NCP(L) stores the step at which
		//C     cluster L is last updated plus M.
		//C

		AN2[L] = AA / (AA + 1.0);									//	AN2(L) = AA / (AA + ONE)
		AN1[L] = DBL_MAX;											//	AN1(L) = BIG

		if (AA > 1.0)												//	IF (AA .GT. ONE) AN1(L) = AA / (AA - ONE)
		{
			AN1[L] = AA / (AA - 1.0);
		}

		ITRAN[L] = 1;												//	ITRAN(L) = 1
		NCP[L] = 0;													//	NCP(L) = -1
	}																//  120 CONTINUE

	INDX = 0;														//      INDX = 0

	for (IJ = 1; IJ <= ITER; IJ++)									//      DO 140 IJ = 1, ITER
	{
		//C
		//C     In this stage, there is only one pass through the data.   Each
		//C     point is re-allocated, if necessary, to the cluster that will
		//C     induce the maximum reduction in within-cluster sum of squares.
		//C

		//	CALL OPTRA(A, M, N, C, K, IC1, IC2, NC, AN1, AN2, NCP, D,
		//     *        ITRAN, LIVE, INDX)

		OPTRA(IC1, IC2, NC, NCP, ITRAN, LIVE, A, D, C, AN1, AN2, M, N, K, INDX);

		//C
		//C     Stop if no transfer took place in the last M optimal transfer
		//C     steps.
		//C

		if (INDX == M)												//	IF (INDX .EQ. M) GO TO 150
			break;

		//C
		//C     Each point is tested in turn to see if it should be re-allocated
		//C     to the cluster to which it is most likely to be transferred,
		//C     IC2(I), from its present cluster, IC1(I).   Loop through the
		//C     data until no further change is to take place.
		//C

		//	CALL QTRAN(A, M, N, C, K, IC1, IC2, NC, AN1, AN2, NCP, D,
		//     *       ITRAN, INDX)

		QTRAN(IC1, IC2, NC, NCP, ITRAN, A, D, C, AN1, AN2, M, N, K, INDX);

		//C
		//C     If there are only two clusters, there is no need to re-enter the
		//C     optimal transfer stage.
		//C

		if (K == 2)													//	IF (K .EQ. 2) GO TO 150
			break;
	
		//C
		//C     NCP has to be set to 0 before entering OPTRA.
		//C

		for (L = 0; L < K; L++)										//	DO 130 L = 1, K
		{
			NCP[L] = -1;											//  130   NCP(L) = 0
		}

	}																//  140 CONTINUE

	//C
	//C     Since the specified number of iterations has been exceeded, set
	//C     IFAULT = 2.   This may indicate unforeseen looping.
	//C

	if (IJ > ITER)
		IFAULT = 2;													//      IFAULT = 2

	//C
	//C     Compute within-cluster sum of squares for each cluster.
	//C

	for (L = 0; L < K; L++)											//  150 DO 160 L = 1, K
	{
		WSS[L] = 0.0;												//	WSS(L) = ZERO

		for (J = 0; J < N; J++)										//	DO 160 J = 1, N
		{
			C[L][J] = 0.0;											//	  C(L,J) = ZERO
		}															//  160 CONTINUE
	}

	for (I = 0; I < M; I++)											//      DO 170 I = 1, M
	{
		II = IC1[I];												//	II = IC1(I)

		for (J = 0; J < N; J++)										//	DO 170 J = 1, N
		{
			C[II][J] = C[II][J] + A[I][J];							//	  C(II,J) = C(II,J) + A(I,J)
		}
	}																//  170 CONTINUE

	for (J = 0; J < N; J++)											//      DO 190 J = 1, N
	{
		for (L = 0; L < K; L++)										//	DO 180 L = 1, K
		{
			C[L][J] = C[L][J] / NC[L];								//  180   C(L,J) = C(L,J) / FLOAT(NC(L))
		}

		for (I = 0; I < M; I++)										//	DO 190 I = 1, M
		{
			II = IC1[1];											//	  II = IC1(I)
			DA = A[I][J] - C[II][J];								//	  DA = A(I,J) - C(II,J)
			WSS[II] = WSS[II] + DA*DA;								//	  WSS(II) = WSS(II) + DA*DA
		}
	}																//  190 CONTINUE

	//C

//	copy(IC1.begin(), IC1.end(), ostream_iterator<int>(cout, "\n"));
//	copy(IC2.begin(), IC2.end(), ostream_iterator<int>(cout, "\n"));
//	copy(NC.begin(), NC.end(), ostream_iterator<int>(cout, "\n"));
//	copy(NCP.begin(), NCP.end(), ostream_iterator<int>(cout, "\n"));
//	copy(ITRAN.begin(), ITRAN.end(), ostream_iterator<int>(cout, "\n"));
//	copy(LIVE.begin(), LIVE.end(), ostream_iterator<int>(cout, "\n"));
//	copy(D.begin(), D.end(), ostream_iterator<double>(cout, "\n"));
//	copy(AN1.begin(), AN1.end(), ostream_iterator<double>(cout, "\n"));
//	copy(AN2.begin(), AN2.end(), ostream_iterator<double>(cout, "\n"));
//	copy(WSS.begin(), WSS.end(), ostream_iterator<double>(cout, "\n"));
//	copy(DT.begin(), DT.end(), ostream_iterator<double>(cout, "\n"));

Error:
	vector<unsigned int> sets(IC1, IC1 + M);

	delete []IC1;
	delete []IC2;
	delete []NC;
	delete []NCP;
	delete []ITRAN;
	delete []LIVE;
	delete []D;
	delete []AN1;
	delete []AN2;
	delete []WSS;
	delete []DT;

	return pair<vector<vector<double> >, vector<unsigned int> >(C, sets);	//      RETURN
}																	//      END

		//C
		//C

void OPTRA(	unsigned int *IC1,
			int *IC2,
			int *NC,
			int *NCP,
			int *ITRAN,
			int *LIVE,
			vector<vector<double> > &A,
			double *D,
			vector<vector<double> > &C,
			double *AN1,
			double *AN2,
			int M,
			int N,
			int K,
			int &INDX )
{
	//      SUBROUTINE OPTRA(A, M, N, C, K, IC1, IC2, NC, AN1, AN2, NCP, D,
	//     *      ITRAN, LIVE, INDX)
	//C
	//C     ALGORITHM AS 136.1  APPL. STATIST. (1979) VOL.28, NO.1
	//C
	//C     This is the optimal transfer stage.
	//C
	//C     Each point is re-allocated, if necessary, to the cluster that
	//C     will induce a maximum reduction in the within-cluster sum of
	//C     squares.
	//C
	//      INTEGER IC1(M), IC2(M), NC(K), NCP(K), ITRAN(K), LIVE(K)
	//      REAL    A(M,N), D(M), C(K,N), AN1(K), AN2(K), ZERO, ONE
	//C
	//C     Define BIG to be a very large positive number.
	//C
	//      DATA BIG /1.0E30/, ZERO /0.0/, ONE/1.0/
	//C
	//C     If cluster L is updated in the last quick-transfer stage, it
	//C     belongs to the live set throughout this stage.   Otherwise, at
	//C     each step, it is not in the live set if it has not been updated
	//C     in the last M optimal transfer steps.
	//C

	int I;
	int J;
	int L;
	int L1;
	int L2;
	int LL;
	double AL1;
	double AL2;
	double ALT;
	double ALW;
	double DA;
	double DB;
	double DC;
	double DD;
	double DE;
	double DF;
	double R2;
	double RR;

	for (L = 0; L < K; L++)											//      DO 10 L = 1, K
	{
		if (ITRAN[L] == 1)											//	IF (ITRAN(L) .EQ. 1) LIVE(L) = M + 1
			LIVE[L] = M + 1;
	}																//   10 CONTINUE

	for (I = 0; I < M; I++)											//      DO 100 I = 1, M
	{
		INDX = INDX + 1;											//	INDX = INDX + 1
		L1 = IC1[I];												//	L1 = IC1(I)
		L2 = IC2[I];												//	L2 = IC2(I)
		LL = L2;													//	LL = L2

		//C
		//C     If point I is the only member of cluster L1, no transfer.
		//C

		if (NC[L1] == 1)											//	IF (NC(L1) .EQ. 1) GO TO 90
			goto ninety;

		//C
		//C     If L1 has not yet been updated in this stage, no need to
		//C     re-compute D(I).
		//C

		if (NCP[L1] != -1)											//	IF (NCP(L1) .EQ. 0) GO TO 30
		{
			DE = 0;													//	DE = ZERO

			for (J = 0; J < N; J++)									//	DO 20 J = 1, N
			{
				DF = A[I][J] - C[L1][J];							//	  DF = A(I,J) - C(L1,J)
				DE = DE + DF*DF;									//	  DE = DE + DF*DF
			}														//   20   CONTINUE

			D[I] = DE * AN1[L1];									//	D(I) = DE * AN1(L1)
		}

	//C
	//C     Find the cluster with minimum R2.
	//C

		DA = 0.0;													//   30   DA = ZERO

		for (J = 0; J < N; J++)										//	DO 40 J = 1, N
		{
			DB = A[I][J] - C[L2][J];								//	  DB = A(I,J) - C(L2,J)
			DA = DA + DB*DB;										//	  DA = DA + DB*DB
		}															//   40   CONTINUE

		R2 = DA * AN2[L2];											//	R2 = DA * AN2(L2)

		for (L = 0; L < K; L++)										//	DO 60 L = 1, K
		{

			//C
			//C     If I >= LIVE(L1), then L1 is not in the live set.   If this is
			//C     true, we only need to consider clusters that are in the live set
			//C     for possible transfer of point I.   Otherwise, we need to consider
			//C     all possible clusters.
			//C

			if ((I >= LIVE[L1] && I >= LIVE[L]) || L == L1 || L == LL)
																	//	  IF (I .GE. LIVE(L1) .AND. I .GE. LIVE(L) .OR. L .EQ. L1 .OR.
																	//     *        L .EQ. LL) GO TO 60
				goto sixty;

			RR = R2 / AN2[L];										//	  RR = R2 / AN2(L)
			DC = 0.0;												//	  DC = ZERO

			for (J = 0; J < N; J++)									//	  DO 50 J = 1, N
			{
                DD = A[I][J] - C[L][J];								//	    DD = A(I,J) - C(L,J)
				DC = DC + DD*DD;									//	    DC = DC + DD*DD

				if (DC >= RR)										//	    IF (DC .GE. RR) GO TO 60
					goto sixty;
			}														//   50     CONTINUE

			R2 = DC * AN2[L];										//	  R2 = DC * AN2(L)
			L2 = L;													//	  L2 = L

sixty:																//   60     CONTINUE
			;
		}

		if (R2 >= D[I])												//	  IF (R2 .LT. D(I)) GO TO 70
		{
			//C
			//C     If no transfer is necessary, L2 is the new IC2(I).
			//C

			IC2[I] = L2;											//	  IC2(I) = L2
			goto ninety;											//	  GO TO 90
		}

		//C
		//C     Update cluster centres, LIVE, NCP, AN1 & AN2 for clusters L1 and
		//C     L2, and update IC1(I) & IC2(I).
		//C

		INDX = 0;													//   70     INDX = 0
		LIVE[L1] = M + I;											//	  LIVE(L1) = M + I
		LIVE[L2] = M + I;											//	  LIVE(L2) = M + I
		NCP[L1] = I;												//	  NCP(L1) = I
		NCP[L2] = I;												//	  NCP(L2) = I
		AL1 = NC[L1];												//	  AL1 = NC(L1)
		ALW = AL1 - 1.0;											//	  ALW = AL1 - ONE
		AL2 = NC[L2];												//	  AL2 = NC(L2)
		ALT = AL2 + 1.0;											//	  ALT = AL2 + ONE

		for (J = 0; J < N; J++)										//	  DO 80 J = 1, N
		{
			C[L1][J] = (C[L1][J] * AL1 - A[I][J]) / ALW;			//	    C(L1,J) = (C(L1,J) * AL1 - A(I,J)) / ALW
			C[L2][J] = (C[L2][J] * AL2 + A[I][J]) / ALT;			//	    C(L2,J) = (C(L2,J) * AL2 + A(I,J)) / ALT
		}															//   80     CONTINUE

		NC[L1] = NC[L1] - 1;										//	  NC(L1) = NC(L1) - 1
		NC[L2] = NC[L2] + 1;										//	  NC(L2) = NC(L2) + 1
		AN2[L1] = ALW / AL1;										//	  AN2(L1) = ALW / AL1
		AN1[L1] = DBL_MAX;											//	  AN1(L1) = BIG

		if (ALW > 1.0)												//	  IF (ALW .GT. ONE) AN1(L1) = ALW / (ALW - ONE)
			AN1[L1] = ALW / (ALW - 1.0);

		AN1[L2] = ALT / AL2;										//	  AN1(L2) = ALT / AL2
		AN2[L2] = ALT / (ALT + 1.0);								//	  AN2(L2) = ALT / (ALT + ONE)
		IC1[I] = L2;												//	  IC1(I) = L2
		IC2[I] = L1;												//	  IC2(I) = L1

ninety:																//   90   CONTINUE
		if (INDX == M)												//	IF (INDX .EQ. M) RETURN
			return;
	}																//  100 CONTINUE

	for (L = 0; L < K; L++)											//      DO 110 L = 1, K
	{
		//C
		//C     ITRAN(L) = 0 before entering QTRAN.   Also, LIVE(L) has to be
		//C     decreased by M before re-entering OPTRA.
		//C

		ITRAN[L] = 0;												//	ITRAN(L) = 0
		LIVE[L] = LIVE[L] - M;										//	LIVE(L) = LIVE(L) - M
	}																//  110 CONTINUE

	//C

	return;															//      RETURN
}																	//      END


void QTRAN(	unsigned int *IC1,
			int *IC2,
			int *NC,
			int *NCP,
			int *ITRAN,
			vector<vector<double> > &A,
			double *D,
			vector<vector<double> > &C,
			double *AN1,
			double *AN2,
			int M,
			int N,
			int K,
			int &INDX )
{
	//C
	//C
	//      SUBROUTINE QTRAN(A, M, N, C, K, IC1, IC2, NC, AN1, AN2, NCP, D,
	//     *    ITRAN, INDX)
	//C
	//C     ALGORITHM AS 136.2  APPL. STATIST. (1979) VOL.28, NO.1
	//C
	//C     This is the quick transfer stage.
	//C     IC1(I) is the cluster which point I belongs to.
	//C     IC2(I) is the cluster which point I is most likely to be
	//C         transferred to.
	//C     For each point I, IC1(I) & IC2(I) are switched, if necessary, to
	//C     reduce within-cluster sum of squares.  The cluster centres are
	//C     updated after each step.
	//C
	//      INTEGER IC1(M), IC2(M), NC(K), NCP(K), ITRAN(K)
	//      REAL    A(M,N), D(M), C(K,N), AN1(K), AN2(K), ZERO, ONE
	//C
	//C     Define BIG to be a very large positive number
	//C
	//      DATA BIG /1.0E30/, ZERO /0.0/, ONE /1.0/
	//C
	//C     In the optimal transfer stage, NCP(L) indicates the step at which
	//C     cluster L is last updated.   In the quick transfer stage, NCP(L)
	//C     is equal to the step at which cluster L is last updated plus M.
	//C

	int I;
	int J;
	int ICOUN;
	int ISTEP;
	int L1;
	int L2;
	double AL1;
	double AL2;
	double ALT;
	double ALW;
	double DA;
	double DB;
	double DD;
	double DE;
	double R2;

	ICOUN = 0;														//      ICOUN = 0
	ISTEP = 0;														//      ISTEP = 0

ten:
	for (I = 0; I < M; I++)											//   10 DO 70 I = 1, M
	{
		ICOUN = ICOUN + 1;											//	ICOUN = ICOUN + 1
		ISTEP = ISTEP + 1;											//	ISTEP = ISTEP + 1
		L1 = IC1[I];												//	L1 = IC1(I)
		L2 = IC2[I];												//	L2 = IC2(I)

		//C
		//C     If point I is the only member of cluster L1, no transfer.
		//C

		if (NC[L1] != 1)											//	IF (NC(L1) .EQ. 1) GO TO 60
		{

			//C
			//C     If ISTEP > NCP(L1), no need to re-compute distance from point I to
			//C     cluster L1.   Note that if cluster L1 is last updated exactly M
			//C     steps ago, we still need to compute the distance from point I to
			//C     cluster L1.
			//C

			if (ISTEP <= NCP[L1])									//	IF (ISTEP .GT. NCP(L1)) GO TO 30
			{
				DA = 0.0;											//	DA = ZERO

				for (J = 0; J < N; J++)								//	DO 20 J = 1, N
				{
					DB = A[I][J] - C[L1][J];						//	  DB = A(I,J) - C(L1,J)
					DA = DA + DB*DB;								//	  DA = DA + DB*DB
				}													//   20   CONTINUE

				D[I] = DA * AN1[L1];								//	D(I) = DA * AN1(L1)
			}

			//C
			//C     If ISTEP >= both NCP(L1) & NCP(L2) there will be no transfer of
			//C     point I at this step.
			//C

			if (ISTEP < NCP[L1] || ISTEP < NCP[L2])					//   30   IF (ISTEP .GE. NCP(L1) .AND. ISTEP .GE. NCP(L2)) GO TO 60
			{

				R2 = D[I] / AN2[L2];								//	R2 = D(I) / AN2(L2)
				DD = 0.0;											//	DD = ZERO

				for (J = 0; J < N; J++)								//	DO 40 J = 1, N
				{
					DE = A[I][J] - C[L2][J];						//	  DE = A(I,J) - C(L2,J)
					DD = DD + DE*DE;								//	  DD = DD + DE*DE

					if (DD >= R2)									//	  IF (DD .GE. R2) GO TO 60
						goto sixty2;
				}													//   40   CONTINUE

				//C
				//C     Update cluster centres, NCP, NC, ITRAN, AN1 & AN2 for clusters
				//C     L1 & L2.   Also update IC1(I) & IC2(I).   Note that if any
				//C     updating occurs in this stage, INDX is set back to 0.
				//C

				ICOUN = 0;											//	ICOUN = 0
				INDX = 0;											//	INDX = 0
				ITRAN[L1] = 1;										//	ITRAN(L1) = 1
				ITRAN[L2] = 1;										//	ITRAN(L2) = 1
				NCP[L1] = ISTEP + M;								//	NCP(L1) = ISTEP + M
				NCP[L2] = ISTEP + M;								//	NCP(L2) = ISTEP + M
				AL1 = NC[L1];										//	AL1 = NC(L1)
				ALW = AL1 - 1.0;									//	ALW = AL1 - ONE
				AL2 = NC[L2];										//	AL2 = NC(L2)
				ALT = AL2 + 1.0;									//	ALT = AL2 + ONE

				for (J = 0; J < N; J++)								//	DO 50 J = 1, N
				{
					C[L1][J] = (C[L1][J] * AL1 - A[I][J]) / ALW;	//	  C(L1,J) = (C(L1,J) * AL1 - A(I,J)) / ALW
					C[L2][J] = (C[L2][J] * AL2 + A[I][J]) / ALT;	//	  C(L2,J) = (C(L2,J) * AL2 + A(I,J)) / ALT
				}													//   50   CONTINUE

				NC[L1] = NC[L1] - 1;								//	NC(L1) = NC(L1) - 1
				NC[L2] = NC[L2] + 1;								//	NC(L2) = NC(L2) + 1
				AN2[L1] = ALW / AL1;								//	AN2(L1) = ALW / AL1
				AN1[L1] = DBL_MAX;									//	AN1(L1) = BIG

				if (ALW > 1.0)										//	IF (ALW .GT. ONE) AN1(L1) = ALW / (ALW - ONE)
					AN1[L1] = ALW / (ALW - 1.0);

				AN1[L2] = ALT / AL2;								//	AN1(L2) = ALT / AL2
				AN2[L2] = ALT / (ALT + 1.0);						//	AN2(L2) = ALT / (ALT + ONE)
				IC1[I] = L2;										//	IC1(I) = L2
				IC2[I] = L1;										//	IC2(I) = L1
			}
		}

		//C
		//C     If no re-allocation took place in the last M steps, return.
		//C
sixty2:
		if (ICOUN == M)												//   60   IF (ICOUN .EQ. M) RETURN
			return;
	}																//   70 CONTINUE

	goto ten;														//      GO TO 10
}																	//      END