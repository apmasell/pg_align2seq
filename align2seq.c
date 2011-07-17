/*
 * align2seq
 *
 * Functions to perform sequence alignment directly in PostgreSQL
 * Andre Masella <andre@masella.no-ip.org>
 */

#include <postgres.h>
#include <fmgr.h>
#include <c.h>
#include <funcapi.h>

#include <string.h>
#include <math.h>

#ifdef PG_MODULE_MAGIC
PG_MODULE_MAGIC;
#endif

int4 n2id (char c);
int4 p2id (char c);
int4 align_and_score (text* seq1, text* seq2, int4 (*scoring)(int4,int4), int4 (*sanitise)(char));
Datum align_n (PG_FUNCTION_ARGS);
Datum align_p (PG_FUNCTION_ARGS);
Datum n2p (PG_FUNCTION_ARGS);
int4 dnascoring(int4 a, int4 b);
int4 blosum62scoring(int4 a, int4 b);

/* Convert a nucleotide character to an arbitrary id number. */
int4 n2id (char c) {
	switch (c) {
		case 'U':
		case 'u':
		case 'T':
		case 't':
			return 1;
		case 'C':
		case 'c':
			return 2;
		case 'A':
		case 'a':
			return 3;
		case 'G':
		case 'g':
			return 4;
		default:
			return 0;
	}
}

/* Convert a single-letter amino acid character to an arbitrary id number. */
int4 p2id (char c) {
	switch (c) {
		case 'C':
		case 'c':
			return 1;
		case 'S':
		case 's':
			return 2;
		case 'T':
		case 't':
			return 3;
		case 'P':
		case 'p':
			return 4;
		case 'A':
		case 'a':
			return 5;
		case 'G':
		case 'g':
			return 6;
		case 'N':
		case 'n':
			return 7;
		case 'D':
		case 'd':
			return 8;
		case 'E':
		case 'e':
			return 9;
		case 'Q':
		case 'q':
			return 10;
		case 'H':
		case 'h':
			return 11;
		case 'R':
		case 'r':
			return 12;
		case 'K':
		case 'k':
			return 13;
		case 'M':
		case 'm':
			return 14;
		case 'I':
		case 'i':
			return 15;
		case 'L':
		case 'l':
			return 16;
		case 'V':
		case 'v':
			return 17;
		case 'F':
		case 'f':
			return 18;
		case 'Y':
		case 'y':
			return 19;
		case 'W':
		case 'w':
			return 20;
		default:
			return 0;
	}
}

/* Align two sequences using the Smith-Waterman algorithm and return the best score. */
#define MATRIX(i, j) (matrix[(i)*(VARSIZE(seq2) - VARHDRSZ) + (j)])
#define GAP_EXTEND_PENALTY 1
#define GAP_OPEN_PENALTY 3
int4 align_and_score (text* seq1, text* seq2, int4 (*scoring)(int4, int4), int4 (*sanitise)(char)) {
	int4* matrix;
	int4 i;
	int4 j;
	/* Best score seen so far. */
	int4 highest = 0;

	matrix = palloc((VARSIZE(seq1) - VARHDRSZ) * (VARSIZE(seq2) - VARHDRSZ) * sizeof(int4));

	/* Initialise first row and column of scoreing matrix. */
	for (i = 0; i < (VARSIZE(seq1) - VARHDRSZ); i++) {
		MATRIX(i, 0) = scoring(sanitise(VARDATA(seq1)[i]), sanitise(VARDATA(seq2)[0]));
	}

	for (j = 1; j < (VARSIZE(seq2) - VARHDRSZ); j++) {
		MATRIX(0, j) = scoring(sanitise(VARDATA(seq1)[0]), sanitise(VARDATA(seq2)[j]));
	}

	/* Fill in the rest of the matrix. */
	for (i = 1; i < (VARSIZE(seq1) - VARHDRSZ); i++) {
		for (j = 1; j < (VARSIZE(seq2) - VARHDRSZ); j++) {
			int4 gap;
			int4 bestchoice;

			/* Calculate the score of aligning the two characters. */
			bestchoice = abs(MATRIX(i - 1, j - 1)) + scoring(sanitise(VARDATA(seq1)[i]), sanitise(VARDATA(seq2)[j]));
			if (bestchoice < 0) bestchoice = 0;

			/* Try inserting or extending a gap in one sequence. */
			if (MATRIX(i - 1, j) < 0) {
				gap = MATRIX(i - 1, j) + GAP_EXTEND_PENALTY;
			} else {
				gap = -MATRIX(i - 1, j) + GAP_OPEN_PENALTY;
			}
			if (gap < 0 && abs(gap) > abs(bestchoice)) bestchoice = gap;

			/* Try inserting or extending a gap in the other sequence. */
			if (MATRIX(i, j - 1) < 0) {
				gap = MATRIX(i, j - 1) + GAP_EXTEND_PENALTY;
			} else {
				gap = -MATRIX(i, j - 1) + GAP_OPEN_PENALTY;
			}
			if (gap < 0 && abs(gap) > abs(bestchoice)) bestchoice = gap;

			/* Save this score. */
			MATRIX(i, j) = bestchoice;
			if (abs(bestchoice) > highest) highest = abs(bestchoice);
		}
	}

	pfree(matrix);
	return highest;
}

/* How to score nucleotide alignments. */
int4 dnascoringmatrix[5][5] = {
	/*X*/ { -1, -1, -1, -1, -1 },
	/*U*/ { -1,  3,  1,  1,  1 },
	/*C*/ { -1,  1,  3,  1,  1 },
	/*A*/ { -1,  1,  1,  3,  1 },
	/*G*/ { -1,  1,  1,  1,  3 }
};

int4 dnascoring(int4 a, int4 b) {
	return dnascoringmatrix[a][b];
}

/* How to score protein alignments. */
int4 blosum62scoringmatrix[21][21] = {
/*X*/ { -20, -20, -20, -20, -20, -20, -20, -20, -20, -20, -20, -20, -20, -20, -20, -20, -20, -20, -20, -20, -20},
/*C*/ { -20,   9,  -1,  -1,  -3,   0,  -3,  -3,  -3,  -4,  -3,  -3,  -3,  -3,  -1,  -1,  -1,  -1,  -2,  -2,  -2},
/*S*/ { -20,  -1,   4,   1,  -1,   1,   0,   1,   0,   0,   0,  -1,  -1,   0,  -1,  -2,  -2,  -2,  -2,  -2,  -3},
/*T*/ { -20,  -1,   1,   4,   1,  -1,   1,   0,   1,   0,   0,   0,  -1,   0,  -1,  -2,  -2,  -2,  -2,  -2,  -3},
/*P*/ { -20,  -3,  -1,   1,   7,  -1,  -2,  -1,  -1,  -1,  -1,  -2,  -2,  -1,  -2,  -3,  -3,  -2,  -4,  -3,  -4},
/*A*/ { -20,   0,   1,  -1,  -1,   4,   0,  -1,  -2,  -1,  -1,  -2,  -1,  -1,  -1,  -1,  -1,  -2,  -2,  -2,  -3},
/*G*/ { -20,  -3,   0,   1,  -2,   0,   6,  -2,  -1,  -2,  -2,  -2,  -2,  -2,  -3,  -4,  -4,   0,  -3,  -3,  -2},
/*N*/ { -20,  -3,   1,   0,  -2,  -2,   0,   6,   1,   0,   0,  -1,   0,   0,  -2,  -3,  -3,  -3,  -3,  -2,  -4},
/*D*/ { -20,  -3,   0,   1,  -1,  -2,  -1,   1,   6,   2,   0,  -1,  -2,  -1,  -3,  -3,  -4,  -3,  -3,  -3,  -4},
/*E*/ { -20,  -4,   0,   0,  -1,  -1,  -2,   0,   2,   5,   2,   0,   0,   1,  -2,  -3,  -3,  -3,  -3,  -2,  -3},
/*Q*/ { -20,  -3,   0,   0,  -1,  -1,  -2,   0,   0,   2,   5,   0,   1,   1,   0,  -3,  -2,  -2,  -3,  -1,  -2},
/*H*/ { -20,  -3,  -1,   0,  -2,  -2,  -2,   1,   1,   0,   0,   8,   0,  -1,  -2,  -3,  -3,  -2,  -1,   2,  -2},
/*R*/ { -20,  -3,  -1,  -1,  -2,  -1,  -2,   0,  -2,   0,   1,   0,   5,   2,  -1,  -3,  -2,  -3,  -3,  -2,  -3},
/*K*/ { -20,  -3,   0,   0,  -1,  -1,  -2,   0,  -1,   1,   1,  -1,   2,   5,  -1,  -3,  -2,  -3,  -3,  -2,  -3},
/*M*/ { -20,  -1,  -1,  -1,  -2,  -1,  -3,  -2,  -3,  -2,   0,  -2,  -1,  -1,   5,   1,   2,  -2,   0,  -1,  -1},
/*I*/ { -20,  -1,  -2,  -2,  -3,  -1,  -4,  -3,  -3,  -3,  -3,  -3,  -3,  -3,   1,   4,   2,   1,   0,  -1,  -3},
/*L*/ { -20,  -1,  -2,  -2,  -3,  -1,  -4,  -3,  -4,  -3,  -2,  -3,  -2,  -2,   2,   2,   4,   3,   0,  -1,  -2},
/*V*/ { -20,  -1,  -2,  -2,  -2,   0,  -3,  -3,  -3,  -2,  -2,  -3,  -3,  -2,   1,   3,   1,   4,  -1,  -1,  -3},
/*F*/ { -20,  -2,  -2,  -2,  -4,  -2,  -3,  -3,  -3,  -3,  -3,  -1,  -3,  -3,   0,   0,   0,  -1,   6,   3,   1},
/*Y*/ { -20,  -2,  -2,  -2,  -3,  -2,  -3,  -2,  -3,  -2,  -1,   2,  -2,  -2,  -1,  -1,  -1,  -1,   3,   7,   2},
/*W*/ { -20,  -2,  -3,  -3,  -4,  -3,  -2,  -4,  -4,  -3,  -2,  -2,  -3,  -3,  -1,  -3,  -2,  -3,   1,   2,  11}};

int4 blosum62scoring(int4 a, int4 b) {
	return blosum62scoringmatrix[a][b];
}

/* Codon translation table. X = unknown, _ = stop */
char translation[5][5][5] = {
/*X*/
	{
		/*X*/ { 'X', 'X', 'X', 'X', 'X' }, 
		/*U*/ { 'X', 'X', 'X', 'X', 'X' }, 
		/*C*/ { 'X', 'X', 'X', 'X', 'X' }, 
		/*A*/ { 'X', 'X', 'X', 'X', 'X' }, 
		/*G*/ { 'X', 'X', 'X', 'X', 'X' }, 
	},

/*U*/
	{
		/*X*/ { 'X', 'X', 'X', 'X', 'X' }, 
		/*U*/ { 'X', 'F', 'F', 'L', 'L' }, 
		/*C*/ { 'S', 'S', 'S', 'S', 'S' }, 
		/*A*/ { 'X', 'Y', 'Y', '_', '_' }, 
		/*G*/ { 'X', 'C', 'C', '_', 'W' }, 
	},

/*C*/
	{
		/*X*/ { 'X', 'X', 'X', 'X', 'X' }, 
		/*U*/ { 'L', 'L', 'L', 'L', 'L' }, 
		/*C*/ { 'P', 'P', 'P', 'P', 'P' }, 
		/*A*/ { 'X', 'H', 'H', 'Q', 'Q' }, 
		/*G*/ { 'R', 'R', 'R', 'R', 'R' }, 
	},

/*A*/
	{
		/*X*/ { 'X', 'X', 'X', 'X', 'X' }, 
		/*U*/ { 'X', 'I', 'I', 'I', 'M' }, 
		/*C*/ { 'T', 'T', 'T', 'T', 'T' }, 
		/*A*/ { 'X', 'N', 'N', 'K', 'K' }, 
		/*G*/ { 'X', 'S', 'S', 'R', 'R' }, 
	},

/*G*/
	{
		/*X*/ { 'X', 'X', 'X', 'X', 'X' }, 
		/*U*/ { 'V', 'V', 'V', 'V', 'V' }, 
		/*C*/ { 'A', 'A', 'A', 'A', 'A' }, 
		/*A*/ { 'X', 'D', 'D', 'E', 'E' }, 
		/*G*/ { 'G', 'G', 'G', 'G', 'G' }, 
	}
};

/* Exported functions. */
PG_FUNCTION_INFO_V1(align_n);
Datum align_n (PG_FUNCTION_ARGS) {
	text* seq1 = PG_GETARG_TEXT_P(0);
	text* seq2 = PG_GETARG_TEXT_P(1);
	
	PG_RETURN_INT32(align_and_score(seq1, seq2, &dnascoring, &n2id));
}

PG_FUNCTION_INFO_V1(align_p);
Datum align_p (PG_FUNCTION_ARGS) {
	text* seq1 = PG_GETARG_TEXT_P(0);
	text* seq2 = PG_GETARG_TEXT_P(1);
	
	PG_RETURN_INT32(align_and_score(seq1, seq2, &blosum62scoring, &p2id));
}

PG_FUNCTION_INFO_V1(n2p);
Datum n2p (PG_FUNCTION_ARGS) {
	text* nucleotides = PG_GETARG_TEXT_P(0);
	int4 skip = PG_GETARG_INT32(1);
	text* aminoacids;
	int4 i;
	int4 len = (VARSIZE(nucleotides) - skip - VARHDRSZ) / 3;

	if (len <= 0) {
		aminoacids = palloc(1 + VARHDRSZ);
	} else { 
		aminoacids = palloc(len + 1 + VARHDRSZ);

		for (i = 0; i < len; i++) {
			VARDATA(aminoacids)[i] = translation[n2id(VARDATA(nucleotides)[3*i+skip])][n2id(VARDATA(nucleotides)[3*i+1+skip])][n2id(VARDATA(nucleotides)[3*i+2+skip])];
		}
	}
	VARDATA(aminoacids)[len] = 0;
	SET_VARSIZE(aminoacids, len + VARHDRSZ);
	PG_RETURN_TEXT_P(aminoacids);
}
