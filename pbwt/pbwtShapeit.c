#include "pbwt.h"

static void countHelp(int *x, int start, int end, int *cc, int **u, int *f, int *g) {
  for ( int k = start; k < end; ++k ) {
    if ((*f) < (*g)) {
      (*f) = x[k] ? cc[k] + ((*f) - u[k][(*f)]) : u[k][(*f)];
      (*g) = x[k] ? cc[k] + ((*g) - u[k][(*g)]) : u[k][(*g)]; 
    } else {
      (*f) = 0; (*g) = 0;
      return ;
    }
  }
}

void pbwtMatchCount (PBWT *p, int L) /* reporting the match number for each segment */ 
{
  if (!p || !p->yz) die ("option -longWithin called without a PBWT") ;
  if (L < 0) die ("L %d for longWithin must be >= 0", L) ;

  uchar **reference = pbwtHaplotypes (p) ; /* haplotypes for reference  (M * N)  */
  int *x;                 /* use for current query */
  PbwtCursor *up = pbwtCursorCreate (p, TRUE, TRUE) ;
  int **a, **d, **u ;   /* stored indexes */
  int **f1, **g1 ;     /* start of match, and pbwt interval as in algorithm 5 */
  int **f2, **g2 ;    /* next versions of the above, e' etc in algorithm 5 */
  int i, j, k, N = p->N, M = p->M ;
  int s, seg_num = 1; /* for the segment number and current segment */
  int *pos; /* for the 1 position for ref */
  int num_1 = 0; /* for the num of 1 */

  /* build indexes */
  a = myalloc (N+1,int*) ; for (i = 0 ; i < N+1 ; ++i) a[i] = myalloc (p->M, int) ;
  d = myalloc (N+1,int*) ; for (i = 0 ; i < N+1 ; ++i) d[i] = myalloc (p->M+1, int) ;
  u = myalloc (N,int*) ; for (i = 0 ; i < N ; ++i) u[i] = myalloc (p->M+1, int) ;
  x = myalloc (N, int*) ; 
  pos = myalloc (M, int*) ;
  f1 = myalloc (8, int*);
  g1 = myalloc (8, int*);
  f2 = myalloc (64, int*);
  g2 = myalloc (64, int*);
  int *cc = myalloc (p->N, int) ;

  for (k = 0 ; k < N ; ++k)
    { memcpy (a[k], up->a, M*sizeof(int)) ;
      memcpy (d[k], up->d, (M+1)*sizeof(int)) ;
      cc[k] = up->c ;
      pbwtCursorCalculateU (up) ;
      memcpy (u[k], up->u, (M+1)*sizeof(int)) ;
      pbwtCursorForwardsReadAD (up, k) ;
    }
  memcpy (a[k], up->a, M*sizeof(int)) ;
  memcpy (d[k], up->d, (M+1)*sizeof(int)) ;
  pbwtCursorDestroy (up) ;

  fprintf (stderr, "Made indices: ") ; timeUpdate () ;
 
  
  struct timeval tstart, tend;
  gettimeofday( &tstart, NULL );

  int t;  //multi_time
  int TIMES = 1;
  for (t = 0; t < TIMES; ++t) {
    // L = rand()%(M/2); 
    
    seg_num = 1;
    num_1 = 0;
    for ( i = 0; i < 8; ++i)  { free(f1[i]); free(g1[i]); }
    for ( i = 0; i < 64; ++i) { free(f2[i]); free(g2[i]); }

    for ( i = 0; i < N; ++i)
      x[i] = reference[L * 2][i] + reference[L * 2 + 1][i];
    
    for ( i = 0, j = 0; i < N; ++i) {
      if (x[i] == 1) {
        ++j;
        pos[num_1++] = i;
        if (j == 3) {
          ++seg_num;
          j = 0;
        }
      }
      x[i] = x[i] / 2;  //change 0->0, 1->0, 2->1; 
    }
    for ( i = 0; i < 8; ++i) { 
      f1[i] = myalloc(seg_num, int*);
      g1[i] = myalloc(seg_num, int*);
    }
    for ( i = 0; i < 64; ++i) { 
      f2[i] = myalloc(seg_num, int*);
      g2[i] = myalloc(seg_num, int*);
    }
    //initial f1,g1
    for (i = 0; i < 8; ++i)
      for (j = 0; j < seg_num; ++j) {
        f1[i][j] = 0; g1[i][j] = M;
      }
    //one segment count
    // last segment may not 3 1s;
    int start = 0, end;
    for ( i = 0; i < seg_num - 1; ++i) {
      end = pos[i * 3 + 2] + 1;
      //000
      countHelp(x, start, end, cc, u, &f1[0][i], &g1[0][i]);
      //001
      x[pos[i * 3 + 2]] = 1;
      countHelp(x, start, end, cc, u, &f1[1][i], &g1[1][i]);
      //011
      x[pos[i * 3 + 1]] = 1;
      countHelp(x, start, end, cc, u, &f1[3][i], &g1[3][i]);    
      //111
      x[pos[i * 3 + 0]] = 1;
      countHelp(x, start, end, cc, u, &f1[7][i], &g1[7][i]); 
      //110
      x[pos[i * 3 + 2]] = 0;
      countHelp(x, start, end, cc, u, &f1[6][i], &g1[6][i]); 
      //100
      x[pos[i * 3 + 1]] = 0;
      countHelp(x, start, end, cc, u, &f1[4][i], &g1[4][i]); 
      //101
      x[pos[i * 3 + 2]] = 1;
      countHelp(x, start, end, cc, u, &f1[5][i], &g1[5][i]); 
      //010
      x[pos[i * 3 + 0]] = 0;
      x[pos[i * 3 + 1]] = 1;
      x[pos[i * 3 + 2]] = 0;
      countHelp(x, start, end, cc, u, &f1[2][i], &g1[2][i]); 
      start = end;
    }

   //  print_one_seg();
    printf("the first segment\n");
    for (i = 0; i < 8; ++i)
      printf("%d\n", g1[i][0] - f1[i][0]);
    printf("\n\n");


    //initial f2, g2
    for (i = 0; i < 64; ++i)
      for (j = 0; j < seg_num; ++j) {
        f2[i][j] = f1[i/8][j]; 
        g2[i][j] = g1[i/8][j];
      }
    start = pos[2] + 1;  
    for ( i = 1; i < seg_num - 1; ++i) {
      end = pos[i * 3 + 2] + 1;
      for ( j = 0; j < 8; ++j){
      x[pos[i * 3 + 0]] = 0; //initial
      x[pos[i * 3 + 1]] = 0;
      x[pos[i * 3 + 2]] = 0;
      //000
      countHelp(x, start, end, cc, u, &f2[0 + j * 8][i - 1], &g2[0 + j * 8][i - 1]);
      //001
      x[pos[i * 3 + 2]] = 1;
      countHelp(x, start, end, cc, u, &f2[1 + j * 8][i - 1], &g2[1 + j * 8][i - 1]);
      //011
      x[pos[i * 3 + 1]] = 1;
      countHelp(x, start, end, cc, u, &f2[3 + j * 8][i - 1], &g2[3 + j * 8][i - 1]);    
      //111
      x[pos[i * 3 + 0]] = 1;
      countHelp(x, start, end, cc, u, &f2[7 + j * 8][i - 1], &g2[7 + j * 8][i - 1]); 
      //110
      x[pos[i * 3 + 2]] = 0;
      countHelp(x, start, end, cc, u, &f2[6 + j * 8][i - 1], &g2[6 + j * 8][i - 1]); 
      //100
      x[pos[i * 3 + 1]] = 0;
      countHelp(x, start, end, cc, u, &f2[4 + j * 8][i - 1], &g2[4 + j * 8][i - 1]); 
      //101
      x[pos[i * 3 + 2]] = 1;
      countHelp(x, start, end, cc, u, &f2[5 + j * 8][i - 1], &g2[5 + j * 8][i - 1]); 
      //010
      x[pos[i * 3 + 0]] = 0;
      x[pos[i * 3 + 1]] = 1;
      x[pos[i * 3 + 2]] = 0;
      countHelp(x, start, end, cc, u, &f2[2 + j * 8][i - 1], &g2[2 + j * 8][i - 1]); 
      }
      start = end;
    }

   for ( i = 0; i < seg_num - 2; ++i) {
      printf("segment num\t%d\n", i + 2);
      printf("index\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",0,1,2,3,4,5,6,7);
      for (j = 0; j < 8; ++j)
        printf("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",j,
          g2[8*j+0][i] - f2[8*j+0][i],
          g2[8*j+1][i] - f2[8*j+1][i],
          g2[8*j+2][i] - f2[8*j+2][i],
          g2[8*j+3][i] - f2[8*j+3][i],
          g2[8*j+4][i] - f2[8*j+4][i],
          g2[8*j+5][i] - f2[8*j+5][i],
          g2[8*j+6][i] - f2[8*j+6][i],
          g2[8*j+7][i] - f2[8*j+7][i]);
      printf("\n\n");
    } 
  }

  gettimeofday( &tend, NULL );
  int timeuse = 1000000 * ( tend.tv_sec - tstart.tv_sec ) + tend.tv_usec -tstart.tv_usec;
  timeuse /= TIMES;
  printf("time: %d us\n", timeuse);

  /* cleanup */
  free (cc) ;
  for (j = 0 ; j < p->M ; ++j) free(reference[j]) ; free (reference) ;
  free(pos);
  for ( i = 0; i < 8; ++i)  { free(f1[i]); free(g1[i]); }
  for ( i = 0; i < 64; ++i) { free(f2[i]); free(g2[i]); }
  free(f1); free(g1); free(f2); free(g2);
  for (j = 0 ; j < N ; ++j) free(a[j]) ; free (a) ;
  for (j = 0 ; j < N ; ++j) free(d[j]) ; free (d) ;
  for (j = 0 ; j < N ; ++j) free(u[j]) ; free (u) ;
}

/******************* end of file *******************/