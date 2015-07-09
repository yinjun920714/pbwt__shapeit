#include "pbwt.h"


static void countHelp(uchar *x, int start, int end, int *cc, int **u, int *f, int *g) {
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

static void setSeq(uchar *dir, int *pos, int seq, int index) {
  if( index > 7 || index < 0){
    printf("error\n");
    return;
  }

  dir[pos[seq * 3]] = index/4 > 0 ? 1 : 0; 
  index %= 4;
  dir[pos[seq * 3 + 1]] = index/2 > 0 ? 1 : 0;
  index %= 2;
  dir[pos[seq * 3 + 2]] = index > 0 ? 1 : 0;

  return;
}

static int compare(uchar **reference, int L, uchar *shape1, uchar *shape2, int N) {
    int total1, total2;
    int *fragment1, *fragment2;
    fragment1 = myalloc(N, int);
    fragment2 = myalloc(N, int);
    total1 = 0; total2 = 0;
    int count1 = 1;
    int count2 = 1;
    for (int i = 0; i < N; ++i) {
      if (reference[L * 2][i] == shape1[i])
        count1++;
      else {
        fragment1[total1++] = count1;
        count1 = 1;
      }

      if (reference[L * 2 + 1][i] == shape1[i])
        count2++;
      else {
        fragment2[total2++] = count2;
        count2 = 1;
      }
    }
    /*
    for (int i = 0; i < N; ++i){
      printf("%u", reference[L * 2][i]);
      printf("%u", reference[L * 2 + 1][i]);
      printf("%u", shape1[i]);
      printf("%u", shape2[i]);
      printf("\n");
    } 
    printf("\n");
    */

    // printf("%*.*s\n", N, N, shape1);
    free (fragment1) ;
    free (fragment2) ;
    return (total1 > total2 ? total2 : total1);
/*
    std::map<int, int> map1; //shape1 compare 2*L;
    std::map<int, int> map2; //shape1 compare 2*L + 1;
    std::map<int, int> ::iterator it;

    int count1, count2;
    count1 = 1; count2 = 1;
    for (int i = 0; i < N; ++i) {
      if (reference[L * 2][i] == shape1[i])
        count1++;
      else {
        it = map1.find(count1);
        if (it != map1.end())
          map1[count1]++;
        else
          map1.insert(std::pair<int, int>(count1, 1));
        count1 = 1;
      }

      if (reference[L * 2][i] == shape1[i])
        count2++;
      else {
        it = map2.find(count2);
        if (it != map2.end())
          map2[count2]++;
        else
          map2.insert(std::pair<int, int>(count2, 1));
        count2 = 1;
      }
    }
    fprintf (stderr, "map1 size: %d\n", map1.size()) ;
    fprintf (stderr, "map2 size: %d\n", map2.size()) ;

    map1.clear();
    map2.clear();
*/
}
static void Normalized(double **data, int s) {
  double total;
  for ( int i = 0; i < 8; ++i)
    total += data[i][s];
  for ( int i = 0; i< 8; ++i)
    data[i][s] = data[i][s]/total;
}

void pbwtMatchCount (PBWT *p, int L) /* reporting the match number for each segment */ 
{
  if (!p || !p->yz) die ("option -longWithin called without a PBWT") ;
  if (L < 0) die ("L %d for longWithin must be >= 0", L) ;

  uchar **reference = pbwtHaplotypes (p) ; /* haplotypes for reference  (M * N)  */
  uchar *x;                 /* use for current query */
  PbwtCursor *up = pbwtCursorCreate (p, TRUE, TRUE) ;
  int **a, **d, **u ;   /* stored indexes */
  int **f1, **g1 ;     /* start of match, and pbwt interval as in algorithm 5 */
  int **f2, **g2 ;    /* next versions of the above, e' etc in algorithm 5 */
  int i, j, k, N = p->N, M = p->M ;
  int s, seg_num = 1; /* for the segment number and current segment */
  int *pos; /* for the 1 position for ref */
  int num_1 = 0; /* for the num of 1 */
  uchar *shape1;  /* for the shape seq1 */
  uchar *shape2;  /* for the shape seq2 */

  /* build indexes */
  a = myalloc (N+1,int*) ; for (i = 0 ; i < N+1 ; ++i) a[i] = myalloc (p->M, int) ;
  d = myalloc (N+1,int*) ; for (i = 0 ; i < N+1 ; ++i) d[i] = myalloc (p->M+1, int) ;
  u = myalloc (N,int*) ; for (i = 0 ; i < N ; ++i) u[i] = myalloc (p->M+1, int) ;
  x = myalloc (N, uchar*) ; 
  pos = myalloc (M, int*) ;
  f1 = myalloc (8, int*);
  g1 = myalloc (8, int*);
  f2 = myalloc (64, int*);
  g2 = myalloc (64, int*);
  int *cc = myalloc (p->N, int) ;

  shape1 = myalloc (N, uchar);
  shape2 = myalloc (N, uchar);

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

  fprintf (stderr, "Made indices: \n") ; timeUpdate () ;
 
  
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

    memcpy (shape1, x, N*sizeof(uchar)) ;
    memcpy (shape2, x, N*sizeof(uchar)) ;
    fprintf (stderr, "\n\nInitial frag_num: \t\t\t\t%d\n", compare(reference, L, shape1, shape2, N));

    
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
/*
    //  print_one_seg();
    printf("the first segment\n");
    for (i = 0; i < 8; ++i)
      printf("%d\n", g1[i][0] - f1[i][0]);
    printf("\n\n");
*/

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
/* print 
   for ( i = 0; i < seg_num - 2; ++i) {
      printf("segment num\t%d\n", i + 2);
      printf("index\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",0,1,2,3,4,5,6,7);
      for (j = 0; j < 8; ++j)
        printf("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", j,
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
*/
  }

  gettimeofday( &tend, NULL );
  int timeuse = 1000000 * ( tend.tv_sec - tstart.tv_sec ) + tend.tv_usec -tstart.tv_usec;
  timeuse /= TIMES;
  printf("time: %d us\n", timeuse);
  //Shape it
  MostLikelySampling(g1, f1, g2, f2, pos, seg_num, shape1, shape2) ;
  fprintf (stderr, "\n\nAfter ML Sampling frag_num: \t\t\t\t%d\n", compare(reference, L, shape1, shape2, N));

  globalOptimalSampling(g1, f1, g2, f2, pos, seg_num, shape1, shape2) ;
  fprintf (stderr, "\n\nAfter global optimal Sampling frag_num :\t\t\t\t%d\n", compare(reference, L, shape1, shape2, N));

/*  
  for (int i = 0; i < N; ++i){
    printf("%u", reference[L * 2][i]);
    printf("%u", reference[L * 2 + 1][i]);
    printf("%u", shape1[i]);
    printf("%u", shape2[i]);
    printf("\n");
  } 
  printf("\n");
*/
  /* cleanup */
  free (cc) ;
  free (shape1) ; free (shape2) ;
  for (j = 0 ; j < p->M ; ++j) free(reference[j]) ; free (reference) ;
  free(pos);
  for ( i = 0; i < 8; ++i)  { free(f1[i]); free(g1[i]); }
  for ( i = 0; i < 64; ++i) { free(f2[i]); free(g2[i]); }
  free(f1); free(g1); free(f2); free(g2);
  for (j = 0 ; j < N ; ++j) free(a[j]) ; free (a) ;
  for (j = 0 ; j < N ; ++j) free(d[j]) ; free (d) ;
  for (j = 0 ; j < N ; ++j) free(u[j]) ; free (u) ;
}

void Sampling(int **g1, int **f1, int **g2, int **f2, int *pos, int seg_num, uchar *shape1, uchar *shape2){


}
void MostLikelySampling(int **g1, int **f1, int **g2, int **f2, int *pos, int seg_num, uchar *shape1, uchar *shape2){
  //MostLike
  int max_idx = 0;
  int prev1,prev2;
  int max_num = 0; 
  int i, s;
  for ( i = 0; i  < 4; ++i) {
    if ((g1[7-i][0] - f1[7-i][0]) * (g1[i][0] - f1[i][0]) > max_num) {
      max_idx = i;
      max_num = (g1[7-i][0] - f1[7-i][0]) * (g1[i][0] - f1[i][0]);
    }
  }
  prev1 = max_idx;
  prev2 = 7 - prev1;
  //set shape1, shape2 to be prev1, prev2;
  setSeq(shape1, pos, 0, prev1);
  setSeq(shape2, pos, 0, prev2);

  for ( s = 0; s < seg_num - 2; s++) {
    max_idx = 0;
    max_num = 0; 
    for ( i = 0; i < 8; ++i) {
      if ((g2[prev1 * 8 + i][s] - f2[prev1 * 8 + i][s]) 
        * (g2[prev2 * 8 + 7 - i][s] - f2[prev2 * 8 + 7 - i][s]) > max_num) {
        max_idx = i;
        max_num =  (g2[prev1 * 8 + i][s] - f2[prev1 * 8 + i][s]) 
                  * (g2[prev2 * 8 + 7 - i][s] - f2[prev2 * 8 + 7 - i][s]);  
      }
    }
    prev1 = max_idx;
    prev2 = 7 - prev1;
    //set shape1, shape2 to be prev1, prev2;
    setSeq(shape1, pos, s + 1, prev1);
    setSeq(shape2, pos, s + 1, prev2);
  }
}

void globalOptimalSampling(int **g1, int **f1, int **g2, int **f2, int *pos, int seg_num, uchar *shape1, uchar *shape2){
  int i,j,s;
  double **data;
  data = myalloc(8, double* );
  for ( i = 0; i < 8; ++i ) {
    data[i] = myalloc(seg_num, double*); }
  int **phis;
  phis = myalloc(8, int* );
  for ( i = 0; i < 8; ++i ) {
    phis[i] = myalloc(seg_num, int*); }

  int total = 0;
  for( i = 0; i < 8; ++i) {
    total += ( g1[i][0] - f1[i][0] ); }
  for( i = 0; i < 8; ++i) {
 //   data[i][0] = ( total == 0 ? 0 : (double)((g1[i][0] - f1[i][0]) * (g1[7 - i][0] - f1[7 - i][0])) / (total * total) ); }
    data[i][0] = ( total == 0 ? 0 : (double) (g1[i][0] - f1[i][0]) / total ); }

  Normalized(data, 0);

  int *totalCon;
  int prev;
  totalCon = myalloc(8, int);
  for( s = 0; s < seg_num - 2; ++s) {
    for ( i = 0; i < 8; ++i)
        totalCon[i] = (g1[i][s] - f1[i][s]);
    for ( i = 0; i < 8; ++i) {
      int maxIdx = 0;
      double maxVal = (totalCon[0] * totalCon[7] == 0 ? 0
                      : data[0][s] * data[0][s] * ((double)(g2[i][s]-f2[i][s])/totalCon[0]) * ((double)(g2[63 - i][s]-f2[63 - i][s])/totalCon[7]));
      for ( j = 1; j < 8; ++j) {
        double val = (totalCon[j] * totalCon[7 - j] == 0 ? 0 
                      : data[j][s] * data[j][s] * ((double)(g2[j * 8 + i][s]-f2[j * 8 + i][s])/totalCon[j]) * ((double)(g2[(7 - j) * 8 + 7 - i][s]-f2[(7 - j) * 8 + 7 - i][s])/totalCon[7 - j]));
        if (val > maxVal) { maxIdx = j; maxVal = val;} }
      data[i][s + 1] = maxVal;
      phis[i][s + 1] = maxIdx;    //from 1 - seg_Num - 2
    }
    Normalized(data, s + 1);
 
    for( i = 0; i < 8; ++i) 
      printf ("data[%d][%d] %f\t", i, s+1, data[i][s+1]);
    printf ("\n");

  }
/*
  for( s = 0; s < seg_num - 2; ++s) {
    for ( i = 0; i < 8; ++i)
        totalCon[i] = (g1[i][s] - f1[i][s]);
    for ( i = 0; i < 8; ++i) {
      int maxIdx = 0;
      double maxVal = (totalCon[0] * totalCon[7] == 0 ? 0 :
                      data[0][s] * ((double)(g2[i][s]-f2[i][s])/totalCon[0]) 
                      * data[7][s] * ((double)(g2[63 - i][s]-f2[63 - i][s])/totalCon[7]));
      for ( j = 1; j < 8; ++j) {
        double val = (totalCon[j] * totalCon[7 - j] == 0 ? 0 
                      : data[j][s] * ((double)(g2[j * 8 + i][s]-f2[j * 8 + i][s])/totalCon[j]) 
                      * data[7 - j][s] * ((double)(g2[(7 - j) * 8 + 7 - i][s]-f2[(7 - j) * 8 + 7 - i][s])/totalCon[7 - j]));
        if (val > maxVal) { maxIdx = j; maxVal = val;} }
      data[i][s + 1] = maxVal;
      phis[i][s + 1] = maxIdx;    //from 1 - seg_Num - 2
    }
    Normalized(data, s + 1);
 
    for( i = 0; i < 8; ++i) 
      printf ("data[%d][%d] %f\t", i, s+1, data[i][s+1]);
    printf ("\n");

  }
*/

  free(totalCon);

  // backtrack path
  int *path;
  path = myalloc(seg_num, int);
  double maxData = data[0][seg_num - 2];
  path[seg_num - 2] = 0;

  for ( i = 1; i < 8; ++i) {
    if ( maxData < data[i][seg_num - 2]) {
      maxData = data[i][seg_num - 2];
      path[seg_num - 2] = i; 
    }
  }

  for ( s = seg_num - 3; s >= 0; --s) {
    path[s] = phis[path[s + 1]][s + 1];
  }

  for ( s = 0; s < seg_num - 1; ++s) {
    setSeq(shape1, pos, s, path[s]);
    setSeq(shape2, pos, s, 7 - path[s]);
  }
  
  free(path);
  for (i = 0 ; i < 8 ; ++i) free(data[i]) ; free (data) ;
  for (i = 0 ; i < 8 ; ++i) free(phis[i]) ; free (phis) ;
}

/******************* end of file *******************/