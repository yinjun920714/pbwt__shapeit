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

static int countHelp2(uchar *x, int start, int end, int *cc, int **u, int ini_f, int ini_g, int *f, int *g) {
  (*f) = ini_f;
  (*g) = ini_g;
  for ( int k = start; k < end; ++k ) {
    if ((*f) < (*g)) {
      (*f) = x[k] ? cc[k] + ((*f) - u[k][(*f)]) : u[k][(*f)];
      (*g) = x[k] ? cc[k] + ((*g) - u[k][(*g)]) : u[k][(*g)]; 
    } else {
      (*f) = 0; (*g) = 0;
      return 0;
    } 
  }
  
  return (*f) < (*g) ? 1 : 0;
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

static void setSeq2(uchar *dir, int *pos, int start, int end, int index) {
  for (int bits = end - start; bits >= 0; bits--) {
    dir[pos[start++]] = index / (1<<bits) > 0 ? 1 : 0; 
    index %= (1<<bits);
  }
  return;
}

static int compare(uchar **origin, uchar *shape1, uchar *shape2, int N) {
    int count = 0;
    uchar *cur;
    cur = shape1;
    for ( int i = 0; i < N; ++i ) {
      if(cur[i] != origin[0][i] ) {
        count++;
        cur = ( cur == shape1 ? shape2 : shape1 ); 
      }  
    }
    return count;
}

static void Normalized(double **data, int s) {
  double total;
  for ( int i = 0; i < 8; ++i)
    total += data[i][s];
  for ( int i = 0; i< 8; ++i)
    data[i][s] = data[i][s]/total;
}

static void addWeight(double **data, int s, double w) {
  for ( int i = 0; i < 8; ++i)
    data[i][s] = w / 8 + (1 - w) * data[i][s];
}

static void Normalized2(double **data, int s, int size) {
  double total;
  for ( int i = 0; i < size; ++i)
    total += data[i][s];
  for ( int i = 0; i< size; ++i)
    data[i][s] = data[i][s]/total;
}

static void addWeight2(double **data, int s, double w, int size) {
  if (size == 0) return;
  for ( int i = 0; i < size; ++i)
    data[i][s] = w / size + (1 - w) * data[i][s];
}

//must used after normalized data.
static int randomChoose(double **data, double *p, int s){
  p[0] = 0;
  for (int i = 1; i < 8; ++i)
    p[i] = data[i - 1][s] + p[i - 1];

  double num = rand() * 1.0;
  for (int i = 7; i > 0; --i)
    if (num >= RAND_MAX * p[i])
      return i;
  return 0;
}

void pbwtMatchCount1 (PBWT *p, FILE *fp) /* reporting the match number for each segment */ 
{
  if (!p || !p->yz) die ("option -longWithin called without a PBWT") ;
  //if (L < 0) die ("L %d for longWithin must be >= 0", L) ;
  //uchar **reference = pbwtHaplotypes (p) ; /* haplotypes for reference  (M * N)  */
  uchar **reference;
  reference = myalloc(p->M, uchar*) ; for (int i = 0; i < p->M; ++i) reference[i] = myalloc(p->N, uchar*);
  uchar ch;
  for (int j = 0; j < p->N; ++j)
    for (int i = 0; i < p->M; ++i)
      { 
        ch = fgetc(fp);
        while(ch == ' ' || ch == '\n')
          ch = fgetc(fp);  
        reference[i][j] = ch -'0'; 
      }
  /*
  for (int j = 0; j < p->N; ++j) {
    for (int i = 0; i < p->M; ++i)
      printf("%u ", reference[i][j]);
    printf("\n");
  }
  */
  uchar **origin;
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
  double w = 1.0 / M;
  
  uchar **newHap = myalloc(M, uchar*) ; for (i = 0; i < M; ++i) newHap[i] = myalloc(N, uchar*);  
  
  /* build indexes */
  a = myalloc (N+1,int*) ; for (i = 0 ; i < N+1 ; ++i) a[i] = myalloc (p->M, int) ;
  d = myalloc (N+1,int*) ; for (i = 0 ; i < N+1 ; ++i) d[i] = myalloc (p->M+1, int) ;
  u = myalloc (N,int*) ; for (i = 0 ; i < N ; ++i) u[i] = myalloc (p->M+1, int) ;
  x = myalloc (N, uchar*) ; 
  pos = myalloc (N, int*) ;
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
  

  origin = myalloc (2, uchar*); for (i = 0; i < 2; ++i) origin[i] = myalloc (p->N, uchar*);

  //  for (j = 0 ; j < p->M ; ++j) free(reference[j]) ; free (reference) ;
  for (j = 0 ; j < N ; ++j) free(a[j]) ; free (a) ;
  for (j = 0 ; j < N ; ++j) free(d[j]) ; free (d) ;
  
  fprintf (stderr, "Made indices: \n") ; timeUpdate () ;

  /* for time cal
  struct timeval tstart, tend;
  gettimeofday( &tstart, NULL );
  */

  int t;  //multi_time
  int TIMES = M/2;
  int L;
  for (t = 0; t < TIMES; ++t) {
    // for time repeat
    //L = rand()%(M/2); 
    L = t;
    seg_num = 0;
    num_1 = 0;

    if (t != 0){
      for ( i = 0; i < 8; ++i)  { free(f1[i]); free(g1[i]); }
      for ( i = 0; i < 64; ++i) { free(f2[i]); free(g2[i]); }
    }
    memcpy (origin[0], reference[2*L], N*sizeof(uchar));
    memcpy (origin[1], reference[2*L + 1], N*sizeof(uchar));
    
    for ( i = 0; i < N; ++i)
      x[i] = origin[0][i] + origin[1][i];
    
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
      for ( j = 0; j < 8; ++j) {
        setSeq(x, pos, i, j);
        countHelp(x, start, end, cc, u, &f1[j][i], &g1[j][i]);
      } 
      start = end;
    }

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
        for (int idx = 0; idx < 8; ++idx) {
          setSeq(x, pos, i, idx);
          countHelp(x, start, end, cc, u, &f2[idx + j * 8][i - 1], &g2[idx + j * 8][i - 1]);
        }
      }
      start = end;
    }
   
    //minus the f1,g1 by 1 for the origin sequence
    for ( i = 0; i < seg_num - 1; ++i) {
      int index = origin[0][pos[i * 3]] * 4 + origin[0][pos[i * 3 + 1]] * 2 + origin[0][pos[i * 3 + 2]];
      f1[index][i]++;      //for origin[0];
      f1[7 - index][i]++;  //for origin[1];
    }

    for ( i = 0; i < seg_num - 2; ++i) {
      int index = origin[0][pos[i * 3]] * 32 + origin[0][pos[i * 3 + 1]] * 16 + origin[0][pos[i * 3 + 2]] * 8  //previous block
                  + origin[0][pos[i * 3 + 3]] * 4 + origin[0][pos[i * 3 + 4]] * 2 + origin[0][pos[i * 3 + 5]]; //current block
      f2[index][i]++;      //for origin[0];
      f2[63 - index][i]++;  //for origin[1];
    }

    //print f1, g1
    /*
    for (i = 0; i < seg_num - 1; ++i) {
      printf("segment num\t%d\n", i+1);
      for (j = 0; j < 8; ++j)
        printf("%d\t", g1[j][i] - f1[j][i]);
      printf("\n");
    }
    // print the count;
    // print_one_seg;
    printf("the first segment\n");
    for (i = 0; i < 8; ++i)
      printf("%d\n", g1[i][0] - f1[i][0]);
    printf("\n\n");
    // print following 
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
    
  


  //fprintf (stderr, "countTheMatch\n");
  //timeUpdate () ;
  // gettimeofday( &tend, NULL );
  // int timeuse = 1000000 * ( tend.tv_sec - tstart.tv_sec ) + tend.tv_usec -tstart.tv_usec;
  // timeuse /= TIMES;
  // printf("time: %d us\n", timeuse);
 
  //Shape it
  //fprintf (stderr, "MostLikelySampling\n");
  //MostLikelySampling(g1, f1, g2, f2, pos, seg_num, shape1, shape2, w) ;
  //fprintf (stderr, "After ML Sampling frag_num: \t\t\t\t%d\n", compare(origin,  shape1, shape2, N));
  //timeUpdate () ;
  



  //fprintf (stderr, "globalSamplming\n");
  viterbiSampling1(g1, f1, g2, f2, pos, seg_num, shape1, shape2, w) ;
  memcpy (newHap[2 * t], shape1, N*sizeof(uchar));
  memcpy (newHap[2 * t + 1], shape2, N*sizeof(uchar));
  //fprintf (stderr, "After global optimal Sampling frag_num :\t\t\t\t%d\n", compare(origin, shape1, shape2, N));
  }
  
  for ( j = 0; j < N; ++j) {
    for ( i = 0; i < M; ++i)
      printf("%u ", newHap[i][j]);
    printf("\n");
  }
  
  /* cleanup */
  free (cc) ;
  for (j = 0 ; j < p->M ; ++j) free(reference[j]) ; free (reference) ;
  for (j = 0 ; j < p->M ; ++j) free(newHap[j]) ; free (newHap) ;
  free (shape1) ; free (shape2) ;
  free(pos);
  for (j = 0 ; j < 2 ; ++j) free(origin[j]) ; free (origin) ;
  for ( i = 0; i < 8; ++i)  { free(f1[i]); free(g1[i]); }
  for ( i = 0; i < 64; ++i) { free(f2[i]); free(g2[i]); }
  free(f1); free(g1); free(f2); free(g2);
  for (j = 0 ; j < N ; ++j) free(u[j]) ; free (u) ;
}

void pbwtMatchCount2 (PBWT *p, FILE *fp) /* reporting the match number for each segment */ 
{
  if (!p || !p->yz) die ("option -longWithin called without a PBWT") ;
  //if (L < 0) die ("L %d for longWithin must be >= 0", L) ;
  //uchar **reference = pbwtHaplotypes (p) ; /* haplotypes for reference  (M * N)  */
  uchar **reference;
  reference = myalloc(p->M, uchar*) ; for (int i = 0; i < p->M; ++i) reference[i] = myalloc(p->N, uchar*);
  uchar ch;
  for (int j = 0; j < p->N; ++j)
    for (int i = 0; i < p->M; ++i)
      { 
        ch = fgetc(fp);
        while(ch == ' ' || ch == '\n')
          ch = fgetc(fp);  
        reference[i][j] = ch -'0'; 
      }
  /*
  for (int j = 0; j < p->N; ++j) {
    for (int i = 0; i < p->M; ++i)
      printf("%u ", reference[i][j]);
    printf("\n");
  }
  */
  uchar **origin;
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
  int **seg;      /* store the seg info */
  double w = 1.0 / M;
  
  uchar **newHap = myalloc(M, uchar*) ; for (i = 0; i < M; ++i) newHap[i] = myalloc(N, uchar*);  
  
  /* build indexes */
  a = myalloc (N+1,int*) ; for (i = 0 ; i < N+1 ; ++i) a[i] = myalloc (p->M, int) ;
  d = myalloc (N+1,int*) ; for (i = 0 ; i < N+1 ; ++i) d[i] = myalloc (p->M+1, int) ;
  u = myalloc (N,int*) ; for (i = 0 ; i < N ; ++i) u[i] = myalloc (p->M+1, int) ;
  x = myalloc (N, uchar*) ; 
  pos = myalloc (N, int*) ;
  f1 = myalloc (8, int*);
  g1 = myalloc (8, int*);
  f2 = myalloc (64, int*);
  g2 = myalloc (64, int*);
  int *cc = myalloc (p->N, int) ;
  seg = myalloc (11, int *) ; for (i = 0; i < 11; ++i) seg[i] = myalloc (N/3 + 1, int);

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
  

  origin = myalloc (2, uchar*); for (i = 0; i < 2; ++i) origin[i] = myalloc (p->N, uchar*);

  //  for (j = 0 ; j < p->M ; ++j) free(reference[j]) ; free (reference) ;
  for (j = 0 ; j < N ; ++j) free(a[j]) ; free (a) ;
  for (j = 0 ; j < N ; ++j) free(d[j]) ; free (d) ;
  
  fprintf (stderr, "Made indices: \n") ; timeUpdate () ;

  /* for time cal
  struct timeval tstart, tend;
  gettimeofday( &tstart, NULL );
  */

  int t;  //multi_time
  int TIMES = M/2;
  int L;
  for (t = 0; t < TIMES; ++t) {
    // for time repeat
    //L = rand()%(M/2); 
    L = t;
    seg_num = 1;
    num_1 = 0;
    if (t != 0){
      for ( i = 0; i < 8; ++i)  { free(f1[i]); free(g1[i]); }
      for ( i = 0; i < 64; ++i) { free(f2[i]); free(g2[i]); }
    }
    memcpy (origin[0], reference[2*L], N*sizeof(uchar));
    memcpy (origin[1], reference[2*L + 1], N*sizeof(uchar));
    
    for ( i = 0; i < N; ++i)
      x[i] = origin[0][i] + origin[1][i];
    
    for ( i = 0, j = 0; i < N; ++i) {
      if (x[i] == 1) {
        pos[num_1++] = i;
      }
      x[i] = x[i] / 2;  //change 0->0, 1->0, 2->1; 
    }

    memcpy (shape1, x, N*sizeof(uchar)) ;
    memcpy (shape2, x, N*sizeof(uchar)) ;
    
    for ( i = 0; i < 8; ++i) { 
      f1[i] = myalloc(N/3 + 1, int*);
      g1[i] = myalloc(N/3 + 1, int*);
    }
    for ( i = 0; i < 64; ++i) { 
      f2[i] = myalloc(N/3 + 1, int*);
      g2[i] = myalloc(N/3 + 1, int*);
    }
    //one segment count
    int start = 0, end;
    s = 0;
    int count, new_count;
    int tmp_f1[8];
    int tmp_g1[8];
    int idx1, idx2;
    for ( i = 0; i < num_1 - 2;) {
      count = 0;
      seg[8][s] = i;
      seg[9][s] = i + 2;
      end = pos[seg[9][s]] + 1;
      for ( j = 0; j < 8; ++j) {
        setSeq2(x, pos, seg[8][s], seg[9][s], j);
        seg[count][s] = j; f1[count][s] = 0; g1[count][s] = M;
        countHelp(x, start, end, cc, u, &f1[count][s], &g1[count][s]);
        if (g1[count][s] - f1[count][s] > 0)
          count++;
      }

      while(count < 5) {
        if (seg[9][s] - seg[8][s] > 2) break;
        new_count = 0;
        if (seg[9][s] < num_1 - 1) {
          start = pos[seg[9][s]] + 1;
          seg[9][s]++;
          end = pos[seg[9][s]] + 1;
        }
        else 
          break;  
        
        for (int ii = 0; ii < count; ++ii) {
          idx1 = seg[ii][s] << 1;
          idx2 = (seg[ii][s] << 1) + 1;
          
          setSeq2(x, pos, seg[8][s], seg[9][s], idx1);
          if (countHelp2(x, start, end, cc, u, f1[ii][s], g1[ii][s], &tmp_f1[new_count], &tmp_g1[new_count])) {
            seg[new_count++][s] = idx1;
          }

          setSeq2(x, pos, seg[8][s], seg[9][s], idx2);
          if (countHelp2(x, start, end, cc, u, f1[ii][s], g1[ii][s], &tmp_f1[new_count], &tmp_g1[new_count])) {
            seg[new_count++][s] = idx2;  
          }
        }

        for (int ii = 0; ii < new_count; ++ii) {
          f1[ii][s] = tmp_f1[ii];
          g1[ii][s] = tmp_g1[ii];
        }

      count = new_count;
      }
      start = pos[seg[9][s]] + 1;
      i = seg[9][s] + 1;
      seg[10][s] = count;
      s++;
    }

    seg_num = s;

    //initial f2, g2
    for (i = 0; i < 64; ++i)
      for (j = 0; j < seg_num; ++j) {
        f2[i][j] = f1[i/8][j]; 
        g2[i][j] = g1[i/8][j];
      }

    start = pos[seg[8][1]] + 1;  
    for ( s = 1; s < seg_num; ++s) {
      end = pos[seg[9][s]] + 1;
      for ( i = 0; i < seg[10][s - 1]; ++i ) {
        for ( j = 0; j < seg[10][s]; ++j ) {
          setSeq2(x, pos, seg[8][s], seg[9][s], seg[j][s]);
          countHelp(x, start, end, cc, u, &f2[j + i * 8][s - 1], &g2[j + i * 8][s - 1]);
        }
      } 
      start = end;
    }
    
    
    //minus the f1,g1 by 1 for the origin sequence
    int target;
    int cpl_target;
    int index;
    int cpl_index;
    int bits;
    for ( s = 0; s < seg_num; ++s) {
      target = 0;
      bits = seg[9][s] - seg[8][s];
      for (i = 0; i <= bits; ++i)
          target += (origin[0][pos[seg[8][s] + i]] << (bits - i));
      cpl_target = (1 << (bits + 1)) - 1 - target;
      cpl_index = seg[10][s];
      for (i = 0; i < seg[10][s]; ++i) {
          if (seg[i][s] == target)
              index = i;
          else if (seg[i][s] ==  cpl_target)
              cpl_index = i;
      }  
      f1[index][s]++;      //for origin[0];
//something wrong!      
      if (cpl_index < seg[10][s])
         f1[cpl_index][s]++;  //for origin[1];
    }


    int target1, target2;
    int cpl_target1, cpl_target2;
    int index1, index2;
    int cpl_index1, cpl_index2;
    int bits1, bits2;
    
    for ( s = 0; s < seg_num - 1; ++s) {
      target1 = 0;
      bits1 = seg[9][s] - seg[8][s];
      for (i = 0; i <= bits1; ++i)
          target1 += (origin[0][pos[seg[8][s] + i]] << (bits1 - i));
      cpl_target1 = (1 << (bits1 + 1)) - 1 - target1;
      cpl_index1 = seg[10][s];
      for (i = 0; i < seg[10][s]; ++i) {
          if (seg[i][s] == target1)
              index1 = i;
          else if (seg[i][s] ==  cpl_target1)
              cpl_index1 = i;
      }

      target2 = 0;
      bits2 = seg[9][s+1] - seg[8][s+1];
      for (i = 0; i <= bits2; ++i)
          target2 += (origin[0][pos[seg[8][s+1] + i]] << (bits2 - i));
      cpl_target2 = (1 << (bits2 + 1)) - 1 - target2;
      cpl_index2 = seg[10][s+1];
      for (i = 0; i < seg[10][s+1]; ++i) {
          if (seg[i][s+1] == target2)
              index2 = i;
          else if (seg[i][s+1] ==  cpl_target2)
              cpl_index2 = i;
      }

      f2[index1 * 8 + index2][s]++;      //for origin[0];
//something wrong      
      if (cpl_index1 < seg[10][s] && cpl_index2 < seg[10][s+1])
      	   f2[cpl_index1 * 8 + cpl_index2][s]++;  //for origin[1];
    }

    //print f1, g1
    /*
    for (s = 0; s < seg_num; ++s) {
      printf("segment num\t%d\n", s);
      for (j = 0; j < seg[10][s]; ++j)
        printf("%d\t:%d\t", seg[j][s], g1[j][s] - f1[j][s]);
      printf("\n");
    }
    // print the count;
    // print_one_seg;
    printf("the first segment\n");
    for (i = 0; i < 8; ++i)
      printf("%d\n", g1[i][0] - f1[i][0]);
    printf("\n\n");
    // print following 
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
  


  //fprintf (stderr, "countTheMatch\n");
  //timeUpdate () ;
  // gettimeofday( &tend, NULL );
  // int timeuse = 1000000 * ( tend.tv_sec - tstart.tv_sec ) + tend.tv_usec -tstart.tv_usec;
  // timeuse /= TIMES;
  // printf("time: %d us\n", timeuse);
 
  //Shape it
  //fprintf (stderr, "MostLikelySampling\n");
  //MostLikelySampling(g1, f1, g2, f2, pos, seg_num, shape1, shape2, w) ;
  //fprintf (stderr, "After ML Sampling frag_num: \t\t\t\t%d\n", compare(origin,  shape1, shape2, N));
  //timeUpdate () ;
  
  //fprintf (stderr, "globalSamplming\n");
  //viterbiRandomSampling(g1, f1, g2, f2, pos, seg_num, shape1, shape2, w) ;
  //memcpy (newHap[2 * t], shape1, N*sizeof(uchar));
  //memcpy (newHap[2 * t + 1], shape2, N*sizeof(uchar));
  //fprintf (stderr, "After global optimal Sampling frag_num :\t\t\t\t%d\n", compare(origin, shape1, shape2, N));

  viterbiSampling2(seg, g1, f1, g2, f2, pos, seg_num, shape1, shape2, w) ;
  memcpy (newHap[2 * t], shape1, N*sizeof(uchar));
  memcpy (newHap[2 * t + 1], shape2, N*sizeof(uchar));
  }
  
  
  for ( j = 0; j < N; ++j) {
    for ( i = 0; i < M; ++i)
      printf("%u ", newHap[i][j]);
    printf("\n");
  }
  

  /* cleanup */
  free (cc) ;
  for (j = 0 ; j < p->M ; ++j) free(reference[j]) ; free (reference) ;
  for (j = 0 ; j < p->M ; ++j) free(newHap[j]) ; free (newHap) ;
  for (j = 0 ; j < 11; ++j) free(seg[j]) ; free (seg);
  free (shape1) ; free (shape2) ;
  free(pos);
  for (j = 0 ; j < 2 ; ++j) free(origin[j]) ; free (origin) ;
  for ( i = 0; i < 8; ++i)  { free(f1[i]); free(g1[i]); }
  for ( i = 0; i < 64; ++i) { free(f2[i]); free(g2[i]); }
  free(f1); free(g1); free(f2); free(g2);
  for (j = 0 ; j < N ; ++j) free(u[j]) ; free (u) ;
}

void Sampling(int **g1, int **f1, int **g2, int **f2, int *pos, int seg_num, uchar *shape1, uchar *shape2, double w){


}
void MostLikelySampling(int **g1, int **f1, int **g2, int **f2, int *pos, int seg_num, uchar *shape1, uchar *shape2, double w){
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

void viterbiSampling1(int **g1, int **f1, int **g2, int **f2, int *pos, int seg_num, uchar *shape1, uchar *shape2, double w){
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
    data[i][0] = (total == 0 ? 0 : (double)((g1[i][0] - f1[i][0]) * (g1[7 - i][0] - f1[7 - i][0])) / (total * total)); }
  //    data[i][0] = ( total == 0 ? 0 : (double) (g1[i][0] - f1[i][0]) / total ); }
  addWeight(data, 0, w);
  Normalized(data, 0);

  int *totalCon;
  int prev;
  totalCon = myalloc(8, int);
  for( s = 0; s < seg_num - 2; ++s) {
    for ( i = 0; i < 8; ++i)
        totalCon[i] = (g1[i][s] - f1[i][s]);
    for ( i = 0; i < 8; ++i) {
      int maxIdx = 0;
      double maxVal = data[0][s] * ((double)(g2[i][s] - f2[i][s] + 1.0/8) / (totalCon[0] + 1)) 
                      * data[7][s] * ((double)(g2[63 - i][s] - f2[63 - i][s] + 1.0/8) / (totalCon[7] + 1));
      for ( j = 1; j < 8; ++j) {
        double val = data[j][s] * ((double)(g2[j * 8 + i][s] - f2[j * 8 + i][s] + 1.0/8) / (totalCon[j] + 1))
                      * data[7 - j][s] * ((double)(g2[(7 - j) * 8 + 7 - i][s] - f2[(7 - j) * 8 + 7 - i][s] + 1.0/8) / (totalCon[7 - j] + 1));
        if (val > maxVal) { maxIdx = j; maxVal = val;} }
      data[i][s + 1] = maxVal;
      phis[i][s + 1] = maxIdx;    //from 1 - seg_Num - 2
    }
    addWeight(data, s + 1, w);
    Normalized(data, s + 1);
  /*
    for( i = 0; i < 8; ++i) 
      printf ("data[%d][%d] %f\t", i, s+1, data[i][s+1]);
    printf ("\n");
  */
  }

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

void viterbiSampling2(int **seg, int **g1, int **f1, int **g2, int **f2, int *pos, int seg_num, uchar *shape1, uchar *shape2, double w){
  int i, j, cpl_i, cpl_j, s = 0;
  int target;
  double **data;
  data = myalloc(8, double* );
  for ( i = 0; i < 8; ++i ) {
    data[i] = myalloc(seg_num, double*); }
  int **phis;
  phis = myalloc(8, int* );
  for ( i = 0; i < 8; ++i ) {
    phis[i] = myalloc(seg_num, int*); }

  int total = 0;
  for( i = 0; i < seg[10][s]; ++i) {
    total += ( g1[i][0] - f1[i][0] ); }

  for( i = 0; i < seg[10][s]; ++i) {
    target = (1 << (seg[9][s] - seg[8][s] + 1)) - 1 - seg[i][s];    
    for (int cpl_i = 0; cpl_i < seg[10][s]; ++cpl_i) {
      if(seg[cpl_i][s] == target)
        break;  
    }
    data[i][0] = (total == 0 ? 0 : (double)((g1[i][0] - f1[i][0]) * (cpl_i < seg[10][s] ? (g1[cpl_i][0] - f1[cpl_i][0]) / total : 0.001)) / total); 
  }
  //    data[i][0] = ( total == 0 ? 0 : (double) (g1[i][0] - f1[i][0]) / total ); }
  addWeight2(data, 0, w, seg[10][0]);
  Normalized2(data, 0, seg[10][0]);

  int *totalCon;
  int prev;
  int target1, target2;
  totalCon = myalloc(8, int);
  //i, target1 for current block;
  //j, target2 for prev block;
  for( s = 0; s < seg_num - 1; ++s) {
    for ( j = 0; j < seg[10][s]; ++j)
        totalCon[j] = (g1[j][s] - f1[j][s]);
    
    for ( i = 0; i < seg[10][s+1]; ++i) {
      int maxIdx = 0;
      double maxVal = 0;
      target1 = (1 << (seg[9][s+1] - seg[8][s+1] + 1)) - 1 - seg[i][s+1]; 
      for (cpl_i = 0; cpl_i < seg[10][s+1]; ++cpl_i) {
        if(seg[cpl_i][s+1] == target1)
          break;  
      }

      for ( j = 0; j < seg[10][s]; ++j) {
        target2 = (1 << (seg[9][s] - seg[8][s] + 1)) - 1 - seg[j][s]; 
        for (cpl_j = 0; cpl_j < seg[10][s]; ++cpl_j) {
          if(seg[cpl_j][s] == target2)
            break;  
        }
        double val;
        if (cpl_i == seg[10][s+1] || cpl_j == seg[10][s])
          //val = 0.0001;
          val = data[j][s] * ((double)(g2[j * 8 + i][s] - f2[j * 8 + i][s] + 1.0/8) / (totalCon[j] + 1)) * 0.001;
        else
          val = data[j][s] * ((double)(g2[j * 8 + i][s] - f2[j * 8 + i][s] + 1.0/8) / (totalCon[j] + 1))
                      * data[cpl_j][s] * ((double)(g2[cpl_j * 8 + cpl_i][s] - f2[cpl_j * 8 + cpl_i][s] + 1.0/8) / (totalCon[cpl_j] + 1));
        if (val > maxVal) { maxIdx = j; maxVal = val;} }
      data[i][s + 1] = maxVal;
      phis[i][s + 1] = maxIdx;    //from 1 - seg_Num - 2
    }
    addWeight2(data, s + 1, w, seg[10][s+1]);
    Normalized2(data, s + 1, seg[10][s+1]);
  /* 
    for( i = 0; i < seg[10][s+1]; ++i) 
      printf ("data[%d][%d] %f %d\t", i, s+1, data[i][s+1], phis[i][s+1]);
    printf ("\n");
  */
  }

  free(totalCon);

  // backtrack path
  int *path;
  path = myalloc(seg_num, int);
  double maxData = data[0][seg_num - 1];
  path[seg_num - 1] = 0;

  for ( i = 1; i < seg[10][seg_num - 1]; ++i) {
    if ( maxData < data[i][seg_num - 1]) {
      maxData = data[i][seg_num - 1];
      path[seg_num - 1] = i; 
    }
  }

  for ( s = seg_num - 2; s >= 0; --s) {
    path[s] = phis[path[s + 1]][s + 1];
  }
  //debug
  //for ( s = 0; s < seg_num; ++s)
  //  fprintf (stderr, "path[%d]\t %d: actually:  %d\t%d\t%d\t%d\n", s, path[s], seg[path[s]][s], seg[8][s], seg[9][s], pos[seg[8][s]]) ;

  for ( s = 0; s < seg_num; ++s) {
    setSeq2(shape1, pos, seg[8][s], seg[9][s], seg[path[s]][s]);
    setSeq2(shape2, pos, seg[8][s], seg[9][s], (1 << (seg[9][s] - seg[8][s] + 1)) - 1 - seg[path[s]][s]);
  }
  
  free(path);
  for (i = 0 ; i < 8 ; ++i) free(data[i]) ; free (data) ;
  for (i = 0 ; i < 8 ; ++i) free(phis[i]) ; free (phis) ;
}

void viterbiRandomSampling(int **g1, int **f1, int **g2, int **f2, int *pos, int seg_num, uchar *shape1, uchar *shape2, double w){
  int i,j,s;
  double **data;
  data = myalloc(8, double* );
  for ( i = 0; i < 8; ++i ) {
    data[i] = myalloc(seg_num, double*); }

  int *path;
  path = myalloc(seg_num, int);

  double *p;
  p = myalloc(8, double);

  int total = 0;
  for( i = 0; i < 8; ++i) {
    total += ( g1[i][0] - f1[i][0] ); }
  for( i = 0; i < 8; ++i) {
    data[i][0] = (total == 0 ? 0 : (double)((g1[i][0] - f1[i][0]) * (g1[7 - i][0] - f1[7 - i][0])) / (total * total)); }
  addWeight(data, 0, w);
  Normalized(data, 0);
  path[0] = randomChoose(data, p, 0);

  int prev = path[0];
  int totalCon1, totalCon2;
  for( s = 0; s < seg_num - 2; ++s) {
    totalCon1 = (g1[prev][s] - f1[prev][s]);
    totalCon2 = (g1[7 - prev][s] - f1[7 - prev][s]);
    for ( i = 0; i < 8; ++i) {
      data[i][s + 1] = ((double)(g2[prev * 8 + i][s] - f2[prev * 8 + i][s] + 1.0/8) / (totalCon1 + 1))
                      * ((double)(g2[(7 - prev) * 8 + 7 - i][s] - f2[(7 - prev) * 8 + 7 - i][s] + 1.0/8) / (totalCon2 + 1));
    }
    addWeight(data, s + 1, w);
    Normalized(data, s + 1);
    path[s + 1] = randomChoose(data, p, s + 1);
    prev = path[s + 1];
  /*
    for( i = 0; i < 8; ++i) 
      printf ("data[%d][%d] %f\t", i, s+1, data[i][s+1]);
    printf ("\n");
  */
  }

  for ( s = 0; s < seg_num - 1; ++s) {
    setSeq(shape1, pos, s, path[s]);
    setSeq(shape2, pos, s, 7 - path[s]);
  }

  free(p);
  free(path);
  for (i = 0 ; i < 8 ; ++i) free(data[i]) ; free (data) ;
}

/******************* end of file *******************/
