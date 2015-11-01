#include "pbwt.h"
//pbwtShapeIt1 for fixed segment size (size = 3 heterozyogous genotypes)
//pbwtShapeIt2 for extensible segment size(min = 3)
//pbwtShapeIt3 for combining the original conditional probability and heterozyogous only conditional probability

/****************************************************************
Go through sequence "x" from the site position "start" to "end"
find the exactly match sequence index, start in "f" and end in "g" 
****************************************************************/
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

/****************************************************************
Initial, match sequence index, start in "ini_f" and end in "ini_g"
Extend the match sequence from the site position "start" to "end"
find the exactly match sequence index, start in "f" and end in "g" 
****************************************************************/
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

/* Set the heterozyogous sites in "seq"th blocks */
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

/* Set the heterozyogous sites in "seq"th blocks */
static void subSetSeq(uchar *dir, int seq, int index) {
  if( index > 7 || index < 0){
    printf("error\n");
    return;
  }

  dir[seq * 3] = index/4 > 0 ? 1 : 0; 
  index %= 4;
  dir[seq * 3 + 1] = index/2 > 0 ? 1 : 0;
  index %= 2;
  dir[seq * 3 + 2] = index > 0 ? 1 : 0;

  return;
}

/* Set the heterozyogous sites from start to end */
static void setSeq2(uchar *dir, int *pos, int start, int end, int index) {
  for (int bits = end - start; bits >= 0; bits--) {
    dir[pos[start++]] = index / (1<<bits) > 0 ? 1 : 0; 
    index %= (1<<bits);
  }
  return;
}

/* Normalized the data array, s th row */
static void Normalized(double **data, int s) {
  double total;
  for ( int i = 0; i < 8; ++i)
    total += data[i][s];
  for ( int i = 0; i< 8; ++i)
    data[i][s] = data[i][s]/total;
}

/* Add the weight for the data array, s th row */
static void addWeight(double **data, int s, double w) {
  for ( int i = 0; i < 8; ++i)
    data[i][s] = w / 8 + (1 - w) * data[i][s];
}

/****************************************************************
Normalized the data array, s th row
This row has size elements
****************************************************************/
static void Normalized2(double **data, int s, int size) {
  double total;
  for ( int i = 0; i < size; ++i)
    total += data[i][s];
  for ( int i = 0; i< size; ++i)
    data[i][s] = data[i][s]/total;
}

/****************************************************************
Same as the Normalized2, used in randomSampling
Normalized the data array
This row has size elements
****************************************************************/
static void NormalizedRandom(double *data, int size) {
  double total;
  for ( int i = 0; i < size; ++i)
    total += data[i];
  for ( int i = 0; i< size; ++i)
    data[i] = data[i]/total;
}

/****************************************************************
Add the weight for the data array, s th row
This row has size elements
****************************************************************/
static void addWeight2(double **data, int s, double w, int size) {
  if (size == 0) return;
  for ( int i = 0; i < size; ++i)
    data[i][s] = w / size + (1 - w) * data[i][s];
}

/****************************************************************
Same as the addWeight2, used in randomSampling
Add the weight for the data array
This row has size elements
****************************************************************/
static void addWeightRandom(double *data, double w, int size) {
  if (size == 0) return;
  for ( int i = 0; i < size; ++i)
    data[i] = w / size + (1 - w) * data[i];
}

/****************************************************************
must used after normalized data. 
randomChoose a index for s th row 
choosing i probability is proportional to data[i][s]
*****************************************************************/
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

/****************************************************************
Same as the randomChoose, used in randomSampling
must used after normalized data. 
choosing i probability is proportional to data[i]
*****************************************************************/
static int randomChooseRandom(double *data, int size){
  double tmp = data[0];
  double pre = data[0];
  data[0] = 0;
  for (int i = 1; i < size; ++i) {
    tmp = data[i];
    data[i] = pre;
    pre += tmp;
  }

  double num = rand() * 1.0;
  for (int i = size - 1; i > 0; --i)
    if (num >= RAND_MAX * data[i])
      return i;
  return 0;
}

/* create new pbwt index for heterozyogous only */
PBWT *myReadHap(uchar **reference, int *pos, int num_1, int mm) {

  PBWT *p = 0;
  int j ;
  uchar *x ;    /* original, sorted, compressed */
  int *a ;
  Array xArray = arrayCreate (10000, uchar) ;
  PbwtCursor *u ;
  int loop = 0;

  while (loop < num_1) {
    int m = 0;
    while (m < mm) {
      array(xArray,m,uchar) = reference[m][pos[loop]]; 
      m++;
    }

    if (p && m != p->M) die ("length mismatch reading haps line") ;
    
    if (!p)
    { p = pbwtCreate (m, 0) ;
      p->sites = arrayCreate(4096, Site) ;
      array(xArray,p->M,uchar) = Y_SENTINEL ; /* sentinel required for packing */
    }
    ++p->N ;
    if (!p->yz) {  /* first line; p was just made! */
      p->yz = arrayCreate(4096*32, uchar) ;
      u = pbwtCursorCreate (p, TRUE, TRUE) ;
    }
    x = arrp(xArray,0,uchar) ;
    for (j = 0 ; j < p->M ; ++j) u->y[j] = x[u->a[j]] ;
      pbwtCursorWriteForwards (u) ;
    loop++;
  }

  pbwtCursorToAFend (u, p) ;
  arrayDestroy(xArray) ; pbwtCursorDestroy (u) ;
  
  return p;  
}

/****************************************************************
Read the PBWT file and ShapeIt, output result hap into FILE out 
Each Block has exactly three heterozyogous genotypes
****************************************************************/
void pbwtShapeIt1 (PBWT *p, FILE *out) //fix heter number
{
  if (!p || !p->yz) die ("option -longWithin called without a PBWT") ;
  
  uchar **reference = pbwtHaplotypes (p) ; /* haplotypes for reference  (M * N)  */
  uchar **origin;
  uchar *x;                 /* use for current query */
  PbwtCursor *up = pbwtCursorCreate (p, TRUE, TRUE) ;
  int **u ;   /* stored indexes */
  int **f1, **g1 ;      /* one block match, start in index f1 and end in index g1 */
  int **f2, **g2 ;      /* two continuously blocks match, start in index f2, and end in index g2*/
  int i, j, k, N = p->N, M = p->M ;
  int s, seg_num = 1; /* for the segment number and current segment */
  int *pos;           /* record the heterozyogous position */
  int num_1 = 0;      /* for the num of heterozyogous */
  uchar *shape1;      /* for the shapeit seq1 */
  uchar *shape2;      /* for the shapeit seq2 */
  double w = 1.0 / M; /* the small weight for the add weight */
  
  /* build indexes */
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

  /* make pbwt index */
  for (k = 0 ; k < N ; ++k)
    { 
      cc[k] = up->c ;
      pbwtCursorCalculateU (up) ;
      memcpy (u[k], up->u, (M+1)*sizeof(int)) ;
      pbwtCursorForwardsReadAD (up, k) ;
    }
  pbwtCursorDestroy (up) ;
  origin = myalloc (2, uchar*); for (i = 0; i < 2; ++i) origin[i] = myalloc (p->N, uchar*);
  fprintf (stderr, "Made indices: \n") ; timeUpdate () ;

  int t;            //multi_time
  int TIMES = M/2;
  int L;
  for (t = 0; t < TIMES; ++t) {
    L = t;
    seg_num = 1;
    num_1 = 0;

    memcpy (origin[0], reference[2*L], N*sizeof(uchar));
    memcpy (origin[1], reference[2*L + 1], N*sizeof(uchar));
    
    /* calculate the genotype */
    for ( i = 0; i < N; ++i)
      x[i] = origin[0][i] + origin[1][i];
    
    /* find the heterozyogous position and record */
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

    /* initial f1,g1 */
    for (i = 0; i < 8; ++i)
      for (j = 0; j < seg_num; ++j) {
        f1[i][j] = 0; g1[i][j] = M;
      }

    /* one segment match count */
    int start = 0, end;
    for ( i = 0; i < seg_num - 1; ++i) {
      end = pos[i * 3 + 2] + 1;
      for ( j = 0; j < 8; ++j) {
        setSeq(x, pos, i, j);   //set the sequence for 8 states
        countHelp(x, start, end, cc, u, &f1[j][i], &g1[j][i]);   //find the match number for each state
      } 
      start = end;
    }

    //initial f2, g2
    for (i = 0; i < 64; ++i)
      for (j = 0; j < seg_num; ++j) {
        f2[i][j] = f1[i/8][j];    // initial the match number, based on the first block re
        g2[i][j] = g1[i/8][j];
      }
    start = pos[2] + 1;  
    for ( i = 1; i < seg_num - 1; ++i) {
      end = pos[i * 3 + 2] + 1;
      for ( j = 0; j < 8; ++j){
        for (int idx = 0; idx < 8; ++idx) {
          setSeq(x, pos, i, idx);    //set the sequence for 8 * 8 states
          countHelp(x, start, end, cc, u, &f2[idx + j * 8][i - 1], &g2[idx + j * 8][i - 1]);  //find the match number for each state
        }
      }
      start = end;
    }
   
    //minus the match number by 1 for the origin sequence
    for ( i = 0; i < seg_num - 1; ++i) {
      //find the state of this block
      int index = origin[0][pos[i * 3]] * 4 + origin[0][pos[i * 3 + 1]] * 2 + origin[0][pos[i * 3 + 2]];
      //increase the start position by 1, equal to minus the match number by 1
      f1[index][i]++;      //for origin[0];
      f1[7 - index][i]++;  //for origin[1];
    }

    for ( i = 0; i < seg_num - 2; ++i) {
      //find the state of two continue blocks.
      int index = origin[0][pos[i * 3]] * 32 + origin[0][pos[i * 3 + 1]] * 16 + origin[0][pos[i * 3 + 2]] * 8  //previous block
                  + origin[0][pos[i * 3 + 3]] * 4 + origin[0][pos[i * 3 + 4]] * 2 + origin[0][pos[i * 3 + 5]]; //current block
      f2[index][i]++;      //for origin[0];
      f2[63 - index][i]++;  //for origin[1];
    }

  //shapeit for this individual
  viterbiSampling1(g1, f1, g2, f2, pos, seg_num, shape1, shape2, w) ;
  memcpy (reference[2 * t], shape1, N*sizeof(uchar));
  memcpy (reference[2 * t + 1], shape2, N*sizeof(uchar));
     
  /* cleanup */
  for ( i = 0; i < 8; ++i)  { free(f1[i]); free(g1[i]); }
  for ( i = 0; i < 64; ++i) { free(f2[i]); free(g2[i]); }
  }
  
  //output the shapeit result.
  for ( j = 0; j < N; ++j) {
    for ( i = 0; i < M; ++i)
      fprintf(out, "%u ", reference[i][j]);
    fprintf(out, "\n");
  }
  fclose(out);
  
  /* cleanup */
  free (cc) ;
  for (j = 0 ; j < p->M ; ++j) free(reference[j]) ; free (reference) ;
  free (shape1) ; free (shape2) ;
  free(pos);
  for (j = 0 ; j < 2 ; ++j) free(origin[j]) ; free (origin) ;
  free(f1); free(g1); free(f2); free(g2);
  for (j = 0 ; j < N ; ++j) free(u[j]) ; free (u) ;
}

/****************************************************************
Read the PBWT file and ShapeIt, output result hap into FILE out 
Each Block has extensible heterozyogous genotypes number, but the max num is maxGeno
****************************************************************/
void pbwtShapeIt2 (PBWT *p, int maxGeno, FILE *out) //extensibile heter number
{
  if (!p || !p->yz) die ("option -longWithin called without a PBWT") ;
  uchar **reference = pbwtHaplotypes (p) ; /* haplotypes for reference  (M * N)  */
  /***************** pbwt part *******************/
  uchar *x;                 /* use for current query */
  PbwtCursor *up = pbwtCursorCreate (p, TRUE, TRUE) ;
  int **u ;   /* stored indexes */
  int *cc = myalloc (p->N, int) ;
  int i, j, k, N = p->N, M = p->M ;
  /* build indexes */
  u = myalloc (N,int*) ; for (i = 0 ; i < N ; ++i) u[i] = myalloc (p->M+1, int) ;
  x = myalloc (N, uchar*) ; 
  for (k = 0 ; k < N ; ++k)
  { 
    cc[k] = up->c ;
    pbwtCursorCalculateU (up) ;
    memcpy (u[k], up->u, (M+1)*sizeof(int)) ;
    pbwtCursorForwardsReadAD (up, k) ;
  }
  pbwtCursorDestroy (up) ;
  fprintf (stderr, "Made indices: \n") ; timeUpdate () ; 
  /**************************************/

  /***************** my algorithm part ************/
  
  int **f1, **g1 ;     /* one block match, start in index f1 and end in index g1 */
  int **f2, **g2 ;     /* two continuously blocks match, start in index f2, and end in index g2*/
  uchar **origin;
  int s, seg_num = 1;   /* for the segment number and current segment */
  int *pos;             /* record the heterozyogous position */
  int num_1 = 0;        /* for the num of heterozyogous */
  uchar *shape1;        /* for the shape seq1 */
  uchar *shape2;        /* for the shape seq2 */
  int **seg;            /* store the each segment info */
  double w = 1.0 / M;   /* the small weight for the add weight */

  pos = myalloc (N, int*) ;
  f1 = myalloc (8, int*);
  g1 = myalloc (8, int*);
  f2 = myalloc (64, int*);
  g2 = myalloc (64, int*);
  shape1 = myalloc (N, uchar);
  shape2 = myalloc (N, uchar);
  origin = myalloc (2, uchar*); for (i = 0; i < 2; ++i) origin[i] = myalloc (p->N, uchar*);

  int t;  //multi_time
  int TIMES = M/2;
  int L;
  for (t = 0; t < TIMES; ++t) {
    L = t;
    seg_num = 1;
    num_1 = 0;
    
    memcpy (origin[0], reference[2*L], N*sizeof(uchar));
    memcpy (origin[1], reference[2*L + 1], N*sizeof(uchar));
    
    /* calculate the genotype */
    for ( i = 0; i < N; ++i)
      x[i] = origin[0][i] + origin[1][i];
    
    /* find the heterozyogous position and record */
    for ( i = 0, j = 0; i < N; ++i) {
      if (x[i] == 1) {
        pos[num_1++] = i;
      }
      x[i] = x[i] / 2;  //change 0->0, 1->0, 2->1; 
    }

    memcpy (shape1, x, N*sizeof(uchar)) ;
    memcpy (shape2, x, N*sizeof(uchar)) ;
     
    seg = myalloc (11, int *) ; for (i = 0; i < 11; ++i) seg[i] = myalloc (num_1/3 + 1, int);
    /*
    store the info. for each block.
    seg[0]~seg[7]: store the 8 states for this block
    seg[8]: the first heterozyogous position in this block
    seg[9]: the last  heterozyogous position in this block
    seg[8],seg[9] decide how many heterozyogous in this block
    seg[10]: how many state in this block(max 8, min 5)
    */

    for ( i = 0; i < 8; ++i) { 
      f1[i] = myalloc(num_1/3 + 1, int*);
      g1[i] = myalloc(num_1/3 + 1, int*);
    }
    for ( i = 0; i < 64; ++i) { 
      f2[i] = myalloc(num_1/3 + 1, int*);
      g2[i] = myalloc(num_1/3 + 1, int*);
    }
    
    //one segment count
    int start = 0, end;
    s = 0;
    int count, new_count;
    int tmp_f1[8];
    int tmp_g1[8];
    int tmp_seg[8];
    int idx1, idx2;
    for ( i = 0; i < num_1 - 2;) {
      count = 0;   //record the vaild state number(match number > 0 states)
      seg[8][s] = i;   //At begin, block has 3 heterozyogous genotype
      seg[9][s] = i + 2;
      end = pos[seg[9][s]] + 1;
      for ( j = 0; j < 8; ++j) {
        setSeq2(x, pos, seg[8][s], seg[9][s], j);  //set the sequence for 8 states
        seg[count][s] = j; f1[count][s] = 0; g1[count][s] = M;
        countHelp(x, start, end, cc, u, &f1[count][s], &g1[count][s]);  //find the match number for each state
        if (g1[count][s] - f1[count][s] > 0)
          count++;
      }


      while(count < 5) {
        //heterozyogous genotype num exceeds the setting maximum, break.
        if (seg[9][s] - seg[8][s] > (maxGeno - 2)) break;
        new_count = 0;
        //end of the sequence, break;
        if (seg[9][s] < num_1 - 1) {
          start = pos[seg[9][s]] + 1;
          seg[9][s]++;
          end = pos[seg[9][s]] + 1;
        }
        else 
          break;  
        
        //extend one more heterozyogous genotype (each state increase to two states, idx1, idx2)
        for (int ii = 0; ii < count; ++ii) {
          idx1 = seg[ii][s] << 1;
          idx2 = (seg[ii][s] << 1) + 1;
          
          setSeq2(x, pos, seg[8][s], seg[9][s], idx1);
          if (countHelp2(x, start, end, cc, u, f1[ii][s], g1[ii][s], &tmp_f1[new_count], &tmp_g1[new_count])) {
            tmp_seg[new_count++] = idx1;
          }

          setSeq2(x, pos, seg[8][s], seg[9][s], idx2);
          if (countHelp2(x, start, end, cc, u, f1[ii][s], g1[ii][s], &tmp_f1[new_count], &tmp_g1[new_count])) {
            tmp_seg[new_count++] = idx2;
          }
        }

        //record the new states and match nums
        for (int ii = 0; ii < new_count; ++ii) {
          f1[ii][s] = tmp_f1[ii];
          g1[ii][s] = tmp_g1[ii];
          seg[ii][s] = tmp_seg[ii];
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
        f2[i][j] = f1[i/8][j]; // initial the match number, based on the first block result
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
    
    //minus the match number by 1 for the origin sequence (for one segment match number)
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
      if (cpl_index < seg[10][s])
         f1[cpl_index][s]++;  //for origin[1];
    }

    //minus the match number by 1 for the origin sequence (for two segment match number)
    int target1, target2;           //previous and current block encoded num
    int cpl_target1, cpl_target2;   //complementary(cpl) sequence encoded num
    int index1, index2;             //previous and current block index in seg array
    int cpl_index1, cpl_index2;     //complementary(cpl) sequence block index in seg array
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
      if (cpl_index1 < seg[10][s] && cpl_index2 < seg[10][s+1])
      	   f2[cpl_index1 * 8 + cpl_index2][s]++;  //for origin[1];
    }
  
  //shapeit for this individual
  //viterbiSampling2(seg, g1, f1, g2, f2, pos, seg_num, shape1, shape2, w) ;
  randomSampling(seg, g1, f1, g2, f2, pos, seg_num, shape1, shape2, w);
  memcpy (reference[2 * t], shape1, N*sizeof(uchar));
  memcpy (reference[2 * t + 1], shape2, N*sizeof(uchar));
  
  /* cleanup */
  for ( i = 0; i < 8; ++i)  { free(f1[i]); free(g1[i]); }
  for ( i = 0; i < 64; ++i) { free(f2[i]); free(g2[i]); }
  for ( j = 0 ; j < 11; ++j) free(seg[j]) ; 
  }
  
  for ( j = 0; j < N; ++j) {
    for ( i = 0; i < M; ++i)
      fprintf(out, "%u ", reference[i][j]);
    fprintf(out, "\n");
  }
  fclose (out);
  

  /* cleanup */
  free (cc) ;
  for (j = 0 ; j < p->M ; ++j) free(reference[j]) ; free (reference) ;
  free(x), free (shape1) ; free (shape2) ;
  free(pos);
  for (j = 0 ; j < 2 ; ++j) free(origin[j]) ; free (origin) ;
  free (seg); free(f1); free(g1); free(f2); free(g2);
  for (j = 0 ; j < N ; ++j) free(u[j]) ; free (u) ;
}

/****************************************************************
Read the PBWT file and ShapeIt, output result hap into FILE out 
Each Block has exactly three heterozyogous genotypes
Combine the whole site pbwt and  heterozyogous site only sub-pbwt
****************************************************************/
void pbwtShapeIt3 (PBWT *p, int percent, FILE *out) 
{
  if (!p || !p->yz) die ("option -longWithin called without a PBWT") ;
  uchar **reference = pbwtHaplotypes (p) ; /* haplotypes for reference  (M * N)  */
  /***************** pbwt part *******************/
  uchar *x;                 /* use for current query */
  PbwtCursor *up = pbwtCursorCreate (p, TRUE, TRUE) ;
  int **u ;   /* stored indexes */
  int *cc = myalloc (p->N, int) ;
  int i, j, k, N = p->N, M = p->M ;
  /* build indexes */
  u = myalloc (N,int*) ; for (i = 0 ; i < N ; ++i) u[i] = myalloc (p->M+1, int) ;
  x = myalloc (N, uchar*) ; 
  for (k = 0 ; k < N ; ++k)
  { 
    cc[k] = up->c ;
    pbwtCursorCalculateU (up) ;
    memcpy (u[k], up->u, (M+1)*sizeof(int)) ;
    pbwtCursorForwardsReadAD (up, k) ;
  }
  pbwtCursorDestroy (up) ;
  //  fprintf (stderr, "Made indices: \n") ; timeUpdate () ; 
  /**************************************/

  /***************** my algorithm part ************/
  int *pos;               /* record the heterozyogous position */
  int s, seg_num;         /* for the segment number and current segment */
  int num_1;              /* for the num of heterozyogous */
  double w = 1.0 / M;     /* the small weight for the add weight */
  double coefficent = (double)percent / 100;  /* percentage of heterozyogous only in decision make */

  pos = myalloc (N, int*) ;
  uchar *shape1;  /* for the shape seq1 */
  uchar *shape2;  /* for the shape seq2 */
  uchar **origin;
  shape1 = myalloc (N, uchar);
  shape2 = myalloc (N, uchar);
  origin = myalloc (2, uchar*); for (i = 0; i < 2; ++i) origin[i] = myalloc (p->N, uchar*);

  /***************** variable for whole index ************/
  int **f1, **g1 ;     /* one block match, start in index f1 and end in index g1 */
  int **f2, **g2 ;     /* two continuously blocks match, start in index f2, and end in index g2*/
  f1 = myalloc (8, int*);
  g1 = myalloc (8, int*);
  f2 = myalloc (64, int*);
  g2 = myalloc (64, int*);
  /******************************************************/
  /***************** variable for sub index ************/
  int **subf1, **subg1 ;     /* heterozyogous only */
  int **subf2, **subg2 ;     /* heterozyogous only */
  subf1 = myalloc (8, int*);
  subg1 = myalloc (8, int*);
  subf2 = myalloc (64, int*);
  subg2 = myalloc (64, int*);
  /******************************************************/


  int t;  //multi_time
  int TIMES = M/2;
  int L;
  for (t = 0; t < TIMES; ++t) {
    L = t;
    seg_num = 1; 
    num_1 = 0; 
    memcpy (origin[0], reference[2*L], N*sizeof(uchar));
    memcpy (origin[1], reference[2*L + 1], N*sizeof(uchar));
    
    /* calculate the genotype */
    for ( i = 0; i < N; ++i)
      x[i] = origin[0][i] + origin[1][i];
    
    /* find the heterozyogous position and record */
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
  

    /* create new pbwt index  */
    PBWT *subp = myReadHap(reference, pos, num_1, M);
    /**************************/
    /***************** pbwt part *******************/
    PbwtCursor *subup = pbwtCursorCreate (subp, TRUE, TRUE) ;
    int **subu ;   /* stored indexes */
    int *subcc = myalloc (subp->N, int) ;
    int subN = subp->N, subM = subp->M ;
    /* build indexes */
    subu = myalloc (subN,int*) ; for (i = 0 ; i < subN ; ++i) subu[i] = myalloc (subp->M+1, int) ;
    uchar *subx;  /* use for sub index query */
    subx = myalloc (subN, uchar*) ; 
    for (k = 0 ; k < subN ; ++k)
    { 
      subcc[k] = subup->c ;
      pbwtCursorCalculateU (subup) ;
      memcpy (subu[k], subup->u, (subM+1)*sizeof(int)) ;
      pbwtCursorForwardsReadAD (subup, k) ;
    }
    pbwtCursorDestroy (subup) ;    
    /**************************************/

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

    for ( i = 0; i < 8; ++i) { 
      subf1[i] = myalloc(seg_num, int*);
      subg1[i] = myalloc(seg_num, int*);
    }
    for ( i = 0; i < 64; ++i) { 
      subf2[i] = myalloc(seg_num, int*);
      subg2[i] = myalloc(seg_num, int*);
    }

    //initial f1,g1
    for (i = 0; i < 8; ++i)
      for (j = 0; j < seg_num; ++j) {
        f1[i][j] = 0; g1[i][j] = M;
      }

    //initial subf1, subg1
    for (i = 0; i < 8; ++i)
      for (j = 0; j < seg_num; ++j) {
        subf1[i][j] = 0; subg1[i][j] = M;
      }

    //one segment count
    int start = 0, end;
    for ( i = 0; i < seg_num - 1; ++i) {
      end = pos[i * 3 + 2] + 1;
      for ( j = 0; j < 8; ++j) {
        setSeq(x, pos, i, j);     //set the sequence for 8 states
        countHelp(x, start, end, cc, u, &f1[j][i], &g1[j][i]);  //find the match number for each state
        subSetSeq(subx, i, j);   //set the heterozyogous only sequence for 8 states
        countHelp(subx, i * 3, (i + 1) * 3, subcc, subu, &subf1[j][i], &subg1[j][i]); //find the match number for each state
      } 
      start = end;
    }

    //initial f2, g2
    for (i = 0; i < 64; ++i)
      for (j = 0; j < seg_num; ++j) {
        f2[i][j] = f1[i/8][j];  // initial the match number, based on the first block result
        g2[i][j] = g1[i/8][j];
      }
    //initial subf2, subg2
    for (i = 0; i < 64; ++i)
      for (j = 0; j < seg_num; ++j) {
        subf2[i][j] = subf1[i/8][j];   // initial the match number, based on the first block result
        subg2[i][j] = subg1[i/8][j];
      }

    start = pos[2] + 1;  
    for ( i = 1; i < seg_num - 1; ++i) {
      end = pos[i * 3 + 2] + 1;
      for ( j = 0; j < 8; ++j){
        for (int idx = 0; idx < 8; ++idx) {
          setSeq(x, pos, i, idx);
          countHelp(x, start, end, cc, u, &f2[idx + j * 8][i - 1], &g2[idx + j * 8][i - 1]);
          subSetSeq(subx, i, idx);
          countHelp(subx, i * 3, (i + 1) * 3, subcc, subu, &subf2[idx + j * 8][i - 1], &subg2[idx + j * 8][i - 1]);
        }
      }
      start = end;
    }
   
    //minus the match number by 1 for the origin sequence
    for ( i = 0; i < seg_num - 1; ++i) {
      int index = origin[0][pos[i * 3]] * 4 + origin[0][pos[i * 3 + 1]] * 2 + origin[0][pos[i * 3 + 2]];
      //increase the start position by 1, equal to minus the match number by 1
      f1[index][i]++;      //for origin[0];
      f1[7 - index][i]++;  //for origin[1];
      subf1[index][i]++;      //for heterozyogous only sequence
      subf1[7 - index][i]++;  //for heterozyogous only sequence
    }

    for ( i = 0; i < seg_num - 2; ++i) {
      //find the state of two continue blocks.
      int index = origin[0][pos[i * 3]] * 32 + origin[0][pos[i * 3 + 1]] * 16 + origin[0][pos[i * 3 + 2]] * 8  //previous block
                  + origin[0][pos[i * 3 + 3]] * 4 + origin[0][pos[i * 3 + 4]] * 2 + origin[0][pos[i * 3 + 5]]; //current block
      f2[index][i]++;      //for origin[0];
      f2[63 - index][i]++;  //for origin[1];
      subf2[index][i]++;      //for heterozyogous only sequence
      subf2[63 - index][i]++;  //for heterozyogous only sequence
    }
    
    //shapeit for this individual
    viterbiSampling3(g1, f1, g2, f2, subg1, subf1, subg2, subf2, pos, seg_num, shape1, shape2, w, coefficent) ;
    memcpy (reference[2 * t], shape1, N*sizeof(uchar));
    memcpy (reference[2 * t + 1], shape2, N*sizeof(uchar));

    /* cleanup */
    free (subcc); free (subx);
    for (j = 0 ; j < subN ; ++j) free(subu[j]) ; free (subu) ;
    if (subp) pbwtDestroy(subp) ;
    for ( i = 0; i < 8; ++i)  { free(f1[i]); free(g1[i]); }
    for ( i = 0; i < 64; ++i) { free(f2[i]); free(g2[i]); }
    for ( i = 0; i < 8; ++i)  { free(subf1[i]); free(subg1[i]); }
    for ( i = 0; i < 64; ++i) { free(subf2[i]); free(subg2[i]); }
  }

  //output the shapeit result.
  for ( j = 0; j < N; ++j) {
    for ( i = 0; i < M; ++i)
      fprintf(out, "%u ", reference[i][j]);
    fprintf(out, "\n");
  }
  
  fclose (out);
  
  /* cleanup */
  free (cc) ;
  for (j = 0 ; j < p->M ; ++j) free(reference[j]) ; free (reference) ;
  free(x), free (shape1) ; free (shape2) ; free(pos);
  for (j = 0 ; j < 2 ; ++j) free(origin[j]) ; free (origin) ;
  free(f1); free(g1); free(f2); free(g2);
  free(subf1); free(subg1); free(subf2); free(subg2);
  for (j = 0 ; j < N ; ++j) free(u[j]) ; free (u) ;
}

/****************************************************************
Using the most likely algorithm to sampling the sequence
Each block is decided based on the condition probability, condition on previous block.
****************************************************************/
void MostLikelySampling(int **g1, int **f1, int **g2, int **f2, int *pos, int seg_num, uchar *shape1, uchar *shape2, double w){
  //MostLike
  int max_idx = 0;
  int prev1,prev2;
  int max_num = 0; 
  int i, s;
  //initial, first block probablity.
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

  /* go trhough each block
     each time, choose a max condition probability */
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

/****************************************************************
Used in pbwtShapeIt1
Using viterbi algorithm to sampling the sequence
Each block has exactly 3 heterozyogous genotype
****************************************************************/
void viterbiSampling1(int **g1, int **f1, int **g2, int **f2, int *pos, int seg_num, uchar *shape1, uchar *shape2, double w){
  int i,j,s;
  double **data; //condition probability
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
    //first block, marginal probability
    data[i][0] = (total == 0 ? 0 : (double)((g1[i][0] - f1[i][0]) * (g1[7 - i][0] - f1[7 - i][0])) / (total * total)); }
  //adding small weight for 0 condition
  addWeight(data, 0, w);
  //normalized the probability
  Normalized(data, 0);

  int *totalCon;
  int prev;
  totalCon = myalloc(8, int);
  //go through each block
  for( s = 0; s < seg_num - 2; ++s) {
    for ( i = 0; i < 8; ++i)
        totalCon[i] = (g1[i][s] - f1[i][s]);
    
    //calculate the condition probability and record the max path. 
    for ( i = 0; i < 8; ++i) {
      int maxIdx = 0;
      double maxVal = data[0][s] * ((double)(g2[i][s] - f2[i][s] + 1.0/8) / (totalCon[0] + 1)) 
                      * data[7][s] * ((double)(g2[63 - i][s] - f2[63 - i][s] + 1.0/8) / (totalCon[7] + 1));
      for ( j = 1; j < 8; ++j) {
        double val = data[j][s] * ((double)(g2[j * 8 + i][s] - f2[j * 8 + i][s] + 1.0/8) / (totalCon[j] + 1))
                      * data[7 - j][s] * ((double)(g2[(7 - j) * 8 + 7 - i][s] - f2[(7 - j) * 8 + 7 - i][s] + 1.0/8) / (totalCon[7 - j] + 1));
        if (val > maxVal) { maxIdx = j; maxVal = val;} }
      data[i][s + 1] = maxVal;
      phis[i][s + 1] = maxIdx;    
    }
    addWeight(data, s + 1, w);
    Normalized(data, s + 1);
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

/****************************************************************
Used in pbwtShapeIt2
Using viterbi algorithm to sampling the sequence
Each block has flexible heterozyogous genotype
****************************************************************/
void viterbiSampling2(int **seg, int **g1, int **f1, int **g2, int **f2, int *pos, int seg_num, uchar *shape1, uchar *shape2, double w){
  int i, j, cpl_i, cpl_j, s = 0;
  int target;
  double **data; //condition probability
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
    //find the complementary(cpl) sequence index
    for (cpl_i = 0; cpl_i < seg[10][s]; ++cpl_i) {
      if(seg[cpl_i][s] == target)
        break;  
    }
    //first block, marginal probability
    data[i][0] = (total == 0 ? 0 : (double)((g1[i][0] - f1[i][0]) * (cpl_i < seg[10][s] ? (g1[cpl_i][0] - f1[cpl_i][0]) / total : 0.001)) / total); 
  }
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
      //find the complementary(cpl) sequence index of current block
      target1 = (1 << (seg[9][s+1] - seg[8][s+1] + 1)) - 1 - seg[i][s+1]; 
      for (cpl_i = 0; cpl_i < seg[10][s+1]; ++cpl_i) {
        if(seg[cpl_i][s+1] == target1)
          break;  
      }
  
      //find the complementary(cpl) sequence index of previous block
      for ( j = 0; j < seg[10][s]; ++j) {
        target2 = (1 << (seg[9][s] - seg[8][s] + 1)) - 1 - seg[j][s]; 
        for (cpl_j = 0; cpl_j < seg[10][s]; ++cpl_j) {
          if(seg[cpl_j][s] == target2)
            break;  
        }
        double val;
        if (cpl_i == seg[10][s+1] || cpl_j == seg[10][s])
          //didn't have the complementary sequence, give a small weight
          val = data[j][s] * ((double)(g2[j * 8 + i][s] - f2[j * 8 + i][s] + 1.0/8) / (totalCon[j] + 1)) * 0.001;
        else
          //calculate the condition probability for each states
          val = data[j][s] * ((double)(g2[j * 8 + i][s] - f2[j * 8 + i][s] + 1.0/8) / (totalCon[j] + 1))
                      * data[cpl_j][s] * ((double)(g2[cpl_j * 8 + cpl_i][s] - f2[cpl_j * 8 + cpl_i][s] + 1.0/8) / (totalCon[cpl_j] + 1));
        if (val > maxVal) { maxIdx = j; maxVal = val;} }
      data[i][s + 1] = maxVal;
      phis[i][s + 1] = maxIdx; 
    }
    addWeight2(data, s + 1, w, seg[10][s+1]);
    Normalized2(data, s + 1, seg[10][s+1]);

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
    //fprintf (stderr, "s: %d, path[s+1]: %d, phis[path[s+1]][s+1]: %d  \n", s, path[s+1], phis[path[s+1]][s+1]) ;
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

/****************************************************************
Used in pbwtShapeIt3
Using viterbi algorithm to sampling the sequence
Combine the condition probability of complete blocks and heterozyogous only blocks
Each block has exactly 3 heterozyogous genotype
****************************************************************/
void viterbiSampling3(int **g1, int **f1, int **g2, int **f2, int **subg1, int **subf1, int **subg2, int **subf2, 
  int *pos, int seg_num, uchar *shape1, uchar *shape2, double w, double coefficent) {
  int i,j,s;
  double **data; //condition probability
  data = myalloc(8, double* );
  for ( i = 0; i < 8; ++i ) {
    data[i] = myalloc(seg_num, double*); }
  int **phis;
  phis = myalloc(8, int* );
  for ( i = 0; i < 8; ++i ) {
    phis[i] = myalloc(seg_num, int*); }

  int total = 0;
  int subtotal = 0;
  double tmp = 0;
  for( i = 0; i < 8; ++i) {
    total += ( g1[i][0] - f1[i][0] ); 
    subtotal += ( subg1[i][0] - subf1[i][0] ); 
  }

  //first block, marginal probability
  for( i = 0; i < 8; ++i) {
    //origin sequence probability
    data[i][0] = (total == 0 ? 0 : (double)((g1[i][0] - f1[i][0]) * (g1[7 - i][0] - f1[7 - i][0])) / (total * total)); 
    //heterozyogous only probability
    tmp = (subtotal == 0 ? 0 : (double)((subg1[i][0] - subf1[i][0]) * (subg1[7 - i][0] - subf1[7 - i][0])) / (subtotal * subtotal)); 
    //combine them two together
    data[i][0] = (1 - coefficent) * data[i][0] + coefficent * tmp;
  }
  addWeight(data, 0, w);
  Normalized(data, 0);

  int *totalCon;
  int *subTotalCon; 
  int prev;
  totalCon = myalloc(8, int);
  subTotalCon = myalloc(8, int);

  for( s = 0; s < seg_num - 2; ++s) {
    for ( i = 0; i < 8; ++i) {
        totalCon[i] = (g1[i][s] - f1[i][s]);
        subTotalCon[i] = (subg1[i][s] - subf1[i][s]);
    }

    //calculate the condition probability for each states and record the max path. 
    for ( i = 0; i < 8; ++i) {
      int maxIdx = 0;
      //origin sequence probability
      double maxVal = data[0][s] * ((double)(g2[i][s] - f2[i][s] + 1.0/8) / (totalCon[0] + 1)) 
                      * data[7][s] * ((double)(g2[63 - i][s] - f2[63 - i][s] + 1.0/8) / (totalCon[7] + 1));
      //heterozyogous only probability
      tmp = data[0][s] * ((double)(subg2[i][s] - subf2[i][s] + 1.0/8) / (subTotalCon[0] + 1)) 
                      * data[7][s] * ((double)(subg2[63 - i][s] - subf2[63 - i][s] + 1.0/8) / (subTotalCon[7] + 1));
      //combine them two together
      maxVal = (1 - coefficent) * maxVal + coefficent * tmp;

      for ( j = 1; j < 8; ++j) {
        double val = data[j][s] * ((double)(g2[j * 8 + i][s] - f2[j * 8 + i][s] + 1.0/8) / (totalCon[j] + 1))
                      * data[7 - j][s] * ((double)(g2[(7 - j) * 8 + 7 - i][s] - f2[(7 - j) * 8 + 7 - i][s] + 1.0/8) / (totalCon[7 - j] + 1));
        tmp = data[j][s] * ((double)(subg2[j * 8 + i][s] - subf2[j * 8 + i][s] + 1.0/8) / (subTotalCon[j] + 1))
                   * data[7 - j][s] * ((double)(subg2[(7 - j) * 8 + 7 - i][s] - subf2[(7 - j) * 8 + 7 - i][s] + 1.0/8) / (subTotalCon[7 - j] + 1));
        val = (1 - coefficent) * val + coefficent * tmp;
        if (val > maxVal) { maxIdx = j; maxVal = val;} }
      data[i][s + 1] = maxVal;
      phis[i][s + 1] = maxIdx;   
    }
    addWeight(data, s + 1, w);
    Normalized(data, s + 1);
  }

  free(totalCon);
  free(subTotalCon);

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

/****************************************************************
Not Used  Yet 
Using viterbi algorithm to sampling the sequence
Each block is choosed randomly, but the probability is proportional to the conditional probability.
****************************************************************/
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
  }

  for ( s = 0; s < seg_num - 1; ++s) {
    setSeq(shape1, pos, s, path[s]);
    setSeq(shape2, pos, s, 7 - path[s]);
  }

  free(p);
  free(path);
  for (i = 0 ; i < 8 ; ++i) free(data[i]) ; free (data) ;
}

/****************************************************************
Not Used  Yet 
Using Random algorithm to sampling the sequence
Each block is choosed randomly, but the probability is proportional to the conditional probability.
****************************************************************/
void randomSampling(int **seg, int **g1, int **f1, int **g2, int **f2, int *pos, int seg_num, uchar *shape1, uchar *shape2, double w){
  int i, j, cpl_i, cpl_j, s = 0;
  int target;

  double *data; //condition probability
  data = myalloc(8, double);
  
  int *path;
  path = myalloc(seg_num, int);
  
  int total = 0;
  for( i = 0; i < seg[10][s]; ++i) {
    total += ( g1[i][0] - f1[i][0] ); }

  for( i = 0; i < seg[10][s]; ++i) {
    target = (1 << (seg[9][s] - seg[8][s] + 1)) - 1 - seg[i][s];    
    //find the complementary(cpl) sequence index
    for (cpl_i = 0; cpl_i < seg[10][s]; ++cpl_i) {
      if(seg[cpl_i][s] == target)
        break;  
    }
    //first block, marginal probability
    data[i] = (total == 0 ? 0 : (double)((g1[i][0] - f1[i][0]) * (cpl_i < seg[10][s] ? (g1[cpl_i][0] - f1[cpl_i][0]) / total : 0.001)) / total); 
  }
  addWeightRandom(data, w, seg[10][0]);
  NormalizedRandom(data, seg[10][0]);
  path[0] = randomChooseRandom(data, seg[10][0]);

  int prev;
  int target1, target2;
  int totalCon1, totalCon2;
  //i, target1 for current block;
  //j, target2 for prev block;
  //j = path[s];
  for( s = 0; s < seg_num - 1; ++s) {
    j = path[s];
    totalCon1 = g1[j][s] - f1[j][s];
    
    //find the complementary(cpl) sequence index of previous block
    target2 = (1 << (seg[9][s] - seg[8][s] + 1)) - 1 - seg[j][s]; 
    for (cpl_j = 0; cpl_j < seg[10][s]; ++cpl_j) {
      if(seg[cpl_j][s] == target2)
        break;
    }
    if (cpl_j == seg[10][s])
        totalCon2 = g1[cpl_j][s] - f1[cpl_j][s];
    else
        totalCon2 = 0;

    for ( i = 0; i < seg[10][s+1]; ++i) { 
      //find the complementary(cpl) sequence index of current block
      target1 = (1 << (seg[9][s+1] - seg[8][s+1] + 1)) - 1 - seg[i][s+1]; 
      for (cpl_i = 0; cpl_i < seg[10][s+1]; ++cpl_i) {
        if(seg[cpl_i][s+1] == target1)
          break;  
      }
      
      double val;
      if (cpl_i == seg[10][s+1] || cpl_j == seg[10][s])
        //didn't have the complementary sequence, give a small weight
        data[i] = ((double)(g2[j * 8 + i][s] - f2[j * 8 + i][s] + 1.0/8) / (totalCon1 + 1)) * 0.001;
      else
        //calculate the condition probability for each states
        data[i] = ((double)(g2[j * 8 + i][s] - f2[j * 8 + i][s] + 1.0/8) / (totalCon1 + 1))
                    * ((double)(g2[cpl_j * 8 + cpl_i][s] - f2[cpl_j * 8 + cpl_i][s] + 1.0/8) / (totalCon2 + 1));
    }
    addWeightRandom(data, w, seg[10][s+1]);
    NormalizedRandom(data, seg[10][s+1]);
    path[s+1] = randomChooseRandom(data, seg[10][s+1]);
  }
  
  for ( s = 0; s < seg_num; ++s) {
    setSeq2(shape1, pos, seg[8][s], seg[9][s], seg[path[s]][s]);
    setSeq2(shape2, pos, seg[8][s], seg[9][s], (1 << (seg[9][s] - seg[8][s] + 1)) - 1 - seg[path[s]][s]);
  }
  
  free(path);
  free (data);
}
/******************* end of file *******************/
