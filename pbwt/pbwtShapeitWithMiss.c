#include "pbwt.h"
typedef struct Tablestruct {
  uchar* id;
  int* data;
} TableNode ;

typedef struct Tablesstruct {
  TableNode* array;
  int cap;
  int num;
} Tables ;

Tables* tablesCreate(int cap) {
  Tables* tables;
  tables = myalloc(1, Tables);
  tables->array = myalloc(cap, TableNode*);
  tables->cap = cap;
  tables->num = 0;
  return tables;
}

void tablesDestroy(Tables *tables) {
  for (int i = 0; i < tables->num; ++i) {
    free(tables->array[i].id);
    free(tables->array[i].data);
  }

  free(tables->array);
  free(tables);
}

void resizeTables(Tables *tables) {
  Tables *newTables = myalloc(1, Tables);
  newTables->array = myalloc(2 * tables->cap, TableNode);
  memcpy(newTables->array, tables->array, tables->cap * sizeof(TableNode));
  newTables->cap = 2 * tables->cap;
  newTables->num = tables->cap;
  free(tables->array);
  free(tables);
  tables = newTables;
}

void updateTable(int *het, int depth, uchar *seq, Tables *tables, int count) {
  int index = 0; 
  int start = het[5] - depth + 1;
fprintf (stderr, "up 1~~~~~~~~~~  \n");
  for (int i = 0; i < 6; ++i) {
    index << 1;
    index += seq[het[i]]; 
  }
  
fprintf (stderr, "up 2~~~~~~~~~~  \n");
fprintf (stderr, "depth  %d~~~~~~~~~~  \n", depth);
  uchar *target = myalloc(depth - 6, uchar*);
fprintf (stderr, "depth  %d~~~~~~~~~~  \n", depth);
  for (int i = 0, j = 0; j < depth; ++j) {
    
fprintf (stderr, "i =  %d,  j = %d,  \n", i, j);
    if ((j + start) == het[i]){
      i++;
      j++;  
    } else {
      target[j - i] = seq[j]; 
    }
  }

fprintf (stderr, "up 3~~~~~~~~~~  \n");
  int i;
  for (i = 0; i < tables->num; ++i) {
    if (!strcmp(target, tables->array[i].id)) {
      break;
    }
  }

fprintf (stderr, "up 4~~~~~~~~~~  \n");
  if (i != tables->num) {
    tables->array[i].data[index] = count;
    free(target);
  } else {
    if (tables->num == tables->cap) {
fprintf (stderr, "resize!!!!!!!!!!!!!!!!!!!!!!~~~~~~~~~~~~~~~~~ \n");
       resizeTables(tables);
    }

fprintf (stderr, "up 5~~~~~~~~~~  \n");
    tables->array[tables->num].data = myalloc(64, int);
    memset(tables->array[tables->num].data, 0, 64 * sizeof(int));
    tables->array[tables->num].data[index] = count;
    tables->array[tables->num].id = target;
    tables->num++;
  }

fprintf (stderr, "up 6~~~~~~~~~~  \n");
}

void extendMatch(int *het, int cur, int num, int depth, uchar *seq, int *cc, int **u, int f, int g, Tables *tables) {
fprintf (stderr, "start extend cur  %d   \t  table size  %d  \n", cur,  tables->num);
  if (f >= g) {
fprintf (stderr, "finish extend cur  %d   \t  table size  %d  \n", cur,  tables->num);
    return;
}
  //update the tables.
fprintf (stderr, "ex 1~~~~~~~~~~  \n");
  if (num == depth) {
fprintf (stderr, "ex 2~~~~~~~~~~  \n");
    updateTable(het, depth, seq, tables, g - f);
fprintf (stderr, "ex 3~~~~~~~~~~  \n");
fprintf (stderr, "update cur  %d   \t  table size  %d  \n", cur,  tables->num);
  } else {
fprintf (stderr, "3~~~~~~~~~~  \n");
    seq[num] = 1;
fprintf (stderr, "4~~~~~~~~~~  \n");
    extendMatch(het, cur + 1, num + 1, depth, seq, cc, u, cc[cur] + f - u[cur][f], cc[cur] + g - u[cur][g], tables); //1
fprintf (stderr, "5~~~~~~~~~~  \n");
    seq[num] = 0;
fprintf (stderr, "6~~~~~~~~~~  \n");
    extendMatch(het, cur + 1, num + 1, depth, seq, cc, u, u[cur][f], u[cur][g], tables); //0
fprintf (stderr, "finish extend cur  %d   \t  table size  %d  \n", cur,  tables->num);
  }
}

//void pbwtShapeItWithMiss (PBWT *p, int *geno, FILE *out) {
void pbwtShapeItWithMiss (PBWT *p, FILE *out) {
 
  if (!p || !p->yz) die ("option -longWithin called without a PBWT") ;
  
  /********   ref  *****/
  uchar **reference = pbwtHaplotypes (p) ; /* haplotypes for reference  (M * N)  */
  int *geno = myalloc(p->N, int);
  for (int i = 0; i < p->N; ++i)
    geno[i] = reference[0][i] + reference[1][i];
  /*********************/

  uchar *x;                 /* use for current query */
  PbwtCursor *up = pbwtCursorCreate (p, TRUE, TRUE) ;
  int **u ;   /* stored indexes */
  int i, j, k, N = p->N, M = p->M ;
  int num_1 = 0;      /* for the num of heterozyogous */
  int s, seg_num = 1; /* for the segment number and current segment */
  
  /* build indexes */
  u = myalloc (N,int*) ; for (i = 0 ; i < N ; ++i) u[i] = myalloc (p->M+1, int) ;
  x = myalloc (N, uchar*) ; 
  int *cc = myalloc (p->N, int) ;

  /* make pbwt index */
  for (k = 0 ; k < N ; ++k)
    { 
      cc[k] = up->c ;
      pbwtCursorCalculateU (up) ;
      memcpy (u[k], up->u, (M+1)*sizeof(int)) ;
      pbwtCursorForwardsReadAD (up, k) ;
    }
  pbwtCursorDestroy (up) ;

fprintf (stderr, "1~~~~~~~~~~~~~~~~~ \n");
  int *pos;           /* record the heterozyogous position */
  pos = myalloc (N, int*) ;
  /* find the heterozyogous position and record */
  for ( i = 0, j = 0; i < N; ++i) {
    if (geno[i] == 1) {
      ++j;
      pos[num_1++] = i;
      if (j == 3) {
        ++seg_num;
        j = 0;
      }
    }
  }

fprintf (stderr, "2~~~~~~~~~~~~~~~~~ \n");
  int start, depth;
  start = 0;
  int *het = myalloc(6, int);
  int *seq;
  Tables *tables = 0;
  //for (s = 0; s < seg_num - 2; ++s) {
  for (s = 0; s < seg_num - 2; ++s) {
      for (i = 0; i < 6; ++i) {
        het[i] = pos[s * 3 + i];
fprintf (stderr, "%d \t", het[i]);
      }
fprintf (stderr, "~~~~~~~~~~~~~~~~~ \n");

      depth = het[5] - start + 1;
      seq = myalloc(depth, uchar);
      tables = tablesCreate(20);
      extendMatch(het, start, 0, depth, seq, cc, u, 0, M, tables);
      /***********/
      fprintf (stderr, "s = %d \t table num    %d  \n", s, tables->num);
      /***********/
      start = het[2] + 1;
fprintf (stderr, "M  6~~~~~~~~~~~~~~~~~ \n");
      free(seq);
fprintf (stderr, "M  7~~~~~~~~~~~~~~~~~ \n");
      tablesDestroy(tables);
fprintf (stderr, "M  8~~~~~~~~~~~~~~~~~ \n");
  }
  free(het);

  /* cleanup */
  free(x); free(pos); free(cc);
  for (j = 0 ; j < N ; ++j) free(u[j]) ; free (u) ;
  for (j = 0 ; j < M ; ++j) free(reference[j]) ; free (reference) ;
  free(geno);
}
