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
  Tables* tables = 0;
  tables = myalloc(1, Tables);
//fprintf (stderr, "fun: c  1~~~~~~\n");
  tables->array = myalloc(cap, TableNode);
//fprintf (stderr, "fun: c  2~~~~~~\n");
  tables->cap = cap;
  tables->num = 0;
//fprintf (stderr, "fun: c  4~~~~~~~~~~~~~~~~~ \n");
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

void tablesDisplay(Tables *tables) {
  if (!tables || tables->num == 0) {
    fprintf(stderr, "empty tables\n");
    return;
  }
  fprintf(stderr, "tables cap : %d   \t  tables num : %d  \n", tables->cap, tables->num);
  for (int i = 0; i < tables->num; ++i) {
    fprintf(stderr, "\n tables %d,   tables id :  %s \n", i, tables->array[i].id);
    for (int j = 0; j < 8; ++j) {
      for (int k = 0; k < 8; ++k) {
        fprintf(stderr, "%d\t", tables->array[i].data[j * 8 + k]);
      }
      fprintf(stderr, "\n");
    }
  }
}

Tables* resizeTables(Tables *tables) {
//fprintf (stderr, "in Resize Table %d \n", tables);
  Tables *newTables = myalloc(1, Tables);
  newTables->array = myalloc(2 * tables->cap, TableNode);
  //memcpy(newTables->array, tables->array, tables->cap * sizeof(TableNode));
  for (int i = 0; i < tables->cap; ++i) {
  	newTables->array[i].id = tables->array[i].id;
  	newTables->array[i].data = tables->array[i].data;
  }
  newTables->cap = 2 * tables->cap;
  newTables->num = tables->cap;
  free(tables->array);
  free(tables);
  //tablesDestroy(tables);
fprintf (stderr, "return resize Table %d \n", newTables);
  return newTables;
}

void updateTable(int *het, int depth, uchar *seq, Tables **tables, int count) {

  int index = 0; 
  int start = het[5] - depth + 1;
//fprintf (stderr, "fun: update 1  %d ~~~~~~~~~~  \n", count);
//fprintf (stderr, "fun: update 1.5  %s ~~~~~~~~~~  \n", seq);
  for (int i = 0; i < 6; ++i) {
    index <<= 1;
    index += (seq[het[i] - start] - '0'); 
  }
  
//fprintf (stderr, "fun: update 2  index %d  ~~~~~~~~~~  \n", index);
//fprintf (stderr, "fun: update depth  %d~~~~~~~~~~  \n", depth);
  uchar *target = myalloc(depth - 5, uchar);
//fprintf (stderr, "fun: update depth  %d~~~~~~~~~~  \n", depth);
  for (int i = 0, j = 0; j < depth; ++j) {
    if ((j + start) == het[i]){
      i++;  
    } else {
      target[j - i] = seq[j]; 
    }
  }
  target[depth - 6] = '\0';

//fprintf (stderr, "fun: update 3~~~~~~~~~~  \n");
//fprintf (stderr, "fun: num : %d ~~~~~~~~~~  \n", tables->num);
//fprintf (stderr, "fun: %d ~~~~~~~~~~  \n", tables->array[11]);
//fprintf (stderr, "fun: success ~~~~~~~~~~  \n");

  int i;
  for (i = 0; i < (*tables)->num; ++i) {
//fprintf (stderr, "i  :  %d     tables->num  ;  %d\n ", i, tables->num);
//fprintf (stderr, "%s \n", target);
//fprintf (stderr, "%s \n", tables->array[i].id);
    if (!strcmp(target, (*tables)->array[i].id)) {
      break;
    }
  }

//fprintf (stderr, "fun: update 4~~~~~~~~~~  \n");
  if (i != (*tables)->num) {
    (*tables)->array[i].data[index] = count;
    free(target);
  } else {
    if ((*tables)->num == (*tables)->cap) {
fprintf (stderr, "resize!!!!!!!!!!!!!!!!!!!!!!~~~~~~~~~~~~~~~~~ \n");
fprintf (stderr, "begin resizeTable %d \t  %d \n", (*tables), (*tables)->cap);
       (*tables) = resizeTables((*tables));
fprintf (stderr, "after resizeTable %d \t  %d \n", (*tables), (*tables)->cap);
    }

    (*tables)->array[(*tables)->num].data = myalloc(64, int);
    memset((*tables)->array[(*tables)->num].data, 0, 64 * sizeof(int));
    (*tables)->array[(*tables)->num].data[index] = count;
    (*tables)->array[(*tables)->num].id = target;
    //tables->array[tables->num].id = myalloc(depth - 5, uchar);
    //memcpy(tables->array[tables->num].id, target, depth - 5);
    (*tables)->num++;
//fprintf (stderr, "after resizeTable %d, %d \n", tables, tables->num);
  }
//fprintf (stderr, "fun: update 6~~~~~~~~~~  \n");
}

void extendMatch(int *het, int cur, int num, int depth, uchar *seq, int *cc, int **u, int f, int g, Tables *tables) {
//fprintf (stderr, "fun:extend 0~~~~~~~~~~  \n");
//fprintf (stderr, "fun:extend start extend cur  %d   \t  table size  %d  \n", cur,  tables->num);
//fprintf (stderr, "seq in extend : %s \n", seq);
  if (f >= g) {
//fprintf (stderr, "fun:extend finish extend cur  %d   \t  table size  %d  \n", cur,  tables->num);
    return;
}
  //update the tables.
//fprintf (stderr, "fun:extend 1~~~~~~~~~~  \n");
  if (num == depth) {
//fprintf (stderr, "fun:extend 2~~~~~~~~~~  \n");
//fprintf (stderr, "begin updateTable %d, tables->num : %d \n", tables, tables->num);
fprintf (stderr, "1~~begin resizeTable %d \t  %d \n", tables, tables->cap);
    updateTable(het, depth, seq, &tables, g - f);
fprintf (stderr, "1~~after resizeTable %d \t  %d \n", tables, tables->cap);
//fprintf (stderr, "after updateTable %d, tables->num : %d  \n", tables, tables->num);
//fprintf (stderr, "fun:extend 3~~~~~~~~~~  \n");
//fprintf (stderr, "fun:extend update cur  %d   \t  table size  %d  \n", cur,  tables->num);
  } else {
//fprintf (stderr, "fun:extend 3~~~~~~~~~~  \n");
    seq[num] = '1';
//fprintf (stderr, "fun:extend 4~~~~~~~~~~  \n");
    extendMatch(het, cur + 1, num + 1, depth, seq, cc, u, cc[cur] + f - u[cur][f], cc[cur] + g - u[cur][g], tables); //1
//fprintf (stderr, "fun:extend 5~~~~~~~~~~  \n");
    seq[num] = '0';
//fprintf (stderr, "fun:extend 6~~~~~~~~~~  \n");
    extendMatch(het, cur + 1, num + 1, depth, seq, cc, u, u[cur][f], u[cur][g], tables); //0
//fprintf (stderr, "fun:extend finish extend cur  %d   \t  table size  %d  \n", cur,  tables->num);
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
  x = myalloc (N, uchar) ; 
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

//fprintf (stderr, "fun: main 1~~~~~~~~~~~~~~~~~ \n");
  int *pos;           /* record the heterozyogous position */
  pos = myalloc (N, int) ;
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

//fprintf (stderr, "fun: main  2~~~~~~~~~~~~~~~~~ \n");
  int start, depth;
  int *het = myalloc(6, int);


fprintf (stderr, "seg_num  %d \n", seg_num);
  for (s = 0; s < seg_num - 2; ++s) {
  //for (s = 0; s < 100; ++s) {
      Tables *tables = 0;
      uchar *seq;
      if (!s)
	start = 0;
      else
	start = pos[s * 3 - 1] + 1;
      
      for (i = 0; i < 6; ++i) {
        het[i] = pos[s * 3 + i];
      }

      depth = het[5] - start + 1;
      seq = myalloc(depth + 1, uchar);
      memset(seq, '2', depth * sizeof(uchar));
      seq[depth] = '\0';
      tables = tablesCreate(500);
//fprintf (stderr, "fun: main  5~~~~~~~~~~~~~~~~~ \n");
      extendMatch(het, start, 0, depth, seq, cc, u, 0, M, tables);
fprintf (stderr, "display  s = %d,  depth = %d \t table size = %d\n", s, depth, tables->num);
//tablesDisplay(tables);      
//fprintf (stderr, "fun: main  6~~~~~~~~~~~~~~~~~ \n");
      free(seq);
//fprintf (stderr, "fun: main  7~~~~~~~~~~~~~~~~~ \n");
      tablesDestroy(tables);
//fprintf (stderr, "fun: main  8~~~~~~~~~~~~~~~~~ \n");
  }
  free(het);

  /* cleanup */
  free(x); free(pos); free(cc);
  for (j = 0 ; j < N ; ++j) free(u[j]) ; free (u) ;
  for (j = 0 ; j < M ; ++j) free(reference[j]) ; free (reference) ;
  free(geno);
}
