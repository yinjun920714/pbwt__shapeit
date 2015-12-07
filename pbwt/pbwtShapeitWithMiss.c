#include "pbwt.h"
typedef struct TableNodeStruct {
  uchar* id;
  int* data;
} TableNode ;

typedef struct TablesStruct {
  TableNode* array;
  int cap;
  int num;
} Tables ;

typedef struct TreeNodeStruct {
  int left;      // '0'
  int right;     // '1'
} TreeNode ;

typedef struct TreeStruct {
  TreeNode* root;
  int cap;
  int num;
} Tree ;

typedef struct LeafNodeStruct {
  uchar* id;
  int count;
} LeafNode ;

typedef struct LeafStruct {
  LeafNode* head;
  int cap;
  int num;
} Leaf ;



Tables* tablesCreate(int cap) {
  Tables* tables = 0;
  tables = myalloc(1, Tables);
  tables->array = myalloc(cap, TableNode);
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
  Tables *newTables = myalloc(1, Tables);
  newTables->array = myalloc(2 * tables->cap, TableNode);
  for (int i = 0; i < tables->cap; ++i) {
  	newTables->array[i].id = tables->array[i].id;
  	newTables->array[i].data = tables->array[i].data;
  }
  newTables->cap = 2 * tables->cap;
  newTables->num = tables->cap;
  free(tables->array);
  free(tables);
  return newTables;
}

void updateTable(int *het, int depth, uchar *seq, Tables **tables, int count) {

  int index = 0; 
  int start = het[5] - depth + 1;
  for (int i = 0; i < 6; ++i) {
    index <<= 1;
    index += (seq[het[i] - start] - '0'); 
  }
  
  uchar *target = myalloc(depth - 5, uchar);
  for (int i = 0, j = 0; j < depth; ++j) {
    if ((j + start) == het[i]){
      i++;  
    } else {
      target[j - i] = seq[j]; 
    }
  }
  target[depth - 6] = '\0';

  int i;
  for (i = 0; i < (*tables)->num; ++i) {
    if (!strcmp(target, (*tables)->array[i].id)) {
      break;
    }
  }

  if (i != (*tables)->num) {
    (*tables)->array[i].data[index] = count;
    free(target);
  } else {
    if ((*tables)->num == (*tables)->cap) {
       (*tables) = resizeTables((*tables));
    }

    (*tables)->array[(*tables)->num].data = myalloc(64, int);
    memset((*tables)->array[(*tables)->num].data, 0, 64 * sizeof(int));
    (*tables)->array[(*tables)->num].data[index] = count;
    (*tables)->array[(*tables)->num].id = target;
    (*tables)->num++;
  }
}

void extendMatch(int *het, int cur, int num, int depth, uchar *seq, int *cc, int **u, int f, int g, Tables **tables) {
  if (f >= g) {
    return;
}
  if (num == depth) {
    updateTable(het, depth, seq, tables, g - f);
  } else {
    seq[num] = '1';
    extendMatch(het, cur + 1, num + 1, depth, seq, cc, u, cc[cur] + f - u[cur][f], cc[cur] + g - u[cur][g], tables); //1
    seq[num] = '0';
    extendMatch(het, cur + 1, num + 1, depth, seq, cc, u, u[cur][f], u[cur][g], tables); //0
  }
}

void pbwtShapeItWithMiss (PBWT *p, FILE *out) {
 
  if (!p || !p->yz) die ("option -longWithin called without a PBWT") ;
  
  /********   ref  *****/
  uchar **reference = pbwtHaplotypes (p) ; /* haplotypes for reference  (M * N)  */

  /*********************/

  uchar *x;                 /* use for current query */
  PbwtCursor *up = pbwtCursorCreate (p, TRUE, TRUE) ;
  int **u ;   /* stored indexes */
  int i, j, k, N = p->N, M = p->M ;
  int num_1;      /* for the num of heterozyogous */
  int s, seg_num; /* for the segment number and current segment */
  
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
  int time = 150;
  int **geno;
  geno = myalloc(time, int*);
  for (i = 0; i < time; ++i) geno[i] = myalloc (p->N, int);
  
  for (j = 0; j < time; ++j) {
     for (i = 0; i < p->N; ++i) {
        geno[j][i] = reference[j * 2][i] + reference[j * 2 + 1][i];
     }
  }

  //clean up
  pbwtCursorDestroy (up) ;
  for (j = 0 ; j < M ; ++j) free(reference[j]) ; free (reference) ;

  fprintf (stderr, "Made indices: \n") ; timeUpdate ();

  int *pos;           /* record the heterozyogous position */
  pos = myalloc (N, int) ;
for (int t = 0; t < time; ++t) {
  
  num_1 = 0;      /* for the num of heterozyogous */
  s = 0;
  seg_num = 1; /* for the segment number and current segment */
  
  /* find the heterozyogous position and record */
  for ( i = 0, j = 0; i < N; ++i) {
    if (geno[t][i] == 1) {
      ++j;
      pos[num_1++] = i;
      if (j == 3) {
        ++seg_num;
        j = 0;
      }
    }
  }

  int start, depth;
  int *het = myalloc(6, int);


fprintf (stderr, "seg_num  %d \n", seg_num);
  for (s = 0; s < seg_num - 2; ++s) {
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
      extendMatch(het, start, 0, depth, seq, cc, u, 0, M, &tables);
//fprintf (stderr, "display  s = %d,  depth = %d \t table size = %d\n", s, depth, tables->num);
//tablesDisplay(tables);      
      free(seq);
      tablesDestroy(tables);
  }
  free(het);
}
  fprintf (stderr, "finished \n") ; timeUpdate ();
  /* cleanup */
  free(x); free(pos); free(cc);
  for (j = 0 ; j < N ; ++j) free(u[j]) ; free (u) ;
  for (j = 0 ; j < time; ++j) free(geno[j]) ; free (geno);
}
