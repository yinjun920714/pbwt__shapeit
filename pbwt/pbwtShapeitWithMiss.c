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

//for the parent of leaf node, the left will present the index in LeafStruct;
typedef struct TreeNodeStruct {
  int left;      // '0'
  int right;     // '1'
} TreeNode ;

typedef struct TreeStruct {
  TreeNode* array;
  int cap;
  int num;
} Tree ;

typedef struct LeafNodeStruct {
  uchar* id;
  int state;
  int count;
} LeafNode ;

typedef struct LeafStruct {
  LeafNode* array;
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

Tree* treeCreate(int cap) {
  Tree* tree = 0;
  tree = myalloc(1, Tree);
  tree->array = myalloc(cap, TreeNode);
  memset(tree->array, 0, cap * sizeof(TreeNode));
  tree->cap = cap;
  tree->num = 1;
  return tree;
}

Leaf* leafCreate(int cap) {
  Leaf* leaf = 0;
  leaf = myalloc(1, Leaf);
  leaf->array = myalloc(cap, LeafNode);
  leaf->cap = cap;
  leaf->num = 0;
  return leaf;
}

void tablesDestroy(Tables *tables) {

  for (int i = 0; i < tables->num; ++i) {
    free(tables->array[i].id);
    free(tables->array[i].data);
  }

  free(tables->array);
  free(tables);
}

void treeDestroy(Tree *tree) {
  free(tree->array);
  free(tree);
}

void leafDestroy(Leaf *leaf) {

  for (int i = 0; i < leaf->num; ++i) {
    free(leaf->array[i].id);
  }

  free(leaf->array);
  free(leaf);
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

void leafDisplay(Leaf *leaf) {
  if (!leaf || leaf->num == 0) {
    fprintf(stderr, "empty leaf\n");
    return;
  }
  fprintf(stderr, "leaf cap : %d   \t  leaf num : %d  \n", leaf->cap, leaf->num);
  for (int i = 0; i < leaf->num; ++i) {
    fprintf(stderr, "\n leaf %d, %s \t  count : %d \t state : %d\n", i, leaf->array[i].id, leaf->array[i].count, leaf->array[i].state);
  }
}

Tables* resizeTables(Tables *tables) {
  Tables *newTables = myalloc(1, Tables);
  newTables->array = myalloc(2 * tables->cap, TableNode);
  /*
  for (int i = 0; i < tables->cap; ++i) {
  	newTables->array[i].id = tables->array[i].id;
  	newTables->array[i].data = tables->array[i].data;
  }
  */
  memcpy(newTables->array, tables->array, tables->cap * sizeof(TableNode));
  newTables->cap = 2 * tables->cap;
  newTables->num = tables->cap;
  free(tables->array);
  free(tables);
  return newTables;
}

Tree* resizeTree(Tree *tree) {
  Tree *newTree = myalloc(1, Tree);
  newTree->array = myalloc(2 * tree->cap, TreeNode);
  memset(newTree->array, 0, 2 * tree->cap * sizeof(TreeNode));
  memcpy(newTree->array, tree->array, tree->cap * sizeof(TreeNode));
  newTree->cap = 2 * tree->cap;
  newTree->num = tree->num;
  free(tree->array);
  free(tree);
  return newTree;
}

Leaf* resizeLeaf(Leaf *leaf) {
  Leaf *newLeaf = myalloc(1, Leaf);
  newLeaf->array = myalloc(2 * leaf->cap, LeafNode);
  memcpy(newLeaf->array, leaf->array, leaf->cap * sizeof(LeafNode));
  newLeaf->cap = 2 * leaf->cap;
  newLeaf->num = leaf->num;
  free(leaf->array);
  free(leaf);
  return newLeaf;
}

/*
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
*/

void insertTree(int *het, int depth, uchar *seq, Tree **tree, Leaf **leaf, int count) {
    int index = 0; 
    int start = het[5] - depth + 1;
    for (int i = 0; i < 6; ++i) {
      index <<= 1;
      index += (seq[het[i] - start] - '0'); 
    }
    
    // uchar *target = myalloc(depth + 1, uchar);
    uchar *target = myalloc(depth - 5, uchar);
 
    for (int i = 0, j = 0; j < depth; ++j) {
      if ((j + start) == het[i]){
        i++;  
      } else {
        target[j - i] = seq[j]; 
      }
    }

    //memcpy(target, seq, depth * sizeof(uchar));
    target[depth - 6] = '\0';

    int i;
    int cur = 0;
    /*
    for (i = 0; i < depth - 6; ++i) {
      if ((*tree)->num == (*tree)->cap) {
            (*tree) = resizeTree((*tree));
      }

      if (target[i] == '0') {
        //find the left node
        if ((*tree)->array[cur].left == 0) {
          (*tree)->array[cur].left = (*tree)->num;
          (*tree)->num++;
        }
        if ((*tree)->array[cur].left >= (*tree)->cap) fprintf(stderr, "out of scope 1 \t %d \t %d\n", 1, (*tree)->cap);
        cur = (*tree)->array[cur].left;
        } else {
        //find the right node
        if ((*tree)->array[cur].right == 0) {
          (*tree)->array[cur].right = (*tree)->num;
          (*tree)->num++;
        }
        if ((*tree)->array[cur].right>= (*tree)->cap) fprintf(stderr, "out of scope 2 \t %d \t %d \n", 0, (*tree)->cap);
        cur = (*tree)->array[cur].right;
      }  
    }
    */

    if ((*leaf)->num == (*leaf)->cap) {
      (*leaf) = resizeLeaf((*leaf));
    }
    /*
    (*tree)->array[cur].left = (*leaf)->num;
    (*tree)->array[cur].right = (*leaf)->num;
    (*leaf)->num++;
    (*leaf)->array[(*tree)->array[cur].left].id = target;
    (*leaf)->array[(*tree)->array[cur].left].count = count;
    (*leaf)->array[(*tree)->array[cur].left].state = index;
    */
    (*leaf)->array[(*leaf)->num].id = target;
    (*leaf)->array[(*leaf)->num].count = count;
    (*leaf)->array[(*leaf)->num].state = index;
    (*leaf)->num++;
}

//store in tree and leaf
void extendMatch(int *het, int cur, int num, int depth, uchar *seq, int *cc, int **u, int f, int g, Tree **tree, Leaf **leaf) {

  if (f >= g) {
    return;
  } 
  if (num == depth) {
    insertTree(het, depth, seq, tree, leaf, g - f);
  } else {
    seq[num] = '1';
    extendMatch(het, cur + 1, num + 1, depth, seq, cc, u, cc[cur] + f - u[cur][f], cc[cur] + g - u[cur][g], tree, leaf); //1
    seq[num] = '0';
    extendMatch(het, cur + 1, num + 1, depth, seq, cc, u, u[cur][f], u[cur][g], tree, leaf); //0
  }
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

void calculateConditionalTable(uchar *homSeq, Leaf **leaf, double** conditionalTable, int s, double epsilon) {
  //  fprintf (stderr, "homSeq : %s \t leaf_size : %d \n", homSeq, (*leaf)->num);
  for (int i = 0; i < 64; ++i)
	conditionalTable[i][s] = 0.0;  
  
  for (int i = 0; i < (*leaf)->num; ++i) {
    int count = 0;
    int j = 0;
    while (homSeq[j] != '\0') {
    	if (homSeq[j] != (*leaf)->array[i].id[j])
	       count++;
      j++;
    } 
    double factory = 1.0;
    while(count > 0) {
	factory *= epsilon;
        count--;
    }

    conditionalTable[(*leaf)->array[i].state][s] += ((*leaf)->array[i].count * factory); 
  }
}

void viterbiSamplingWithMiss(int** g1, int** f1, double** conditionTables, int *pos, int seg_num, uchar *shape1, uchar *shape2, double w) {
  int i, j, s;
  int **phis;
  double **data; //condition probability
  data = myalloc(8, double* );
  for ( i = 0; i < 8; ++i ) {
    data[i] = myalloc(seg_num, double); }
  phis = myalloc(8, int* );
  for ( i = 0; i < 8; ++i ) {
    phis[i] = myalloc(seg_num, int); }
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
      double maxVal = data[0][s] * (conditionTables[i][s] + 1.0/8) / (totalCon[0] + 1) 
                      * data[7][s] * (conditionTables[63 - i][s] + 1.0/8) / (totalCon[7] + 1);
      for ( j = 1; j < 8; ++j) {
        double val = data[j][s] * (conditionTables[j * 8 + i][s] + 1.0/8) / (totalCon[j] + 1)
                      * data[7 - j][s] * (conditionTables[(7 - j) * 8 + 7 - i][s] + 1.0/8) / (totalCon[7 - j] + 1);
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

void pbwtShapeItWithMiss (PBWT *p, int thousandths, FILE *out) {
 
  if (!p || !p->yz) die ("option -longWithin called without a PBWT") ;
  
  /********   ref  *****/
  uchar **reference = pbwtHaplotypes (p) ; /* haplotypes for reference  (M * N)  */
  uchar **origin;
  /*********************/

  uchar *x;                 /* use for current query */
  PbwtCursor *up = pbwtCursorCreate (p, TRUE, TRUE) ;
  int **u ;   /* stored indexes */
  int **f1, **g1 ;      /* one block match, start in index f1 and end in index g1 */
  double **conditionTable;
  int i, j, k, N = p->N, M = p->M ;
  int num_1;      /* for the num of heterozyogous */
  int s, seg_num; /* for the segment number and current segment */
  uchar *shape1;      /* for the shapeit seq1 */
  uchar *shape2;      /* for the shapeit seq2 */
  uchar *homSeq;
  double w = 1.0 / M; /* the small weight for the add weight */
  double epsilon = (double)thousandths/1000;

  /* build indexes */
  u = myalloc (N,int*) ; for (i = 0 ; i < N ; ++i) u[i] = myalloc (p->M+1, int) ;
  x = myalloc (N, uchar) ; 
  int *cc = myalloc (p->N, int) ;
  f1 = myalloc (8, int*);
  g1 = myalloc (8, int*);
  conditionTable = myalloc (64, double*);
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
  //clean up
  pbwtCursorDestroy (up) ;
  origin = myalloc (2, uchar*); for (i = 0; i < 2; ++i) origin[i] = myalloc (p->N, uchar);
  fprintf (stderr, "Made indices: \n") ; timeUpdate ();

  int time = M/2;
  //  int time = 5;
  int *pos;           /* record the heterozyogous position */
  pos = myalloc (N, int) ;
  for (int t = 0; t < time; ++t) {
    num_1 = 0;      /* for the num of heterozyogous */
    s = 0;
    seg_num = 1; /* for the segment number and current segment */
    
    memcpy (origin[0], reference[2*t], N*sizeof(uchar));
    memcpy (origin[1], reference[2*t + 1], N*sizeof(uchar));
    
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
    
    int depth;
    int *het = myalloc(6, int);
    for ( i = 0; i < 8; ++i) { 
      f1[i] = myalloc(seg_num, int);
      g1[i] = myalloc(seg_num, int);
    }
    for ( i = 0; i < 64; ++i) { 
      conditionTable[i] = myalloc(seg_num, double);
    }

    /* initial f1,g1 */
    for (i = 0; i < 8; ++i) {
      for (j = 0; j < seg_num; ++j) {
        f1[i][j] = 0; g1[i][j] = M;
      }
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


    //fprintf (stderr, "seg_num  %d \n", seg_num);

    for (s = 0; s < seg_num - 2; ++s) {
        Tree *tree = 0;
        Leaf *leaf = 0;
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
        tree = treeCreate(5000);
        leaf = leafCreate(500);
        extendMatch(het, start, 0, depth, seq, cc, u, 0, M, &tree, &leaf);
      	homSeq = myalloc(depth - 5, uchar);
      	for (int ii = 0, jj = 0; jj < depth; ++jj) {
		if ((jj + start) == het[ii]){
       			 ii++;  
          	} else {
			homSeq[jj - ii] = x[jj + start] == 0 ? '0' : '1'; 
      	  	}		
	}
	homSeq[depth - 6] = '\0';
	calculateConditionalTable(homSeq, &leaf, conditionTable, s, epsilon);
	
	free(homSeq);
        //tablesDisplay(tables);      
        //leafDisplay(leaf);      
        
        free(seq);
        treeDestroy(tree);
        leafDestroy(leaf);
    }
    viterbiSamplingWithMiss(g1, f1, conditionTable, pos, seg_num, shape1, shape2, w);
    memcpy (reference[2 * t], shape1, N*sizeof(uchar));
    memcpy (reference[2 * t + 1], shape2, N*sizeof(uchar));
    
    free(het);
    /* cleanup */
    for ( i = 0; i < 8; ++i)  { free(f1[i]); free(g1[i]); }
/*
    //debug    
    for( s = 1000; s < 1010; ++s) {
    	for (int loop = 0; loop < 64; ++loop) {
      		fprintf(out, "%f ", conditionTable[loop][s]);
    	}
	fprintf(out, "\n");		
    }
*/
    for ( i = 0; i < 64; ++i) { free(conditionTable[i]); }

  }
  //output the shapeit result.
  for ( j = 0; j < N; ++j) {
    for ( i = 0; i < M; ++i)
      fprintf(out, "%u ", reference[i][j]);
    fprintf(out, "\n");
  }
  fclose(out);

  /* cleanup */
  free(x); free(pos); free(cc);
  for (j = 0 ; j < N ; ++j) free(u[j]) ; free (u) ;
  free (shape1) ; free (shape2) ;
  free(f1); free(g1); free(conditionTable);
  for (j = 0 ; j < M ; ++j) free(reference[j]) ; free (reference) ;
}
