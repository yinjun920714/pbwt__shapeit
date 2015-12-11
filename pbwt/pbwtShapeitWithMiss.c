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
  newTree->num = tree->cap;
  free(tree->array);
  free(tree);
  return newTree;
}

Leaf* resizeLeaf(Leaf *leaf) {
  Leaf *newLeaf = myalloc(1, Leaf);
  newLeaf->array = myalloc(2 * leaf->cap, LeafNode);
  memcpy(newLeaf->array, leaf->array, leaf->cap * sizeof(LeafNode));
  newLeaf->cap = 2 * leaf->cap;
  newLeaf->num = leaf->cap;
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

//    memcpy(target, seq, depth * sizeof(uchar));
    target[depth - 6] = '\0';

    int i;
    TreeNode *cur;
    cur = (*tree)->array;
    for (i = 0; i < depth - 6; ++i) {
      if ((*tree)->num == (*tree)->cap) {
            (*tree) = resizeTree((*tree));
      }

      if (target[i] == '0') {
        //find the left node
        if (cur->left == 0) {
    //      if ((*tree)->num == (*tree)->cap) {
    //        (*tree) = resizeTree((*tree));
    //      }
          cur->left = (*tree)->num;
          (*tree)->num++;
        }
        cur = &((*tree)->array[cur->left]);
      } else {
        //find the right node
        if (cur->right == 0) {
    //        if ((*tree)->num == (*tree)->cap) {
    //          (*tree) = resizeTree((*tree));
    //      }
          cur->right = (*tree)->num;
          (*tree)->num++;
        }
        cur = &((*tree)->array[cur->right]);
      }  
    }

    if ((*leaf)->num == (*leaf)->cap) {
      (*leaf) = resizeLeaf((*leaf));
    }
    cur->left = (*leaf)->num;
    cur->right = cur->left;
    (*leaf)->num++;
    (*leaf)->array[cur->left].id = target;
    (*leaf)->array[cur->left].count = count;
    (*leaf)->array[cur->left].state = index;
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

double* calculateConditionalTable(uchar *homSeq, Leaf **leaf, double epsilon) {
//  fprintf (stderr, "homSeq : %s \t leaf_size : %d \n", homSeq, (*leaf)->num);
  double *table = myalloc(64, double);
  for (int i = 0; i < 64; ++i)
	table[i] = 0.0;
  
  for (int i = 0; i < (*leaf)->num; ++i) {
    int count = 0;
    int j = 0;
    while (homSeq[j] != '\0') {
    	if (homSeq[j] != (*leaf)->array[i].id[j])
	   count++;
        if (count > 3)
	   break;
        j++;
    } 
    if (count > 3)
	continue;
    double factory = 1.0;
    while(count > 0) {
	factory *= epsilon;
        count--;
    }
    table[(*leaf)->array[i].state] += (*leaf)->array[i].count * factory; 
  }
  return table;
}

void displayTable(double *table) {
  for (int i = 0; i < 8; ++i) {
  	for (int j = 0; j < 8; ++j)
  		fprintf (stderr, "%f \t", table[i * 8 + j]);
  	fprintf (stderr, "\n");
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
  int time = 20;
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
	//debug only
	double* conditionTable;
    	uchar *homSeq = myalloc(depth - 5, uchar);
    	for (int ii = 0, jj = 0; jj < depth; ++jj) {
		if ((jj + start) == het[ii]){
       			ii++;  
            	} else {
       			homSeq[jj - ii] = (geno[t][jj + start] == 2) ? '1' : '0'; 
      	    	}	
	}	
	homSeq[depth - 6] = '\0';
	conditionTable = calculateConditionalTable(homSeq, &leaf, 0.1);
	if (s == 1000) {
		fprintf (stderr, "homSeq : %s \t leaf_size : %d \n", homSeq, leaf->num);
		displayTable(conditionTable);
	}
	free(homSeq);
	free(conditionTable);
     //tablesDisplay(tables);      
     //leafDisplay(leaf);      
        free(seq);
        treeDestroy(tree);
        leafDestroy(leaf);
    }
    free(het);
  } 
  /* cleanup */
  free(x); free(pos); free(cc);
  for (j = 0 ; j < N ; ++j) free(u[j]) ; free (u) ;
  for (j = 0 ; j < time; ++j) free(geno[j]) ; free (geno);
}
