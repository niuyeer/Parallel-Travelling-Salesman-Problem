#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <mpi.h>
#include <memory.h>
//#include <time.h>
#include "time.h"//clock() get time
MPI_Status status;

#define INFINITY ((u_short) 0xFFFF)

#define INIT_CAPACITY 64

typedef unsigned int u_int;
typedef unsigned short u_short;

void printResult();
int p,rank; //进程个个数、当前进程的编号,笛卡尔进程编号
void init();
void clear();

void babSearch(u_int, u_int);
void babSearchLeft(u_int, u_int, u_int);
void babSearchRight(u_int, u_int, u_int);

u_int estimate();
u_int estLeft(u_int row, u_int col);
u_int estRight(u_int row, u_int col);

static u_int n;
static u_short *mat = 0;

static u_int selectedPathLen = 0;

static u_int bound = (u_int)0xFFFFFFFF;

static u_int origin = 0;

static u_int *path;

static struct{
    u_int size;
    bool *row;
    bool *column;
} visited;

typedef struct{
    u_int row;
    u_int column;
} Node;

Node *newNode(u_int row, u_int column);
void deleteNode(Node *node);

inline Node *newNode(u_int row, u_int column){
    Node *node = malloc(sizeof(Node));
    node->row = row;
    node->column = column;
    return node;
}

inline void deleteNode(Node *node){
    free(node);
}

typedef struct{
    Node **nodes;
    u_int size;
    u_int capacity;
} Stack;

Stack stk;

void stack_init(Stack *s){
    s->size = 0;
    s->capacity = INIT_CAPACITY;
    s->nodes = malloc(sizeof(Node *) * INIT_CAPACITY);
}

void stack_destroy(Stack *s){
    free(s->nodes);
    s->nodes = 0;
    s->size = 0;
}

void stack_push(Stack *s, Node *n){
    if(s->size == s->capacity){
        s->capacity *= 2;
        Node **tmp = malloc(sizeof(Node *) * s->capacity);
        memcpy(tmp, s->nodes, sizeof(Node *) * s->size);
        free(s->nodes);
        s->nodes = tmp;
    }
    s->nodes[s->size++] = n;
}

Node* stack_pop(Stack *s){
    Node *ret = 0;
    if(s->size != 0){
        ret = s->nodes[--s->size];
        if(s->size == s->capacity / 3){
            s->capacity /= 2;
            Node **tmp = malloc(sizeof(Node *) * s->capacity);
            memcpy(tmp, s->nodes, sizeof(Node *) * s->size);
            free(s->nodes);
            s->nodes = tmp;
        }
    }
    return ret;
}


void babSearch(u_int currentValue, u_int start){
    if(visited.size == n-1){ //all nodes are visited
        u_int row_l;
        u_int col_l;
        for(row_l = 0; row_l < n; ++row_l){
            if(!visited.row[row_l]) break;
        }
        for(col_l = 0; col_l < n; ++col_l){
            if(!visited.column[col_l]) break;
        }
        stack_push(&stk, newNode(row_l, col_l));
        if(currentValue < bound){
            bound = currentValue;
            for(u_int i = 0; i < stk.size; ++i){
                //path[i] = stk.nodes[i]->row;
                path[i] = stk.nodes[i]->row + 1;
            }
        }
        deleteNode(stack_pop(&stk));
        return;
    }
    
    int ci = -1;
    u_short min = INFINITY;
    for(u_int j = 0; j < n; ++j){
        if(visited.column[j]) continue;
        if(mat[start*n+j] < min){
            min = mat[start*n+j];
            ci = j;
        }
    }
    
    if(ci == -1){ //WHY?!
        return;
    }
    
    u_int el = estLeft(start, ci);
    u_int er = estRight(start, ci);
    
    if(el <= er){ //goto left first
        if(el < bound){
            babSearchLeft(el, start, ci);
        }
        if(er < bound){
            babSearchRight(er, start, ci);
        }
    }else{ //goto right first
        if(er < bound){
            babSearchRight(er, start, ci);
        }
        if(el < bound){
            babSearchLeft(el, start, ci);
        }
    }
}

inline void babSearchLeft(u_int currentValue, u_int row, u_int col){
    stack_push(&stk, newNode(row, col));
    
    selectedPathLen += (u_int)mat[row*n+col];
    
    u_short raw = mat[col*n+col];
    mat[col*n+col] = INFINITY;
    visited.row[row] = true;
    visited.column[col] = true;
    visited.size++;
    
    babSearch(currentValue, col);
    
    mat[col*n+col] = raw;
    visited.row[row] = false;
    visited.column[col] = false;
    visited.size--;
    
    selectedPathLen -= (u_int)mat[row*n+col];
    
    deleteNode(stack_pop(&stk));
}

inline void babSearchRight(u_int currentValue, u_int row, u_int col){
    u_short raw = mat[row*n+col];
    mat[row*n+col] = INFINITY;
    
    babSearch(currentValue, row);
    
    mat[row*n+col] = raw;
}

inline u_int estLeft(u_int row, u_int col){
    selectedPathLen += (u_int)mat[row*n+col];
    
    visited.row[row] = true;
    visited.column[col] = true;
    visited.size++;
    
    u_int ret = estimate() + selectedPathLen;
    
    visited.row[row] = false;
    visited.column[col] = false;
    visited.size--;
    
    selectedPathLen -= (u_int)mat[row*n+col];
    
    return ret;
}

inline u_int estRight(u_int row, u_int col){
    u_int raw = mat[row*n+col];
    mat[row*n+col] = INFINITY;
    
    u_int ret = estimate() + selectedPathLen;
    
    mat[row*n+col] = raw;
    
    return ret;
}


static u_short *minEachRow = 0;
static bool *usedColumn = 0;

u_int estimate(){
    u_int ret = 0;
    memset(minEachRow, 0, sizeof(u_short) * n);
    memset(usedColumn, 0, sizeof(bool) * n);
    for(u_int i = 0; i < n; ++i){
        if(visited.row[i]) continue; //deleted row
        u_short min = INFINITY;
        u_int columnIndex = 0;
        for(u_int j = 0; j < n; ++j){
            if(visited.column[j]) continue; //deleted column
            if(mat[i*n+j] < min){
                min = mat[i*n+j];
                columnIndex = j;
            }
        }
        if(min == INFINITY){
            ret = INFINITY;
            goto L;
        }
        ret += min;
        minEachRow[i] = min;
        usedColumn[columnIndex] = true;
    }
    for(u_int j = 0; j < n; ++j){
        if(visited.column[j]) continue; //deleted column
        if(!usedColumn[j]){
            u_short min = INFINITY;
            for(u_int i = 0; i < n; ++i){
                if(visited.row[i]) continue; //deleted row
                if(mat[i*n+j] == INFINITY) continue;
                u_int r = (u_int)mat[i*n+j] - (u_int)minEachRow[i];
                if(r < min){
                    min = r;
                }
            }
            if(min == INFINITY){
                ret = INFINITY;
                goto L;
            }
            ret += min;
        }
    }
    
L:
    return ret;
}

void init(){
    
    printf("Scan the data...\n");
  //  scanf("%u", &n);
    n=12;
    visited.size = 0;
    visited.row = malloc(sizeof(bool) * n);
    visited.column = malloc(sizeof(bool) * n);
    for(u_int i = 0; i < n; ++i){
        visited.row[i] = false;
        visited.column[i] = false;
    }
    
    path = malloc(sizeof(u_int) * n);
    
    mat = malloc(sizeof(u_short) * n * n);
    u_int i = 0;
    FILE *fp;
    fp=fopen("./data1.txt","r");
    if(fp==NULL)
        printf("Fail to open the file！\n");
    else
    {
        //读入A,就像正常写入一样，把scanf换成fscanf
        int i=0;
        printf("Reading the matrix A:\n");
        while(i < n * n){
              fscanf(fp,"%u",&mat[i]);
           // scanf("%u", &mat[i]);
            
            if(i%n==i/n){
                mat[i] = INFINITY;
            }
            ++i;
        }
       
    }
        int j;
        for(i=0;i<n*n;i++)
        {
            if (i%n==0)
            { printf("\n");}
           printf("%3d ",mat[i]);
         
        }
     printf("\n");
    stack_init(&stk);
    
    minEachRow = malloc(sizeof(u_short) * n);
    usedColumn = malloc(sizeof(bool) * n);
}

void clear(){
    free(visited.row);
    free(visited.column);
    
    free(path);
    
    if(mat != 0){
        free(mat);
        mat = 0;
    }
    
    stack_destroy(&stk);
    
    free(minEachRow);
    free(usedColumn);
}

void printResult(){
    printf("%d\n", bound);
    for(u_int i = 0; i < n; ++i){
        printf("%d ", path[i]);
    }
    printf("%d\n", path[0]);
}

int main(int argc, char *argv[]){
    clock_t starttime,endtime;
   
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      starttime=clock();   //开始时间
    printf("My rank= %d\n",rank);
      init();
    if(rank==0)
    {
      
         MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }
    else
    {
         printf("My rank= %d n=%d\n ",rank,n);
       // MPI_Recv(a, np*np, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
       // MPI_Recv(b, np*np, MPI_INT, 0, 2, MPI_COMM_WORLD, &status);
    }
    babSearch(0, origin);
   // if(rank==0)
   // {
    printResult();
    clear();
  //  }
    endtime=clock();
    printf("Answer \n");
    printf("My rank= %d time=%f \n",rank,endtime-starttime);
    MPI_Finalize();
    // MPI_Comm_free(&comm);
    return 0;
}


