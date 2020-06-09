#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>

#define DIM 101       // 8 liver,  101 circle, 5 iris,   10  glass,
#define NUM 200     // 345 liver, 200 circle, 150 iris, 214 glass,
#define STEP 500
#define SAMPLE 1
#define NL 100  //  1000/10

//===========================

/* Period parameters */  
#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

static unsigned long mt[N]; /* the array for the state vector  */
static int mti=N+1; /* mti==N+1 means mt[N] is not initialized */

void init_genrand(unsigned long s);
void init_by_array(unsigned long init_key[], int key_length);
unsigned long genrand_int32(void);
long genrand_int31(void);
double genrand_real1(void);
double genrand_real2(void);
double genrand_real3(void);
double genrand_real53(void);

//=============================================

typedef struct {
  double x;
  double y;
} Mvector;

typedef struct {
  double x[DIM];
} Dvector;

typedef struct {
  Mvector r;
  Mvector v;
  int nx, ny;
  Dvector o;
} Datoid;

Datoid initDatoid(double obj[], double w, double h, double ls); 
void normalization(void);
void similarDist(void);
double euclidDist( Dvector o1, Dvector o2 );
void distance(void);
void distance_1(int n);
double mapDist( Mvector m1, Mvector m2 );
void shaffle(int mask[], int m );
int update( int n );
Mvector sim_nn_K( int n );
int dis_min(int n);
int sim_min(int n);
void printState(void);
void shiftMask(int l, int i);
double analysis(void);
double G_K_factor(void);
int avoidOverlap(int ix, int iy,int n, double ls);

double width = 1000.0;
double height = 1000.0;
double lattice_constant = 10.0;
double V_max = 100.0;

double w_separation=0.5;
double w_alignment = 0.9;
double w_cohesion=0.9;
double w;
int K_near;

double dmax;
double Td = 30.0;
double Td2 = 22.0;  // 0608

double obj[NUM][DIM];
Datoid datoid[NUM];
double similar[NUM][NUM];
double sd[NUM][NUM];
double dd[NUM][NUM];
double dist[NUM][NUM];
int mask[NUM];
int place[NL][NL];
Mvector null_vec = {0.0, 0.0};
double measure_sim=0.0;
double measure_GK=0.0;
double ms[STEP/100+1];
double mgk[STEP/100+1];
int countD;

FILE *fp, *fp2, *fp3;

int main(void) {
  int i_c, i_s;
  int i, j, dim, n, l;
  int samp, step;
  int ix, iy;
  int pre_x, pre_y;
  int count;
  int countDD;
  int k, m;
  int flag;
  int ll, nn;

  dmax = sqrt( width*width + height*height );
  //  printf("dmax = %lf\n", dmax);

  //  init_genrand((unsigned)time(NULL));
  //  init_genrand(131);

  fp = fopen("./circle200.txt", "r");
  //  if( fp == NULL) printf("file error\n");
  fp2 = fopen("./circle_w_9_5_2T.csv", "w");
  fp3 = fopen("./circle_w_9_5_2Tmap.csv", "w");

  for(i=0; i<NUM; i++) {
    mask[i] = i;
    for(dim=0; dim<DIM; dim++) {
      fscanf(fp,"%lf", &obj[i][dim]);
      //      printf("%lf, ", obj[i][dim]);
    }
    //    printf("\n");
  }
  normalization();

  init_genrand(1311);
  for(i=0; i<NUM; i++)  mask[i] = i;
  measure_sim=0.0;
  measure_GK=0.0;

  for(i=0; i<=STEP/100; i++) {
    ms[i] = 0.0;
  }
  for(i=0; i<=STEP/100; i++) {
    mgk[i] = 0.0;
  }
  for( samp = 0; samp<SAMPLE; samp++) {
    //    printf("sample = %d\n", samp);
    for( i=0; i<NL; i++ ) {
      for( j=0; j<NL; j++ ) {
	place[i][j] = -1;
      }
    }

    for(i=0; i<NUM; i++) {
      datoid[i] = initDatoid( obj[i], width, height, lattice_constant );
      ix = datoid[i].nx;
      iy = datoid[i].ny;
      place[ix][iy] = i;
    }
	
    similarDist();
    for(i=0; i<NUM; i++) {
      for(j=0; j<NUM; j++) {
	if( similar[i][j] < 0.0) printf("similar is negative  ");
	else if( similar[i][j] > 1.0) printf("similar is larger than one  ");
	//	    printf("%lf\n", similar[i][j]);
      }
    }
    distance();
    
    K_near = (int)(NUM/4.0);
    w = 0.7;
    
    for( step = 0; step <= STEP; step++ ) {
      //      printf("%d step\n", step);
      shaffle( mask, NUM );

      count = 0;
      n = mask[count];
      nn = n;
      pre_x = datoid[n].nx;
      pre_y = datoid[n].ny;
      while( count<NUM ) {
	l = update( n );

	if( l >= 0 && l != n ) shiftMask( l, count);

	ix = datoid[n].nx;
	iy = datoid[n].ny;
	place[ix][iy] = n;
	if( count >= NUM-1 ) {
	  if( !( l==n || l<0 )  ) {
	    place[pre_x][pre_y] = l;
	    datoid[l].nx = pre_x;
	    datoid[l].ny = pre_y;
	    datoid[l].r.x = pre_x*lattice_constant;
	    datoid[l].r.y = pre_y*lattice_constant;
	  }
	  else if( l<0 ) {
	    place[pre_x][pre_y] = -1;
	  }
	  else if( l ==  n) {
	    if( !( datoid[n].nx == pre_x && datoid[n].ny == pre_y) ) {
	      place[pre_x][pre_y] = -1;
	    }
	  }
	}
	else {
	  if( l == n || l<0 ) {
	    if( l<0 ) {
	      place[pre_x][pre_y] = -1;
	    }
	    else if( l ==  n) {
	      if( !( datoid[n].nx == pre_x && datoid[n].ny == pre_y) ) {
		place[pre_x][pre_y] = -1;
	      }
	    }
	    n = mask[count+1];
	    pre_x = datoid[n].nx;
	    pre_y = datoid[n].ny;
	  }
	  else {
	    if( !(datoid[n].nx == pre_x && datoid[n].ny == pre_y ) ) {
	      n = l;
	    }
	  }
	}
	count++;	
      }
      if( w > 0.2 )  w *= 0.999999;
      if(K_near > 15 && step%100==0 ) K_near--; 
      

      if( step%100 == 0 ) {
	//	    printState();
	//	    fprintf(fp2, "%d, %e\n", step, analysis() );
	ms[step/100] += analysis();
      }

	  
      if(step%1000 == 0 ) {
	mgk[step/1000] += G_K_factor();
	    //	    measure_GK += G_K_factor();
	    //	    measure_sim = analysis();
	    //	    measure_GK = G_K_factor();
	    //	    fprintf(fp2, "%d %e\n", step, measure_sim);
      } 

      /*	  
      if( step < (int)( 0.2*STEP) ) {
	w_alignment = 0.9 - 0.4*step/(0.2*STEP);
	w_cohesion = 0.9 - 0.4*step/(0.2*STEP);
      }
      else {
	w_alignment = 0.5;
	w_cohesion = 0.5;
      }
      */  
      /*
      countD=0;
      for( k=0; k<NL; k++ ) {
	for( m=0; m<NL; m++ ) {
	  if( place[k][m]>=0 ) countD++; 
	}
      }
      printf("countD = %d\n", countD);
      */

    }// step
	//	measure_GK += G_K_factor();
    printState();
  }//sample loop
      //      measure_sim += analysis();
      //      printf("analysis done\n");
      //      measure_GK += G_K_factor();
      //      measure_sim /= (double)SAMPLE;
      //      measure_GK /= (double)SAMPLE;
      //      fprintf(fp2, "%e %e %e %e\n", w_cohesion, w_separation, measure_sim, measure_GK);
      
  for(i=0; i<=STEP/100; i++) {
    ms[i] /= (double)SAMPLE;
    //	fprintf(fp2, "%d %e\n", w_alignment, w_cohesion, w_separation, i*100, ms[i]);
    fprintf(fp2, "%d %e\n", i*100, ms[i]);
  }
  for(i=0; i<=STEP/1000; i++) {
    mgk[i] /= (double)SAMPLE;
    fprintf(fp2, "%d %e\n", i*1000, mgk[i]);
  }
      
      //      fprintf(fp2, "GK = %e\n", measure_GK/(double)SAMPLE );
      //      printState();
      //      printf("%e %e %e %e\n", w_cohesion, w_separation, measure_sim, measure_GK);
      //    }
    //    fprintf(fp2,"\n");
      //  }

  fclose(fp);
  fclose(fp2);
  fclose(fp3);

  return 0;
}

void shiftMask(int l, int m) {
  int i, tmp;
  for(i=0; i<NUM; i++) {
    if(mask[i] == l && i > m ) {
      tmp = mask[m+1];
      mask[m+1] = l;
      mask[i] = tmp;
    }
  }
}

Datoid initDatoid(double obj[], double w, double h, double ls ) {
  Datoid d;
  int dim;

  do {
    d.nx = (int)( (w/ls - 4.0)*genrand_real1() ) + 2;
    d.ny = (int)( (w/ls - 4.0)*genrand_real1() ) + 2;
  } while( place[d.nx][d.ny] >= 0 );

  d.r.x = ls*d.nx;
  d.r.y = ls*d.ny;

  d.v.x = 2.0*V_max*( genrand_real1() - 0.5 );
  d.v.x = 2.0*V_max*( genrand_real1() - 0.5 );

  for(dim=0; dim<DIM; dim++) {
    d.o.x[dim] = obj[dim];
  }
  return d;
}

void normalization(void) {
  double max[DIM-1];
  double min[DIM-1];
  int dim, i;

  for(dim=0; dim<DIM-1; dim++) {
    max[dim] = obj[0][dim]; 
    min[dim] = obj[0][dim];
  }
  for(i=0; i<NUM; i++) {
    for(dim=0; dim<DIM-1; dim++) {
      if( obj[i][dim] > max[dim] ) max[dim] = obj[i][dim];
      if( obj[i][dim] < min[dim] ) min[dim] = obj[i][dim];
    }
  }
  for(i=0; i<NUM; i++) {
    for(dim=0; dim<DIM-1; dim++) {
      if(min[dim] != max[dim]) {
        obj[i][dim] = ( obj[i][dim] - min[dim] )/( max[dim] - min[dim] );
      }
      else {
	obj[i][dim] = obj[i][dim]/max[dim];
      }
    }
  }
}

void similarDist(void) {
  int i, j;
  double max=0.0;

  for(i=0; i<NUM; i++) {
    similar[i][i] = 1.0;
    for(j=i+1; j<NUM; j++) {
      similar[i][j] = euclidDist( datoid[i].o, datoid[j].o );
      //     similar[i][j] = hammingDist( datoid[i].o, datoid[j].o );
      similar[j][i] = similar[i][j];
      if( max < similar[i][j] ) max = similar[i][j];
    }
  }
  for(i=0; i<NUM; i++) {
    for(j=i+1; j<NUM; j++) {
      //      similar[i][j] = 1.0 - similar[i][j]/sqrt(DIM-1.0); 
      //      similar[j][i] = 1.0 - similar[j][i]/sqrt(DIM-1.0);
      similar[i][j] = 1.0 - similar[i][j]/max;
      similar[j][i] = 1.0 - similar[j][i]/max;
    }
  }
}

double euclidDist( Dvector o1, Dvector o2 ) {
  int dim;
  double d = 0.0;

  for(dim=0; dim<DIM-1; dim++) {
    d += ( o1.x[dim] - o2.x[dim])*( o1.x[dim] - o2.x[dim]);
  }
  d = sqrt(d);

  return d;
}

void distance(void) {
  int i, j;
  
  for(i=0; i<NUM; i++) {
    dist[i][i] = 0.0;
    for(j=i+1; j<NUM; j++) {
      dist[i][j] = mapDist( datoid[i].r, datoid[j].r );
      dist[j][i] = dist[i][j];
    }
  }
  for(i=0; i<NUM; i++) {
    for(j=0; j<NUM; j++) {
      sd[i][j] = similar[i][j]*dist[i][j] + (1 - similar[i][j])*dmax;
      dd[i][j] = ( 1 - similar[i][j])*dist[i][j] + similar[i][j]*dmax;
    }
  }
}

double mapDist( Mvector m1, Mvector m2) {
  int dim;
  double d = 0.0;

  d = ( m1.x - m2.x )*( m1.x - m2.x )
    + ( m1.y - m2.y )*( m1.y - m2.y );
  
  d = sqrt(d);

  return d;
}

void shaffle(int mask[], int mm) {
  int i, m;
  int tmp;

  for(i=0; i<mm; i++) {
    m = (int)( NUM*genrand_real3() );
    tmp = mask[i];
    mask[i] = mask[m];
    mask[m] = tmp;
  }
}

int update(int n) {
  int n_sim;
  int n_dis;
  int l = -1;
  double X, Y;
  //  double preX, preY;
  Mvector V_a, V_c, V_r, C_sim, V_r2 = {0, 0};
  int nx, ny, i;
  int Vcount=0;

  n_sim = sim_min(n);
  if( n_sim < 0 || n_sim>=NUM) printf("n_sim out of range\n");
  V_a.x = w_alignment*( datoid[n_sim].v.x - datoid[n].v.x );
  V_a.y = w_alignment*( datoid[n_sim].v.y - datoid[n].v.y );

  C_sim = sim_nn_K( n );
  V_c.x = ( similar[n_sim][n]*genrand_real1()*w_cohesion )*( C_sim.x - datoid[n].r.x );
  V_c.y = ( similar[n_sim][n]*genrand_real1()*w_cohesion )*( C_sim.y - datoid[n].r.y );

  n_dis = dis_min(n);
  if( n_dis < 0 || n_dis >= NUM) printf("n_dis out of range\n");
  if( dist[n_dis][n] < Td ) {
    V_r.x = ( 1 - similar[n_sim][n] )*w_separation*( datoid[n].r.x - datoid[n_dis].r.x );
    V_r.y = ( 1 - similar[n_sim][n] )*w_separation*( datoid[n].r.y - datoid[n_dis].r.y );
  }
  else {
    V_r = null_vec;
  }

  // from here 0608
  for( i=0; i<NUM; i++) {
    if( i != n ) {
      if( distance[i][n] < Td2 ) {
        V_r2.x -= datoid[i].r.x;
        V_r2.y -= datoid[i].r.y;
        Vcount++;
      }
    }
  }
  V_r2.x /= (double)Vcount;
  V_r2.y /= (double)Vcount;
  V_r2.x += datoid[n].r.x;
  V_r2.y += datoid[n].r.y;
  V_r2.x *= w_separation;
  V_r2.y *= w_separation;
  // 0608

  datoid[n].v.x *= w;
  datoid[n].v.y *= w;

  datoid[n].v.x += V_a.x + V_c.x + V_r.x + V_r2.x; // 0608
  datoid[n].v.y += V_a.y + V_c.y + V_r.y + V_r2.y; // 0608

  datoid[n].r.x += datoid[n].v.x;
  datoid[n].r.y += datoid[n].v.y;
    
  if( datoid[n].r.x > (width - lattice_constant) ) {
    //    datoid[n].r.x = 2*(width-lattice_constant) - datoid[n].r.x;
    datoid[n].r.x = width -1.7*lattice_constant;
    datoid[n].v.x = -datoid[n].v.x;
  }
  else if(datoid[n].r.x < lattice_constant ) {
    //    datoid[n].r.x = 2*lattice_constant - datoid[n].r.x;
    datoid[n].r.x = 0.3*lattice_constant;
    datoid[n].v.x = -datoid[n].v.x;
  }
  if( datoid[n].r.y > (height - lattice_constant) ) {
    //    datoid[n].r.y = 2*(height-lattice_constant) - datoid[n].r.y;
    datoid[n].r.y = height-1.7*lattice_constant;
    datoid[n].v.y = -datoid[n].v.y;
  }
  else if(datoid[n].r.y < lattice_constant ) {
    //    datoid[n].r.y = 2*lattice_constant - datoid[n].r.y;
    datoid[n].r.y = 0.3*lattice_constant;
    datoid[n].v.y = -datoid[n].v.y;
  }

  //  X = datoid[n].r.x;
  //  Y = datoid[n].r.y;
  nx = (int)( datoid[n].r.x/lattice_constant + 0.5 );
  ny = (int)( datoid[n].r.y/lattice_constant + 0.5 );
  datoid[n].nx = nx;
  datoid[n].ny = ny;
  datoid[n].r.x = lattice_constant*nx;
  datoid[n].r.y = lattice_constant*ny;

  //  if(nx < 0 || nx >= NL || ny < 0 || ny >= NL) printf(" (nx, ny) is out of range\n");

  l = place[ nx ][ ny ];

  //   println("width"+width);
  //   println("hight"+height);
    
  distance_1( n );
/*
  if( l >= 0 ) {
    if( l>=NUM) printf("error in l\n");
    distance_1( l );
  }
*/
  return l;
}

int sim_min(int n) {
  int i, index=-1;
  double min = dmax; 

  for(i = 0; i < NUM; i++) {
    if( i != n) {
      if( sd[i][n] < min ) {
        min = sd[i][n];
        index = i;
      }
    }
  }
  if( index < 0 ) printf("1:index < 0");
  return index;
}

int dis_min(int n) {
  int i, index=-1;
  double min = dmax;  // dmax;
  
  for( i = 0; i < NUM; i++) {
    if( i != n) {
      if( dd[i][n] < min ) {
        min = dd[i][n];
        index = i;
      }
    }
  }

  if( index < 0 ) {
    printf("2:index < 0");
    printf("(%d:%lf) ", n, min);
  }
  return index;
}

Mvector sim_nn_K( int n ) {
  int index[NUM];
  double tmp; 
  double d[NUM];
  int i, j, in, count=0;
  Mvector C = {0.0, 0.0};

  //  C = null_vec;
  
  for(i = 0; i<NUM; i++) {
    index[i] = i;
    d[i] = sd[i][n];
  }
  
  for(i=0; i<K_near+1; i++) {
    for(j=1; j<NUM-i; j++) {
      if( d[j] > d[j-1] ) {
	tmp = d[j];
	in = index[j];
	d[j] = d[j-1];
	index[j] = index[j-1];
	d[j-1] = tmp;
	index[j-1] = in;
      }
    }
  }
  i=NUM-1;
  while(count<K_near) {
    if( index[i] != n) {
      C.x += datoid[ index[i] ].r.x;
      C.y += datoid[ index[i] ].r.y;
      count++;
    }
    i--;
  }

  //  if(count==0) printf("error count=0!\n");
  C.x /= count;
  C.y /= count;
  
  return C;
}

void printState(void) {
  int i, j;
  int no;

  fprintf(fp3, "===============================\n");
  for(i=0; i<NL; i++) {
    for(j=0; j<NL; j++) {
      if( place[i][j] >= 0 ) {
	//	no = datoid[ place[i][j] ].o.x[DIM-1];
	no = (int)datoid[ place[i][j] ].o.x[DIM-1];
	fprintf(fp3, "%d", no);
      }
      else fprintf(fp3, " ");
    }
    fprintf(fp3, "\n");
  }
  fprintf(fp3, "===============================\n");
}

void distance_1(int n) {
  int i;
  
  for(i=0; i<NUM; i++) {
    dist[n][i] = mapDist( datoid[i].r, datoid[n].r );
    dist[i][n] = dist[n][i];
  }
  for(i=0; i<NUM; i++) {
    sd[i][n] = similar[i][n]*dist[i][n] + (1 - similar[i][n])*dmax;
    dd[i][n] = ( 1 - similar[i][n])*dist[i][n] + similar[i][n]*dmax;
    sd[n][i] = sd[i][n];
    dd[n][i] = dd[i][n];
  }
}

int avoidOverlap(int ix, int iy,int n, double ls) {
  int flag=0, l=-1, ret = -1;
  double tmpX, tmpY;
  double simn=0, siml=0;
  int lx, ly;
  //  int ix, iy;
   
  lx = datoid[n].nx;
  ly = datoid[n].ny;
  //  ix = (int)(preX/ls);
  //  iy = (int)(preY/ls);
  //  if(ix<0 || ix>=NL || iy<0 || iy>= NL ) printf("(ix, iy) is out of range\n");
  //  if(lx<0 || lx>=NL || ly<0 || ly>= NL ) printf("(lx, ly) is out of range\n");
  if( place[lx][ly] >= 0 ) {
    if( !(lx == ix && ly == iy ) ) {
      l = place[lx][ly];
      flag = 1;
      datoid[l].r.x = ls*ix;
      datoid[l].r.y = ls*iy;
      place[lx][ly] = n;
      place[ix][iy] = l;
      datoid[l].nx = ix;
      datoid[l].ny = iy;
      ret = l;
    }
  }
  else {
    place[ix][iy] = -1;
    place[lx][ly] = n;
  }
  return ret;
}

double analysis(void) {
  double ave_similar=0.0;
  int i, ix, iy, n;
  int count, total=0;
  double tmp;
  
  for(i=1; i<NUM; i++) {
    ix = datoid[i].nx;
    iy = datoid[i].ny;
    if(ix >= 1 && ix <=NL-2 && iy >= 1 && iy <=NL-2 ) {
      tmp = 0.0;
      count = 0;
      n = place[ix+1][iy-1];
      if( n >= 0) {
	tmp += similar[ i ][ n ];
	count++;
      }
      n = place[ix+1][iy];
      if( n >= 0) {
	tmp += similar[ i ][ n ];
	count++;
      }
      n = place[ix+1][iy+1];
      if( n >= 0) {
	tmp += similar[ i ][ n ];
	count++;
      }
      n = place[ix][iy-1];
      if( n >= 0) {
	tmp += similar[ i ][ n ];
 	count++;
      }
      n = place[ix][iy+1];
      if( n >= 0) {
	tmp += similar[ i ][ n ];
	count++;
      }
      n = place[ix-1][iy-1];
      if( n >= 0) {
	tmp += similar[ i ][ n ];
	count++;
      }
      n = place[ix-1][iy];
      if( n >= 0) {
	tmp += similar[ i ][ n ];
	count++;
      }
      n = place[ix-1][iy+1];
      if( n >= 0) {
	tmp += similar[ i ][ n ];
	count++;
      }

      if(count>0) {
	ave_similar += tmp/(double)count;
	total++;
      }
    }
  }
  if(total>0) ave_similar /= (double)total;

  return ave_similar;
  //  fprintf(fp2,"ave_sim = %lf\n", ave_similar);
}

double G_K_factor(void) {
  int i, j, k, l, NN;
  double C=0.0, D=0.0, E=0.0;
  double GK;
  double distD[NUM][NUM];
  double tmp;

  NN = NUM*(NUM-1)/2;

  for(i=0; i<NUM; i++) {
    for(j=0; j<NUM; j++) {
      distD[i][j] = 1.0 - similar[i][j];
    }
  }
  for(i=0; i<NUM; i++) {
    for(j=i+1; j<NUM; j++) {
      for(k=0; k<NUM; k++) {
	for(l=k+1; l<NUM; l++) {
	  tmp = (dist[i][j] - dist[k][l])*(distD[i][j] - distD[k][l]);
	  ( tmp > 0.0 ? C++ : (tmp < 0.0 ?  D++ : E++) ); 
	}
      }
    }
  }
  if( C+D > 0.2 )  GK = (C-D)/(C+D);
  else GK = 0.0;
  
  return GK;
  //  fprintf(fp2, "GK = %e", GK);
}

/* initializes mt[N] with a seed */
void init_genrand(unsigned long s)
{
    mt[0]= s & 0xffffffffUL;
    for (mti=1; mti<N; mti++) {
        mt[mti] = 
	    (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti); 
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        mt[mti] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
}

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
/* slight change for C++, 2004/2/26 */
void init_by_array(unsigned long init_key[], int key_length)
{
    int i, j, k;
    init_genrand(19650218UL);
    i=1; j=0;
    k = (N>key_length ? N : key_length);
    for (; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL))
          + init_key[j] + j; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++; j++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
        if (j>=key_length) j=0;
    }
    for (k=N-1; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL))
          - i; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
    }

    mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */ 
}

/* generates a random number on [0,0xffffffff]-interval */
unsigned long genrand_int32(void)
{
    unsigned long y;
    static unsigned long mag01[2]={0x0UL, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (mti >= N) { /* generate N words at one time */
        int kk;

        if (mti == N+1)   /* if init_genrand() has not been called, */
            init_genrand(5489UL); /* a default initial seed is used */

        for (kk=0;kk<N-M;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (;kk<N-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mti = 0;
    }
  
    y = mt[mti++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}

/* generates a random number on [0,0x7fffffff]-interval */
long genrand_int31(void)
{
    return (long)(genrand_int32()>>1);
}

/* generates a random number on [0,1]-real-interval */
double genrand_real1(void)
{
    return genrand_int32()*(1.0/4294967295.0); 
    /* divided by 2^32-1 */ 
}

/* generates a random number on [0,1)-real-interval */
double genrand_real2(void)
{
    return genrand_int32()*(1.0/4294967296.0); 
    /* divided by 2^32 */
}

/* generates a random number on (0,1)-real-interval */
double genrand_real3(void)
{
    return (((double)genrand_int32()) + 0.5)*(1.0/4294967296.0); 
    /* divided by 2^32 */
}

/* generates a random number on [0,1) with 53-bit resolution*/
double genrand_res53(void) 
{ 
    unsigned long a=genrand_int32()>>5, b=genrand_int32()>>6; 
    return(a*67108864.0+b)*(1.0/9007199254740992.0); 
} 
/* These real versions are due to Isaku Wada, 2002/01/09 added */
