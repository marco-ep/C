#include <FPT.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "D3d_matrix.h"

int numobjects;
int numpoints[10];
double x[10][10000], y[10][10000], z[10][10000];
int numpolys[10000];
int psize[10][10000];
int con[10][10000][200];
double red[10][5000], grn[10][5000], blu[10][5000];
double cx,cy,cz;

double lightx, lighty,lightz;
double rinit, ginit, binit;
double ambient, maxdiff, specpow ; 

typedef
struct {
  int objnum ;
  int polynum ;
  double dist ;
}
  THING;

THING allpolys[10000];

void print_array(int n) {
  int i ;
  for (i = 0 ; i < n ; i++) {
    printf("%d %d %lf\n",allpolys[i].objnum, allpolys[i].polynum, allpolys[i].dist) ;
  }
  printf("\n") ;
}

int compare (const void *p, const void *q) {
  THING *a, *b ;

  a = (THING*)p ;
  b = (THING*)q ;

  if  (((*a).dist) < ((*b).dist)) return -1 ;
  else if (((*a).dist) > ((*b).dist)) return 1 ;
  else return 0;
}


int init_array() {  

  int allpolysN = 0;
  int i;
  int j;

  for (i = 0; i < numobjects; i++) { //COUNT TOTAL NUMBER OF POLYGONS IN SORT and INITIALIZE 
    for (j = 0; j < numpolys[i]; j++) {
      allpolys[allpolysN].objnum = i;
      allpolys[allpolysN].polynum = j;
      allpolys[allpolysN].dist = z[i][con[i][j][0]];
      allpolysN = allpolysN + 1;
    }  
  }

  qsort (allpolys, allpolysN, sizeof(THING), compare);

  return allpolysN;
}


int crossprod(double *v,double *u, double w[3]){

  w[0]= v[1]*u[2] - v[2]*u[1];
  w[1]= -(v[0]*u[2] - v[2]*u[0]);
  w[2]= v[0]*u[1] - v[1]*u[0];

  // nx = vY*vZZ - vZ*vYY
  // ny = -vX*vZZ - vZ*vXX 
  // nz = vX*vYY - vY*vXX
	
  return 1;

} 

int vecs(double *v,double *u,double w[3]){
	
  w[0]= v[0]-u[0];
  w[1]= v[1]-u[1];
  w[2]= v[2]-u[2];

  /*
  if(w[0] == w[1] && w[1] == w[2]){ 
	printf("scalar multiples");
	return -1;
  }
  */	
  return 1;
}

double vecsize(double *v){

  return sqrt((v[0]*v[0])+(v[1]*v[1])+(v[2]*v[2]));

}

double dotprod(double *v,double *u){ //which = cos(beta)

  return u[0]*v[0]+u[1]*v[1]+u[2]*v[2];

}





double xmax(int onum){
  int i;
  double max= x[onum][0];
  for(i=0;i<numpoints[onum];i++){
    if(x[onum][i] > max){
      max = x[onum][i];
    }
  }
  return max;
}

double xmin(int onum){
  int i;
  double min= x[onum][0];
  for(i=0;i<numpoints[onum];i++){
    if(x[onum][i] < min){
      min = x[onum][i];
    }
  }
  return min;
}
double ymax(int onum){
  int i;
  double max= y[onum][0];
  for(i=0;i<numpoints[onum];i++){
    if(y[onum][i] > max){
      max = y[onum][i];
    }
  }
  return max;
}

double ymin(int onum){
  int i;
  double min= y[onum][0];
  for(i=0;i<numpoints[onum];i++){
    if(y[onum][i] < min){
      min = y[onum][i];
    }
  }
  return min;
}

double zmax(int onum){
  int i;
  double max= z[onum][0];
  for(i=0;i<numpoints[onum];i++){
    if(z[onum][i] > max){
      max = z[onum][i];
    }
  }
  return max;
}

double zmin(int onum){
  int i;
  double min= z[onum][0];
  for(i=0;i<numpoints[onum];i++){
    if(z[onum][i] < min){
      min = z[onum][i];
    }
  }
  return min;
}

int lightcolor(THING allpolys,double *shade){

  double intensity;
  double cosA ;
  double Nu[3], Lu[3], Ru[3], Eu[3]; 
  double v[3], u[3];

  //find v and u inside polygon:

  double p[3], q[3], r[3];


  p[0]= x[allpolys.objnum][con[allpolys.objnum][allpolys.polynum][0]]; //index
  p[1]= y[allpolys.objnum][con[allpolys.objnum][allpolys.polynum][0]];
  p[2]= z[allpolys.objnum][con[allpolys.objnum][allpolys.polynum][0]];

  q[0]= x[allpolys.objnum][con[allpolys.objnum][allpolys.polynum][1]]; //index1
  q[1]= y[allpolys.objnum][con[allpolys.objnum][allpolys.polynum][1]];
  q[2]= z[allpolys.objnum][con[allpolys.objnum][allpolys.polynum][1]]; 

  r[0]= x[allpolys.objnum][con[allpolys.objnum][allpolys.polynum][2]]; //index2
  r[1]= y[allpolys.objnum][con[allpolys.objnum][allpolys.polynum][2]];
  r[2]= z[allpolys.objnum][con[allpolys.objnum][allpolys.polynum][2]]; 

  vecs(q,p,v); //v[0] = vX      v[1] = vY     v[2] = vZ
  vecs(r,p,u) ;//u[0] = vXX     u[1] = vYY    u[2] = vZZ

  crossprod(v, u, Nu);

  double len ;

  len = vecsize(Nu) ;
  Nu[0] = Nu[0] / len ;
  Nu[1] = Nu[1] / len ;
  Nu[2] = Nu[2] / len ;

  /*
  Nu[0] = Nu[0] / vecsize(Nu) ; 
  Nu[1] = Nu[1] / vecsize(Nu) ; 
  Nu[2] = Nu[2] / vecsize(Nu) ;
  */

  Lu[0] = lightx - p[0] ;
  Lu[1] = lighty - p[1] ;
  Lu[2] = lightz - p[2]  ;
  len = vecsize(Lu) ;
  Lu[0] = Lu[0]/len ;
  Lu[1] = Lu[1]/len ;
  Lu[2] = Lu[2]/len ;  
  

  Eu[0] = 0 - p[0] ;
  Eu[1] = 0 - p[1] ;
  Eu[2] = 0 - p[2] ;
  len = vecsize(Eu) ;
  Eu[0] = Eu[0]/len ;
  Eu[1] = Eu[1]/len ;
  Eu[2] = Eu[2]/len ;    



  
  if( (dotprod(Nu,Eu)<0 && dotprod(Nu,Lu)>0) || (dotprod(Nu,Eu)>0 && dotprod(Nu,Lu)<0) ){ //intensity = ambient
        shade[0] = ambient/(ambient+maxdiff)*rinit; 
	shade[1] = ambient/(ambient+maxdiff)*ginit; 
	shade[2] = ambient/(ambient+maxdiff)*binit; 
	return 1;
  } // skip the rest and go to the next polygon



  
  
  if(dotprod(Nu,Lu) < 0 && dotprod(Nu,Eu) < 0){
    Nu[0] = -1*Nu[0]; Nu[1] = -1*Nu[1]; Nu[2] = -1*Nu[2];
   } //reverse Nu

  double cosB = dotprod(Nu,Lu) ;
  double h = cosB ;
  h = 2*h ;
  
  Ru[0] = (h*Nu[0] - Lu[0])  ;
  Ru[1] = (h*Nu[1] - Lu[1])  ;
  Ru[2] = (h*Nu[2] - Lu[2])  ;



  cosA = dotprod(Eu,Ru) ;
  intensity = ambient + maxdiff*cosB + (1-ambient-maxdiff)*pow(cosA,specpow);

  if(intensity<ambient+maxdiff){
  shade[0] = (intensity/(ambient+maxdiff))*rinit;    
  shade[1] = (intensity/(ambient+maxdiff))*ginit;
  shade[2] = (intensity/(ambient+maxdiff))*binit;
  }
  else{
    double f = (intensity - (ambient+maxdiff))/(1- (ambient+maxdiff)) ;
    shade[0] = rinit + f*(1-rinit) ;
    shade[1] = ginit + f*(1-ginit) ;
    shade[2] = binit + f*(1-binit) ;
  
  }
  return 1;
}



void draw_file(int onum){

  int i,j,k;
  double x_current[10000], y_current[10000], z_current[10000];
  double light[3];

  int allpolysN = init_array();

  for(k=allpolysN; k>=0; k--){
    onum = allpolys[k].objnum;
    i = allpolys[k].polynum;	
	

    //for(i=0; i<numpolys[onum];i++){
    for(j=0; j<psize[onum][i]; j++){
      // finding the original values of x, y and z points
      x_current[j] = x[onum][con[onum][i][j]];
      y_current[j] = y[onum][con[onum][i][j]];
      z_current[j] = z[onum][con[onum][i][j]];			 

      // find the distance each point has to projection screen
      y_current[j] = y_current[j] / z_current[j];
      x_current[j] = x_current[j] / z_current[j];
      z_current[j] = 1;  
	
      //translate the H,H perspective to actual graphics window
      x_current[j] = (400/tan(20*M_PI/180))*x_current[j] + 400;
      y_current[j] = (400/tan(20*M_PI/180))*y_current[j] + 400;
    }


    //Instead of changing the color for each object here, you should call the light function to change each polygons color!!!
    //give it a polygon

    lightcolor(allpolys[k], light);
    G_rgb(light[0],light[1],light[2]);
    //    G_rgb(1,0,0) ;    
    G_fill_polygon(x_current,y_current,psize[onum][i]);
    G_rgb(0,0,0) ;
    G_polygon(x_current, y_current, psize[onum][i]);
  }
}

void center(int onum){

  double m[4][4];
  double minv[4][4];

  D3d_make_identity(m);
  D3d_make_identity(minv);

  double xBig =  xmax(onum);
  double xSmall = xmin(onum);
  double yBig =  ymax(onum);
  double ySmall = ymin(onum);
  double zBig = zmax(onum);
  double zSmall = zmin(onum);


  double xCenter = (xBig + xSmall) / 2 ;
  double yCenter = (yBig + ySmall) / 2;
  double zCenter = (zBig + zSmall) / 2;

  D3d_translate(m,minv,-xCenter, -yCenter, -zCenter);

  double sf ;
  if(xBig-xSmall > yBig-ySmall && xBig-xSmall > zBig-zSmall){
    sf = 1/(xBig-xSmall); 
  } else if(yBig-ySmall > xBig-xSmall && yBig-ySmall > zBig-zSmall) {  
    sf = 1/(yBig-ySmall) ; 
  } else {
    sf = 1/(zBig-zSmall);
  }

  D3d_scale(m,minv,sf,sf,sf); 
  D3d_translate(m,minv,0,0,20);
  D3d_mat_mult_points(x[onum],y[onum],z[onum],m,x[onum],y[onum],z[onum],numpoints[onum]);
}


void objCenter(int onum){

  double m[4][4];
  double minv[4][4];
  double xSum, ySum, zSum;
  int i;

  D3d_make_identity(m);
  D3d_make_identity(minv);

  xSum = ySum = zSum = 0 ;

  for(i=0; i<numpoints[onum]; i++){
    xSum += x[onum][i];
    ySum += y[onum][i];
    zSum += z[onum][i];
  }

  cx = xSum / numpoints[onum];
  cy = ySum / numpoints[onum];
  cz = zSum / numpoints[onum];

} 

int main(int argc, char **argv){

  FILE *f ;
  int a,b,c,d,e,g,i,j,onum ;

  double m[4][4];
  double minv[4][4];
  onum = 0;
  i = 0;


  /*
  //let user input values for LIGHT etc
  printf("input x,y,z coordinates for the light source \n");
  printf("x: \n");
  scanf("%lf", &lightx);
  printf("y: \n");
  scanf("%lf", &lighty);
  printf("z: \n");
  scanf("%lf", &lightz);

  //input values for R,G,B
  printf("enter RGB values for color all in one line with a space inbetween\n");
  scanf("%lf %lf %lf", &rinit, &ginit, &binit);

  printf("input values for Ambient, MaxDiffuse and Specpow \n");
  printf("Ambient: \n");
  scanf("%lf", &ambient);
  printf("MaxDiffuse: \n");
  scanf("%lf", &maxdiff);
  printf("Specpow: \n");
  scanf("%lf", &specpow);
  */

  lightx = 100 ;
  lighty = 200 ;
  lightz = 0 ;
  rinit = 0.62 ;
  ginit = 0.37 ;
  binit = 0.19 ;

  ambient = 0.2 ;
  maxdiff = 0.5 ;
  specpow = 50 ;


  //initialize Graphics

  G_init_graphics(800,800);
  G_rgb(0,0,0);
  G_clear();	

  numobjects = argc-1;
    
  //attempt to open file
  for(i=0; i<numobjects; i++){

    f= fopen(argv[i+1], "r");
    if(f == NULL) {
      printf("can't open file\n");
      exit(0);
    }

    //read file
	
    fscanf(f,"%d",&numpoints[i]);

    for(a=0; a<numpoints[i];a++){
      fscanf(f, "%lf %lf %lf", &x[i][a], &y[i][a], &z[i][a]);
    }

    fscanf(f,"%d",&numpolys[i]);

    for(b=0; b<numpolys[i]; b++){
      fscanf(f, "%d", &psize[i][b]);	
      for(c=0; c < psize[i][b]; c++){
	fscanf(f, "%d", &con[i][b][c]);	
      }
    }
    center(i);
  }


  //draw stage

  char v = '0' ;
  int go = 0; 

     
  //Q TO QUIT

  while(v != 'q'){

    G_rgb(0,0,0);
    D3d_make_identity(m);
    D3d_make_identity(minv);
    G_clear();
    objCenter(go);
    
    if( v == 'x'){   //ROTATE X
      D3d_translate(m,minv, -cx, -cy,-cz);
      D3d_rotate_x(m,minv, .05);
      D3d_translate(m,minv, cx, cy, cz);}
    else if( v == 'y'){   // ROTATE Y
      D3d_translate(m,minv, -cx, -cy,-cz);
      D3d_rotate_y(m,minv, .05);
      D3d_translate(m,minv, cx, cy, cz);}
    else if( v == 'z'){   // ROTATE Z         
      D3d_translate(m,minv, -cx, -cy,-cz);
      D3d_rotate_z(m,minv, .05);
      D3d_translate(m,minv, cx, cy, cz);}

    
    else if( v == 't'){ //FORWARD X
      D3d_translate(m,minv,1,0,0);}
    else if( v == 'w'){ // FORWARD Y
      D3d_translate(m,minv,0,1,0);}
        else if( v == 'd'){ // FORWARD Z
      D3d_translate(m,minv,0,0,-1);}

    
    else if( v == 'g'){ // BACKWARD X
      D3d_translate(m,minv,-1,0,0);}
    else if( v == 's'){ // BACKWARD Y
      D3d_translate(m,minv,0,-1,0);}
        else if( v == 'e'){ // BACKWARD Z
      D3d_translate(m,minv,0,0,1);  
      z[go][numpoints[go]] += c;
    }

    
    else{
      go = v - 48; //ASCII code
    }

    D3d_mat_mult_points(x[go],y[go],z[go],m,x[go],y[go],z[go],numpoints[go]);
    draw_file(go);
    v = G_wait_key();
  }
}
