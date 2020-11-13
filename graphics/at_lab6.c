#include <FPT.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "D3d_matrix.h"

int numobjects;
int numpoints[10];
double x[10][10000], y[10][10000], z[10][10000];
int numpolys[10000];
int psize[10][1000];
int con[10][5000][2000]; 
double red[10][5000], grn[10][5000], blu[10][5000];
double cx,cy,cz;

typedef
struct {
  int objnum ;
  int polynum ;
  double dist ;
}
THING ;

THING allpolys[20000];
int allpolysN = 0;
int init_array() {
  
  int i;
  int j;

  allpolysN = 0;
  for (i = 0; i < numobjects; i++) {
    for (j = 0; j < numpolys[i]; j++) {
      allpolys[allpolysN].objnum = i;
      allpolys[allpolysN].polynum = j;
      allpolys[allpolysN].dist = z[i][con[i][j][0]];
      allpolysN++ ;
    }  
  }
}

void print_array(int n) {
  int i ;
  for (i = 0 ; i < n ; i++) {
    printf("allpolys[%d] objnum = %d polynum = %d z = %lf\n",i, allpolys[i].objnum, allpolys[i].polynum, allpolys[i].dist) ;
  }
  printf("\n") ;
}

int compare(const void *p, const void *q){
  THING *a, *b ;

  a = (THING*)p ;
  b = (THING*)q ;

  if (((*a).dist) < ((*b).dist)) return -1 ;
  else if (((*a).dist) > ((*b).dist)) return 1 ;
  else return 0 ;
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



void draw_file(int total){

  G_rgb(0,0,0) ;
  G_clear() ;
  int r, g ; 

  double x_current[1000], y_current[1000], z_current[1000] ;
  int onumCurr, polynumCurr; 

for(int i=total-1; i>=0; i--){  

  onumCurr = allpolys[i].objnum ;
  polynumCurr = allpolys[i].polynum ;
  r = 0 ; g = 0 ;
  
  for(int j=0; j<psize[onumCurr][polynumCurr]; j++){
    // finding the original values of x, y and z points
    x_current[j] = x[onumCurr][con[onumCurr][polynumCurr][j]];
    y_current[j] = y[onumCurr][con[onumCurr][polynumCurr][j]];
    z_current[j] = z[onumCurr][con[onumCurr][polynumCurr][j]];
    
    // find the distance each point has to projection screen
    y_current[j] = y_current[j] / z_current[j];
    x_current[j] = x_current[j] / z_current[j];
    z_current[j] = 1;  
    
    //translate the H,H perspective to actual graphics window
    x_current[j] = (400/tan(20*M_PI/180))*x_current[j] + 400;
    y_current[j] = (400/tan(20*M_PI/180))*y_current[j] + 400;
    
    //you could do the two steps above in one step...
  } // end for j
  if(onumCurr==0){
    r = 1 ;
  }
  else{
    g = 1 ;
  }
  G_rgb(r, g, 0) ;
  G_fill_polygon(x_current, y_current, psize[onumCurr][polynumCurr]) ;

  G_rgb(0,0,0) ;
  G_polygon(x_current, y_current, psize[onumCurr][polynumCurr]) ;  
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
        sf = 1/(xBig-xSmall);  // it's 1 instead of 400 bc draw is centering it at 400,400
    } else if(yBig-ySmall > xBig-xSmall && yBig-ySmall > zBig-zSmall) {  
      sf = 1/(yBig-ySmall) ; // fits it into 1 by 1 box
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
    numobjects = argc - 1 ;

    //initialize graphics
    G_init_graphics(800,800);
    G_rgb(0,0,0);
    G_clear();	

    int total = 0 ;
    //attempt to open file 
        for(i=0; i<argc-1; i++){
      
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
      total = total + numpolys[i] ;
      for(b=0; b<numpolys[i]; b++){
	fscanf(f, "%d", &psize[i][b]);	
	for(c=0; c < psize[i][b]; c++){
	  fscanf(f, "%d", &con[i][b][c]);
	  //printf("%lf\n", z[i][con[i][b][c]]) ;
	}
      }
       center(i);
    }

	init_array() ;	
	printf("total = %d  allpolysN = %d\n",total, allpolysN) ;
	qsort (allpolys, total, sizeof(THING), compare) ;
      //print_array(total) ;
	draw_file(total) ;
      
    //draw stage
    
    int v = G_wait_key() ;
    int go = 0; //go is the object number (onum) when you press keys 0 1 2 so on == onum
    
    while(v != 'q'){
      D3d_make_identity(m);
      D3d_make_identity(minv);
      objCenter(go);
       //find the bounding box
      if( v == 'x'){  //rotate x
	D3d_translate(m,minv, -cx, -cy,-cz);
	D3d_rotate_x(m,minv, .05);
	D3d_translate(m,minv, cx, cy, cz);
       }
       else if( v == 'y'){   // rotate y
	 D3d_translate(m,minv, -cx, -cy,-cz);
	 D3d_rotate_y(m,minv, .05);
	 D3d_translate(m,minv, cx, cy, cz);}
       else if( v == 'z'){   // rotate z
	 D3d_translate(m,minv, -cx, -cy,-cz);
	 D3d_rotate_z(m,minv, .05);
	 D3d_translate(m,minv, cx, cy, cz);}
      
       else if( v == 't'){ // forward x
	 D3d_translate(m,minv,.5,0,0);}
       else if( v == 'w'){ // forward y
	 D3d_translate(m,minv,0,.5,0);
       }
       else if( v == 'd'){ // forward in z direction
	 D3d_translate(m,minv,0,0,-.5);
       }

       else if( v == 'g'){ // backward in x direction
	 D3d_translate(m,minv,-.5,0,0);}
       else if( v == 's'){ // backward in y direction
	 D3d_translate(m,minv,0,-.5,0);}
       else if( v == 'e'){ //backward in z direction
	 D3d_translate(m,minv,0,0,.5);
       }      
      
       else{
	 go = v - 48 ;
       }
       
       D3d_mat_mult_points(x[go],y[go],z[go],m,x[go],y[go],z[go],numpoints[go]);

       init_array() ;	
      qsort (allpolys, total, sizeof(THING), compare) ;
       draw_file(total);
       
       v = G_wait_key();
    }
}

