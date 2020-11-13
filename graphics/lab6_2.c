#include <FPT.h>
#include <D3d_matrix.h>


typedef
struct {
  int objnum;
  int polynum;
  double dist;
}
  THING;

FILE *f;
int numpoints, numpolys,d;
double x[10][20000], y[10][20000], z[10][20000];
double newx[20000], newy[20000], newz[20000];
int connect[10][10000][10];
double xmax,ymax,zmax, xmin,ymin,zmin, xcenter,ycenter,zcenter;
THING t[10000];


void draw(double H){
  double m[4][4], minv[4][4];
  int i,j,k;

  G_rgb(0,0,0);
  G_clear();

  //for each THING in t[], starting from end
  for(i=d-1; i>=0; i--){

    // draw t[i]

    //arrays of x,y,z of just this poly
    for(j=0; j<connect[t[i].objnum][t[i].polynum][0]; j++){
      newx[j]=x[t[i].objnum][connect[t[i].objnum][t[i].polynum][j+1]];
      newy[j]=y[t[i].objnum][connect[t[i].objnum][t[i].polynum][j+1]];
      newz[j]=z[t[i].objnum][connect[t[i].objnum][t[i].polynum][j+1]];
    }


    for(j=0; j<connect[t[i].objnum][t[i].polynum][0]; j++){
      newx[j]=newx[j]/newz[j];
      newy[j]=newy[j]/newz[j];
      newz[j]=1;
    }

    D3d_make_identity(m);
    D3d_make_identity(minv);
    D3d_scale(m,minv,400/H,400/H,0);
    D3d_translate(m,minv,400,400,0);
    D3d_mat_mult_points(newx,newy,newz,m,newx,newy,newz,numpoints);


    if(t[i].objnum==0){G_rgb(1,0,0);}
    else if(t[i].objnum==1){G_rgb(0,1,0);}
    else if(t[i].objnum==2){G_rgb(0,0,1);}

    G_fill_polygon(newx,newy,connect[t[i].objnum][t[i].polynum][0]);

    G_rgb(0,0,0);
    G_polygon(newx,newy,connect[t[i].objnum][t[i].polynum][0]);
    
  }
  
  }

/////////////////////////////////////////////////////////////////////////////
void init_array(int objnum, int d)
{
  int i,j;
  double m[4][4], minv[4][4];
  
  //for each poly in object
  for(i=0; i<numpolys; i++){
    t[d].objnum=objnum;
    t[d].polynum=i;

    for(j=0; j<connect[objnum][i][0]; j++){
      newx[j]=x[objnum][connect[objnum][i][j+1]];
      newy[j]=y[objnum][connect[objnum][i][j+1]];
      newz[j]=z[objnum][connect[objnum][i][j+1]];
      }
    D3d_make_identity(m);
    D3d_make_identity(minv);
    double dz; dz=fabs(3-zmin);
    D3d_translate(m,minv,0,0,dz);
    D3d_mat_mult_points(newx,newy,newz,m,newx,newy,newz,numpoints);

    t[d].dist=newz[0];
    
    /*for(j=0; j<connect[objnum][i][0]; j++){
      z[objnum][connect[objnum][i][j+1]]=newz[j];
      }*/
    
    d=d+1;
  }

  /*for(j=0; j<numpoints; j++){
    printf("%lf %lf %lf\n",x[objnum][j],y[objnum][j],z[objnum][j]);
    }*/
  
    D3d_make_identity(m);
    D3d_make_identity(minv);
    double dz; dz=fabs(3-zmin);
    D3d_translate(m,minv,0,0,dz);
    D3d_mat_mult_points(x[objnum],y[objnum],z[objnum],m,x[objnum],y[objnum],z[objnum],numpoints);
    
}

/////////////////////////////////////////////////////////////////////////////
void print_array(THING *x,int n)
{
  int i ;
  for (i = 0 ; i < n ; i++) {
    printf("%d %d %lf\n",x[i].objnum, x[i].polynum, x[i].dist) ;
  }
  printf("\n") ;
}

/////////////////////////////////////////////////////////////////////////////
int compare (const void *p, const void *q)
{
  THING *a, *b ;

  a = (THING*)p ;
  b = (THING*)q ;

  if  (((*a).dist) < ((*b).dist)) return -1 ;
  else if (((*a).dist) > ((*b).dist)) return 1 ;
  else return 0 ;
}

/////////////////////////////////////////////////////////////////////////////
void translate(double H, int objnum){
  int q,i;
  int a,c;
  a=1;
  c=1;
  double m[4][4],minv[4][4];
  q=G_wait_key();

  
  if(q==120){
    while(a==1){
      D3d_make_identity(m);
      D3d_make_identity(minv);
      D3d_translate(m,minv,c,0,0);
      D3d_mat_mult_points(x[objnum],y[objnum],z[objnum],m,x[objnum],y[objnum],z[objnum],numpoints);
      //change .dist value for each t entry of this obj
      for(i=0; i<d; i++){
	if(t[i].objnum==objnum){
	  t[i].dist=z[objnum][connect[objnum][t[i].polynum][1]];
	}
      }
      qsort(t,d,sizeof(THING),compare);
      draw(H);
      q=G_wait_key();
      if(q==99){c=c*(-1);}
      else if(q==116||q==114||q==32||q==121||q==122){break;}
    }
  }
  if(q==121){
    while(a==1){
      D3d_make_identity(m);
      D3d_make_identity(minv);
      D3d_translate(m,minv,0,c,0);
      D3d_mat_mult_points(x[objnum],y[objnum],z[objnum],m,x[objnum],y[objnum],z[objnum],numpoints);
      for(i=0; i<d; i++){
	if(t[i].objnum==objnum){
	  t[i].dist=z[objnum][connect[objnum][t[i].polynum][1]];
	}
      }
      qsort(t,d,sizeof(THING),compare);
      draw(H);
      q=G_wait_key();
      if(q==99){c=c*(-1);}
      else if(q==116||q==114||q==32||q==120||q==122){break;}
    }
  }
  if(q==122){
    while(a==1){
      D3d_make_identity(m);
      D3d_make_identity(minv);
      D3d_translate(m,minv,0,0,c);
      D3d_mat_mult_points(x[objnum],y[objnum],z[objnum],m,x[objnum],y[objnum],z[objnum],numpoints);
      for(i=0; i<d; i++){
	if(t[i].objnum==objnum){
	  t[i].dist=z[objnum][connect[objnum][t[i].polynum][1]];
	}
      }
      qsort(t,d,sizeof(THING),compare);
      draw(H);
      q=G_wait_key();
      if(q==99){c=c*(-1);}
      else if(q==116||q==114||q==32||q==120||q==121){break;}
    }
  }
}

//////////////////////////////////////////////////////////////////////////////
void rotate(double H,int objnum){
  int q;
  int a,i;
  double c;
  a=1;
  c=2*(M_PI/180);
  double m[4][4],minv[4][4];

  xmax=0; ymax=0; zmax=0;
  xmin=x[objnum][0]; ymin=y[objnum][0]; zmin=z[objnum][0];
  for(i=0; i<numpoints; i++){
    if(x[objnum][i]>xmax){xmax=x[objnum][i];}
    if(y[objnum][i]>ymax){ymax=y[objnum][i];}
    if(z[objnum][i]>zmax){zmax=z[objnum][i];}

    if(x[objnum][i]<xmin){xmin=x[objnum][i];}
    if(y[objnum][i]<ymin){ymin=y[objnum][i];}
    if(z[objnum][i]<zmin){zmin=z[objnum][i];}
  }
  xcenter=(xmax+xmin)/2;
  ycenter=(ymax+ymin)/2;
  zcenter=(zmax+zmin)/2;


  
  q=G_wait_key();
  if(q==120){
    while(q!=116&&q!=114&&q!=32&&q!=121&&q!=122){
      D3d_make_identity(m);
      D3d_make_identity(minv);
      D3d_translate(m,minv,-xcenter,-ycenter,-zcenter);
      D3d_rotate_x(m,minv,c);
      D3d_mat_mult_points(x[objnum],y[objnum],z[objnum],m,x[objnum],y[objnum],z[objnum],numpoints);

      D3d_make_identity(m);
      D3d_make_identity(minv);
      D3d_translate(m,minv,xcenter,ycenter,zcenter);
      D3d_mat_mult_points(x[objnum],y[objnum],z[objnum],m,x[objnum],y[objnum],z[objnum],numpoints);
      for(i=0; i<d; i++){
	if(t[i].objnum==objnum){
	  t[i].dist=z[objnum][connect[objnum][t[i].polynum][1]];
	}
      }
      qsort(t,d,sizeof(THING),compare);
      draw(H);
      q=G_wait_key();
      if(q==99){c=c*(-1);}
      else if(q==116||q==114||q==32||q==121||q==122){break;}
    }
  }
  if(q==121){
    while(q!=116&&q!=114&&q!=32&&q!=120&&q!=122){
      D3d_make_identity(m);
      D3d_make_identity(minv);
      D3d_translate(m,minv,-xcenter,-ycenter,-zcenter);
      D3d_rotate_y(m,minv,c);
      D3d_mat_mult_points(x[objnum],y[objnum],z[objnum],m,x[objnum],y[objnum],z[objnum],numpoints);

      D3d_make_identity(m);
      D3d_make_identity(minv);
      D3d_translate(m,minv,xcenter,ycenter,zcenter);
      D3d_mat_mult_points(x[objnum],y[objnum],z[objnum],m,x[objnum],y[objnum],z[objnum],numpoints);
      for(i=0; i<d; i++){
	if(t[i].objnum==objnum){
	  t[i].dist=z[objnum][connect[objnum][t[i].polynum][1]];
	}
      }
      qsort(t,d,sizeof(THING),compare);
      draw(H);
      q=G_wait_key();
      if(q==99){c=c*(-1);}
      else if(q==116||q==114||q==32||q==120||q==122){break;}
    }
  }
  if(q==122){
    while(q!=116&&q!=114&&q!=32&&q!=120&&q!=121){
      D3d_make_identity(m);
      D3d_make_identity(minv);
      D3d_rotate_z(m,minv,c);
      D3d_mat_mult_points(x[objnum],y[objnum],z[objnum],m,x[objnum],y[objnum],z[objnum],numpoints);
      for(i=0; i<d; i++){
	if(t[i].objnum==objnum){
	  t[i].dist=z[objnum][connect[objnum][t[i].polynum][1]];
	}
      }
      qsort(t,d,sizeof(THING),compare);
      draw(H);
      q=G_wait_key();
      if(q==99){c=c*(-1);}
      else if(q==116||q==114||q==32||q==120||q==121){break;}
    }
  }
 }

////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
  int r,s,i,j,a,q,obj;
  double halfangle,H;
  halfangle=45*(M_PI/180);
  H=tan(halfangle);
  a=1;

  
  G_init_graphics(800,800);
  G_rgb(0,0,0);
  G_clear();

  d=0;

  for(r=0; r<10000; r++){
  //fill arrays for each object
  for(s=1; s<argc; s++){
    f=fopen(argv[s],"r");
  
    if(f==NULL){
      printf("can't open file %s\n", argv[s]);
      exit(0);
    }

    //get number of points
    fscanf(f,"%d",&numpoints);

    //fill x, y, and z arrays
    for(i=0; i<numpoints; i++){
      fscanf(f,"%lf %lf %lf", &x[s-1][i],&y[s-1][i],&z[s-1][i]);
    }

    //get number of polygons
    fscanf(f,"%d",&numpolys);

    //fill connectivity array
    for(i=0; i<numpolys; i++){
      fscanf(f,"%d",&connect[s-1][i][0]);
      for(j=1; j<=connect[s-1][i][0]; j++){
	fscanf(f,"%d",&connect[s-1][i][j]);
      }
    }

    zmax=0;
    zmin=z[s-1][0];
    for(i=0; i<numpoints; i++){
      if(z[s-1][i]>zmax){zmax=z[s-1][i];}
      if(z[s-1][i]<zmin){zmin=z[s-1][i];}
    }
    zcenter=(zmax+zmin)/2;
 
    init_array(s-1,d);
    d=d+numpolys;
    
  }

  //print_array(t,d);
  qsort(t,d,sizeof(THING),compare);
  //print_array(t,d);
 
  //draw
  //G_rgb(1,0,0);
  draw(H);

  while(a==1){
    obj=G_wait_key();
    q=G_wait_key();
    if(q==112){
      print_array(t,d);
      break;
    }
    if(q==116){
      translate(H,obj-48);
    }
    if(q==114){
      rotate(H,obj-48);
      }
    if(q==32){break;}
    }
  q=G_wait_key();
  if(q==32){break;}
}


  
//G_wait_key();

}

