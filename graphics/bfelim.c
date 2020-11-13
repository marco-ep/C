  #include <FPT.h>
  #include <D3d_matrix.h>

double xcenter, ycenter,zcenter, xmax, ymax,zmax, xmin, ymin, zmin;
void scale(double x[],double y[],double z[],int n){
    int i;

    xcenter=0; ycenter=0;zcenter=0;
    xmax=0; ymax=0; zmax=0;
    xmin=x[0]; ymin=y[0]; zmin=z[0];
    
    for(i=0; i<n; i++){
      if(x[i]>xmax){xmax=x[i];}
      if(y[i]>ymax){ymax=y[i];}
      if(z[i]>zmax){zmax=z[i];}
      if(x[i]<xmin){xmin=x[i];}
      if(y[i]<ymin){ymin=y[i];}
      if(z[i]<zmin){zmin=z[i];}
    } 
    xcenter=(xmax+xmin)/2;
    ycenter=(ymax+ymin)/2;
    zcenter=(zmax+zmin)/2;   
  }

void draw(int psize[],int connect[][20],double x[], double y[], double z[],
	  int numpolys, int numpoints,double H){
  double m[4][4], minv[4][4];
  int i,j;
  double P[3];
  double E[3],A[3],B[3];
  double dotval;
  
  G_rgb(1,0,0);
  double tempx[10000],tempy[10000], tempz[10000];
    D3d_make_identity(m);
    D3d_make_identity(minv);
    D3d_translate(m, minv, -xcenter, -ycenter, -zcenter);
    
  for(i=0; i<numpolys; i++){
    for(j=0; j<psize[i]; j++){
	tempx[j]=x[connect[i][j]];
	tempy[j]=y[connect[i][j]];
	tempz[j]=z[connect[i][j]];

	/*tempx[j]=tempx[j]/tempz[j];
	tempy[j]=tempy[j]/tempz[j];
	tempz[j]=1;

	tempx[j] =(400/H)*tempx[j]+400;
        tempy[j] =(400/H)*tempy[j]+400;*/
    }
     D3d_make_identity(m);
    D3d_make_identity(minv);
    double dz; dz=fabs(3-zmin);
    D3d_translate(m,minv,0,0,dz);
    D3d_mat_mult_points(tempx,tempy,tempz,m, tempx,tempy,tempz, numpoints);
    
    //
    A[0]=tempx[1]-tempx[0];
    A[1]=tempy[1]-tempy[0];
    A[2]=tempz[1]-tempz[0];
    A[3]=1;

    B[0]=tempx[2]-tempx[0];
    B[1]=tempy[2]-tempy[0];
    B[2]=tempz[2]-tempz[0];
    B[3]=1;

    D3d_x_product(P,A,B);

    E[0]=-tempx[0];
    E[1]=-tempy[0];
    E[2]=-tempz[0];
    E[3]=1;
    
    dotval=((E[0]*P[0])+(E[1]*P[1])+(E[2]*P[2]));
    //
    
    for(j=0; j<psize[i]; j++){
      tempx[j]=tempx[j]/tempz[j];
      tempy[j]=tempy[j]/tempz[j];
      tempz[j]=1;
      }
    D3d_scale(m,minv,400/H,400/H,0);
    D3d_translate(m,minv, 400,400,0);
    D3d_mat_mult_points(tempx,tempy,tempz,m, tempx,tempy,tempz, psize[i]);
    
    if(dotval*-1>0){
    G_polygon(tempx,tempy,psize[i]);
    }
  
  }
}

void center(double x[], double y[], double z[], int numpoints){
  double m[4][4], minv[4][4];
  double xs=0,ys=0, zs=0;
  int i;
  D3d_make_identity(m);
    D3d_make_identity(minv);
    for(i=0;i<numpoints;i++){xs+=x[i]; ys+=y[i];zs+=z[i];}
    xcenter=xs/numpoints;
    ycenter=ys/numpoints;
    zcenter=zs/numpoints;
  
}

  int main(int argc, char **argv)
{
 
  int numpoints, numpolys;
  double x[5000],y[5000], z[5000];
  int psize[4000];
  int connect[4000][20];
  int g,i,j,q;
 double halfangle =45*M_PI/180;
double H= tan(halfangle);
  g=0;
  G_init_graphics(800,800);
 
  if(argc > 1){
    for(int g=1;g<argc;g++){
       G_rgb(0,0,0);
  G_clear();
    FILE *fp;
    fp=fopen(argv[g],"r");
    if(fp==NULL){
      printf("can't open file\n");
      exit(0);
    }

    fscanf(fp,"%d",&numpoints);
    for(i=0; i<numpoints; i++){
      fscanf(fp,"%lf %lf %lf", &x[i],&y[i],&z[i]);
    }
    fscanf(fp,"%d",&numpolys);
    for(i=0; i<numpolys; i++){
      fscanf(fp,"%d",&psize[i]);    
      for(j=0; j<psize[i]; j++){
	fscanf(fp,"%d",&connect[i][j]); 
      }
    }

    
  //draw
    scale(x,y,z,numpoints);
    draw(psize,connect,x,y,z,numpolys,numpoints,H);
    q=G_wait_key();
    //rotate
    double e[4][4], einv[4][4];
    double thing = 3*(M_PI)/180;
    double thingr = thing;
    while(q!='q'){
      
      D3d_make_identity(e);
      D3d_make_identity(einv);
      center(x,y,z,numpoints);
      G_rgb(0,0,0);
      G_clear();
      if(q=='c'){thingr = -thingr;}
      if(q=='x'){
	D3d_translate(e, einv, -xcenter,-ycenter,-zcenter);
	D3d_rotate_x(e, einv, thingr);
	D3d_translate(e, einv, xcenter,ycenter,zcenter);
      }
      if(q=='z'){
	D3d_translate(e, einv, -xcenter,-ycenter,-zcenter);
	D3d_rotate_z(e, einv, thingr);
	D3d_translate(e, einv, xcenter,ycenter,zcenter);
      }
      if(q=='y'){
	D3d_translate(e, einv, -xcenter,-ycenter,-zcenter);
	D3d_rotate_y(e, einv, thingr);
	D3d_translate(e, einv, xcenter,ycenter,zcenter);
      }
      if(q=='e'){
	D3d_translate(e, einv, thing, 0, 0);
      }
      if(q=='d'){
	D3d_translate(e, einv, -thing, 0, 0);
      }
      if(q=='w'){
	D3d_translate(e, einv, 0, 0, thing);
      }
      if(q=='s'){
	D3d_translate(e, einv, 0, 0, -thing);
      }
      if(q=='7'){
	D3d_translate(e, einv, 0,thing,0);
      }
      if(q=='u'){
	D3d_translate(e, einv, 0,-thing,0);
      }
      
     
        D3d_mat_mult_points(x,y,z,e, x,y,z,numpoints); //important
	draw(psize,connect,x,y,z,numpolys, numpoints,H);
	
      q=G_wait_key();
    }
  }
 }
}
