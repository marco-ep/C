 #include <FPT.h>
  #include <D3d_matrix.h>

typedef
struct {
  int objnum ;
  int polynum ;
  double dist ;
}
  THING ;

  THING allp[10000];
  int numpoints[10], numpolys[10];
  double x[10][5000],y[10][5000], z[10][5000];
  int psize[10][4000];
  int connect[10][4000][20];
  double xcenter, ycenter,zcenter, xmax, ymax,zmax, xmin, ymin, zmin;
  int allpn = 0;
  int arn;

  int init_array() {
  int i;
  int j;

  allpn=0;
  for (i = 0; i < arn; i++) {
    for (j = 0; j < numpolys[i]; j++) {
      allp[allpn].objnum = i;
      allp[allpn].polynum = j;
      allp[allpn].dist = z[i][connect[i][j][0]];
      allpn++ ;
    }  
  }
}

void scale(int oc){
  int i;
    xcenter=0; ycenter=0;zcenter=0;
    xmax=0; ymax=0; zmax=0;
    xmin=x[oc][0]; ymin=y[oc][0]; zmin=z[oc][0];
    
    for(i=0; i<numpoints[oc]; i++){
      if(x[oc][i]>xmax){xmax=x[oc][i];}
      if(y[oc][i]>ymax){ymax=y[oc][i];}
      if(z[oc][i]>zmax){zmax=z[oc][i];}
      if(x[oc][i]<xmin){xmin=x[oc][i];}
      if(y[oc][i]<ymin){ymin=y[oc][i];}
      if(z[oc][i]<zmin){zmin=z[oc][i];}
    } 
    xcenter=(xmax+xmin)/2;
    ycenter=(ymax+ymin)/2;
    zcenter=(zmax+zmin)/2;   
  
}
int compare (const void *p, const void *q)
{
  THING *a, *b ;

  a = (THING*)p ;
  b = (THING*)q ;

  if  (((*a).dist) < ((*b).dist)) return -1 ;
  else if (((*a).dist) > ((*b).dist)) return 1 ;
  else return 0 ;
}

void draw(double H,int t){
  double m[4][4], minv[4][4];
  int i,j,r,g;
  double tempx[10000],tempy[10000], tempz[10000];
  int Ocurr, Pcurr; 
   
  for(int i=t-1; i>=0; i--){
    r=0; g=0;
    D3d_make_identity(m);
    D3d_make_identity(minv);
    D3d_translate(m, minv, -xcenter, -ycenter, -zcenter);
    Ocurr = allp[i].objnum ;
    Pcurr = allp[i].polynum ;
    
    for(j=0; j<psize[Ocurr][Pcurr]; j++){
	tempx[j]=x[Ocurr][connect[Ocurr][Pcurr][j]];
	tempy[j]=y[Ocurr][connect[Ocurr][Pcurr][j]];
	tempz[j]=z[Ocurr][connect[Ocurr][Pcurr][j]];

     }
     D3d_make_identity(m);
     D3d_make_identity(minv);
     double dz; dz=fabs(3-zmin);
     D3d_translate(m,minv,0,0,dz);
     D3d_mat_mult_points(tempx,tempy,tempz,m, tempx,tempy,tempz, numpoints[Ocurr]);
    
     for(j=0; j<psize[Ocurr][Pcurr]; j++){
         tempx[j]=tempx[j]/tempz[j];
	 tempy[j]=tempy[j]/tempz[j];
	 tempz[j]=1;
      }
    D3d_scale(m,minv,400/H,400/H,0);
    D3d_translate(m,minv, 400,400,0);
    D3d_mat_mult_points(tempx,tempy,tempz,m, tempx,tempy,tempz, psize[Ocurr][Pcurr]);
    
    if(Ocurr==0){r = 1;
    }else{g = 1 ;}
  
    G_rgb(r,g,0);
    G_fill_polygon(tempx,tempy,psize[Ocurr][Pcurr]);
    G_rgb(0,0,0) ;
    G_polygon(tempx,tempy,psize[Ocurr][Pcurr]) ; 
  
  }
}

void center(int oc){//modify with 2d arrays and oc
  double m[4][4], minv[4][4];
  double xs=0,ys=0, zs=0;
  int i;
  D3d_make_identity(m);
    D3d_make_identity(minv);
    for(i=0;i<numpoints[oc];i++){xs+=x[oc][i]; ys+=y[oc][i];zs+=z[oc][i];}
    xcenter=xs/numpoints[oc];
    ycenter=ys/numpoints[oc];
    zcenter=zs/numpoints[oc];
  
}

  int main(int argc, char **argv)
{
  int g,i,j,q,t;
  double halfangle =45*M_PI/180;
  double H= tan(halfangle);
  G_init_graphics(800,800);
  G_rgb(0,0,0);
  G_clear();
  arn=argc;
  if(argc > 1){
    for(g=1;g<argc;g++){
        int onum = g-1;
	FILE *fp;
	fp=fopen(argv[g],"r");
	if(fp==NULL){
	  printf("can't open file\n");
	  exit(0);
	}

	fscanf(fp,"%d",&numpoints[onum]);
	for(i=0; i<numpoints[onum]; i++){
	    fscanf(fp,"%lf %lf %lf", &x[onum][i],&y[onum][i],&z[onum][i]);
	}
	fscanf(fp,"%d",&numpolys[onum]);
	for(i=0; i<numpolys[onum]; i++){
	    fscanf(fp,"%d",&psize[onum][i]);    
	    for(j=0; j<psize[onum][i]; j++){
	        fscanf(fp,"%d",&connect[onum][i][j]); 
	    }
	}
    t+=numpolys[onum];
    }
    
    init_array() ;
    qsort(allp,t,sizeof(THING),compare);
    scale(0);
    scale(1);
    draw(H,t);
  
    q=G_wait_key();
    //rotate
    double e[4][4], einv[4][4];
    double thing = 3*(M_PI)/180;
    double thingr = thing;
    int onum;
    while(q!='q'){
      //if(q >= '0' && q <= '9'){
      //	g=q-'0';
      //  }else{
      if(q == '1'){onum = 0;}
      if(q == '2'){onum=1;}
      if(q == '3'){onum=2;}
      
      D3d_make_identity(e);
      D3d_make_identity(einv);
      center(onum); //onlygive the current objects stuff not both
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
  
      D3d_mat_mult_points(x[onum],y[onum],z[onum],e, x[onum],y[onum],z[onum],numpoints[onum]); //important
      init_array() ;	
      qsort (allp, t, sizeof(THING), compare) ;
      draw(H,t);
		
      q=G_wait_key();
    }
   
 }
}
