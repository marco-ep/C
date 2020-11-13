 #include <FPT.h>
  #include <D3d_matrix.h>

typedef
struct {
  int objnum ;
  int polynum ;
  double dist ;
}
  THING ;

  THING all[10000];
  int numpoints[10], numpolys[10];
  double x[10][5000],y[10][5000], z[10][5000];
  int psize[10][4000];
  int connect[10][4000][20];
  double xcenter, ycenter,zcenter, xmax, ymax,zmax, xmin, ymin, zmin;
  int alln = 0;
double xL, yL, zL ;
// light model variables
double ambient, diffuse, specular ;
double colorkey ;
// inherent color values
double ir[10], ig[10], ib[10] ;
  int init_array() {
  int i;
  int j;

  alln=0;
  // for (i = 0; i < 1; i++) {
    for (j = 0; j < numpolys[i]; j++) {
      all[alln].objnum = i;
      all[alln].polynum = j;
      all[alln].dist = z[i][connect[i][j][0]];
      alln++ ;
    }  
    // }
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
double light_model(int OC, int PC) {
  double diff, spec ;
  double intensity ;
  // vectors for light
  double l[3], r[3], n[3], e[3] ;
  double ml, mr, mn, me ;
  // vectors to get normal
  double a[3], b[3] ;
  double tempdot1, tempdot2, tempdot3 ;

  diff = diffuse ;
  spec = specular ;
  
  // get normal unit vector and magnitude
  a[0] = x[OC][connect[OC][PC][1]] - x[OC][connect[OC][PC][0]] ;
  a[1] = y[OC][connect[OC][PC][1]] - y[OC][connect[OC][PC][0]] ;
  a[2] = z[OC][connect[OC][PC][1]] - z[OC][connect[OC][PC][0]] ;
  b[0] = x[OC][connect[OC][PC][2]] - x[OC][connect[OC][PC][0]] ;
  b[1] = y[OC][connect[OC][PC][2]] - y[OC][connect[OC][PC][0]] ;
  b[2] = z[OC][connect[OC][PC][2]] - z[OC][connect[OC][PC][0]] ;
  D3d_x_product(n, a, b) ;
  mn = sqrt(pow(n[0], 2) + pow(n[1], 2) + pow(n[2], 2)) ;
  n[0] /= mn ; n[1] /= mn ; n[2] /= mn ;
  
  // get light unit vector and magnitude
  l[0] = xL - x[OC][connect[OC][PC][0]] ;
  l[1] = yL - y[OC][connect[OC][PC][0]] ;
  l[2] = zL - z[OC][connect[OC][PC][0]] ;
  ml = sqrt(pow(l[0], 2) + pow(l[1], 2) + pow(l[2], 2)) ;
  l[0] /= ml ; l[1] /= ml ; l[2] /= ml ;
  
  // get eye unit vector and magnitude
  e[0] = -x[OC][connect[OC][PC][0]] ;
  e[1] = -y[OC][connect[OC][PC][0]] ;
  e[2] = -z[OC][connect[OC][PC][0]] ;
  me = sqrt(pow(e[0], 2) + pow(e[1], 2) + pow(e[2], 2)) ;
  e[0] /= me ; e[1] /= me ; e[2] /= me ;
  
  // dot product of n * l
  tempdot1 = (n[0]*l[0]) + (n[1]*l[1]) + (n[2]*l[2]) ;
  // dot product of n * e
  tempdot3 = (n[0]*e[0]) + (n[1]*e[1]) + (n[2]*e[2]) ;

  // if normal is on the opposite plane of eye
  if((tempdot1 < 0) && (tempdot3 < 0)) {
    n[0] *= -1 ; n[1] *= -1 ; n[2] *= -1 ;
    tempdot1 = (n[0]*l[0]) + (n[1]*l[1]) + (n[2]*l[2]) ;
    tempdot3 = (n[0]*e[0]) + (n[1]*e[1]) + (n[2]*e[2]) ;
  }
  
  // get reflection unit vector and magnitude
  r[0] = (2 * tempdot1) * n[0] - l[0] ;
  r[1] = (2 * tempdot1) * n[1] - l[1] ;
  r[2] = (2 * tempdot1) * n[2] - l[2] ;
  mr = sqrt(pow(r[0], 2) + pow(r[1], 2) + pow(r[2], 2)) ;
  r[0] /= mr ; r[1] /= mr ; r[2] /= mr ;
  
  // dot product of e * r
  tempdot2 = (e[0]*r[0]) + (e[1]*r[1]) + (e[2]*r[2]) ;
  
  diff *= tempdot1 ;
  spec *= pow(tempdot2, 50) ;

  // if plane is inbetween eye and light
  if((tempdot1 > 0) && (tempdot3 < 0)) { intensity = ambient ; }
  else { intensity = ambient + diff + spec ; }
  return intensity ;
}
void draw(double H,int t){
  double m[4][4], minv[4][4];
  int i,j,r,g;
  double tempx[10000],tempy[10000], tempz[10000];
  int OC, PC;
  double intensity, inpro ;
  double rr, rg, rb ;
   
  for(int i=t-1; i>=0; i--){
    r=0; g=0;
    D3d_make_identity(m);
    D3d_make_identity(minv);
    D3d_translate(m, minv, -xcenter, -ycenter, -zcenter);
    OC = all[i].objnum ;
    PC = all[i].polynum ;
    
    for(j=0; j<psize[OC][PC]; j++){
	tempx[j]=x[OC][connect[OC][PC][j]];
	tempy[j]=y[OC][connect[OC][PC][j]];
	tempz[j]=z[OC][connect[OC][PC][j]];

     }
     D3d_make_identity(m);
     D3d_make_identity(minv);
     double dz; dz=fabs(3-zmin);
     D3d_translate(m,minv,0,0,dz);
     D3d_mat_mult_points(tempx,tempy,tempz,m, tempx,tempy,tempz, numpoints[OC]);
    
     for(j=0; j<psize[OC][PC]; j++){
         tempx[j]=tempx[j]/tempz[j];
	 tempy[j]=tempy[j]/tempz[j];
	 tempz[j]=1;
      }
    D3d_scale(m,minv,400/H,400/H,0);
    D3d_translate(m,minv, 400,400,0);
    D3d_mat_mult_points(tempx,tempy,tempz,m, tempx,tempy,tempz, psize[OC][PC]);
    
    //if(onumCurr==0){r = 1;
    //}else{g = 1 ;}
  intensity = light_model(OC, PC) ;
   rr = ir[0] ;
  rg = ig[0] ;
  rb = ib[0];
  if(intensity != colorkey) {
      inpro = intensity/colorkey ;
	rr *= inpro ;
	rg *= inpro ;
	rb *= inpro ;
  }
    G_rgb(rr,rg,rb);
    G_fill_polygon(tempx,tempy,psize[OC][PC]);
    //G_rgb(0,0,0) ;
    //G_polygon(tempx,tempy,psize[OC][PC]) ; 
  
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

  printf("Enter coordinates of light source, separate with spaces:\n") ;
  scanf("%lf %lf %lf", &xL, &yL, &zL) ;
  printf("\n") ;
  printf("Enter ambient value:\n") ;
  scanf("%lf", &ambient) ;
  printf("\n") ;
  printf("Enter diffuse value:\n") ;
  scanf("%lf", &diffuse) ;  
  specular = 1 - ambient - diffuse ;
  colorkey = ambient + diffuse ;

  //for(i = 1 ; i < argc; i++) {
    printf("\n") ;
    printf("Enter inherent color values for object %d with spaces:\n", i) ;
    scanf("%lf %lf %lf", &ir[0], &ig[0], &ib[0]) ;
    //}
  
  printf("\n") ;
  printf("Light source is at (%lf, %lf, %lf)\n", xL, yL, zL) ;
  printf("Ambient:  %lf\n", ambient) ;
  printf("Diffuse:  %lf\n", diffuse) ;
  printf("Specular: %lf\n", specular) ;
  
  G_init_graphics(800,800);
  G_rgb(0,0,0);
  G_clear();
  
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
    qsort(all,t,sizeof(THING),compare);
    // scale(0);
    //scale(1);
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
      qsort (all, t, sizeof(THING), compare) ;
      draw(H,t);
		
      q=G_wait_key();
    }
   
 }
}
