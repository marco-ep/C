#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main(){
  double co[1000] ;
  double fac[1000], ffac[1000], lfac[1000] ;
  int deg, i, j, k ;
  int fi, li, maxj, maxk ;
  double first, last, n, sum, tdeg ;

  printf("enter degree \n") ;
  scanf("%d", &deg) ;

  i = 0 ;
  printf("enter coefficients \n") ;
  while (i <= deg) {
    scanf("%lf", &co[i]) ;
    i++ ;
  }
  last = fabs(co[i - 1]);
  first = fabs(co[0]);
  
  i = 0;
  n = 1;
  //first factors
  while (n <= first) {
    if (fmod(first, n) == 0) {
      ffac[i] = n;
      i++;
    }
    
    n++;
  }
  fi = i - 1;

  i = 0;
  n = 1;
  //last factors
  while (n <= last) {
    if (fmod(last, n) == 0) {
      lfac[i] = n;
      i++;
    }
    
    n++;
  }
  li = i - 1;

  j = 0;
  k = 0;
  //actual potential factors
  while (j <= li) {
    i = 0;

    while (i <= fi) {
      fac[k] = lfac[j]/ffac[i];
      i++;
      k++;
    }
    
    j++;
  }

  maxk = k - 1;

  i = 0;
  k = 0;
  //positive
  while (k <= maxk) {
    sum = 0;
    tdeg = deg;
    i = 0;

    while (tdeg >= 0) {
      sum = sum + co[i]*pow((fac[k]), tdeg);

      i++;
      tdeg--;
    }

    if ((sum >= -0.00001) && (sum <= 0.00001)) {
      printf("root: %lf\n", fac[k]);
    }
    
    k++;
  }

  i = 0;
  k = 0;
  //negative
  while (k <= maxk) {
    sum = 0;
    tdeg = deg;
    i = 0;

    while (tdeg >= 0) {
      sum = sum + co[i]*pow((-fac[k]), tdeg);

      i++;
      tdeg--;
    }

    if ((sum >= -0.00001) && (sum <= 0.00001)) {
      printf("root: %lf\n", -fac[k]);
    }
    
    k++;
  }


}
