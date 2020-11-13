#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main()
{

double n, s, sign, i ;

printf("enter a multiple of four\n") ;
scanf("%lf", &n) ;

s = 0 ;
sign = 1 ;
i = 4 ;

while (i <= n) {
	s = s + sign * (i - 3) * (i - 2) / ((i - 1) + i) ;
	sign = -sign ;
	i = i + 4 ;
}
 printf("the answer is %lf\n", s) ;
}