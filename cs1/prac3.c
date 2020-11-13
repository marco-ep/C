#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main()
{

double n, i, a, sign, den, tot ;

printf("enter positive integer\n") ;
scanf("%lf", &n) ;

i = 1 ;
sign = 1 ;
den = 0 ;
tot = 0 ;

while (i <= n) {
den = den + 1/i ;
a = sign * (i/den) ;
tot = tot + a ;

sign = -sign ;
	i++ ;
}
printf("%lf\n", tot) ;
}