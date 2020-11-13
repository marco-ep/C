#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main()
{

double n, f, a, b, i ;

printf("enter a number\n") ;
scanf("%lf", &n) ; 

if ((n == 1 || n == 2)) {
	f = 1 ;
} else {
	i = 3 ;
	a = 1 ; b = 1 ;
while (i <= n ) {
f = a - b ;

a = b ; 
b = f ;
	i++ ;
}
}
printf("%lf\n", f) ;
}