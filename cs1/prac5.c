#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main()
{

double n, i, c, total ;
printf("enter whole number greater than 1\n") ;
scanf("%lf", &n) ;

i = 2 ;
total = 0 ;

while (i <= n) {
if (fmod(i, 2) == 0) {
	c = sqrt(i) ;
} else {
	c = sqrt(1/i) ;
}

total = total + c ;
	i++ ;
}

printf("%lf\n", total) ;
}