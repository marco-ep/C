#include <math.h>
#include <stdio.h>
#include <stdlib.h>
int main()
{

  double pennies, dollars, quarters, dimes, nickels ;
  printf("How many pennies do you have?") ;
  scanf("%lf", &pennies) ;
  printf("%lf", pennies);
  dollars=floor(pennies/100);
  pennies=fmod(pennies,100);
  quarters=floor(pennies/25);
  pennies=fmod(pennies,25);
  dimes=floor(pennies/10);
  pennies=fmod(pennies,10);
  nickels=floor(pennies/5);
  pennies=fmod(pennies,5);

  printf("%lf dollars and %lf quarters and %lf dimesand %lf nickels and %lf pennies\n", dollars,quarters,dimes,nickels, pennies) ;


}
