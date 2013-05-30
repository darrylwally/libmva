/*! \file mvconstants.c
  \brief Contains various constants used in mvlib.
  \author Darryl Wallace - wallacdj@gmail.com, darryl@wallynet.ca
  \date 2012

  THIS CODE IS COPYRIGHT DARRYL WALLACE (c) 2012.
  UNAUTHORIZED DISTRIBUTION, USE, OR VIEWING THE CONTENTS OF THIS FILE
  IS STRICTLY PROHIBITED.
  IF YOU HAVE RECEIVED THIS CODE WITHOUT AUTHORIZATION IT MUST BE DESTROYED
  IMMEDIATELY
  */

#include "mvconstants.h"

double mvNaN()
{
    unsigned long long __nand = 0x7ff8000000000000;
    return *(double *) &__nand;
}

float mvNaNf()
{
    unsigned int __nanf = 0x7fc00000;
    return *(float *)&__nanf;
}

double mvInf()
{
    unsigned long long __infd = 0x7ff0000000000000;
    return *(double *)&__infd;
}

double mvNegInf()
{
    unsigned long long __infd = 0xfff0000000000000;
    return *(double *)&__infd;
}

float mvInff()
{
    unsigned int __inff = 0x7f800000;
    return *(float *)&__inff;
}

float mvNegInff()
{
    unsigned int __inff = 0xff800000;
    return *(float *)&__inff;
}
