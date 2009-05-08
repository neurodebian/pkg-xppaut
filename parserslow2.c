#include <stdlib.h> 
/*  lots of changes - good luck */
/* NOTE!!!  I am changing the parser a great deal here
   I am making parameters 200-300
   kernels are 1000
   variables are 1100 - 2100 for now
   shift/delay adds 3000 to them
   currently 600 is the maximum to get it up to eg 700
   delete all .o files  change xpplim.h by adding 100 to MAXODE,MAXODE1
   add 50 to MAXPRIMEVAR
   better make worksize bigger - I am not sure how much cant hurt to
   be too big - 300000 = 2.4MB memory which is not too much these days
   change  is_uvar - change 17 to 18   
           is_var  -   "        "

The following are the currently used indices
0-99  functions of 1 variable like sin, cos, etc
100-199  fixed functions of 2 variables like + -, etc
200-399  parameters


600-699   networks
700-799   lookup tables
800-899   dummy arguments for functions
900-999   special stuff
1000-1099 volterra kernels
10000-11000  variables

2400-2499  user defined functions

3200-3300  shifted constants
20000-21000  shifted variables

NOTE - I will shortly make user-defined functions 2400 instead of 400
so that I have the block from 200-599 for parameter

*/
#ifndef WCTYPE
#include <ctype.h>
#else
#include <wctype.h>
#endif

#include <math.h>
/* #include <malloc.h> */
#include <stdio.h>
#include <string.h>

#include "parserslow.h"
#include "xpplim.h"
#include "getvar.h"
#define MAXEXPLEN 1024
#define THOUS 10000
#define DOUB_EPS 2.23E-15 
#define POP stack[--stack_pointer]
double zippy;
#define PUSH(a) zippy=(a); stack[stack_pointer++]=zippy;



#define DFNORMAL 1
#define DFFP 2
#define DFSTAB 3

#define COM(a) my_symb[toklist[(a)]].com
int             ERROUT;
extern int DelayFlag;
int NDELAYS=0;
/*double pow2(); */
double get_delay();
double delay_stab_eval();
double lookup(),network_value();
double atof(),poidev();
double ndrand48();
double ker_val();
double hom_bcs();
double BoxMuller;
int BoxMullerFlag=0;
int RandSeed=12345678;
extern int newseed;
extern int del_stab_flag;

double CurrentIndex=0;
int SumIndex=1;
/*************************
  RPN COMPILER           *
**************************/



/*****************************
*      PARSER.C              *
*
*
*     parses any algebraic expression
*     and converts to an integer array
*     to be interpreted by the rpe_val
*     function.
*
*     the main data structure is a contiguous
*     list of symbols with their priorities
*     and their symbol value
*
*     on the first pass, the expression is converted to
*     a list of integers without any checking except for
*     valid symbols and for numbers
*
*     next this list of integers is converted to an RPN expression
*     for evaluation.
*
*
*  6/95  stuff added to add names to namelist without compilation
*************************************************************/

/* INIT_RPN    */
/*
int matherr(e)
     struct exception *e;
{
char *c;
	switch (e->type) {
	    case DOMAIN:    c = "domain error"; break;
	    case SING  :    c = "argument singularity"; break;
	    case OVERFLOW:  c = "overflow range"; break;
	    case UNDERFLOW: c = "underflow range"; break;
	    default:		c = "(unknown error"; break;
	} 
	fprintf(stderr, "math exception : %d\n", e->type);
	fprintf(stderr, "    name : %s\n", e->name);
	fprintf(stderr, "    arg 1: %e\n", e->arg1);
	fprintf(stderr, "    arg 2: %e\n", e->arg2);
	fprintf(stderr, "    ret  : %e\n", e->retval);
	return 1;
}
*/


init_rpn()
{
	int             i;
	ERROUT = 1;
	NCON = 0;
	NFUN = 0;
	NVAR = 0;
	NKernel=0;
    
	MaxPoints=4000;
	NSYM = STDSYM;
	two_args();
	one_arg();
	add_con("PI", M_PI);
        add_con("I'",0.0);
	/*  add_con("c___1",0.0);
        add_con("c___2",0.0);
        add_con("c___3",0.0); */
	SumIndex=NCON-1;
	init_table();
	if (newseed==1) RandSeed=time(0);
	nsrand48(RandSeed);
}

 /*  FREE_UFUNS   */

free_ufuns()
{
 int i;
 for(i=0;i<NFUN;i++)
 {
  free(ufun[i]);
  free(ufun_def[i]);
 }
}
duplicate_name(junk)
     char *junk;
{
  int i;
  find_name(junk,&i);
  if(i>=0){
    if(ERROUT)printf("%s is a duplicate name\n",junk);
    return(1);
  }
  return(0);
}
  
/*  ADD_CONSTANT   */

add_constant(junk)
char *junk;
{
 int len;
 char string[100];
 
 if(duplicate_name(junk)==1)return(1);
 if(NCON>=397)
 {
  if(ERROUT)printf("too many constants !!\n");
  return(1);
 }
 convert(junk,string);
 len=strlen(string);
  if(len<1){
   printf("Empty parameter - remove spaces\n");
   return 1;
 }
 if(len>MXLEN)len=MXLEN;
 strncpy(my_symb[NSYM].name,string,len);
 my_symb[NSYM].name[len]='\0';
 my_symb[NSYM].len=len;
 my_symb[NSYM].pri=10;
 my_symb[NSYM].arg=0;
 my_symb[NSYM].com=200+NCON-1;
 NSYM++;
 return(0);
}


get_var_index(name)
     char *name;
{
  int in=-1;
  int type,com;
  find_name(name,&type);
  if(type<0)return -1;
  com=my_symb[type].com;
  if(is_uvar(com))
  {
      return(com-10000);
  }
  return(-1);
}


/* GET_TYPE   */


get_type(index)
int index;
{
return(my_symb[index].com);
}

/*   ADD_CON      */

add_con(name,value)
char *name;
double value;
{
 if(NCON>=397)
 {
  if(ERROUT)printf("too many constants !!\n");
  return(1);
 }
 constants[NCON]=value;
 NCON++;
 return(add_constant(name));
}

add_kernel(name,mu,expr)
     char *name,*expr;
     double mu;
{
  char string[100];
  char bexpr[MAXEXPLEN],kexpr[MAXEXPLEN];
  int len,i,in=-1;
  if(duplicate_name(name)==1)return(1);
  if(NKernel==MAXKER){
    printf("Too many kernels..\n");
    return(1);
  }
  if(mu<0||mu>=1.0){
    printf(" mu must lie in [0,1.0) \n");
    return(1);
  }
  convert(name,string);
  len=strlen(string);
  if(len>MXLEN)len=MXLEN;
  strncpy(my_symb[NSYM].name,string,len);
  my_symb[NSYM].name[len]='\0';
  my_symb[NSYM].len=len;
  my_symb[NSYM].pri=10;
  my_symb[NSYM].arg=0;
  my_symb[NSYM].com=1000+NKernel;
  kernel[NKernel].k_n1=0.0;
  kernel[NKernel].mu=mu;
  kernel[NKernel].k_n=0.0;
  kernel[NKernel].k_n1=0.0;
  kernel[NKernel].flag=0;
  for(i=0;i<strlen(expr);i++)
    if(expr[i]=='#')in=i;
  if(in==0||in==(strlen(expr)-1)){
    printf("Illegal use of convolution...\n");
    return(1);
  }
  if(in>0){
    kernel[NKernel].flag=CONV;
    kernel[NKernel].expr=(char *)malloc(strlen(expr)+2-in);
    kernel[NKernel].kerexpr=(char *)malloc(in+1);
    for(i=0;i<in;i++)kernel[NKernel].kerexpr[i]=expr[i];
    kernel[NKernel].kerexpr[in]=0;
    for(i=in+1;i<strlen(expr);i++)kernel[NKernel].expr[i-in-1]=expr[i];
    kernel[NKernel].expr[strlen(expr)-in-1]=0;
    printf("Convolving %s with %s\n",
	   kernel[NKernel].kerexpr,kernel[NKernel].expr);
  }
  else {
    kernel[NKernel].expr=(char *)malloc(strlen(expr)+2);
    strcpy(kernel[NKernel].expr,expr);
  }
  NSYM++;
  NKernel++;
  return(0);
}
  


/*  ADD_VAR          */

add_var(junk, value)
char *junk;
double value;
{
 char string[100];
 int len;
 /*   printf(" variable - %s \n",junk); */
 if(duplicate_name(junk)==1)return(1);
 if(NVAR>=MAXODE1)
 {
  if(ERROUT)printf("too many variables !!\n");
  return(1);
 }
 convert(junk,string);
 len=strlen(string);
 if(len>MXLEN)len=MXLEN;
 strncpy(my_symb[NSYM].name,string,len);
 my_symb[NSYM].name[len]='\0';
 my_symb[NSYM].len=len;
 my_symb[NSYM].pri=10;
 my_symb[NSYM].arg=0;
 my_symb[NSYM].com=10000+NVAR;
 NSYM++;
 variables[NVAR]=value;
 NVAR++;
 return(0);
}

/* ADD_EXPR   */

add_expr(expr,command, length )
char *expr;
int *command, *length;
{
 char dest[1024];
 int my_token[1024];
 int err,i;
 convert(expr,dest);
/* printf(" Making token ...\n");  */
 err=make_toks(dest,my_token); 

/*  i=0;
  while(1){
  printf(" %d %d \n",i,my_token[i]);
  if(my_token[i]==ENDTOK)break;
  i++;
} */ 
 if(err!=0)return(1);
 err = alg_to_rpn(my_token,command);
 if(err!=0)return(1);
  i=0;
   while(command[i]!=ENDEXP)i++;
   *length=i+1;
 /*  for(i=0;i<*length;i++)printf("%d \n",command[i]);  */
   return(0);
}

add_net_name(index,name)
     int index;
     char *name;
{
  char string[50];
  int len=strlen(name);
  printf(" Adding net %s %d \n",name,index);
  if(duplicate_name(name)==1)return(1);  
  convert(name,string);
  if(len>MXLEN)len=MXLEN;
  strncpy(my_symb[NSYM].name,string,len);
  my_symb[NSYM].name[len]='\0';
  my_symb[NSYM].len=len;
  my_symb[NSYM].pri=10;
  my_symb[NSYM].arg=1;
  my_symb[NSYM].com=600+index;
  NSYM++;
  return(0);
   
}

add_vect_name(index,name)
     int index;
     char *name;
{
  char string[50];
  int len=strlen(name);
  printf(" Adding vector %s %d \n",name,index);
  if(duplicate_name(name)==1)return(1);  
  convert(name,string);
  if(len>MXLEN)len=MXLEN;
  strncpy(my_symb[NSYM].name,string,len);
  my_symb[NSYM].name[len]='\0';
  my_symb[NSYM].len=len;
  my_symb[NSYM].pri=10;
  my_symb[NSYM].arg=1;
  my_symb[NSYM].com=VECT_ROOT+index;
  NSYM++;
  return(0);
   
}


/* ADD LOOKUP TABLE   */

add_2d_table(name,file)
     char *name,*file;
{
 printf(" TWO D NOT HERE YET \n");
 return(1);
}

add_file_table(index,file)
     char *file;
     int index;
{
 
  if(load_table(file,index)==0)
    {
      if(ERROUT)printf("Problem with creating table !!\n");
       return(1);
    }
 
    return(0);
}

add_table_name(index,name)
     char *name;
     int index;
{
     char string[50];
     int len=strlen(name);
     if(duplicate_name(name)==1)return(1);  
     convert(name,string);
     if(len>MXLEN)len=MXLEN;
     strncpy(my_symb[NSYM].name,string,len);
     my_symb[NSYM].name[len]='\0';
     my_symb[NSYM].len=len;
     my_symb[NSYM].pri=10;
     my_symb[NSYM].arg=1;
     my_symb[NSYM].com=700+index;
     set_table_name(name,index);
     NSYM++;
     return(0);
   }
/* ADD LOOKUP TABLE   */


add_form_table(index,nn,xlo,xhi,formula)
     char *formula;
     double xlo,xhi;
     int nn;
     int index;
{
 
 
  if(create_fun_table(nn,xlo,xhi,formula,index)==0)
    {
      if(ERROUT)printf("Problem with creating table !!\n");
       return(1);
    }
    return(0);
}
    

set_old_arg_names(narg)
     int narg;
{
  int i;
  for(i=0;i<narg;i++){
    sprintf(my_symb[FIRST_ARG+i].name,"ARG%d",i+1);
    my_symb[FIRST_ARG+i].len=4;
  }
}

set_new_arg_names(narg,args)
     char args[10][11];
     int narg;
{
  int i;
  for(i=0;i<narg;i++){
    strcpy(my_symb[FIRST_ARG+i].name,args[i]);
    my_symb[FIRST_ARG+i].len=strlen(args[i]);
 }
}


/* NEW ADD_FUN for new form_ode code  */

add_ufun_name(name,index,narg)
     char *name;
     int index,narg;
{
  char string[50];
  int len=strlen(name);
 if(duplicate_name(name)==1)return(1);
 if(index>=MAXUFUN)
 {
  if(ERROUT)printf("too many functions !!\n");
  return(1);
 }
  printf(" Added user fun %s \n",name);
  convert(name,string);
  if(len>MXLEN)len=MXLEN;
  strncpy(my_symb[NSYM].name,string,len);
  my_symb[NSYM].name[len]='\0';
  my_symb[NSYM].len=len;
  my_symb[NSYM].pri=10;
  my_symb[NSYM].arg=narg;
  my_symb[NSYM].com=2400+index;
  NSYM++;
  strcpy(ufun_names[index],name);
  return (0);
}
fixup_endfun(u,l,narg)
     int *u;
     int l,narg;
{
 u[l-1]=ENDFUN;
 u[l]=narg;
 u[l+1]=ENDEXP;
}

add_ufun_new(index,narg,rhs,args)
     char *rhs,args[MAXARG][11];
     int narg,index;
{
  
  int i,l;
  int end;
   if(narg>MAXARG){
    printf("Maximal arguments exceeded \n");
    return(1);
  }
  if((ufun[index]=(int *)malloc(1024))==NULL)
    {
      if(ERROUT)printf("not enough memory!!\n");
      return(1);
    }
  if((ufun_def[index]=(char *)malloc(MAXEXPLEN))==NULL)
    {
      if(ERROUT)printf("not enough memory!!\n");
      return(1);
    }
  ufun_arg[index].narg=narg;
  for(i=0;i<narg;i++)
    strcpy(ufun_arg[index].args[i],args[i]);
  set_new_arg_names(narg,args);
  if(add_expr(rhs,ufun[index],&end)==0)
    {
      
      ufun[index][end-1]=ENDFUN;
      ufun[index][end]=narg;
      ufun[index][end+1]=ENDEXP;
      strcpy(ufun_def[index],rhs);
      l=strlen(ufun_def[index]);
      ufun_def[index][l]=0;
      narg_fun[index]=narg;
      set_old_arg_names(narg);
      return(0);
    } 
  
  set_old_arg_names(narg);
  if(ERROUT)printf(" ERROR IN FUNCTION DEFINITION\n");
  return(1);
}

/* ADD_UFUN   */

add_ufun(junk,expr,narg)
char *junk, *expr;
int narg;
{
 char string[50];
 int i,l;
 int end;
 int len=strlen(junk);

 if(duplicate_name(junk)==1)return(1);
 if(NFUN>=MAXUFUN)
 {
  if(ERROUT)printf("too many functions !!\n");
  return(1);
 }
 if((ufun[NFUN]=(int *)malloc(1024))==NULL)
 {
  if(ERROUT)printf("not enough memory!!\n");
  return(1);
 }
 if((ufun_def[NFUN]=(char *)malloc(MAXEXPLEN))==NULL)
 {
  if(ERROUT)printf("not enough memory!!\n");
  return(1);
 }

 convert(junk,string);
 if(add_expr(expr,ufun[NFUN],&end)==0)
 {
  if(len>MXLEN)len=MXLEN;
  strncpy(my_symb[NSYM].name,string,len);
  my_symb[NSYM].name[len]='\0';
  my_symb[NSYM].len=len;
  my_symb[NSYM].pri=10;
  my_symb[NSYM].arg=narg;
  my_symb[NSYM].com=2400+NFUN;
  NSYM++;
  ufun[NFUN][end-1]=ENDFUN;
  ufun[NFUN][end]=narg;
  ufun[NFUN][end+1]=ENDEXP;
  strcpy(ufun_def[NFUN],expr);
  l=strlen(ufun_def[NFUN]);
  ufun_def[NFUN][l-1]=0;
  strcpy(ufun_names[NFUN],junk);
  narg_fun[NFUN]=narg;
  for(i=0;i<narg;i++){
    sprintf(ufun_arg[NFUN].args[i],"ARG%d",i+1);
  }
  NFUN++;
  return(0);
 }
       if(ERROUT)printf(" ERROR IN FUNCTION DEFINITION\n");
       return(1);
}


check_num(tok,value)
int *tok;
double value;
{
 int bob,in,m,i;
 for(i=0;i<NSYM;i++){
	
	if(strncmp(my_symb[i].name,"NUM##",5)==0){
	bob=my_symb[i].com;
	in=bob%100;
	m=bob/100;
	if(m==10)in+=100;
	if(constants[in]==value){
	*tok=i;
	return(1);
	}
	}
 }
 return(0);
}
	



/* is_ufun         */

is_ufun( x)
int x;
{
 if((x/100)==UFUN)return(1);
 else return(0);
}

/* IS_UCON        */

is_ucon( x)
int x;
{
 int xx=x/100;
 if(xx>1&&xx<6)return(1);
 else return(0);
}

/* IS_UVAR       */

is_uvar(x)
int x;
{
 int y=x/100;
 if((y>=100)&&(y<200))return(1);
 else return(0);
}

isvar(y)
     int y;
{
  if((y>=100)&&(y<200))return(1);
  return 0;
}
iscnst(y)
     int y;
{
  if(y>1&&y<6)return(1);
  return 0;
}
isker(y)
     int y;
{
  if(y==10)return(1);
  return 0;
}

is_kernel(x)
int x;
{
  if((x/100)==10)return(1);
  else return(0);
}
is_lookup(x)
int x;
{
 if((x/100)==7)return(1);
 else return(0);
}

find_lookup(name)
     char *name;
{
 int index,com;
 find_name(name,&index);
  if(index==-1)return(-1);
  com=my_symb[index].com;
  if(is_lookup(com))return(com%100);
  return(-1);
}
 

/* FIND_NAME    */

find_name(string, index)
 char *string;
 int *index;
{
  char junk[100];
  int i,len;
  convert(string,junk);
  len=strlen(junk);
  for(i=0;i<NSYM;i++)
  {
   if(len==my_symb[i].len)
    if(strncmp(my_symb[i].name,junk,len)==0)break;
  }
   if(i<NSYM)
    *index=i;
   else *index=-1;
}


get_param_index(name)
     char *name;
{
 int type,com;
  find_name(name,&type);
  if(type<0)return(-1);
  com=my_symb[type].com;
  if(is_ucon(com))
  {
      return(com-200);

  }
    return(-1);
}

/* GET_VAL   */

get_val(name,value)
char *name;
double *value;
{
  int type,com;
  *value=0.0;
  find_name(name,&type);
  if(type<0)return(0);
  com=my_symb[type].com;
  if(is_ucon(com))
  {
   *value=constants[com-200];
   return(1);
  }
  if(is_uvar(com))
  {
      *value=variables[com-10000];
   return(1);
  }
  return(0);
}

/* SET_VAL         */

set_val(name, value)
char *name;
double value;
{
  int type,com;
  find_name(name,&type);
  if(type<0)return(0);
  com=my_symb[type].com;
  if(is_ucon(com))
  {
         constants[com-200]=value;
   
    return(1);
  }
  if(is_uvar(com))
  {

      variables[com-10000]=value;
    return(1);
  }
  return(0);
}



set_ivar(i, value)
int i;
double value;
{
 SETVAR(i,value);
}

double get_ivar(i)
int i;
{       	 return(GETVAR(i));
}




alg_to_rpn(toklist,command)
int *toklist,*command;
{
  int tokstak[500],comptr=0,tokptr=0,lstptr=0,temp;
  int ncomma=0,delflag=0,deltok;
  int loopstk[100];
  int lptr=0;
  int nif=0,nthen=0,nelse=0;
  int newtok,oldtok,zip;
  int i,my_com,my_arg,jmp;
  char ch;
  int jj;
  tokstak[0]=STARTTOK;
  tokptr=1;
  oldtok=STARTTOK;
  while(1)
         {
 getnew:
          newtok=toklist[lstptr++];
   /*    for(zip=0;zip<tokptr;zip++)
	    printf("%d %d\n",zip,tokstak[zip]);  */
/*        check for delay symbol             */
          if(newtok==DELSYM)
	  {
           temp=my_symb[toklist[lstptr+1]].com;
   /* !! */   if(is_uvar(temp))
	   {
	    my_symb[LASTTOK].com=temp+10000; /* create a temporary sybol */
            NDELAYS++;
           toklist[lstptr+1]=LASTTOK;
	  	
	    my_symb[LASTTOK].pri=10;
	 
	   	    }
	   else 
	   {
		printf("Illegal use of DELAY \n");
		return(1);
           }

          
	 }

/*        check for delshft symbol             */
          if(newtok==DELSHFTSYM)
	  {
           temp=my_symb[toklist[lstptr+1]].com;
   /* !! */   if(is_uvar(temp))
	   {
	    my_symb[LASTTOK].com=temp+10000; /* create a temporary sybol */
            NDELAYS++;
           toklist[lstptr+1]=LASTTOK;
	  	
	    my_symb[LASTTOK].pri=10;
	 
	   	    }
	   else 
	   {
		printf("Illegal use of DELAY Shift \n");
		return(1);
           }

          
	 }



/* check for shift  */
	  if(newtok==SHIFTSYM||newtok==ISHIFTSYM)
	  {
           temp=my_symb[toklist[lstptr+1]].com;
/* !! */	   if(is_uvar(temp) || is_ucon(temp))
	   {
             if(is_uvar(temp))my_symb[LASTTOK].com=temp+10000;
	        if(is_ucon(temp))my_symb[LASTTOK].com=temp+3000;
/* create a temporary sybol */
         
           toklist[lstptr+1]=LASTTOK;
	  	
	    my_symb[LASTTOK].pri=10;
	 
	   	    }
	   else 
	   {
		printf("Illegal use of SHIFT \n");
		return(1);
           }

          
       	  }


  
           

 next:
          if((newtok==ENDTOK)&&(oldtok==STARTTOK))break;
         
          if(newtok==LPAREN)
           {
             tokstak[tokptr]=LPAREN;
             tokptr++;
             oldtok=LPAREN;
             goto getnew;
            }
           if(newtok==RPAREN)
           {
            switch(oldtok)
                  {
                     case LPAREN:
                                 tokptr--;
                                 oldtok=tokstak[tokptr-1];
                                 goto getnew;
                     case COMMA:
                                 tokptr--;
                                 ncomma++;
                                 oldtok=tokstak[tokptr-1];
                                 goto next;
                  }
           }
           if((newtok==COMMA)&&(oldtok==COMMA))
           {
            tokstak[tokptr]=COMMA;
            tokptr++;
            goto getnew;
           }
           if(my_symb[oldtok%THOUS].pri>=my_symb[newtok%THOUS].pri)
           {
            command[comptr]=my_symb[oldtok%THOUS].com;
/*             printf("com(%d)=%d\n",comptr,command[comptr]); */
	    if((my_symb[oldtok%THOUS].arg==2)&&
	       (my_symb[oldtok%THOUS].com/100==1))
	      ncomma--;
            my_com=command[comptr];
	                comptr++;
 /*   New code   3/95      */
	   if(my_com==NUMSYM){
/*             printf("tp=%d ",tokptr);  */
	     tokptr--;
/*     printf(" ts[%d]=%d ",tokptr,tokstak[tokptr]); */
	     command[comptr]=tokstak[tokptr-1];
/*	     printf("xcom(%d)=%d\n",comptr,command[comptr]);  */
	     comptr++;
	     tokptr--;
	     command[comptr]=tokstak[tokptr-1];
/*	     printf("xcom(%d)=%d\n",comptr,command[comptr]);  */
	     comptr++;
	   }
 /*   end new code    3/95    */
           if(my_com==SUMSYM){
	     loopstk[lptr]=comptr;
             comptr++;
             lptr++;
             ncomma-=1;
	   }
           if(my_com==ENDSUM){
	     lptr--;
             jmp=comptr-loopstk[lptr]-1;
             command[loopstk[lptr]]=jmp;
	   }
	   if(my_com==MYIF){
         	     loopstk[lptr]=comptr; /* add some space for jump */
                     comptr++;
		     lptr++;
		     nif++;
                }    
	   if(my_com==MYTHEN){ 
		              /* First resolve the if jump */
			lptr--;
			jmp=comptr-loopstk[lptr];  /* -1 is old */
			command[loopstk[lptr]]=jmp;
			   /* Then set up for the then jump */
			loopstk[lptr]=comptr;
			lptr++;
			comptr++;
			nthen++;
			}
	   if(my_com==MYELSE){
			     lptr--;
			     jmp=comptr-loopstk[lptr]-1;
			     command[loopstk[lptr]]=jmp;
			     nelse++;
			     }


			      
			      
             if(my_com==ENDDELAY||my_com==ENDSHIFT||my_com==ENDISHIFT){
        
	     ncomma-=1;
                }
             if(my_com==ENDDELSHFT)
	       ncomma-=2;  
           /*  if(my_com==CONV||my_com==DCONV){
       	      ncomma-=1;
	     }  */

   
            /*    CHECK FOR USER FUNCTION       */
            if(is_ufun(my_com))
            {
             my_arg=my_symb[oldtok].arg;
                         command[comptr]=my_arg;
             comptr++;
             ncomma=ncomma+1-my_arg;
            }
           /*      USER FUNCTION OKAY          */
            tokptr--;
            oldtok=tokstak[tokptr-1];
            goto next;
          }
  /*    NEW code       3/95     */
	  if(newtok==NUMTOK){
	    tokstak[tokptr++]=toklist[lstptr++];
	    tokstak[tokptr++]=toklist[lstptr++];
	  }
 /*  end  3/95     */
          tokstak[tokptr]=newtok;
          oldtok=newtok;
          tokptr++;
          goto getnew;
       }
        if(ncomma!=0){
        printf("Illegal number of arguments\n");
	return(1);
        }
	if((nif!=nelse)||(nif!=nthen)){
	  printf("If statement missing ELSE or THEN \n");
	  return(1);
	    }
        command[comptr]=my_symb[ENDTOK].com;

       /* pr_command(command); */
        return(0);
    }

pr_command(command)
     int *command;
{
 int i=0;
 int token;
 while(1){
  token=command[i];
  printf("%d %d \n",i,token);
  if(token==ENDEXP)return;
   i++;
  }
}




show_where(string,index)
     char *string;
     int index;
{
  char junk[MAXEXPLEN];
  int i;
  for(i=0;i<index;i++)junk[i]=' ';
  junk[index]='^';
  junk[index+1]=0;
  printf("%s\n%s\n",string,junk);
}

function_sym(token) /* functions should have ( after them  */
{
  int com=my_symb[token].com;
  int i1=com/100;

    if(i1==0&&!unary_sym(token))return(1); /* single variable functions */
  if(i1==1&&!binary_sym(token))return(1); /* two-variable function */
  if(i1==UFUN||i1==7||i1==6||i1==5)return(1);
  if(token==DELSHFTSYM||token==DELSYM||token==SHIFTSYM||token==ISHIFTSYM||com==MYIF||com==MYTHEN||com==MYELSE
     ||com==SUMSYM||com==ENDSUM)return(1);
  return(0);
}

unary_sym(token)
{

  if(token==9||token==55)return(1);
  return(0);
}

binary_sym(token)
     int token;
{
  if(token>2&&token<9)return(1);
  if(token>43&&token<51)return(1);
  if(token==54)return(1);
  return(0);
}

pure_number(token)
     int token;
{
  int com=my_symb[token].com;
  int i1=com/100;
/* !! */  if(token==NUMTOK||isvar(i1)||iscnst(i1)||isker(i1)||i1==8||token==INDX)
    return(1);
  return(0);
}


gives_number(token)
     int token;
{
  int com=my_symb[token].com;
  int i1=com/100;
  if(token==INDX)return(1);
  if(token==NUMTOK)return(1);
  if(i1==0&&!unary_sym(token))return(1); /* single variable functions */
  if(i1==1&&!binary_sym(token))return(1); /* two-variable function */
  /* !! */ if(i1==8||isvar(i1)||iscnst(i1)||i1==7||i1==6||i1==5||isker(i1)||i1==UFUN)return(1);
  if(com==MYIF||token==DELSHFTSYM||token==DELSYM||token==SHIFTSYM||token==ISHIFTSYM||com==SUMSYM)return(1);
  return(0);
}


check_syntax(oldtoken,newtoken) /* 1 is BAD!   */
     int oldtoken,newtoken;
{
  int com1=my_symb[oldtoken].com,com2=my_symb[newtoken].com;

/* if the first symbol or (  or binary symbol then must be unary symbol or 
   something that returns a number or another (   
*/



  if(unary_sym(oldtoken)||oldtoken==COMMA||oldtoken==STARTTOK
     ||oldtoken==LPAREN||binary_sym(oldtoken))
   {
     if(unary_sym(newtoken)||gives_number(newtoken)||newtoken==LPAREN)return(0);

     return(1);
   }

/* if this is a regular function, then better have ( 
*/
 
 if(function_sym(oldtoken)){
   if(newtoken==LPAREN)return(0);

   return(1);
 }

/* if we have a constant or variable or ) or kernel then better
   have binary symbol or "then" or "else" as next symbol
*/
   
 if(pure_number(oldtoken)){
   if(binary_sym(newtoken)||newtoken==RPAREN
      ||newtoken==COMMA||newtoken==ENDTOK)
     return(0);

   return(1);
 }

 if(oldtoken==RPAREN){
   if(binary_sym(newtoken)||newtoken==RPAREN
      ||newtoken==COMMA||newtoken==ENDTOK)return(0);
   if(com2==MYELSE||com2==MYTHEN||com2==ENDSUM)return(0);

   return(1);
 }

  printf("Bad token %d \n",oldtoken);
  return(1);
    
}


/******************************
*    PARSER                   *
******************************/




 make_toks(dest,my_token)
 char *dest;
 int *my_token;
 {
 char num[40],junk[10];
 double value;
  int old_tok=STARTTOK,tok_in=0;
 int index=0,err,token,nparen=0,lastindex=0;
 union    /*  WARNING  -- ASSUMES 32 bit int  and 64 bit double  */
   {
     struct {
       int int1;
       int int2;
     } pieces;
     struct {
       double z;
     } num;
   } encoder;

 while(dest[index]!='\0')
  {
   lastindex=index;
   find_tok(dest,&index,&token);
   if((token==MINUS)&&
   ((old_tok==STARTTOK)||(old_tok==COMMA)||(old_tok==LPAREN)))
  token=NEGATE;
  if(token==LPAREN)++nparen;
  if(token==RPAREN)--nparen;
  if(token==NSYM)
    {
      if(do_num(dest,num,&value,&index)){
	show_where(dest,index);
	return(1);
      }
/*    new code        3/95      */
      encoder.num.z=value;
      my_token[tok_in++]=NUMTOK;
      my_token[tok_in++]=encoder.pieces.int1;
      my_token[tok_in++]=encoder.pieces.int2;
      if(check_syntax(old_tok,NUMTOK)==1){
	 printf("Illegal syntax \n");
	 show_where(dest,lastindex);
	 return(1);
       }
      old_tok=NUMTOK;
 
    }
   
   else
     {
       my_token[tok_in++]=token;
       if(check_syntax(old_tok,token)==1){
	 printf("Illegal syntax (Ref:%d %d) \n",old_tok,token);
	 show_where(dest,lastindex);
         tokeninfo(old_tok);
         tokeninfo(token);
	 return(1);
       }

       old_tok=token;
     }
 }

my_token[tok_in++]=ENDTOK;
if(check_syntax(old_tok,ENDTOK)==1){
  printf("Premature end of expression \n");
  show_where(dest,lastindex);
  return(1);
}
if(nparen!=0)
{
 if(ERROUT)printf(" parentheses don't match\n");
 return(1);
}
return(0);

}

tokeninfo(tok)
int tok;
{
 printf(" %s %d %d %d %d \n",
	my_symb[tok].name,my_symb[tok].len,my_symb[tok].com,
        my_symb[tok].arg,my_symb[tok].pri);
}


do_num(source,num,value,ind)
char *source,*num;
double *value;
int *ind;
{
 int j=0,i=*ind,error=0;
 int ndec=0,nexp=0,ndig=0;
 char ch,oldch;
 oldch='\0';
 *value=0.0;
 while(1)
 {
  ch=source[i];
  if(((ch=='+')||(ch=='-'))&&(oldch!='E'))break;
  if((ch=='*')||(ch=='^')||(ch=='/')||(ch==',')||(ch==')')||(ch=='\0')
              || (ch=='|') || (ch=='>') || (ch=='<') || (ch=='&')
               || (ch=='='))break;
  if((ch=='E')||(ch=='.')||(ch=='+')||(ch=='-')||isdigit(ch))
  {
   if(isdigit(ch))ndig++;
   switch(ch)
             {
              case 'E':
                       nexp++;
                       if((nexp==2)||(ndig==0))goto err;break;
              case '.':
                       ndec++;
                       if((ndec==2)||(nexp==1))goto err;break;

             }
   num[j]=ch;
   j++;
   i++;
   oldch=ch;
  }
  else
  {
err:
    num[j]=ch;
    j++;
    error=1;
    break;
  }
  }
  num[j]='\0';
  if(error==0)*value=atof(num);
  else
  if(ERROUT)printf(" illegal expression: %s\n",num);
  *ind=i;
  return(error);
}



convert(source,dest)
char *source,*dest;
{
 char ch;
 int i=0,j=0;
 while (1)
 {
  ch=source[i];
  if(!isspace(ch))dest[j++]=ch;
  i++;
  if(ch=='\0')break;
 }
 strupr(dest);
}



find_tok(source,index,tok)
char *source;
int *index,*tok;
{
 int i=*index,maxlen=0,symlen;
 int k,j,my_tok,match;
 my_tok=NSYM;
 for(k=0;k<NSYM;k++)
 {
  symlen=my_symb[k].len;
  if(symlen<=maxlen)continue;

   match=1;
   for(j=0;j<symlen;j++)
   {
    if(source[i+j]!=my_symb[k].name[j])
     {
      match=0;
      break;
     }
   }
   if(match!=0)
    {
     my_tok=k;
     maxlen=symlen;
    }
 }
   *index=*index+maxlen;
   *tok=my_tok;
}

double pmod(x,y)
     double x,y;
{
  double z=fmod(x,y);
  if(z<0)z+=y;
  return(z);
}

two_args()
{
 fun2[4]=atan2;
 fun2[5]=pow;
 fun2[6]=max;
 fun2[7]=min;
/*  fun2[8]=fmod;  */
 fun2[8]=pmod; /* This always gives an answer in [0,y) for mod(x,y) */
 fun2[9]=dand;
 fun2[10]=dor;
 fun2[11]=dgt;
 fun2[12]=dlt;
 fun2[13]=deq;
 fun2[14]=dge;
 fun2[15]=dle;
 fun2[16]=dne;
 fun2[17]=normal;
 fun2[18]=bessel_j;
 fun2[19]=bessel_y;
 
 




}
/*   These are the Bessel Functions; if you dont have them then
     return some sort of dummy value or else write a program 
     to compute them
*/

double bessel_j(x,y)
     double x,y;
{
 int n=(int)x;
 return(jn(n,y));
}

double bessel_y(x,y)
     double x,y;
{
 int n=(int)x;
 return(yn(n,y));
}

/*
double pow2(z,w)
double z,w;
{
 return(pow(z,w));

 double sign=1.0;
 if(floor(w)==w){
 if(z<0.0)sign=-1.0;
 if(fabs(fmod(w,2.0))<1.e-10)sign=1.0;
  return(sign*pow(fabs(z),w));
 }
 else
 return(pow(fabs(z),w));
 
}
 
 */



/*********************************************
          FANCY DELAY HERE                   *-------------------------<<<
*********************************************/

char *com_name(com)
int com;
{
    int i;
    for( i=0;i<NSYM;i++)
	if( my_symb[i].com == com ) break;
    if( i < NSYM )
	return my_symb[i].name;
    else
	return "";
}
double do_shift(shift,variable)
double shift,variable;
{
  int it, in;
  int i=(int)(variable),ish=(int)shift;

/* printf( "shifting %d (%s) by %d to %d (%s)\n", 
 *	(int)variable, com_name((int)variable), (int)shift, i, com_name(i) );
 */

  if(i<0) return(0.0);
   it=i/100;
  switch(it){
  case 2: 
  case 3:
  case 4:
  case 5:
        in=i-200;
        in+=ish;
	if(in>NCON)
	  return 0.0;
	else
	  return constants[in]; 
	break;
  case 100:
  case 101:
  case 102:
  case 103:
  case 104:
  case 105:
  case 106:
  case 107:
  case 108:
  case 109:
    case 110:
  case 111:
  case 112:
  case 113:
  case 114:
  case 115:
  case 116:
  case 117:
  case 118:
  case 119:
	in=i-10000;
 
        in+=ish;
	if(in>MAXODE)
	  return 0.0;
	else 
	  return variables[in];  
  default:
    printf("This can't happen: Invalid symbol index for SHIFT\n");
    return 0.0;
  }
}
double do_ishift(shift,variable)
double shift,variable;
{
  int it, in;
  int i=(int)(variable),ish=(int)shift;

/* printf( "shifting %d (%s) by %d to %d (%s)\n", 
 *	(int)variable, com_name((int)variable), (int)shift, i, com_name(i) );
 */
 return variable+shift;
 
}

double do_delay_shift(delay,shift,variable)
     double delay,shift,variable;
{
 int it, in;
  int i=(int)(variable),ish=(int)shift;
  if(i<0) return(0.0);
  in=i-10000+ish;

  if(in>MAXODE)
    return 0.0;
 
  if(del_stab_flag>0){
    if(DelayFlag&&delay>0.0)
      return(get_delay(in-1,delay));
    return(variables[in]);
  }
 
  return(delay_stab_eval(delay,in));
  
 
 


}
double do_delay(delay,i)
double delay,i;
{
  double z;
  int variable, it, in;

    variable = i-10000;

  if(del_stab_flag>0){
    if(DelayFlag&&delay>0.0)
      return(get_delay(variable-1,delay));
    return(variables[variable]);
  }
 
  return(delay_stab_eval(delay,(int)variable));
  
}
/*
double Exp(z)
double z;
{
 if(z>700)return(1.01423e+304);
 return(exp(z));
}
double Ln(z)
double z;
{
 if(z<1e-320)return(-736.82724);
 return(log(z));
}
double Log10(z)
double z;
{
 if(z<1e-320)return(-320.);
 return(log10(z));
}
*/

one_arg()
{
 fun1[0]=sin;
 fun1[1]=cos;
 fun1[2]=tan;
 fun1[3]=asin;
 fun1[4]=acos;
 fun1[5]=atan;
 fun1[6]=sinh;
 fun1[7]=tanh;
 fun1[8]=cosh;
 fun1[9]=fabs;
 fun1[10]=exp;
 fun1[11]=log;
 fun1[12]=log10;
 fun1[13]=sqrt;
 fun1[14]=neg;
 fun1[15]=recip;
 fun1[16]=heaviside;
 fun1[17]=signum;
 fun1[18]=floor;
 fun1[19]=rndom;
 fun1[20]=dnot;
 fun1[21]=erf;
 fun1[22]=erfc;
 fun1[23]=hom_bcs;
  fun1[24]=poidev;
}


double normal(mean,std)
     double mean,std;
{
 double fac,r,v1,v2;
 if(BoxMullerFlag==0){ 
   do {
     v1=2.0*ndrand48()-1.0;
     v2=2.0*ndrand48()-1.0;
     r=v1*v1+v2*v2;
   } while(r>=1.0);
   fac=sqrt(-2.0*log(r)/r);
   BoxMuller=v1*fac;
   BoxMullerFlag=1;
   return(v2*fac*std+mean);
 }
 else {
   BoxMullerFlag=0;
   return(BoxMuller*std+mean);
 }
}


double max(x,y)
double x,y;
{
 return(((x>y)?x:y));
}

double min(x,y)
double x,y;
{
 return(((x<y)?x:y));
}

double neg(z)
double z;
{
 return(-z);
}

double recip(z)
double z;
{
 return(1.00/z);
}

double heaviside(z)
double z;
{
 float w=1.0;
 if(z<0)w=0.0;
 return(w);
}

double rndom( z)
double z;
{
 /* return (z*(double)rand()/32767.00); */
  return(z*ndrand48());
}

double signum(z)
double z;

{
  if(z<0.0)return(-1.0);
  if(z>0.0)return(1.0);
  return(0.0);
}


/*  logical stuff  */


double dnot(x)
double x;
{
 return((double)(x==0.0));
}
double dand(x,y)
double x,y;
{
 return((double)(x&&y));
}
double dor(x,y)
double x,y;
{
 return((double)(x||y));
}
double dge(x,y)
double x,y;
{
 return((double)(x>=y));
}
double dle(x,y)
double x,y;
{
 return((double)(x<=y));
}
double deq(x,y)
double x,y;
{
 return((double)(x==y));
}
double dne(x,y)
double x,y;
{
 return((double)(x!=y));
}
double dgt(x,y)
double x,y;
{
 return((double)(x>y));
}
double dlt(x,y)
double x,y;
{
 return((double)(x<y));
}


/*              end of logical stuff    */





 double evaluate(equat)
 int *equat;
 {
  uptr=0;
  stack_pointer=0;
  return(eval_rpn(equat));
 }


 double eval_rpn(equat)
 int *equat;
 {
   int i,it,in,j,*tmpeq;
  int sumcount,ii,new[500],is;
  int isum,aind;
  int low,high,ind,iloop,ijmp;
  double temp,temx,temy,temz;
  double sum;
  union    /*  WARNING  -- ASSUMES 32 bit int  and 64 bit double  */
   {
     struct {
       int int1;
       int int2;
     } pieces;
     struct {
       double z;
     } num;
   } encoder;



  while((i=*equat++)!=ENDEXP)
  {
 
   switch(i)
   {
   case NUMSYM:
     encoder.pieces.int2=*equat++;
     encoder.pieces.int1=*equat++;
     PUSH(encoder.num.z);
     break;
   case ENDFUN:
   		 i=*equat++;
		
    		 uptr-=i;
		
   		 break;

   case MYIF:
		temx=POP;
		ijmp=*equat++;
		if(temx==0.0)equat+=ijmp;
                break;
   case MYTHEN:
	       ijmp=*equat++;
	       equat+=ijmp;
		break;
   case MYELSE:
		break;
  
   case ENDDELSHFT:
     temx=POP;
     temy=POP;
     temz=POP;
     PUSH(do_delay_shift(temx,temy,temz));
     break;
   case ENDDELAY:
		    temx=POP;
		    temy=POP;
                   
                   PUSH(do_delay(temx,temy));
		   break;

   case ENDSHIFT:
                 temx=POP;
                 temy=POP;
                 PUSH(do_shift(temx,temy));
                 break;
   case ENDISHIFT:
                 temx=POP;
                 temy=POP;
                 PUSH(do_ishift(temx,temy));
                 break;
   case SUMSYM:
              temx=POP;
              high=(int)temx;
              temx=POP;
              low=(int)temx;
              ijmp=*equat++;
              sum=0.0;
              if(low<=high){
		for(is=low;is<=high;is++){
		  tmpeq=equat;
		  constants[SumIndex]=(double)is;
		  sum+=eval_rpn(tmpeq);
		}
	      }
             equat+=ijmp;
             PUSH(sum);
             break;

   case ENDSUM:
            return(POP);
   case INDXCOM:
     PUSH(CurrentIndex);
     break;
   default:
   {
   it=i/100;
   in=i%100;
   switch(it)
    {
     case 0: PUSH(fun1[in](POP));
            break;
     case 1:
	    {
	     if(in==0){temx=POP;temy=POP;PUSH(temx+temy);goto bye;}
	     if(in==2){temx=POP;temy=POP;PUSH(temx*temy);goto bye;}
	     if(in==1){temx=POP;temy=POP;PUSH(temy-temx);goto bye;}
	     if(in==3){temx=POP;if(temx==0.0)temx=DOUB_EPS;
		      temy=POP;PUSH(temy/temx);goto bye;}
             temx=POP;
             temy=POP;
	     PUSH(fun2[in](temy,temx));break;
	    }
     case 2: 
     case 3: 
     case 4:
     case 5:
              PUSH(constants[i-200]);break;
     case 6:  PUSH(network_value(POP,in));break;
     case 7: PUSH(lookup(POP,in));break;
            
     case 8:
            PUSH(ustack[uptr-1-in]); break;
     case 10: PUSH(ker_val(in));break;
    case 100:
    case 101:
     case 102: 
     case 103:
     case 104:
    case 105:
    case 106:
    case 107:
    case 108:
    case 109:
      case 110:
  case 111:
  case 112:
  case 113:
  case 114:
  case 115:
  case 116:
  case 117:
  case 118:
  case 119:
             PUSH(variables[i-10000]); break;
       
     /* indexes for shift and delay operators... */
     case 32:
     case 33:
    case 34:
    case 35:
         PUSH((double)(i-3000)); break;
     case 200:
     case 201:
     case 202:
     case 203:
    case 204:
    case 205:
    case 206:
    case 207:
    case 208:
    case 209:
    case 210:
     case 211:
     case 212:
     case 213:
    case 214:
    case 215:
    case 216:
    case 217:
    case 218:
    case 219:
             PUSH((double)(i-10000)); break;

     case UFUN: i=*equat++;
         
            for(j=0;j<i;j++)
            {
            ustack[uptr]=POP;
	 
	    uptr++;
            }
            PUSH(eval_rpn(ufun[in])); 
break;
    }
bye: j=0;
   }
  }
  }
   return(POP);

}


/*  STRING STUFF  */
#ifndef STRUPR
strupr(s)
char *s;
{
 int i=0;
 while(s[i])
 {
  if(islower(s[i]))s[i]-=32;
  i++;
  }
}


strlwr(s)
char *s;
{
 int i=0;
 while(s[i])
 {
  if(isupper(s[i]))s[i]+=32;
  i++;
  }
}
#endif




