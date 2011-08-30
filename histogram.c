#include <stdlib.h> 
#include <stdio.h>
#include <stdio.h>
#include <math.h>

#include "xpplim.h"

#include "fftn.h"

double evaluate();
double ndrand48();

void mycor();
float *get_data_col();
extern int DCURY,MAXSTOR;
extern char uvar_names[MAXODE][12];
typedef struct {
  int nbins,type,col,col2;
  double xlo,xhi;
  char cond[80];
} HIST_INFO;

HIST_INFO hist_inf = {10,0,0,0,1,1,""};

extern int NCON,NSYM,NCON_START,NSYM_START;

extern float **storage;
extern int storind;
int hist_len,four_len;
float *my_hist[MAXODE+1];
float *my_four[MAXODE+1];
int HIST_HERE,FOUR_HERE;

extern int NEQ,NODE,NMarkov,FIX_VAR;

extern char *no_hint[],*info_message;

four_back()
{
 if(FOUR_HERE){
   set_browser_data(my_four,1);
   /*   my_browser.data=my_four;
	my_browser.col0=1; */
   refresh_browser(four_len);
 }
}

hist_back()
{
 if(HIST_HERE){
   set_browser_data(my_hist,1);
   /*
   my_browser.data=my_hist;
   my_browser.col0=1; */
   refresh_browser(hist_len);
 }
}

new_four(nmodes,col)
     int nmodes,col;
{
  int i;
  int length=nmodes+1;
  float *bob;
  if(FOUR_HERE){
   data_back();
   free(my_four[0]);
   free(my_four[1]);
   free(my_four[2]);
   FOUR_HERE=0;
 }
  four_len=nmodes;
  my_four[0]=(float *)malloc(sizeof(float)*length);
  my_four[1]=(float *)malloc(sizeof(float)*length);
  my_four[2]=(float *)malloc(sizeof(float)*length);
  if(my_four[2]==NULL){
   free(my_four[1]);
   free(my_four[2]);
   err_msg("Cant allocate enough memory...");
   return;
 }
 FOUR_HERE=1;
 for(i=3;i<=NEQ;i++)my_four[i]=storage[i];
 for(i=0;i<length;i++)my_four[0][i]=i;
 /*  sft(my_browser.data[col],my_four[1],my_four[2],length,storind);
  */
 bob=get_data_col(col);
    fft(bob,my_four[1],my_four[2],nmodes,storind);
 four_back();
  ping();
}

new_hist(nbins,zlo,zhi,col,col2,condition,which)
     int nbins;
     int col,col2,which;
     double zlo,zhi;
     char *condition;
     
{
  int i,j,n=2,index;
  int command[256];
  int cond=0,flag=1;
  double z,y;
  double dz;
  int length=nbins+1;
  int count=0;
  if(length>=MAXSTOR)
    length=MAXSTOR-1;
  dz=(zhi-zlo)/(double)(length-1);
  if(HIST_HERE){
    data_back();
    free(my_hist[0]);
    free(my_hist[1]);
    HIST_HERE=0;
  }
  hist_len=length;
  my_hist[0]=(float *)malloc(sizeof(float)*length);
  my_hist[1]=(float *)malloc(sizeof(float)*length);
  if(my_hist[1]==NULL){
    free(my_hist[0]);
    err_msg("Cannot allocate enough...");
    return;
  }
  HIST_HERE=1;
  for(i=2;i<=NEQ;i++)my_hist[i]=storage[i];
  for(i=0;i<length;i++){
    my_hist[0][i]=(float)(zlo+dz*i);
    my_hist[1][i]=0.0;
  }
  if(which==0){
    if(strlen(condition)==0)cond=0;
    else
      {
	if(add_expr(condition,command,&i)){
	  err_msg("Bad condition. Ignoring...");
	  
	}
	else {
	  cond=1;
	}
      }
    /* printf(" cond=%d \n condition=%s \n,node=%d\n", 
       cond,condition,NODE);  */
    for(i=0;i<storind;i++)
      {
	flag=1;
	if(cond){
	  for(j=0;j<NODE+1;j++)set_ivar(j,(double)storage[j][i]);
	  for(j=0;j<NMarkov;j++)
	    set_ivar(j+NODE+1+FIX_VAR,(double)storage[j+NODE+1][i]);
	  z=evaluate(command);
	  if(fabs(z)>0.0)flag=1;
	  else flag=0;
	}
	z=(storage[col][i]-zlo)/dz;
	index=(int)z;
	if(index>=0&&index<length&&flag==1){
	  my_hist[1][index]+=1.0;
	  count++;
	}
      }
    NCON=NCON_START;
    NSYM=NSYM_START;
    hist_back();
    ping();
    return;
  }
  if(which==1){
    for(i=0;i<storind;i++){
      for(j=0;j<storind;j++){
	y=storage[col][i]-storage[col][j];
	z=(y-zlo)/dz;
	index=(int)z;
	if(index>=0&&index<length)
	  my_hist[1][index]+=1.0;
      }
    }
    hist_back();
    ping();
    return;
  }
  if(which==2){
    mycor(storage[col],storage[col2],storind,zlo,zhi,nbins,my_hist[1],1);
    hist_back();
    ping();
    return;
  }

}

    
  

column_mean()
{
 int i;
 char bob[100];
 double sum,sum2,ss;
 double mean,sdev;
 if(storind<=1){
   err_msg("Need at least 2 data points!");
   return;
 }
 if(get_col_info(&hist_inf.col)==0)return;
 sum=0.0;
 sum2=0.0;
 for(i=0;i<storind;i++){
   ss=storage[hist_inf.col][i];
   sum+=ss;
   sum2+=(ss*ss);
 }
 mean=sum/(double)storind;
 sdev=sqrt(sum2/(double)storind-mean*mean);
 sprintf(bob,"Mean=%g Std. Dev. = %g ",mean,sdev);
 err_msg(bob);
}
get_col_info(col)
 int *col;
{
 char variable[20];
 if(*col==0)
   strcpy(variable,"t");
 else
   strcpy(variable,uvar_names[*col-1]);
 new_string("Variable ",variable);
 find_variable(variable,col);
 if(*col<0){
   err_msg("No such variable...");
   return(0);
 }
 return(1);
}

compute_power()
{
  int i;
  double s,c;
  float *datx,*daty;
  compute_fourier();
  if((NEQ<2)||(storind<=1))return;
  datx=get_data_col(1);
  daty=get_data_col(2);

  for(i=0;i<four_len;i++){
    c=datx[i];
    s=daty[i];
    datx[i]=sqrt(s*s+c*c);
    daty[i]=atan2(s,c);
  }
 
}

compute_fourier()
{
  int nmodes=10,col=0;
  if(NEQ<2){
    err_msg("Need at least three data columns");
    return;
  }
  /* new_int("Number of modes ",&nmodes); */
  if(storind<=1){
    err_msg("No data!");
    return;
  }
  if(get_col_info(&col)==1)
    nmodes=storind/2-1;
    new_four(nmodes,col);
}
 
compute_correl()
{

  new_int("Number of bins ",&hist_inf.nbins);
  new_float("Low ",&hist_inf.xlo);
  new_float("Hi ",&hist_inf.xhi);
  if(get_col_info(&hist_inf.col)==0)return;
  if(get_col_info(&hist_inf.col2)==0)return;
  new_hist(hist_inf.nbins,hist_inf.xlo,
	   hist_inf.xhi,hist_inf.col,hist_inf.col2,hist_inf.cond,2);
}
compute_stacor()
{
  new_int("Number of bins ",&hist_inf.nbins);
  new_float("Low ",&hist_inf.xlo);
  new_float("Hi ",&hist_inf.xhi);
  if(get_col_info(&hist_inf.col)==0)return;
   new_hist(hist_inf.nbins,hist_inf.xlo,
	   hist_inf.xhi,hist_inf.col,0,hist_inf.cond,1);
}

void mycor(float *x,float *y, int n,  double zlo, double zhi, int nbins, float *z, int flag)
{
  int i,j;
  int k,count=0;
  float sum,avx=0.0,avy=0.0;
  double dz=(zhi-zlo)/(double)nbins,jz;
  if(flag){
    for(i=0;i<n;i++){
      avx+=x[i];
      avy+=y[i];
    }
    avx=avx/(float)n;
    avy=avy/(float)n;
  }
  for(j=0;j<=nbins;j++){
    sum=0.0;
    count=0;
    jz=dz*j+zlo;
    for(i=0;i<n;i++){
      k=i+(int)jz;
      if((k>=0)&&(k<n)){
	count++;
	sum+=(x[i]-avx)*(y[k]-avy);
      }
    }
    if(count>0)
      sum=sum/count; 
    z[j]=sum;
  }
}

compute_hist()
{
  
  new_int("Number of bins ",&hist_inf.nbins);
  new_float("Low ",&hist_inf.xlo);
  new_float("Hi ",&hist_inf.xhi);
  if(get_col_info(&hist_inf.col)==0)return;
  new_string("Condition ",hist_inf.cond);
  new_hist(hist_inf.nbins,hist_inf.xlo,
	   hist_inf.xhi,hist_inf.col,0,hist_inf.cond,0);
}
  
  

sft(data,ct,st,nmodes,grid)
int grid,nmodes;
float *data,*ct,*st;
{
 int i,j;
 double sums,sumc;
 double tpi=6.28318530717959;
 double dx,xi,x;
 dx=tpi/(grid);
 for(j=0;j<nmodes;j++){
   sums=0.0;
   sumc=0.0;
   xi=j*dx;
   for(i=0;i<grid;i++){
    x=i*xi;
     sumc+=(cos(x)*data[i]);
     sums+=(sin(x)*data[i]);
   }
   if(j==0){
   ct[j]=sumc/(float)grid;
   st[j]=sums/(float)grid;
 }
   else{
     ct[j]=2.*sumc/(float)grid;
   st[j]=2.*sums/(float)grid;
   }
 }
}


fft(data,ct,st,nmodes,length)
     float *data,*ct,*st;
     int nmodes,length;
{
  double *im,*re;
  int dim[2],i;
  dim[0]=length;
  re=(double *)malloc(length*sizeof(double));
  im=(double *)malloc(length*sizeof(double));
  for(i=0;i<length;i++){
    im[i]=0.0;
    re[i]=data[i];
  }

   fftn(1,dim,re,im,1,-1);
   ct[0]=re[0];
   st[0]=0.0;
   for(i=1;i<nmodes;i++){
     ct[i]=re[i]*2.0;
     st[i]=im[i]*2.0;
   }
   free(im);
   free(re);
}

   
  








