#ifndef _load_eqn_h_
#define _load_eqn_h_

#include <stdio.h>

/*
Options are set accroding to an order of precendence

command line < mfile < .xpprc < default.opt

Add any options here that you might want to track.
*/
typedef struct {
   int BIG_FONT_NAME;
   int SMALL_FONT_NAME;
   int BACKGROUND;
   int IXPLT;
   int IYPLT;
   int IZPLT;
   int AXES;
   int NMESH;
   int METHOD;
   int TIMEPLOT;
   int MAXSTOR;
   int TEND;
   int DT;
   int T0;
   int TRANS;
   int BOUND;
   int TOLER;
   int DELAY;
   int XLO;
   int XHI;
   int YLO;
   int YHI;  
   int UserBlack;
   int UserWhite;
   int UserMainWinColor;
   int UserDrawWinColor;
   int UserGradients;
   int UserBGBitmap;
   int UserMinWidth;
   int UserMinHeight;
   int YNullColor;
   int XNullColor;
   int StableManifoldColor;
   int UnstableManifoldColor;
   int START_LINE_TYPE;
   int RandSeed;
   int PaperWhite;
   int COLORMAP;
   int NPLOT;
   int DLL_LIB;
   int DLL_FUN;
   int XP;
   int YP;
   int ZP;
   int NOUT;
   int VMAXPTS;
   int TOR_PER;
   int JAC_EPS;
   int NEWT_TOL;
   int NEWT_ITER;
   int FOLD;
   int DTMIN;
   int DTMAX;
   int ATOL;
   int TOL;
   int BANDUP;
   int BANDLO;
   int PHI;
   int THETA;
   int XMIN;
   int XMAX;
   int YMIN;
   int YMAX;
   int ZMIN;
   int ZMAX;
   int POIVAR;
   int OUTPUT;
   int POISGN;  
   int POIEXT;
   int POISTOP;
   int STOCH;
   int POIPLN;
   int POIMAP;
   int RANGEOVER;
   int RANGESTEP;
   int RANGELOW;
   int RANGEHIGH;
   int RANGERESET;
   int RANGEOLDIC;
   int RANGE;
   int NTST;
   int NMAX;
   int NPR;
   int NCOL;
   int DSMIN;
   int DSMAX;
   int DS;
   int PARMAX;
   int NORMMIN;
   int NORMMAX;
   int EPSL;
   int EPSU;
   int EPSS;
   int RUNNOW;
   int SEC;
   int UEC;
   int SPC;
   int UPC;
   int AUTOEVAL;
   int AUTOXMAX;
   int AUTOYMAX;
   int AUTOXMIN;
   int AUTOYMIN;
   int AUTOVAR;
   int PS_FONT;
   int PS_LW;   
   int PS_FSIZE;
   int PS_COLOR;
   int FOREVER;
   int BVP_TOL;
   int BVP_EPS;
   int BVP_MAXIT;
   int BVP_FLAG;
   int SOS;
   int FFT;
   int HIST;
   int PSFlag;
   int ATOLER;
   int MaxEulIter;
   int EulTol;
   int EVEC_ITER;
   int EVEC_ERR;
   int NULL_ERR;
   int NEWT_ERR;
   int NULL_HERE;
 
  
  } OptionsSet;


void dump_torus(FILE *fp, int f);
void load_eqn(void);
void set_X_vals(void);
void set_all_vals(void);
void read_defaults(FILE *fp);
void fil_flt(FILE *fpt, double *val);
void fil_int(FILE *fpt, int *val);
void add_intern_set(char *name, char *does);
void extract_action(char *ptr);
void extract_internset(int j);
void do_intern_set(char *name1, char *value);
int msc(char *s1, char *s2);
void set_internopts(OptionsSet *mask);
void set_internopts_xpprc_and_comline(void);
void split_apart(char *bob, char *name, char *value);
void check_for_xpprc(void);
void stor_internopts(char *s1);
void set_option(char *s1, char *s2,int force,OptionsSet *mask);



#endif
