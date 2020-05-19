
/*SAS code used for data analysis described in the publication:

Li Y, Schaubel DE and He K. Matching methods for obtaining survival 
functions to estimate the effect of a time-dependent treatment. 
Statistics in Biosciences. (2013). 1-22.

In observational studies of survival time featuring a binary time-dependent 
treatment, the hazard ratio (an instantaneous measure) is often used to represent 
the treatment effect. However, investigators are often more interested in the 
difference in survival functions. We propose semiparametric methods to estimate 
the causal effect of treatment among the treated with respect to survival probability. 
The objective is to compare post-treatment survival with the survival function that 
would have been observed in the absence of treatment. For each patient, we compute 
a prognostic score (based on the pre-treatment death hazard) and a propensity score
(based on the treatment hazard). Each treated patient is then matched with an alive, 
uncensored and not-yet-treated patient with similar prognostic and/or propensity scores. 
The experience of each treated and matched patient is weighted using a variant of Inverse 
Probability of Censoring Weighting to account for the impact of censoring. We propose 
estimators of the treatment-specific survival functions (and their difference), computed 
through weighted Nelson-Aalen estimators. Closed-form variance estimators are proposed which 
take into consideration the potential replication of subjects across matched sets. The proposed 
methods are evaluated through simulation, then applied to estimate the effect of kidney 
transplantation on survival among end-stage renal disease patients using data from 
a national organ failure registry.
*/ 

libname DE98SMA "d:\DATA\ESRD98\SURVIVAL\MATCHED\APPLICATION\";
libname DE98Y "d:\DATA\ESRD98\YUN\";


OPTIONS NOCENTER ls=80 PS=108 nofmterr notes MergeNoBy=ERROR nodate;


proc contents data=DE98Y.rrti_tx_25_apr_2009; run;

proc freq data=DE98Y.rrti_tx_25_apr_2009; tables yr_esrd; run;
proc freq data=DE98Y.rrti_tx_25_apr_2009; tables province; run;


%let d_end_obs_pd='31-DEC-98'd;

%let Z = age50by5 female minority yr90 diabetes PD num_comorbid;
%let ZC = age50by5 female minority diabetes PD num_comorbid;

%let D_tol = log(1.05);
%let TX_tol = log(1.05);

%let t_tx_max=1095;
%let t_D_max=1825;
%let month_max=60;

%let m_max=1;



data rrti_;
	set DE98Y.rrti_tx_25_apr_2009;  
if not(province="ON ") then delete;
drop d_censor;
TX=0;
DEAD0=0;
dead1=0;
censored=0;
if (d_tx ne . ) then do;
 X0=D_TX-D_RRTI; X0=max(X0,1);
 if (don_liv=0) then do;
  TX=1;
  t_tx=X0;
  if (d_death ne . ) then do;
   dead1=1; X1=d_death-d_tx; 
  end; else X1=&d_end_obs_pd-d_tX;
  X1=max(X1,1);
 end; 
end; else do;
 if (d_death ne . ) then do;
  dead0=1; X0=d_death-d_rrti; 
 end; else X0=&d_end_obs_pd-d_rrti;
 X0=max(X0,1);
end;

if (tx=1) then X=X0+X1; else X=X0;
censored=(d_death= . ) + (don_liv=1); censored=min(censored,1);

age50by5=(age_rrti-50)/5;
female=(sex="F");

yr90=yr_esrd-90;

GLOM=(PRIMDIAG='GLOM');
POLY=(PRIMDIAG='POLY');
DIAB=(PRIMDIAG='DIAB');
VASC=(PRIMDIAG='VASC');
OTH=1-(glom+poly+vasc+diab);

ASIAN=(RACE='ASI');
BLACK=(RACE='BLA');
NATIVE=(RACE='NAT');
INDIAN=(RACE='IND');
OTHER=(RACE='OTH');
minority = asian+black+native+indian+other;

one=1;
keep idno X X0 X1 TX t_tx dead0 dead1 censored &Z one;
run;
proc sort data=rrti_; by idno; run;

proc freq data=rrti_;
 tables one / out=out_n;
run;

data out_n1; 
 set out_n;
 call symput('n_total', count); 
run;


data rrti;
 set rrti_;
id_num=_N_; drop idno;
run;

proc phreg data=RRTI;
 model X0*dead0(0) = &Z / ties=breslow;
 output out=xbeta_d0 (keep=id_num xbeta_d0) xbeta=xbeta_d0;
 id id_num;
run;

proc phreg data=RRTI;
 model X0*TX(0) = &Z / ties=breslow;
 output out=xbeta_tx (keep=id_num xbeta_tx) xbeta=xbeta_tx;
 id id_num;
run;

data covset0;
 age50by5=0; female=0; minority=0; diabetes=0; PD=0; num_comorbid=0;
run;
proc phreg data=RRTI;
 model X*censored(0) = &ZC / ties=breslow;
 baseline out=baseoutC (keep=X Lam0C ) covariates=covset0 cumhaz=Lam0C / method=ch nomean;
 output out=xbeta_C (keep=id_num xbeta_C) xbeta=xbeta_C;
 id id_num;
run;

data C1;
 do X=1 to (&t_tx_max+&t_D_max);
  output;
 end;
run;

data baseoutC1;
 merge baseoutC C1 (in=in_C); by X;
if (in_C=1);
if (Lam0C ne . ) then Lam0C_=Lam0C;
retain Lam0C_;
if (Lam0C= . ) then Lam0C=Lam0C_;
drop Lam0C_;
run;

data tx1;
 set rrti;
if (TX=1) and (t_tx <= &t_tx_max);
keep id_num t_tx dead1 X1;
run;

proc sort data=tx1; by id_num; run;
proc sort data=xbeta_d0; by id_num; run;
proc sort data=xbeta_tx; by id_num; run;

data tx2;
 merge tx1 (in=in_) xbeta_d0 xbeta_tx; by id_num; 
if (in_=0) then delete;
run;

data tx3;
 set tx2;
match_id=_n_;
one=1;
run;

proc means data=tx3;
 var match_id;
 output out=outmean1 max(match_id)=m_id_max;
run;

data outmean2; 
 set outmean1;
 call symput('m_id_max', m_id_max); 
run;

proc sort data=rrti; by id_num; run;
data dial1;
 merge rrti (in=in_) xbeta_d0 xbeta_tx; by id_num; 
if (in_=0) then delete;
one=1;
keep id_num one X0 dead0 xbeta_d0 xbeta_tx;
run;


%macro matching1 (j_tot);

options nonotes;

 %do j=1 %to &j_tot;

data match1;
 merge tx3 (where=(match_id=&j) 
            rename=(xbeta_tx=xbeta_tx_m xbeta_D0=xbeta_D0_m id_num=id_m)
            keep=one id_num match_id t_tx xbeta_tx xbeta_D0)
	   dial1; by one;
if (id_num ne id_m);
if (X0 > t_tx);
if (abs(xbeta_tx-xbeta_tx_m) <= &tx_tol);
if (abs(xbeta_d0-xbeta_d0_m) <= &d_tol);
Xm=X0-t_tx;
deadm=dead0*(Xm <= &t_D_max);
Xm=min(Xm,&t_D_max);
U=ranuni(0);
keep id_num match_id Xm deadm U X0 t_tx;
run;

%if (&j=1) %then %do;
 data DE98SMA.matched; set match1; run;
%end; %else %do;
 proc append base=DE98SMA.matched data=match1; run;
%end;

%end;

options notes;

proc freq data=DE98SMA.matched;
 tables match_id / out=outfreq1;
run;

%mend matching1;

%matching1(j_tot=&m_id_max);

proc sort data=DE98SMA.matched; by match_id U; run;

data matched1;
 set DE98SMA.matched; by match_id;
if first.match_id then num=0;
 retain num;
 num=num+1;
run;

data matched2;
 set matched1;
if (num <= &m_max);
run;

proc freq data=matched2;
 tables match_id / out=outfreq2 (keep=match_id count rename=(count=m_j));
run;

data outfreq3;
 set outfreq2;
match_id_new=_n_;
run;

data matched3;
 merge outfreq3 (in=in_) matched2; by match_id;
if (in_=1);
match_id=match_id_new;
drop match_id_new;
run;

proc sort data=tx3; by match_id; run;
data tx4;
 merge tx3 outfreq3 (in=in_ keep=match_id match_id_new); by match_id;
if (in_=1);
match_id=match_id_new;
drop match_id_new;
run;

data outfreq4;
 set outfreq3;
match_id=match_id_new;
drop match_id_new;
run;

proc freq data=matched3 noprint;
 where (deadm=1);
 tables Xm / out=outdead_0 (keep=Xm);
run;

%macro match_expand0(j_tot);

options nonotes;

%do j=1 %to &j_tot;

data outdead_0a;
if (_n_=1) then set outfreq4 (where=(match_id=&j));
set outdead_0;
do num=1 to m_j;
 output;
end;
drop m_j;
run;

data matched3_j;
 set matched3;
if (match_id=&j);
run;
proc sort data=matched3_j; by match_id num; run;
proc sort data=outdead_0a; by match_id num; run;

data outdead_0b;
 merge outdead_0a 
       matched3_j (in=in_j keep=id_num m_j t_tx match_id num Xm rename=(Xm=Xm_)); 
 by match_id num; 
if (in_j=1);
if (0 <= Xm <= Xm_);
X=t_tx+Xm;
keep match_id num m_j id_num t_tx Xm X;
run;

proc sort data=outdead_0b; by X; run;
proc sort data=baseoutC1; by X; run;

data outdead_0c;
 merge outdead_0b (in=in_b) baseoutC1; by X;
if (in_b=1);
run;

proc sort data=outdead_0c; by id_num; run;
proc sort data=xbeta_C; by id_num; run;

data outdead_0d;
 merge outdead_0c (in=in_) xbeta_C; by id_num;
if (in_=1);
Wi_X=exp(Lam0C*exp(xbeta_C));
drop Lam0C xbeta_C;
run;

proc sort data=outdead_0d; by match_id num Xm; run;
proc sort data=matched3; by match_id num Xm; run;

data outdead_0e;
 merge outdead_0d (in=in_) matched3 (keep=deadm match_id num Xm); by match_id num Xm; 
if (in_=1);
if (deadm= . ) then deadm=0;
t2=Xm;
t1=t2-0.1;
TX=0;
Wi_X=Wi_X/m_j;
dNw=deadm*Wi_X;
run;

%if (&j=1) %then %do;
 data DE98SMA.match_expand_0; set outdead_0e; run; 
%end; %else %do;
 proc append base=DE98SMA.match_expand_0 data=outdead_0e; run;
%end;

%end;

options notes;

%mend match_expand0;

%match_expand0(j_tot=&m_id_max);

proc contents data=DE98SMA.match_expand_0; run;
proc means data=DE98SMA.match_expand_0; run;

proc sort data=DE98SMA.match_expand_0; by t2; run;
proc means data=DE98SMA.match_expand_0 noprint; 
 var dNw Wi_X;
 output out=out_NY_0 mean(dNw)=dNw_bar mean(Wi_X)=Yw_bar sum(Wi_X)=Yw;
 by t2; 
run;

data out_NY_0a;
 set out_NY_0;
dLam_0=dNw_bar/Yw_bar;
TX=0;
keep t2 dLam_0 dNw_bar Yw_bar Yw TX;
run;

data match_expand_0a;
 merge DE98SMA.match_expand_0 out_NY_0a; by t2;
dPhi_0_t=Wi_X*(deadm-dLam_0)/Yw;
keep id_num t2 dPhi_0_t;
run;

proc freq data=match_expand_0a noprint; 
 tables id_num / out=out_id_0 (keep=id_num);
run;

data out_id_0_expand;
 set out_id_0;
 do t2=0 to 1800 by 30;
  output;
 end;
run;

proc sort data=match_expand_0a; by id_num t2; run;
proc sort data=out_id_0_expand; by id_num t2; run;

data match_expand_0b;
 merge match_expand_0a out_id_0_expand; by id_num t2; 
if (dPhi_0_t= . ) then dPhi_0_t=0;
run;

data match_expand_0c;
 set match_expand_0b; by id_num;
if first.id_num then Phi_0_t=0;
retain Phi_0_t;
Phi_0_t=Phi_0_t+dPhi_0_t;
drop dPhi_0_t;
run;

data match_expand_0d;
 set match_expand_0c; 
if mod(t2,30)=0;
Phi_0_t_sq=Phi_0_t**2;
run;

proc sort data=match_expand_0d; by t2; run;
proc means data=match_expand_0d noprint; 
 var Phi_0_t_sq;
 output out=out_phi_0 sum(Phi_0_t_sq)=V0;
 by t2; 
run;

data out_phi_0a;
 set out_phi_0;
SD_0_t=sqrt(V0);
keep t2 SD_0_t;
run;

data covset0; TX=0; run;
proc phreg data=DE98SMA.match_expand_0 noprint;  
 model t2*deadm(0) =  / entry=t1 ties=breslow;
 strata TX;
 weight Wi_X;
 baseline out=outS0 survival=S0 cumhaz=Lam0 covariates=covset0 / nomean method=CH;
run;

data month1; 
 do t2=0 to 1800 by 30;
  output;
 end;
run;

data outS0a;
 merge outS0 month1; by t2;
if (_N_=1) then do; 
 Lam0_=0;  S0_=1;
end;
retain Lam0_ S0_;
if (Lam0= . ) then Lam0=Lam0_;
if (S0= . ) then S0=S0_;
Lam0_=Lam0;
S0_=S0;
if (TX= . ) then TX=0;
run;

data outS0b;
 set outS0a;
if mod(t2,30)=0;
run;


proc sort data=match_expand_0d; by t2; run;
data match_expand_0e;
 merge match_expand_0d outS0b; by t2;
Sphi_0_t=-S0*phi_0_t;
run;

proc sort data=out_phi_0a; by t2; run;
proc sort data=outS0b; by t2; run;

data outS0c;
 merge outS0b out_phi_0a; by t2; 
Lam0_low=Lam0-1.96*SD_0_t;
Lam0_up=Lam0+1.96*SD_0_t;
S0=exp(-Lam0);
S0_low=exp(-Lam0_up);
S0_up=exp(-Lam0_low);
drop S0_ Lam0_;
run;

proc print data=outS0c;
 var t2 S0 S0_low S0_up Lam0;
run;

data outfreq4tx; set outfreq4; m_j=1; run;

data tx5;
 set tx4;
deadm=dead1*(Xm <= &t_D_max);
Xm=min(X1,&t_D_max);
m_j=1;
num=1;
keep id_num match_id num m_j Xm deadm t_tx;
run;

proc freq data=tx5 noprint;
 where (deadm=1);
 tables Xm / out=outdead_1 (keep=Xm);
run;


%macro match_expand1(j_tot);

options nonotes;

%do j=1 %to &j_tot;

data outdead_1a;
if (_n_=1) then set outfreq4tx (where=(match_id=&j));
set outdead_1;
do num=1 to m_j;
 output;
end;
drop m_j;
run;

data tx5_j;
 set tx5;
if (match_id=&j);
run;
proc sort data=tx5_j; by match_id num; run;
proc sort data=outdead_1a; by match_id num; run;

data outdead_1b;
 merge outdead_1a 
       tx5_j (in=in_j keep=id_num m_j t_tx match_id num Xm rename=(Xm=Xm_)); 
 by match_id num; 
if (in_j=1);
if (0 <= Xm <= Xm_);
X=t_tx+Xm;
keep match_id num m_j id_num t_tx Xm X;
run;

proc sort data=outdead_1b; by X; run;
proc sort data=baseoutC1; by X; run;

data outdead_1c;
 merge outdead_1b (in=in_b) baseoutC1; by X;
if (in_b=1);
run;

proc sort data=outdead_1c; by id_num; run;
proc sort data=xbeta_C; by id_num; run;

data outdead_1d;
 merge outdead_1c (in=in_) xbeta_C; by id_num;
if (in_=1);
Wi_X=exp(Lam0C*exp(xbeta_C));
drop Lam0C xbeta_C;
run;

proc sort data=outdead_1d; by match_id num Xm; run;
proc sort data=tx5; by match_id num Xm; run;

data outdead_1e;
 merge outdead_1d (in=in_) tx5 (keep=deadm match_id num Xm); by match_id num Xm; 
if (in_=1);
if (deadm= . ) then deadm=0;
t2=Xm;
t1=t2-0.1;
dNw=deadm*Wi_X;
TX=1;
run;

%if (&j=1) %then %do;
 data DE98SMA.match_expand_1; set outdead_1e; run;
%end; %else %do;
 proc append base=DE98SMA.match_expand_1 data=outdead_1e; run;
%end;

%end;

options notes;

%mend match_expand1;

%match_expand1(j_tot=&m_id_max);

proc sort data=DE98SMA.match_expand_1; by t2; run;
proc means data=DE98SMA.match_expand_1 noprint; 
 var dNw Wi_X;
 output out=out_NY_1 mean(dNw)=dNw_bar mean(Wi_X)=Yw_bar sum(Wi_X)=Yw;
 by t2; 
run;

data out_NY_1a;
 set out_NY_1;
dLam_1=dNw_bar/Yw_bar;
TX=1;
keep t2 dLam_1 dNw_bar Yw_bar Yw TX;
run;

data match_expand_1a;
 merge DE98SMA.match_expand_1 out_NY_1a; by t2;
dPhi_1_t=Wi_X*(deadm-dLam_1)/Yw;
keep id_num t2 dPhi_1_t;
run;

proc freq data=match_expand_1a noprint; 
 tables id_num / out=out_id_1 (keep=id_num);
run;

data out_id_1_expand;
 set out_id_1;
 do t2=0 to 1800 by 30;
  output;
 end;
run;

proc sort data=match_expand_1a; by id_num t2; run;
proc sort data=out_id_1_expand; by id_num t2; run;

data match_expand_1b;
 merge match_expand_1a out_id_1_expand; by id_num t2; 
if (dPhi_1_t= . ) then dPhi_1_t=0;
run;

data match_expand_1c;
 set match_expand_1b; by id_num;
if first.id_num then Phi_1_t=0;
retain Phi_1_t;
Phi_1_t=Phi_1_t+dPhi_1_t;
drop dPhi_1_t;
run;

data match_expand_1d;
 set match_expand_1c; 
if mod(t2,30)=0;
Phi_1_t_sq=Phi_1_t**2;
run;

proc sort data=match_expand_1d; by t2; run;
proc means data=match_expand_1d noprint; 
 var Phi_1_t_sq;
 output out=out_phi_1 sum(Phi_1_t_sq)=V1;
 by t2; 
run;

data out_phi_1a;
 set out_phi_1;
SD_1_t=sqrt(V1);
keep t2 SD_1_t;
run;
proc print; run;

data covset1; TX=1; run;
proc phreg data=DE98SMA.match_expand_1 noprint;  
 model t2*deadm(0) =  / entry=t1 ties=breslow;
 strata TX;
 weight Wi_X;
 baseline out=outS1 survival=S1 cumhaz=Lam1 covariates=covset1 / nomean method=CH;
run;

data month1; 
 do t2=0 to 1800 by 30;
  output;
 end;
run;

data outS1a;
 merge outS1 month1; by t2;
if (_N_=1) then do; 
 Lam1_=0;  S1_=1;
end;
retain Lam1_ S1_;
if (Lam1= . ) then Lam1=Lam1_;
if (S1= . ) then S1=S1_;
Lam1_=Lam1;
S1_=S1;
if (TX= . ) then TX=1;
run;

data outS1b;
 set outS1a;
if mod(t2,30)=0;
run;

proc sort data=match_expand_1d; by t2; run;
data match_expand_1e;
 merge match_expand_1d outS1b; by t2;
Sphi_1_t=-S1*phi_1_t;
run;

proc sort data=out_phi_1a; by t2; run;
proc sort data=outS1b; by t2; run;

data outS1c;
 merge outS1b out_phi_1a; by t2; 
Lam1_low=Lam1-1.96*SD_1_t;
Lam1_up=Lam1+1.96*SD_1_t;
S1=exp(-Lam1);
S1_low=exp(-Lam1_up);
S1_up=exp(-Lam1_low);
drop S1_ Lam1_;
run;

proc print data=outS1c;
 var t2 S1 S1_low S1_up Lam1;
run;

proc sort data=match_expand_0e; by id_num t2; run;
proc sort data=match_expand_1e; by id_num t2; run;

data match_expand_e_;
 merge match_expand_0e match_expand_1e; by id_num t2; 
if (Sphi_0_t= . ) then Sphi_0_t=0;
if (Sphi_1_t= . ) then Sphi_1_t=0;
Sphi_diff_t=Sphi_1_t-Sphi_0_t;
Sphi_diff_t_sq=Sphi_diff_t**2;
run;

proc sort data=match_expand_e_; by t2; run;
proc means data=match_expand_e_ noprint; 
 var Sphi_diff_t_sq;
 output out=out_phi sum(Sphi_diff_t_sq)=V_diff;
 by t2; 
run;

data out_phi_a;
 set out_phi;
SD_diff_t=sqrt(V_diff);
keep t2 SD_diff_t;
run;
proc print; run;

data diff1;
 merge out_phi_a outS0b outS1b; by t2;
S_diff=S1-S0;
S_diff_low=S_diff-1.96*SD_diff_t;
S_diff_up=S_diff+1.96*SD_diff_t;

if (Lam0 > 0) then CHR_t=Lam1/Lam0; else CHR_t=1; 
run;

proc print data=diff1;
 var t2 S_diff S_diff_low S_diff_up CHR_t;
run;

