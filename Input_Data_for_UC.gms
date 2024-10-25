***************************************************************
*** SETS
***************************************************************
set t time periods                              /t1*t24/;
set i generators                                /i1*i11/;
set b generator blocks                          /b1*b3/;
set j start up cost intervals                   /j1*j8/;
set column auxiliary                            /column1/;
set tr transmission data (1-HP 2-MOV 3-MS 4-YD) /tr1*tr4/;
set c cost function coefficients                /cons,lin,quad/;
set seasons seasons                             /winter,spring,summer,autumn/;
set days days                                   /d1*d7/;
set res renewable generation technologies       /Wind,PV/;

***************************************************************
*** OPTIONS
***************************************************************

scalar transmission_option /4/;
parameter ws_penalty       /0/;
parameter installed_wind,installed_solar;

$onEps

*values are in MW
installed_wind=7;
installed_solar=5.8;

***************************************************************
*** PARAMETERS
***************************************************************


*GENERATOR DATA


table g_max(i,b) generator block generation limit
*$call =xls2gms r=generators_c!be2:bh14 i=Input.xlsx o=block_max.inc
$include block_max.inc
;

table k_option(i,b,tr) slope of each generator cost curve block
*$call =xls2gms r=generators_c!aw2:bc290 i=Input.xlsx o=k.inc
$include k.inc
;
parameter k(i,b);
k(i,b)=sum(tr$(ord(tr)=transmission_option),k_option(i,b,tr));

table suc_sw_option(i,j,tr) generator stepwise start-up cost
*$call =xls2gms r=generators_c!bj2:bp770 i=Input.xlsx o=start_up_sw.inc
$include start_up_sw.inc
;
parameter suc_sw(i,j);
suc_sw(i,j)=sum(tr$(ord(tr)=transmission_option),suc_sw_option(i,j,tr));

table suc_sl(i,j) generator stepwise start-up hourly blocks
*$call =xls2gms r=generators_c!br2:bz14 i=Input.xlsx o=start_up_sl.inc
$include start_up_sl.inc
;

table aux2(i,column)
*$call =xls2gms r=generators_c!d2:e14 i=Input.xlsx o=aux2.inc
$include aux2.inc
;
parameter count_off_init(i) number of time periods each generator has been off;
count_off_init(i)=aux2(i,'column1');

table aux3(i,column)
*$call =xls2gms r=generators_c!g2:h14 i=Input.xlsx o=aux3.inc
$include aux3.inc
;
parameter count_on_init(i) number of time periods each generator has been on;
count_on_init(i)=sum(column,aux3(i,column));

table aux4(i,tr)
*$call =xls2gms r=generators_c!j2:N14 i=Input.xlsx o=aux4.inc
$include aux4.inc
;
parameter a(i) fixed operating cost of each generator;
a(i)=sum(tr$(ord(tr)=transmission_option),aux4(i,tr));

table aux5(i,tr)
*$call =xls2gms r=generators_c!p2:t14 i=Input.xlsx o=aux5.inc
$include aux5.inc
;
parameter ramp_up(i) generator ramp-up limit;
ramp_up(i)=sum(tr$(ord(tr)=transmission_option),aux5(i,tr));

table aux6(i,tr)
*$call =xls2gms r=generators_c!v2:z14 i=Input.xlsx o=aux6.inc
$include aux6.inc
;
parameter ramp_down(i) generator ramp-down limit;
ramp_down(i)=sum(tr$(ord(tr)=transmission_option),aux6(i,tr));

table aux7(i,tr)
*$call =xls2gms r=generators_c!ab2:af14 i=Input.xlsx o=aux7.inc
$include aux7.inc
;
parameter g_down(i) generator minimum down time;
g_down(i)=sum(tr$(ord(tr)=transmission_option),aux7(i,tr));

table aux8(i,tr)
*$call =xls2gms r=generators_c!ah2:al14 i=Input.xlsx o=aux8.inc
$include aux8.inc
;
parameter g_up(i) generator minimum up time;
g_up(i)=sum(tr$(ord(tr)=transmission_option),aux8(i,tr));

table aux9(i,tr)
*$call =xls2gms r=generators_c!an2:ar14 i=Input.xlsx o=aux9.inc
$include aux9.inc
;
parameter g_min(i) generator minimum output;
g_min(i)=sum(tr$(ord(tr)=transmission_option),aux9(i,tr));

table aux10(i,column)
*$call =xls2gms r=generators_c!at2:au14 i=Input.xlsx o=aux10.inc
$include aux10.inc
;
parameter g_0(i) generator generation at t=0;
g_0(i)=sum(column,aux10(i,column));

parameter onoff_t0(i) on-off status at t=0;
onoff_t0(i)$(count_on_init(i) gt 0) = 1;

count_on_init(i)=sum(column,aux3(i,column));
g_0(i)=sum(column,aux10(i,column));
onoff_t0(i)$(count_on_init(i) gt 0) = 1;
count_off_init(i)=sum(column,aux2(i,column));

parameter L_up_min(i) used for minimum up time constraints;
L_up_min(i) = min(card(t), (g_up(i)-count_on_init(i))*onoff_t0(i));

parameter L_down_min(i) used for minimum up time constraints;
L_down_min(i) = min(card(t), (g_down(i)-count_off_init(i))*(1-onoff_t0(i)));

table cost_coefs(i,c) cost function coefficients
*$call =xls2gms r=generators_c!cb2:ce14 i=Input.xlsx o=cost_coefs.inc
$include cost_coefs.inc
;
table auxinertia(i,column)
*$call =xls2gms r=dynamics!a2:b14 i=Input.xlsx o=auxinertia.inc
$include auxinertia.inc
;
parameter inertia(i) inertia of generators (seconds);
inertia(i)=sum(column,auxinertia(i,column));

table auxmbase(i,column)
*$call =xls2gms r=dynamics!d2:e14 i=Input.xlsx o=auxmbase.inc
$include auxmbase.inc
;
parameter mbase(i) mbase of generators (mw);
mbase(i)=sum(column,auxmbase(i,column));

table auxk(i,column)
*$call =xls2gms r=dynamics!g2:h14 i=Input.xlsx o=auxk.inc
$include auxk.inc
;
parameter kg(i) k of generators (R^-1);
kg(i)=sum(column,auxk(i,column));
scalar M number of hours a unit can be on or off /2600/;


*RES DATA
table RES_det(t,days,seasons) available res generation
*$call =xls2gms r=RES!a2:g170 i=Input.xlsx o=RES_seasons.inc
$include RES_seasons.inc
;
table wind_det_aux(t,days,seasons,column)
*$call =xls2gms r=RES!i2:N674 i=Input.xlsx o=aux_wind_seasons.inc
$include aux_wind_seasons.inc
;
parameter wind_det(t,days,seasons);
wind_det(t,days,seasons)=sum(column,wind_det_aux(t,days,seasons,column));

table solar_det_aux(t,days,seasons,column)
*$call =xls2gms r=RES!p2:u674 i=Input.xlsx o=solar_seasons.inc
$include solar_seasons.inc
;
parameter solar_det(t,days,seasons);
solar_det(t,days,seasons)=sum(column,solar_det_aux(t,days,seasons,column));


*DEMAND DATA
table d(t,days,seasons) hourly demand
*$call =xls2gms r=load!a2:g170 i=Input.xlsx o=d_seasons.inc
$include d_seasons.inc
;


