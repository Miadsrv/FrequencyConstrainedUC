
***************************************************************
*** INPUT
***************************************************************
$include Input_Data_for_UC.gms
alias (t,tt);
alias (i,ii,iii);

***************************************************************
*** Analytical FCUC sets and parameters
***************************************************************
set       x1r  x1 segments      /x1r1*x1r6/;
set       x2r  x2 segments      /x2r1*x2r6/;
set       pr   pr segments      /pr1*pr6/;
set       leaf number of leaves /l1*l3/;
set       p2r  p2 segments      /p2r1*p2r5/;

parameter ax1(x1r) linearization intervals of x1^2  /x1r1 5, x1r2 9, x1r3 13, x1r4 17, x1r5 21, x1r6 26/;
parameter ax2(x2r) linearization intervals of x2^2  /x2r1 3, x2r2 7, x2r3 11, x2r4 15, x2r5 19, x2r6 22/;
parameter ap (pr)  linearization intervals of p^2   /pr1 0, pr2 2.5, pr3 5, pr4 7.5, pr5 10, pr6 12.5/;
parameter a2p(p2r) getting from p2                  /p2r1 10, p2r2 40, p2r3 70, p2r4 100,  p2r5 130/;

scalar    alpha     normalizing parameter                  /0.2/;
scalar    betta     normalizing parameter                  /0.2/;
scalar    qss_max   maximum qss frequency range            /0.5/;
scalar    max_rocof maximum rocof allowed                  /2.5/;
scalar    MB        big positive number                    /2000/;
scalar    MP        big positive number                    /500/;
scalar    MM        big positive number                    /10/;
scalar    sBase     sum of all s bases of units            /127.57/;
scalar    fn        nominal frequency                      /50/;
scalar    Td        time constant of generators            /8.26/;
scalar    Ta        approximated time for reserve delivery /4.5/;
scalar    Da        load damping factor                    /0.01/;


***************************************************************
*** SETTINGS
***************************************************************
*$set modelname  ACFCUC
*$set season     summer
*$set day        d4
*$set deltaF_max 2.5
*$set cost       20


***************************************************************
*** VARIABLES
***************************************************************
variable          obj                     objective function variable

variable          c_aux        (t)        auxilliary variable

variable          co_gen       (t,i)      generation cost in each time period

variable          co_su        (t,i)      start-up cost in each time period

positive variable g            (t,i)      generator outputs

positive variable g_lin        (t,i,b)    generator block outputs

binary variable   suc          (t,i,j)    start up cost

binary variable   x            (t,i)      binary variable equal to 1 if generator is producing and 0 otherwise

binary variable   y            (t,i)      binary variable equal to 1 if generator is start-up and 0 otherwise
 
binary variable   z            (t,i)      binary variable equal to 1 if generator is shut-down and 0 otherwise

positive variable inertia_AO   (t,i)      inertia

positive variable reserve_AO   (t,i)      reserve after outage

positive variable aux_pk       (t,i,i)    auxilary variable for linearization ok p*k

binary Variable   yx1          (t,i,x1r)  x1 segment selector

binary Variable   yx2          (t,i,x2r)  x2 segment selector

binary Variable   yp           (t,i,pr)   p segment selector

positive variable lambdaX1     (t,i,x1r)  lambda multiplier for x1

positive variable lambdaX2     (t,i,x2r)  lambda multiplier for x2

positive variable lambdaP      (t,i,pr)   lambda multiplier for p

variable          aux_1        (t,i)      aux_1

variable          aux_2        (t,i)      aux_2

variable          ndr          (t,i)      frequency nadir

binary variable   ul           (leaf,t,i) binary variables corresponding to leaves

positive variable k_AO         (t,i)      weighted gain after outage of each generator

positive variable p_loss       (t,i)      lost generation

positive variable re           (t,i)      available reserve after outage

variable          lr1          (t,i)      logistic regression polynomial1 for the first classification

variable          lr2          (t,i)      logistic regression polynomial2 for the second classification

variable          linr         (leaf,t,i) linear regression polynomial of second leaf

variable          r            (leaf,t,i) aux variables to linearize linear regression

positive variable ufls         (t,i)      estimated ufls amount of each outage in each time step

binary variable   y2p          (t,i,p2r)  p2 segment selector

positive variable lambda2P     (t,i,p2r)  lambda multiplier for p2

positive variable lambda2P_crit(t,i,p2r)  lambda p crit

positive variable P_crit       (t,i)      p cirt

positive variable someError1              error of approximation

positive variable someError2              error of approximation

variable          p2_crit      (t,i)      square of p crit

binary variable   u            (t,i)      sign of p-p_crit

positive variable p_ufls       (t,i)      amount of ufls

binary variable   x_matrix     (t,i,ii)   2D commitment binary variable

positive variable pos_aux      (t,i)      positive auxilary variable

positive variable reserve      (t,i)      reserve

variable abtest(t,i);
variable ab(t,i);
variable aa(t,i);
variable ac(t,i);



***************************************************************
*** EQUATION DECLARATION
***************************************************************
equations
cost                                  objective function
cost_aux                   (t)        auxilliary equation
bin_set1                   (t,i)      setting start-up binary variables
bin_set10                  (t,i)      setting start-up binary variables
bin_set2                   (t,i)      setting start-up binary variables
gen_sum                    (t,i)      summing the generation of blocks per generator
gen_min                    (t,i)      genertor minimum output
cost_gen                   (t,i)      generation cost summation
cost_su                    (t,i)      start-up cost summation
block_output               (t,i,b)    limiting the output of each generator block
min_updown_1               (t,i)      minimum updown time constraint 1
min_updown_2               (t,i)      minimum updown time constraint 2
min_updown_3               (t,i)      minimum updown time constraint 3
ramp_limit_min             (t,i)      ramp-down limit
ramp_limit_max             (t,i)      ramp-up limit
ramp_limit_min_1           (i)        ramp-down limit for the first time period
ramp_limit_max_1           (i)        ramp-up limit for the first time period
start_up_cost1             (t,i,j)    stairwise linear cost function - equation 1
start_up_cost2             (t,i)      stairwise linear cost function - equation 2
power_balance              (t)        power balance for each bus
reserve_requirement_thermal(t,i)      reserve requirement
inertia_eq                 (t,i)      inertia after outage to be used in rocof constraint
rocof                      (t,i)      rocof constraint
simp                       (t,i)      simplifying similar generators
simp2                      (t,i)      simplifying similar generators
pk_eq1                     (t,i,ii)   indicator constraint
pk_eq2                     (t,i,ii)   indicator constraint
req_power_eq               (t,i,ii)   new reserve requirement

*Analytical PFCUC constraints
qss                        (t,i)      quasi-steady-state limit
aux_1_eq                   (t,i)      calculating aux_1
aux_2_eq                   (t,i)      calculating aux_2
variable_change_1          (t,i)      variable change
variable_change_2          (t,i)      variable change
p_eq                       (t,i)      calculating p_eq
lambdaX1_sum               (t,i)      constraints for linearization of square term
lambdaX2_sum               (t,i)      constraints for linearization of square term
lambdaP_sum                (t,i)      constraints for linearization of square term
yx1_sum                    (t,i)      constraints for linearization of square term
yx2_sum                    (t,i)      constraints for linearization of square term
yp_sum                     (t,i)      constraints for linearization of square term
lambdaX1_lim               (t,i,x1r)  constraints for linearization of square term
lambdaX2_lim               (t,i,x2r)  constraints for linearization of square term
lambdaP_lim                (t,i,pr)   constraints for linearization of square term
available_reserve          (t,i)      available reserve after outage
nadir                      (t,i)      resulting nadir constraint

*Machine learning based PFCUC
nadirml3p5                 (t,i)      nadir estimation with logistic regression function for threshold set at 3.5 Hz
nadirml3p4                 (t,i)      nadir estimation with logistic regression function for threshold set at 3.4 Hz
nadirml3p3                 (t,i)      nadir estimation with logistic regression function for threshold set at 3.3 Hz
nadirml3p2                 (t,i)      nadir estimation with logistic regression function for threshold set at 3.2 Hz
nadirml3p1                 (t,i)      nadir estimation with logistic regression function for threshold set at 3.1 Hz
nadirml3                   (t,i)      nadir estimation with logistic regression function for threshold set at 3 Hz
nadirml2p9                 (t,i)      nadir estimation with logistic regression function for threshold set at 2.9 Hz
nadirml2p8                 (t,i)      nadir estimation with logistic regression function for threshold set at 2.8 Hz
nadirml2p7                 (t,i)      nadir estimation with logistic regression function for threshold set at 2.7 Hz
nadirml2p6                 (t,i)      nadir estimation with logistic regression function for threshold set at 2.6 Hz
nadirml2p5                 (t,i)      nadir estimation with logistic regression function for threshold set at 2.5 Hz
nadirml2p4                 (t,i)      nadir estimation with logistic regression function for threshold set at 2.4 Hz
nadirml2p3                 (t,i)      nadir estimation with logistic regression function for threshold set at 2.3 Hz

*Machine learning based CFCUC
req_power_eq_UFLS          (t,i,ii)   new reserve requirement modified with ufls
**Logistic regression contraints
leaf_binary_limit                     assigning to only one leaf
branch1_1                             branch1 constraint1
branch1_2                             branch1 constraint2
branch2_1                             branch2 constraint1
branch2_2                             branch2 constraint2
k_eq                       (t,i)      calculating weighted gain after outage of unit i
loss_eq                    (t,i)      calculating power loss after outage of unit i
avai_reserve               (t,i)      calculating available reserve after outage of unit i
lr_eq1                     (t,i)      logistic regression1 polynomial
lr_eq2                     (t,i)      logistic regression2 polynomial

**Linear regression constraints 
linr_pol1                  (t,i)      linear regression polynomial of leaf 1
linr_pol2                  (t,i)      linear regression polynomial of leaf 2
linr_pol3                  (t,i)      linear regression polynomial of leaf 3
linr_eq1                   (leaf,t,i) resulting constraint for linearizing linear regression constraint
linr_eq2                   (leaf,t,i) resulting constraint for linearizing linear regression constraint
linr_eq3                   (leaf,t,i) resulting constraint for linearizing linear regression constraint
linr_eq4                   (leaf,t,i) resulting constraint for linearizing linear regression constraint
linr_eq5                   (t,i)      resulting constraint for linearizing linear regression constraint
cost_aux_ufls              (t)        new auxilary cost that contains ufls cost

*Analytical CFCUC constraints
cost_aux_aufls             (t)        new auxilary cost for analytical ufls estimation that contains ufls cost
lambda2P_sum               (t,i)      weighting the coefficients for power variable
y2p_sum                    (t,i)      finding the right spine for power variable
lambda2P_lim               (t,i,p2r)  forcing power variable to be presented only by two lambdas
hourly_reserve             (t,i)      the amount of reserve of unit i in hour t
P_crit_eq1                 (t,i)      p crit equation 1
ufls_pos1                  (t,i)
ufls_pos2                  (t,i)
x2D_eq1                    (t,i,i)
x2D_eq2                    (t,i,i)
x2D_eq3                    (t,i,i)
x2D_eq4                    (t,i,i)
x2D_eq5                    (t,i,i)
P_crit_quad_eq             (t,i)
req_power_eq_pcrit         (t,i,i)
ufls_eq                    (t,i)
; 


***************************************************************
*** EQUATIONS
***************************************************************
cost..
         obj =e= sum(t,c_aux(t));

cost_aux(t)..
         c_aux(t) =e= sum(i,co_gen(t,i)+co_su(t,i));
         
cost_aux_ufls(t)..
         c_aux(t) =e= sum(i,co_gen(t,i)+co_su(t,i)+(ufls(t,i)*%cost%*2/(24*365)));
         
*cost_aux_aufls(t)..
*         c_aux(t) =e= sum(i,co_gen(t,i)+co_su(t,i)+(p_ufls(t,i)*%cost%*2/(24*365)));

bin_set1(t,i)$(ord(t) gt 1)..
         y(t,i) - z(t,i) =e= x(t,i) - x(t-1,i);

bin_set10(t,i)$(ord(t) = 1)..
         y(t,i) - z(t,i) =e= x(t,i) - onoff_t0(i);

bin_set2(t,i)..
         y(t,i) + z(t,i) =l= 1;

cost_gen(t,i)..
         co_gen(t,i) =e= a(i)*x(t,i) + sum(b,g_lin(t,i,b)*k(i,b));
         
cost_su(t,i)..
         co_su(t,i) =e= sum(j,suc_sw(i,j)*suc(t,i,j));

gen_sum(t,i)..
         g(t,i) =e= sum(b,g_lin(t,i,b));

gen_min(t,i)..
         g(t,i) =g= g_min(i)*x(t,i);

block_output(t,i,b)..
         g_lin(t,i,b) =l= g_max(i,b)*x(t,i);
         
min_updown_1(t,i)$(L_up_min(i)+L_down_min(i) gt 0 and ord(t) le L_up_min(i)+L_down_min(i))..
         x(t,i) =e= onoff_t0(i);

min_updown_2(t,i)..
         sum(tt$(ord(tt) ge ord(t)-g_up(i)+1 and ord(tt) le ord(t)),y(tt,i)) =l= x(t,i);

min_updown_3(t,i)..
         sum(tt$(ord(tt) ge ord(t)-g_down(i)+1 and ord(tt) le ord(t)),z(tt,i)) =l= 1-x(t,i);

ramp_limit_min(t,i)$(ord(t) gt 1)..
         -ramp_down(i) =l= g(t,i) - g(t-1,i);

ramp_limit_max(t,i)$(ord(t) gt 1)..
         ramp_up(i) =g= g(t,i) - g(t-1,i);

ramp_limit_min_1(i)..
         -ramp_down(i) =l= g('t1',i) - g_0(i);

ramp_limit_max_1(i)..
         ramp_up(i) =g= g('t1',i) - g_0(i);

start_up_cost1(t,i,j)..
         suc(t,i,j) =l= sum(tt$(ord(tt) lt ord(t) and ord(tt) ge suc_sl(i,j) and ord(tt) le suc_sl(i,j+1)-1),z(t-ord(j),i))+
         1$(ord(j) lt card(j) and count_off_init(i)+ord(t)-1 ge suc_sl(i,j) and count_off_init(i)+ord(t)-1 lt suc_sl(i,j+1))+
         1$(ord(j) = card(j) and count_off_init(i)+ord(t)-1 ge suc_sl(i,j));

start_up_cost2(t,i)..
         sum(j,suc(t,i,j)) =e= y(t,i);

power_balance(t)..
         sum(i,g(t,i)) + wind_det(t,"%day%","%season%")+solar_det(t,"%day%","%season%") =e= d(t,"%day%","%season%");
        
reserve_requirement_thermal(t,i)..
         sum(ii$(ord(ii)<>ord(i)),sum(b,g_max(ii,b)*x(t,ii)-g_lin(t,ii,b)))=g=g(t,i);
                 
inertia_eq(t,i)..
         inertia_AO(t,i)=e=sum(ii$(ord(ii) <> ord(i)),inertia(ii)*mbase(ii)*x(t,ii));
        
rocof(t,i)$(ord(i) <> 12)..
         2*inertia_AO(t,i)*max_rocof/fn=g=g(t,i);
        
simp(t,i)$ (ord(i) ge 1 and ord(i) le 3) ..
         x(t,i) =g= x(t,i+1);
        
simp2(t,i)$ (ord(i) ge 8 and ord(i) le 10) ..
         x(t,i) =g= x(t,i+1);
        
req_power_eq(t,i,ii)..
         +kg(i)*mbase(i)*g(t,ii)/Td
         =l=sum(iii$(ord(iii) <> ord(ii)),aux_pk(t,i,iii))+MB*(1-x(t,i));

req_power_eq_UFLS(t,i,ii)..
         +kg(i)*mbase(i)*(g(t,ii)-ufls(t,ii))/Td
         =l=sum(iii$(ord(iii) <> ord(ii)),aux_pk(t,i,iii))+MB*(1-x(t,i));
         
req_power_eq_pcrit(t,i,ii)..
         +kg(i)*mbase(i)*p_crit(t,ii)/Td
         =l=sum(iii$(ord(iii) <> ord(ii)),aux_pk(t,i,iii))+MB*(1-x(t,i));
         
pk_eq1(t,i,ii)..
        aux_pk(t,i,ii)=e=0;

pk_eq2(t,i,ii)..
        aux_pk(t,i,ii)=e=(sum(b,g_max(i,b))-g(t,i))*kg(ii)*mbase(ii)/Td;
       
aux_1_eq(t,i)..
         aux_1(t,i)=e=sum(x1r,ax1(x1r)*lambdaX1(t,i,x1r));
        
aux_2_eq(t,i)..
         aux_2(t,i)=e=sum(x2r,ax2(x2r)*lambdaX2(t,i,x2r));
        
p_eq(t,i)..
         g(t,i)=e=sum(pr,ap(pr)*lambdaP(t,i,pr));
        
lambdaX1_sum(t,i)..
         sum(x1r,lambdaX1(t,i,x1r))=e=1;

lambdaX2_sum(t,i)..
         sum(x2r,lambdaX2(t,i,x2r))=e=1;
        
lambdaP_sum(t,i)..
         sum(pr,lambdaP(t,i,pr))=e=x(t,i);
        
yx1_sum(t,i)..
         sum(x1r$(ord(x1r)>1),yx1(t,i,x1r))=e=1;
        
yx2_sum(t,i)..
         sum(x2r$(ord(x2r)>1),yx2(t,i,x2r))=e=1;
        
yp_sum(t,i)..
         sum(pr$(ord(pr)>1),yp(t,i,pr))=e=x(t,i);
        
lambdaX1_lim(t,i,x1r)..
         lambdaX1(t,i,x1r)=l=yx1(t,i,x1r)$(ord(x1r)>1)+yx1(t,i,x1r+1)$(ord(x1r)<card(x1r));
        
lambdaX2_lim(t,i,x2r)..
         lambdaX2(t,i,x2r)=l=yx2(t,i,x2r)$(ord(x2r)>1)+yx2(t,i,x2r+1)$(ord(x2r)<card(x2r));
        
lambdaP_lim(t,i,pr)..
         lambdaP(t,i,pr)=l=yp(t,i,pr)$(ord(pr)>1)+yp(t,i,pr+1)$(ord(pr)<card(pr));
        
nadir(t,i)..
         (sum(x1r,ax1(x1r)*ax1(x1r)*lambdaX1(t,i,x1r))-sum(x2r,ax2(x2r)*ax2(x2r)*lambdaX2(t,i,x2r)))/(alpha*betta)-(fn*Ta/(4*%deltaF_max%))*sum(pr,ap(pr)*ap(pr)*lambdaP(t,i,pr))
         +Da*d(t,"%day%","%season%")*Ta*g(t,i)*fn/4=g=0;         
         
nadirml3p5(t,i)..
         1.78257655+(-0.76970702*inertia_AO(t,i))+(0.08930268*k_AO(t,i))+(-6.51773368*g(t,i))+(5.51897807*re(t,i))=g=0;
         
nadirml3p4(t,i)..
         1.236353+(-0.68281566*inertia_AO(t,i))+(0.08152437*k_AO(t,i))+(-5.94906082*g(t,i))+(4.78202908*re(t,i))=g=0;

nadirml3p3(t,i)..
         0.00068905+(-0.5684104*inertia_AO(t,i))+(0.07126977*k_AO(t,i))+(-4.91391179*g(t,i))+(3.5711342*re(t,i))=g=0;
         
nadirml3p2(t,i)..
         -0.01823947+(-0.51462321*inertia_AO(t,i))+(0.06750551*k_AO(t,i))+(-4.79446214*g(t,i))+(3.17206649*re(t,i))=g=0;
         
nadirml3p1(t,i)..
         -1.35732672+(-0.48787124*inertia_AO(t,i))+(0.06737791*k_AO(t,i))+(-4.49305877*g(t,i))+(2.64221023*re(t,i))=g=0;
         
nadirml3(t,i)..
         -1.97598845+(-0.456884*inertia_AO(t,i))+(0.06643285*k_AO(t,i))+(-4.33310143*g(t,i))+(2.20196743*re(t,i))=g=0;
         
nadirml2p9(t,i)..
         -2.0623524+(-0.41219794*inertia_AO(t,i))+(0.06580827*k_AO(t,i))+(-4.67333309*g(t,i))+(2.00825938*re(t,i))=g=0;
         
nadirml2p8(t,i)..
         -2.72904767+(-0.40840966*inertia_AO(t,i))+(0.06817966*k_AO(t,i))+(-4.74264797*g(t,i))+(1.72647821*re(t,i))=g=0;
         
nadirml2p7(t,i)..
         -2.46151347+(-0.38630606*inertia_AO(t,i))+(0.06921278*k_AO(t,i))+(-5.1811914*g(t,i))+(1.61108176*re(t,i))=g=0;
         
nadirml2p6(t,i)..
         -2.95403123+(-0.35680933*inertia_AO(t,i))+(0.0676894*k_AO(t,i))+(-5.15784091*g(t,i))+(1.3487068*re(t,i))=g=0;
         
nadirml2p5(t,i)..
         -2.3277676+(-0.34759675*inertia_AO(t,i))+(0.0669719*k_AO(t,i))+(-5.33628921*g(t,i))+(1.26857062*re(t,i))=g=0;
         
nadirml2p4(t,i)..
         -1.44664794+(-0.24919179*inertia_AO(t,i))+(0.05812038*k_AO(t,i))+(-5.63710384*g(t,i))+(1.15780389*re(t,i))=g=0;
         
nadirml2p3(t,i)..
         -0.11243555+(-0.23318268*inertia_AO(t,i))+(0.05465937*k_AO(t,i))+(-5.6795615*g(t,i))+(1.10199865*re(t,i))=g=0;
                
available_reserve(t,i)..
         reserve_AO(t,i)=e=sum(ii$(ord(ii)<>ord(i)),sum(b,g_max(ii,b)*x(t,ii)-g_lin(t,ii,b)));
         
variable_change_1(t,i)..
         (aux_1(t,i)+aux_2(t,i))/alpha=e=inertia_AO(t,i);

variable_change_2(t,i)..
         (aux_1(t,i)-aux_2(t,i))/betta=e=reserve_AO(t,i);
         
qss(t,i)..
         (aux_1(t,i)-aux_2(t,i))/betta=g=g(t,i)-Da*d(t,"%day%","%season%")*qss_max;
        
k_eq(t,i)..
         k_AO(t,i)=e=sum(ii$(ord(ii) <> ord(i)),kg(ii)*x(t,ii)*mbase(ii));
        
loss_eq(t,i)..
         p_loss(t,i)=e=sum(ii$(ord(ii) <> ord(i)),g(t,ii));
        
avai_reserve(t,i)..
         re(t,i)=e=sum(ii,sum(b,g_max(ii,b))*x(t,ii)-g(t,ii))-sum(b,g_max(i,b))*x(t,i)+g(t,i);
        
lr_eq1(t,i)..
         lr1(t,i)=e=0.34105642+(0.07195032*inertia_AO(t,i))+(-0.04314308*k_AO(t,i))+(8.90363545*g(t,i))+(-0.10818889*re(t,i));
        
lr_eq2(t,i)..
         lr2(t,i)=e=0.74227943+(0.05564296*inertia_AO(t,i))+(-0.0563435*k_AO(t,i))+(6.49898422*g(t,i))+(-0.08096212*re(t,i));
        
leaf_binary_limit(t,i)..
         sum(leaf,ul(leaf,t,i)) =e= 1;

branch1_1(t,i)..
         lr1(t,i)+(-MP*ul('l2',t,i)) =g= -MP;
        
branch1_2(t,i)..
         lr1(t,i)+(MP*ul('l3',t,i)) =l= MP;
        
branch2_1(t,i)..
         lr2(t,i)+(-MP*ul('l1',t,i)) =g= -MP;
        
branch2_2(t,i)..
         lr2(t,i)+(MP*ul('l2',t,i)) =l= MP;
        
linr_pol3(t,i)..
         linr('l3',t,i)=e=0;
        
linr_pol2(t,i)..
         linr('l2',t,i)=e=0.0466082814874369+(0.0185344405*inertia_AO(t,i))+(-0.0000516254898*k_AO(t,i))+(0.103736850*g(t,i))+(-0.0492858305*re(t,i));
        
linr_pol1(t,i)..
         linr('l1',t,i)=e=0.30637117880673426+(0.0310096*inertia_AO(t,i))+(-0.00258244*k_AO(t,i))+(0.85428141*g(t,i))+(-0.16353837*re(t,i));
        
linr_eq1(leaf,t,i)..
         linr(leaf,t,i)-(MM*(1-ul(leaf,t,i)))=l=r(leaf,t,i);
        
linr_eq2(leaf,t,i)..
         linr(leaf,t,i)-(-MM*(1-ul(leaf,t,i)))=g=r(leaf,t,i);
        
linr_eq3(leaf,t,i)..
         MM*ul(leaf,t,i)=g=r(leaf,t,i);
        
linr_eq4(leaf,t,i)..
         -MM*ul(leaf,t,i)=l=r(leaf,t,i);
        
linr_eq5(t,i)..
         ufls(t,i)=e=sum(leaf,r(leaf,t,i));
                 
hourly_reserve(t,i)..
         reserve(t,i)=e=sum(b,g_max(i,b)*x(t,i))-g(t,i);
         
lambda2P_sum(t,i)..
         sum(p2r,lambda2P_crit(t,i,p2r))=e=1;
         
y2p_sum(t,i)..
         sum(p2r$(ord(p2r)>1),y2p(t,i,p2r))=e=1;
         
lambda2P_lim(t,i,p2r)..
         lambda2P_crit(t,i,p2r)=l=y2p(t,i,p2r)$(ord(p2r)>1)+y2p(t,i,p2r+1)$(ord(p2r)<card(p2r));
        
P_crit_eq1(t,i)..
         P_crit(t,i)=e=sum(p2r,sqrt(a2p(p2r))*lambda2P_crit(t,i,p2r));
        
ufls_pos1(t,i)..
         ufls(t,i)=l=150*u(t,i);
        
ufls_pos2(t,i)..
         ufls(t,i)=l=g(t,i)-P_crit(t,i)+150*(1-u(t,i));
         
ufls_eq(t,i)..
         ufls(t,i)=g=g(t,i)-P_crit(t,i);
         
x2D_eq1(t,i,ii)$(ord(i)>ord(ii))..
         x_matrix(t,i,ii)=l=x(t,i);
        
x2D_eq2(t,i,ii)$(ord(i)>ord(ii))..
         x_matrix(t,i,ii)=l=x(t,ii);

x2D_eq3(t,i,ii)$(ord(i)>ord(ii))..
         x_matrix(t,i,ii)=g=x(t,i)+x(t,ii)-1;
        
x2D_eq4(t,i,ii)$(ord(i)=ord(ii))..
         x_matrix(t,i,ii)=e=x(t,i);
        
x2D_eq5(t,i,ii)$(ord(i)<ord(ii))..
         x_matrix(t,i,ii)=e=x_matrix(t,ii,i);
         
P_crit_quad_eq(t,i)..
        sum(p2r,a2p(p2r)*lambda2P_crit(t,i,p2r))=e=
        (%deltaF_max%*%deltaF_max%/(fn*fn))*2*sum((ii,iii)$(ord(ii) <> ord(i) and
        ord(iii) <> ord(i)),x_matrix(t,ii,iii)*inertia(ii)*mbase(ii)*kg(iii)*mbase(iii))/Td;
        


***************************************************************
*** MODEL DEFINITIONS
***************************************************************
model UC simple UC
/
cost
cost_aux
bin_set1
bin_set10
bin_set2
gen_sum
gen_min
cost_gen
cost_su
block_output
min_updown_1
min_updown_2
min_updown_3
ramp_limit_min
ramp_limit_max
ramp_limit_min_1
ramp_limit_max_1
start_up_cost1
start_up_cost2
power_balance
reserve_requirement_thermal
simp
simp2
/;

model RCUC RoCoF constrained UC
/
UC
inertia_eq
rocof
/;

model ReCUC Reserve constrained UC
*******Check indicator constraints********
/
UC
-reserve_requirement_thermal
pk_eq1
pk_eq2
req_power_eq
hourly_reserve
/;

model APFCUC analytical preventive FCUC
/
UC
inertia_eq
rocof
qss
aux_1_eq
aux_2_eq
variable_change_1
variable_change_2
p_eq
lambdaX1_sum
lambdaX2_sum
lambdaP_sum
yx1_sum
yx2_sum
yp_sum
lambdaX1_lim
lambdaX2_lim
lambdaP_lim
available_reserve
nadir
/;

model MPFCUC machine learning based preventive FCUC
*******Choose only one nadir constraint based of the preferred threshold*******
/
UC
inertia_eq
rocof
qss
avai_reserve
k_eq
loss_eq
*nadirml3p5           
*nadirml3p4              
*nadirml3p3            
*nadirml3p2         
*nadirml3p1
*nadirml3    
*nadirml2p9
*nadirml2p8
*nadirml2p7
*nadirml2p6
*nadirml2p5
*nadirml2p4              
*nadirml2p3             
/;

model MCFCUC machine learning based corrective FCUC
*******Check indicator constraints********
/
UC
-reserve_requirement_thermal
-cost_aux
cost_aux_ufls
inertia_eq
rocof
pk_eq1
pk_eq2
req_power_eq_UFLS
leaf_binary_limit
branch1_1
branch1_2
branch2_1
branch2_2
k_eq
loss_eq
avai_reserve
lr_eq1
lr_eq2
linr_pol1
linr_pol2
linr_pol3
linr_eq1
linr_eq2
linr_eq3
linr_eq4
linr_eq5
/;

model ACFCUC analytical corrective FCUC
*******Check indicator constraints********
/
UC
-reserve_requirement_thermal
-cost_aux
cost_aux_ufls
inertia_eq
rocof
hourly_reserve,
lambda2P_sum,
y2p_sum,
lambda2P_lim,
k_eq
P_crit_eq1,
x2D_eq1,
x2D_eq2,
x2D_eq3,
x2D_eq4,
x2D_eq5,
P_crit_quad_eq,
req_power_eq_pcrit,
pk_eq1,
pk_eq2,
ufls_pos1,
ufls_pos2,
ufls_eq,
/;


***************************************************************
*** SOLVE
***************************************************************
option reslim = 100000;
option Savepoint=1;
option optcr=0.001;
%modelname%.optfile = 1;
%modelname%.limrow =0;
%modelname%.limcol =0;

$onecho > gurobi.opt
mipfocus 1
threads -1
$offecho
Option
mip = Gurobi;


****** Use the following section for models "ReCUC, ACFCUC, MCFCUC" ******
*$ontext
file fgrb Gurobi Option file / gurobi.opt /; 
loop((t,i,ii),
  put fgrb    
        'indic ' pk_eq1.tn(t,i,ii) '$' x.tn(t,ii) '0' /
        'indic ' pk_eq2.tn(t,i,ii) '$' x.tn(t,ii) '1' /
    );
putclose fgrb;
*$offtext
****** Use the upper section for models "ReCUC, ACFCUC, MCFCUC" **********


solve %modelname% using mip minimizing obj;

parameter time_elapsed, solver_time, var_num, dvar_num, eq_num, total_ufls;
*, nadirtest(t,i);
*,frequency_nadir_approx(t,i), frequency_nadir_exact(t,i);

time_elapsed = %modelname%.etSolver;
solver_time = %modelname%.resUsd;
var_num = %modelname%.numVar;
dvar_num = %modelname%.numDVar;
eq_num = %modelname%.numEqu;
*nadirtest(t,i)=3.6296717+(-0.71166133*inertia_AO.l(t,i))+(0.08045439*k_AO.l(t,i))+(-6.46111846*g.l(t,i))+(5.55345005*re.l(t,i))
*nadirtest(t,i)=-0.07355444+(-0.35844117*inertia_AO.l(t,i))+(0.05426479*k_AO.l(t,i))+(-4.26063478*g.l(t,i))+(2.1794269*re.l(t,i));
*frequency_nadir_approx(t,i)=(sum(pr,ap(pr)*ap(pr)*lambdaP.l(t,i,pr))*50*Ta)/(-4*((sum(x1r,ax1(x1r)*ax1(x1r)*lambdaX1.l(t,i,x1r))-sum(x2r,ax2(x2r)*ax2(x2r)*lambdaX2.l(t,i,x2r)))/(alpha*betta)+0.01*d(t,"%day%","%season%")*Ta*g.l(t,i)*50*0.25));
*frequency_nadir_exact(t,i)=(50*Ta*g.l(t,i)*g.l(t,i)/(4*inertia_AO.l(t,i)*re.l(t,i)+0.01*d(t,"%day%","%season%")*Ta*g.l(t,i)*50))*x.l(t,i);


$onEpsToZero
Execute_Unload '%modelname%_%cost%_%season%_%day%.gdx'