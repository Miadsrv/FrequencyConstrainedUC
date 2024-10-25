# FrequencyConstrainedUC
## _One model to rule them all_
This model is capable of solving four different frequency constrained UC formulations.
The main case study is an island power sytem.
For each formulation please read the corrsponding section and modify accordingly.

- .inc files are included in a seperate folder to avoid loading time for every run.
- Some outputs are shared in outputs folder.
- In the settings section, uncomment the following to choose the model, saeason, day, frequency nadir threshold, and UFLS cost (see the example below).
```sh
$set modelname  ACFCUC
$set season     summer
$set day        d4
$set deltaF_max 2.5
$set cost       20
```

## Analyical Preventive FCUC (APFCUC)

- Choose the desired frequency nadir threshold in the setting section.
- You may need to adjust piece-wise linearization segments for different data.

## Data-driven Preventive FCUC (MPFCUC)

- Choose the desired frequency nadir threshold by uncommenting only one nadirml constraint when defining model (see below).

```sh
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
```

- The dataset and the python code that is used to train each constraint is included in "training" folder.


## Analytical Corrective FCUC (ACFCUC)

- Choose the desired UFLS cost in the setting section.
- Uncomment the indicator constraint section for this model (see below).

```sh
*$ontext
file fgrb Gurobi Option file / gurobi.opt /; 
loop((t,i,ii),
  put fgrb    
        'indic ' pk_eq1.tn(t,i,ii) '$' x.tn(t,ii) '0' /
        'indic ' pk_eq2.tn(t,i,ii) '$' x.tn(t,ii) '1' /
    );
putclose fgrb;
*$offtext
```

- You may need to adjust piece-wise linearization segments for different data.
- Frequency nadir threshold used for this model is 2.5 Hz.


## Data-driven Corrective FCUC (MCFCUC)

- Choose the desired UFLS cost in the setting section.
- Uncomment the indicator constraint section for this model (see below).
```sh
*$ontext
file fgrb Gurobi Option file / gurobi.opt /; 
loop((t,i,ii),
  put fgrb    
        'indic ' pk_eq1.tn(t,i,ii) '$' x.tn(t,ii) '0' /
        'indic ' pk_eq2.tn(t,i,ii) '$' x.tn(t,ii) '1' /
    );
putclose fgrb;
*$offtext
```
- The dataset and the python code that is used to train each constraint is included in "training" folder.

## Contact
> Feel free to contact msarvarizadeh@comillas.edu with any sort of questions or feedback!





