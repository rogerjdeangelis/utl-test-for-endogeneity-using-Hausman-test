# utl-test-for-endogeneity-using-Hausman-test
Test for endogeneity using Hausman test
    Test for endogeneity using Hausman test

    SAS and WPS interface with R (SAS/IML/R WPS/PROC R)

    github
    https://tinyurl.com/y3mwhu9g
    https://github.com/rogerjdeangelis/utl-test-for-endogeneity-using-Hausman-test

    For documentation of this example see
    http://staff.utia.cas.cz/barunik/files/AE/Seminar7/sem_06.html

    SAS Forum
    https://communities.sas.com/t5/SAS-Procedures/please/m-p/534261

    as a side not
    See application that uses 'Vietnam Draft Lottery'
    for instrumental variable.

    INPUT (for building two models)
    ===============================

    options validvarname=upcase;
    libname sd1 "d:/sd1";
    data sd1.have;
     call streaminit(1234);
     do rec=1 to 1500;
        z1=rand("norm",2,1  );
        z2=rand("norm",1.5,1);
        z3=rand("norm",4,1  );
        u1=rand("norm");
        u2=rand("norm");
        x3=rand("norm");
        e1=rand("norm");
        e2=rand("norm");
        output;
     end;
     drop rec;
    run;quit;

    SD1.HAVE total obs=1,500

    Obs     Z1      Z2       Z3      U1       U2       X3       E1       E2

       1    2.87     2.31    2.96    -1.06    -1.75     0.57    -1.46     0.26
       2    1.35     1.92    3.94    -1.83    -0.07     0.45    -0.60    -0.30
       3    2.33     2.09    4.94    -0.56    -0.95     0.78    -0.63     0.10
    ...
    1498    2.14     1.03    4.78     0.49     0.80     0.72    -0.98    -0.90
    1499    1.78     0.02    4.50    -0.07    -0.30     0.13     1.81     1.31
    1500    1.43     1.46    3.94     1.09     0.49    -1.13    -0.17    -0.45

    EXAMPLE OUTPUT
    --------------

     WORK.WANT total obs=3

      V1

      Testing Ho: difference in coefficients not systematic
      Hausman statistics 198.9751 compared with Chi2(4)=9.488
      p-value (Prob>chi2) 0.0000


    PROCESS  (see link above for detailed documentation)
    =====================================================

    %utlfkil(d:/xpt/hausman.xpt); * just in case;

    %utl_submit_r64('
    library(dplyr);
    library(gmm);
    library(systemfit);
    library(AER);
    library(SASxport);
    library(haven);
    have<-read_sas("d:/sd1/have.sas7bdat");
    head(have);
    z1=have$Z1;
    z2=have$Z2;
    z3=have$Z3;
    u1=have$U1;
    u2=have$U2;
    x3=have$X3;
    e1=have$E1;
    e2=have$E2;
    x1 = .5*z1 - z2 + z3 + u1 + e1;
    x2 = 1.5*z1 - 2*z2 + u2 + e2;
    y = 1 + 2*x1 - x2 + x3 + 2*u1 - 1.5*u2;
    X=cbind(1,x1,x2,x3);
    colnames(X)=c("const","beta_1","beta_2","beta_3");
    betaOLS=solve(t(X)%*%X)%*%t(X)%*%y;
    betaOLS;
    X=cbind(1,x1,x2,x3);
    Z=cbind(1,x3,z1,z2);
    colnames(X)=c("const","beta_1","beta_2","beta_3");
    betaIV=solve(t(Z)%*%X)   %*%   t(Z)%*%y;
    betaIV;
    betaOLSlm=lm(y~x1+x2+x3);
    summary(betaOLSlm);
    X=cbind(1,x1,x2,x3);
    Z=cbind(1,x3,z1,z2,z3);
    colnames(X)=c("const","beta_1","beta_2","beta_3");
    beta2SLS=solve(t(X)%*%Z %*%solve(t(Z)%*%Z) %*%t(Z)%*%X)   %*%   t(X)%*%Z %*%solve(t(Z)%*%Z) %*%t(Z)%*%y;
    beta2SLS;
    betaIVreg.4inst=ivreg(y~x1+x2+x3 | x3+z1+z2+z3);
    summary(betaIVreg.4inst);
    dataSim=cbind(y,1,x1,x2,x3,z1,z2,z3);
    iv.moments=function(param,dataSim){
        y=dataSim[ ,1];
        x=dataSim[ ,2:5];
        z=dataSim[ ,c(2,5:8)];
        z*as.vector(y-x%*%param)
    };
    start.vals=c(rep(1,4));
    names(start.vals)=c("const","beta_1","beta_2","beta_3");
    gmm.model=gmm(iv.moments,dataSim,t0=start.vals);
    summary(gmm.model);
    cf_diff=coef(betaIVreg.4inst)-coef(betaOLSlm);
    vc_diff=vcov(betaIVreg.4inst)-vcov(betaOLSlm);
    H=t(cf_diff)%*%solve(vc_diff)%*%cf_diff;
    hausman_p=pchisq(H,df=rankMatrix(vc_diff),lower.tail=FALSE);
    sprintf("Testing Ho: difference in coefficients not systematic");
    sprintf("Hausman statistics %3.4f compared with Chi2(%1.0f)=9.488",H,rankMatrix(vc_diff));
    sprintf("p-value (Prob>chi2) %1.4f",hausman_p);
    lyn1<-sprintf("Testing Ho: difference in coefficients not systematic");
    lyn2<-sprintf("Hausman statistics %3.4f compared with Chi2(%1.0f)=9.488",H,rankMatrix(vc_diff));
    lyn3<-sprintf("p-value (Prob>chi2) %1.4f",hausman_p);
    lyn123<-as.data.frame(rbind(lyn1,lyn2,lyn3), stringsAsFactors = FALSE);
    write.xport(lyn123,file="d:/xpt/hausman.xpt");
    ');


    libname xpt xport "d:/xpt/hausman.xpt";

    proc contents data=xpt._all_;
    run;quit;

    data want;
     set
         xpt.lyn123;
    run;quit;

    proc print data=want;
    run;quit;

    libname xpt clear;

    *_
    | | ___   __ _
    | |/ _ \ / _` |
    | | (_) | (_| |
    |_|\___/ \__, |
             |___/
    ;

         Z1    Z2    Z3     U1      U2      X3     E1      E2
      <dbl> <dbl> <dbl>  <dbl>   <dbl>   <dbl>  <dbl>   <dbl>
    1  2.87  2.31  2.96 -1.06  -1.75    0.567  -1.46   0.260
    2  1.35  1.92  3.94 -1.83  -0.0747  0.454  -0.605 -0.298
    3  2.33  2.09  4.94 -0.555 -0.952   0.782  -0.634  0.0953
    4  2.20  1.49  4.26  0.575  0.654   2.58    0.143  0.396
    5  1.90  1.65  3.92 -1.31  -0.408  -0.0552  0.522  1.69
    6  1.98  2.12  3.79  0.487  1.06    1.06    0.303  0.347
                [,1]
    const  -1.434011
    beta_1  2.710927
    beta_2 -1.413894
    beta_3  1.063335
                 [,1]
    const   2.2899580
    beta_1  1.6437399
    beta_2 -0.8182181
    beta_3  1.0850504

    Call:
    lm(formula = y ~ x1 + x2 + x3)

    Residuals:
        Min      1Q  Median      3Q     Max
    -7.2105 -1.3666 -0.0398  1.4508  6.6195

    Coefficients:
                Estimate Std. Error t value Pr(>|t|)
    (Intercept) -1.43401    0.11017  -13.02   <2e-16 ***
    x1           2.71093    0.02781   97.47   <2e-16 ***
    x2          -1.41389    0.02055  -68.79   <2e-16 ***
    x3           1.06334    0.05135   20.71   <2e-16 ***
    ---
    Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    Residual standard error: 2.014 on 1496 degrees of freedom
    Multiple R-squared:  0.8761,	Adjusted R-squared:  0.8759
    F-statistic:  3526 on 3 and 1496 DF,  p-value: < 2.2e-16

                 [,1]
    const   1.1947593
    beta_1  1.9575101
    beta_2 -0.9629899
    beta_3  1.0782099

    Call:
    ivreg(formula = y ~ x1 + x2 + x3 | x3 + z1 + z2 + z3)

    Residuals:
         Min       1Q   Median       3Q      Max
    -8.21013 -1.78874 -0.02489  1.68339  8.62093

                Estimate Std. Error t value Pr(>|t|)
    (Intercept)  1.19476    0.22790   5.243 1.81e-07 ***
    x1           1.95751    0.06259  31.277  < 2e-16 ***
    x2          -0.96299    0.03893 -24.735  < 2e-16 ***
    x3           1.07821    0.06410  16.820  < 2e-16 ***
    ---
    Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    Residual standard error: 2.514 on 1496 degrees of freedom
    Multiple R-Squared: 0.807,	Adjusted R-squared: 0.8066
    Wald test: 428.8 on 3 and 1496 DF,  p-value: < 2.2e-16


    Call:
    gmm(g = iv.moments, x = dataSim, t0 = start.vals)


    Method:  twoStep

    Kernel:  Quadratic Spectral(with bw =  0.22204 )

    Coefficients:
            Estimate      Std. Error    t value       Pr(>|t|)
    const     1.1716e+00    2.1857e-01    5.3604e+00    8.3048e-08
    beta_1    1.9646e+00    6.0788e-02    3.2319e+01   3.8081e-229
    beta_2   -9.6555e-01    3.7482e-02   -2.5760e+01   2.4622e-146
    beta_3    1.0791e+00    5.9801e-02    1.8046e+01    8.5296e-73

    J-Test: degrees of freedom is 1
                    J-test   P-value
    Test E(g)=0:    1.34725  0.24576

    Initial values of the coefficients
         const     beta_1     beta_2     beta_3
     1.4787864  1.8826937 -0.9277454  1.0893160

    #############
    Information related to the numerical optimization
    Convergence code =  0
    Function eval. =  233
    Gradian eval. =  NA
    [1] "Testing Ho: difference in coefficients not systematic"
    [1] "Hausman statistics 198.9751 compared with Chi2(4)=9.488"
    [1] "p-value (Prob>chi2) 0.0000"
    >
