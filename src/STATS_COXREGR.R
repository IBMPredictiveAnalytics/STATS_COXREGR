#/***********************************************************************
# * Licensed Materials - Property of IBM 
# *
# * IBM SPSS Products: Statistics Common
# *
# * (C) Copyright IBM Corp. 2014
# *
# * US Government Users Restricted Rights - Use, duplication or disclosure
# * restricted by GSA ADP Schedule Contract with IBM Corp. 
# ************************************************************************/

# author__ = "SPSS, JKP"
# version__ = "1.0.1"

# History
# 15-jan-2014 Original Version


helptext="STATS COXREG CENSORTYPE = LEFT or RIGHT* or COUNTING
    EVENT=event variable name
    TIME1 = censor time variable
    TIME2 = second censor time variable for interval censoring
    INDEP = list of independent variables
    STRATA = list of stratification variables
    CLUSTER = clustering variable
/TIME X1=varlist X2=varlist ... X5=varlist
      T1='formula' T2='formula' ... T5='formula'
/OPTIONS 
    TIES=EFRON* or BRESLOW or EXACT
    MISSING = LISTWISE* or INCLUDE or FAIL
    MAXITER=number PCONVERGE=relative tolerance SINGULAR = number
/OUTPUT 
    SURVSPECS = YES or NO*
    PROPORTEST = YES or NO* TESTTRANSFORM=KM* or RANK or IDENTITY
    SURVIVALPLOT = YES or NO PLOTDATA=list of covariate values
/SAVE CASERESULTS=dataset
    RESIDUALS=MARTINGALE DEVIANCE SCORE DFBETA DFBETAS SCHOENFELD SCALEDSCH PARTIAL
    PREDICT=LP RISK EXPECTED
    PREDREF=STRATA or SAMPLE
* indicates default.
CENSORTYPE, EVENT, TIME1, and INDEP are required.

Example: 
stats coxregr censortype=right event=event time1=time
  indep = age quant
  /options ties=breslow
  /save caseresults=residuals residuals=dfbeta
  /output survspecs=yes.

This procedure does not honor weights or split files.

CENSORTYPE specifies the type of censoring.

EVENT specifies the event variable.  The values are interpreted as
  0 = no event (alive), 1 = event occurred (dead).  
  For interval censoring, the values are 0 = right censored,
  1 = event at time, 2 = left censored, 3 = interval censored.

TIME1.  For right-censored data, the values are the follow up (endpoint) time.
  For interval censoring, the values are the interval start times.  
  
TIME2. This is only used with counting process censoring, where it is required.

INDEP specifies the independent variables.  Variables with a
  categorical measurement level - nominal or ordinal, will be automatically 
  converted to factors.
  
STRATA specifies stratification variables.  The same variable
    cannot appear in both the independents and strata fields.

CLUSTER specifies a clustering variable for variance estimation.
This produces an approximate jacknife variance estimate if
the variable values are unique or, if it identifies clusters
of correlated observations, a grouped jacknife.

The TIME subcommand allows for time-varying variables.  Up to
five sets can be specified.  Each X parameter is paired with a T
parameter.  The X parameter lists one or more variables that share
a transformation formula defined by the T variable.  In writing
the formula, use x for the variable, t for time, enum for the
event, and weight, for weighted data, for the weight.  The
T terms must be quoted.  For example
X1 = X Z T1 = 'x + t/12'

Time varying data can also be represented as multiple cases
with one row for each time period using time1 and time2 variables.

MISSING indicates the missing value treatment.  Cases with user missing
    values can be excluded (LISTWISE), or
    the procedure can be stopped if any user missing are encountered (FAIL).
    System missing values always cause the case to be omitted.
    Missing values in the event or time variables are not allowed.
    
MAXITER, PCONVERGE, and SINGULAR specifications control the estimation
    process.  PCONVERGE is the relative tolerance (change) below which
    convergence is determined.

SURVSPECS specifies whether to produce a table of survival time
interpretations for up to the first 50 cases

PROPORTEST specifies whether or not to test the proportional hazards
assumption.  It produces a test table and plot for each variable.

TESTTRANSFORM specifies the time transformation to be performed
before the test is carried out.

SURVIVALPLOT specifies whether or not to produce this plot.
By default, the curve is evaluated at the means of the independent
variables.
PLOTDATA can specify a list of values, one per independent variable,
at which the survival curve is evaluated for the plot.

The SAVE subcommand specifies a dataset to be created with
selected types of residuals and predicted values.
CASERESULS is the dataset name, which must not already be in use.
RESIDUALS specifies one or more residual types to be saved.
The first four types produce a value for each case.  The last
two produce a value for each event.  You cannot save values from
both groups in the same command.  If event residuals are selected,
predicted values are not available.
PREDICT specifies one or more predicted variable types to be saved.
    LP is the linear predictor
    RISK is the reisk score, exp(LP)
    EXPECTED is the expected number of events
Standard errors are provided with each type specified.
PREDREF = STRATA specifies that predicted values are relative
to their stratum.  SAMPLE users the overall sample means.
    
STATS COXREGR /HELP.  prints this information and does nothing else.
"

# override for api to account for extra parameter in V19 and beyond
StartProcedure <- function(procname, omsid) {
    if (substr(spsspkg.GetSPSSVersion(),1, 2) >= 19) {
        spsspkg.StartProcedure(procname, omsid)
    }
    else {
        spsspkg.StartProcedure(omsid)
    }
}

gtxt <- function(...) {
    return(gettext(...,domain="STATS_COXREGR"))
}

gtxtf <- function(...) {
    return(gettextf(...,domain="STATS_COXREGR"))
}

dscheck = function(alldsspecs, warns) {
    # check dataset validation conditions
    
    if (!is.null(alldsspecs)) {
        alldatasets = spssdata.GetDataSetList()
        if ("*" %in% alldatasets) {
            warns$warn(gtxt("The active dataset must have a name if creating new datasets"),
                       dostop=TRUE)
        }
        if (length(intersect(alldsspecs, alldatasets) > 0)) {
            warns$warn(gtxt("One or more specified output dataset names are already in use"),
                       dostop=TRUE)
        }
    }
}

appendl = function(alist, avalue) {
    #return alist with avalue added as a new element at the end
    alist[[length(alist)+1]] = avalue
    return(alist)
}

### MAIN ROUTINE ###
docox = function(censortype, event, time1, time2=NULL, indep=NULL, 
    strata=NULL, cluster=NULL, 
    id=NULL, 
    x1=NULL,t1=NULL, x2=NULL,t2=NULL, x3=NULL,t3=NULL, x4=NULL,t4=NULL, x5=NULL,t5=NULL,
    ties="efron", survspecs=FALSE, proportest=FALSE, testtransform="km",
    caseresults=NULL, resids=NULL, pred=NULL, predref="strata", survivalplot=FALSE,
    plotdata=NULL,
    maxiter=30, pconverge=1e-09, singular=1e-10, missingaction="listwise"
    ) {
    #estimate Cox regression models
    
    setuplocalization("STATS_COXREGR")
    
    # A warnings proc name is associated with the regular output
    # (and the same omsid), because warnings/errors may appear in
    # a separate procedure block following the regular output
    procname=gtxt("Cox Regression")
    warningsprocname = gtxt("Cox Regression: Warnings")
    omsid="STATSCOXREGR"
    warns = Warn(procname=warningsprocname,omsid=omsid)
    strataeffect = "mult"
    allargs = as.list(environment())

    tryCatch(library(survival), error=function(e){
        warns$warn(gtxtf("The R %s package is required but could not be loaded.","survival"),dostop=TRUE)
        }
    )
    if (length(c(indep,x1,x2,x3,x4,x5)) == 0) {
        warns$warn(gtxt("At least one independent variable must be specified"), dostop=TRUE)
    }
    if (censortype == "interval" && is.null(time2)) {
        warns$warn(gtxt("Two time variables are required for interval censoring"), dostop=TRUE)
    }
    if (!is.null(caseresults)) {
        if (is.null(c(resids, pred))) {
            warns$warn(gtxt("An output dataset was specified without any residual or prediction specifications"), 
                dostop=TRUE)
        }
    } else {
        if (!is.null(c(resids, pred))) {
            warns$warn(gtxt("Residuals or predicted values were specified but no dataset name was given"),
                dostop=TRUE)
        }
    }

    # ensure that output dataset names are not in use
    alldsspecs = c(caseresults)
    dscheck(alldsspecs, warns)
    if (length(union(indep, strata)) != length(indep) + length(strata)) {
        warns$warn(gtxt("The same variable cannot appear in both the independents and strata lists"),
            dostop=TRUE)
    }
    alldata = c(event, time1, time2, indep, x1,x2,x3,x4,x5, strata, cluster)
    # if no id variable was specified, the row.label argument will be NULL and is ignored
    dta = spssdata.GetDataFromSPSS(alldata, row.label=id, missingValueToNA=TRUE, factorMode="levels")
    if (any(is.factor(dta[[time1]]), (!is.null(time2) && is.factor(dta[[time2]])), is.factor(dta[[event]]))) {
        warns$warn(gtxt("The event and time variables must have a scale measurement level"), dostop=TRUE)
    }
    for (f in strata) {
        if (!is.factor(dta[[f]])) {
            warns$warn(gtxtf("Strata variables must be categorical: %s", f), dostop=TRUE)
        }
    }
    naaction = list(listwise=na.omit, include=na.pass, fail=na.fail)[[missingaction]]
    if (is.null(time2)) {
        sur = tryCatch(Surv(dta[[time1]], dta[[event]], type=censortype),
            error=function(e) {warns$warn(e$message, dostop=TRUE)
            }
        )
    } else {
        sur = tryCatch(Surv(dta[[time1]], dta[[time2]], event=dta[[event]], type=censortype),
            error=function(e) {warns$warn(e$message, dostop=TRUE)
            }
        )
    }
    frml = buildfrml(indep, strata, cluster, strataeffect, x1,x2,x3,x4,x5,t1,t2,t3,t4,t5)
    ctrl = coxph.control(iter.max=maxiter, eps=pconverge, toler.chol=singular)
    # The tt argument is an empty list if no time functions are specified
    if (!(is.null(pred) || is.null(resid))) {
        x = TRUE   # save X matrix in result
    } else {
        x = FALSE
    }
    # The anova function will require these variables, but when it calls coxph, it neglects
    # to pass them (and won't accept them as extra parameters)
    assign("ties", ties, envir=.GlobalEnv)
    assign("naaction", naaction, envir=.GlobalEnv)
    assign("dta", dta, envir=.GlobalEnv)
    assign("ctrl", ctrl, envir=.GlobalEnv)
    assign("x", x, envir=.GlobalEnv)
    res = tryCatch(coxph(as.formula(frml[[1]]), data=dta, x=x, na.action=naaction, ties=ties,
        control=ctrl, tt=frml[[2]]),
        error = function(e) {
            warns$warn(e$message, dostop=FALSE)
            return(NULL)
        }
    )
    
    # if failed and survival interpretation requested, display that much and stop
    if (is.null(res)) {
        if (survspecs) {
            StartProcedure(gtxt(procname), omsid)
            survtable(sur, dta, event, time1, time2)
        }
        warns$warn(dostop=TRUE)
    }
            
    ressum = summary(res)
    displayresults(allargs, res, ressum, dta, sur, warns)
    
    if (!is.null(caseresults)) {
        docaseresults(caseresults, res, resids, pred, predref, dta, id, warns)
    }
    warns$warn()  # ensure warnings are displayed
    rm(dta, ties, ctrl, x, naaction, envir=.GlobalEnv)
}
survtable = function(sur, dta, event, time1, time2) {
    # Display partial survival interpretation
      # The Surv object actually has two columns and gets special formatting
      # when displayed, but trickery is required for the pivot table formatting

      limit = min(nrow(sur), 50)
      surdf = data.frame(sur[1:limit,])
      attr(surdf, "names") = gtxt("Survival")
      surdf = sapply(surdf, as.character)
      for (v in c(event, time1, time2)) {
          surdf = data.frame(surdf, dta[1:limit, v])
          names(surdf)[length(surdf)] = v
      }
      spsspivottable.Display(surdf, title=gtxtf("Survival Values (First %s Cases)", limit),
        templateName="SURVVALUES",
        caption=gtxt("+ indicates that the event did not occur"),
          isSplit=FALSE
      )
}        

buildfrml = function(indep, strata, cluster, strataeffect,x1,x2,x3,x4,x5,t1,t2,t3,t4,t5, warns) {
    # Return list whose first element is the formula as a string and other elements
    # are strings for ...
    # followed by a list of functions defined from t1,...t5, which may be empty
    # Dependent is always named sur
    # indep is the list of independents
    # cluster and strata are null or names of clustering and strata variables
    # strataeffect is "mult" or "add" indicating multiplicative or additive effect
    # x1,...x5 are lists of independents or NULL for time functions
    # t1,...,t5 are literals defining corresponding time functions
    # warns is the error message object
    
    rhslist = list()
    frmllist = list()
    ttlist = list()
    allx = list(x1,x2,x3,x4,x5)
    allt = list(t1,t2,t3,t4,t5)
    # add time-dependent covariates
    for (i in 1:5) {
        if (!is.null(allx[[i]])) {
            if (is.null(allt[[i]])) {
                warns$warn(gtxt("A list of time variables does not have a corresponding function"),dostop=TRUE)
            }
            for (v in allx[[i]]) {
                rhslist = appendl(rhslist, paste("tt(", v, ")",sep=""))
                ttlist = appendl(ttlist, allt[[i]])
            }
        }
    }
    # add the regular variables
    for (v in indep) {
        rhslist = appendl(rhslist, v)
    }
    if (!is.null(strata)) {
        strata = paste(strata,collapse="+")
        if (strataeffect == "add") {
            rhslist = appendl(rhslist, paste("strata(", strata, ")", sep=""))
        } else {
            rhslist[1] = paste("(", rhslist[[1]], sep="")
            rhslist[length(rhslist)] = paste(rhslist[[length(rhslist)]], ")*strata(", strata, ")", sep="")
        }
    }
    if (!is.null(cluster)) {
        rhslist = appendl(rhslist, paste("cluster(", cluster, ")", sep=""))
    }
    frml = paste("sur~", paste(rhslist, collapse="+"))
    for (tt in seq(along.with=ttlist)) {
        # turn the expression into a function that coxph can call
        ttlist[[tt]] = eval(bquote(function(x, t, enum=NULL, weight=NULL) .(parse(text = ttlist[tt])[[1]])))
    }
    return(list(frml, ttlist))
}

displayresults = function(allargs, res, ressum, dta, sur, warns) {
    StartProcedure(gtxt("Cox Regression"), "STATSCOXREG")
    summarylabels=list(
        gtxt("Event Variable"),
        gtxt("Censoring Type"),
        gtxt("Time"),
        gtxt("Time2"),
        gtxt("Strata"),
        gtxt("Cluster"),
        gtxt("Ties"),
        gtxt("Output Dataset"),
        gtxt("Prediction Reference"),
        gtxt("Missing Values Treatment"),
        gtxt("Number of Cases"),
        gtxt("Missing Value Deletion"),
        gtxt("Number of Events"),
        gtxt("Number of Iterations"),
        gtxt("Log Likelihood (Initial)"),
        gtxt("Log Likelihood (Final)"),
        gtxt("Likelihood Ratio Test"),
        gtxt("LR Degrees of Freedom"),
        gtxt("LR Sig"),
        gtxt("Wald Test"),
        gtxt("Wald Degrees of Freedom"),
        gtxt("Wald Sig"),
        gtxt("Score (Logrank) Test"),
        gtxt("Score Degrees of Freedom"),
        gtxt("Score Sig"),
        gtxt("R Squared"),
        gtxt("Maximum Possible Rsq"),
        gtxt("Concordance"),
        gtxt("Concordance Std. Error")
    )

    loglik = -2 * (res$loglik[1] - res$loglik[2])
    df = ifelse(is.null(res$df), sum(!is.na(res$coefficients)), sum(res$df))
    rsq = 1-exp(-loglik/res$n)
    maxrsq = 1-exp(res$loglik[1]/res$n)
                   
   if (!is.null(res$concordance)) {
       if (is.matrix(res$concordance)) {
           ctemp = colSums(res$concordance)
       }
       else {
           ctemp = res$concordance
       }
       concord = round((ctemp[1] + ctemp[3]/2) / sum(ctemp[1:3]),4)
       concordse = round(ctemp[5] / (2*sum(ctemp[1:3])), 4)
    } else {
       concord = gtxt("--NA--")
       concordse = gtxt("--NA--")
    }
    if (!is.null(res$wald.test)) {
        wald = round(res$wald.test, 4)
        wdf = df
        wsig = round(1 - pchisq(wald, wdf), 4)
    } else {
        wald = gtxt("--NA--")
        wdf = gtxt("--NA--")
        wsig = gtxt("--NA--")
    }
   if (!is.null(res$score)) {
       score = round(res$score, 4)
       sdf = df
       ssig = round(1 - pchisq(score, sdf), 4)
   } else {
       score = gtxt("--NA--")
       sdf = gtxt("--NA--")
       ssig = gtxt("--NA--")
   }  
    if (length(res$na.action) > 0) {
        # get the count of omitted
        nomitted = length(res$na.action)
    } else {
        nomitted = 0
    }
    
    # NOTE: some terms below are from syntax and would need to mapped for
    # translation purposes, e.g., the censortype and ties parameters
    summaryvalues = list(
        allargs[["event"]],
        allargs[["censortype"]],
        allargs[["time1"]],
        ifelse(is.null(allargs[["time2"]]), gtxt("--None--"), allargs[["time2"]]),
        ifelse(is.null(allargs[["strata"]]), gtxt("--None--"), paste(allargs[["strata"]], collapse=" ")),
        ifelse(is.null(allargs[["cluster"]]), gtxt("--None--"), allargs[["cluster"]]),
        res$method,    # ties
        ifelse(is.null(allargs[["caseresults"]]), gtxt("--None--"), allargs[["caseresults"]]),
        ifelse(is.null(allargs[["pred"]]), gtxt("--NA--"), allargs[["predref"]]),
        allargs[["missingaction"]],
        res$n,
        nomitted,
        res$nevent,
        res$iter,
        round(res$loglik[[1]],4), round(res$loglik[[2]],4),
        round(loglik,4), round(df, 4), round(1 - pchisq(loglik, df),4),    # LRT
        wald, wdf, wsig,    #Wald
        score, sdf, ssig,  #score
        round(rsq,4),
        round(maxrsq,4),
        concord,
        concordse
    )

    names(summaryvalues) = summarylabels
    summarydf = data.frame(cbind(summaryvalues))
    colnames(summarydf) = gtxt("Values")
    if ("robust se" %in% dimnames(ressum$coefficients)[[2]]) {
        rolabel = gtxt("Robust\nS.E.")
    } else {
        rolabel = NULL
    }
    spsspivottable.Display(summarydf, title=gtxt("Cox Regression Summary"), 
                           templateName="COXREGSUMMARY",
                           caption=gtxt("Results computed by R survival package"),
                           isSplit=FALSE,
                           format=formatSpec.Count
    )

    # main results - differ for time functions
    if (!is.null(ressum)) {
        resdf = data.frame(ressum$coefficients, ressum$conf.int)
        resdf = resdf[-2]
        names(resdf) = c(gtxt("Coefficient"), gtxt("S.E.\nCoefficient"), 
            rolabel, gtxt("Z"), gtxt("Sig."),
            gtxt("Exp(coef)"), gtxt("Exp(-coef)"), gtxt("Lower\n95% C.I"), gtxt("Upper\n95% C.I"))
        caption = gtxtf("Event variable: %s", allargs[["event"]])
        spsspivottable.Display(resdf, title=gtxt("Cox Regression"), templateName="COXREGRESS",
                               caption=caption, isSplit=FALSE
        )
    } else {
        print("Unable to create pivot table")  #time covariates with pspline
        print(summary(res))
    }

    anytimes = dotimes(list(allargs[['x1']], allargs[['x2']], allargs[['x3']], allargs[['x4']], allargs[['x5']]),
        list(allargs[['t1']], allargs[['t2']], allargs[['t3']], allargs[['t4']], allargs[['t5']])
    )

    ares = tryCatch(anova(res), error=function(e) {return(NULL)}
    )
    if (is.null(ares)) {
        warns$warn(gtxt("The ANOVA table is not available when clustering is used"), dostop=FALSE)
    } else {
        df = data.frame(ares)
        names(df) = c(gtxt("Log Likelihood"), gtxt("Chi-Squared"), gtxt("Degrees of\nFreedom"), gtxt("Sig."))
        spsspivottable.Display(
            df, 
            title = gtxt("Analysis of Deviance"),
            templateName="COXANOVA",
            caption=gtxt("Terms added sequentially")
        )   
    }
    
    if (allargs[["survivalplot"]]) {
        dosurvivalplot(allargs, res, dta, anytimes, warns)
    }
    
    doproptest(allargs, res)

    if (allargs[["survspecs"]]) {
        survtable(sur, dta, allargs[["event"]], allargs[["time1"]], allargs[["time2"]])
    }
    spsspkg.EndProcedure()
}

dotimes = function (vars, funcs) {
    # display pivot table of time transformations, if any
    # vars is the list of variable lists.  Each element is a list
    # funcs is the list of function bodies
    # returns TRUE if there are any time variables else FALSE
    
    vlist = list()
    flist = list()

    for (i in 1:length(vars)) {
        if (!is.null(vars[[i]])) {
            for (v in vars[[i]]) {
                vlist = appendl(vlist, v)
                flist = appendl(flist, funcs[[i]])
            }
        }
    }
    if (length(vlist) > 0) {
        df = data.frame(cbind(flist), row.names=vlist)
        spsspivottable.Display(df, 
            title=gtxt("Time Functions"), 
            templateName="COXREGTIMEFUNCS",
            rowdim=gtxt("Variable"),
            collabels=gtxt("Function"),
            isSplit=FALSE,
            caption=gtxt("Function definitions are interpreted as R syntax"))
    }
    return(length(vlist) > 0)
}

dosurvivalplot <- function (allargs, res, dta, anytimes, warns) {
    # display a survival plot
    # allargs is the command argument list
    # res is the coxph result
    # dta is the data
    # anytimes is TRUE if there are time-varying variables
    
    if (anytimes) {
        warns$warn(gtxt("Survival plot is not available when there are time-varying functions."), dostop=FALSE)
    } else {
    if (!is.null(allargs[["plotdata"]])) {
        # convert indep var values to appropriate type
        # This code assumes that categorical variables have been converted to R factors.
        # The value types for the survival plot are very sensitive to the type
        # The survfit function will convert factor values appropriately but
        # expects strings, not numeric values.
        ddd = spssdictionary.GetDictionaryFromSPSS(allargs[["indep"]])
        tryCatch({
            values = list()
            for (i in 1:length(allargs[["plotdata"]])) {
                if (ddd[["varMeasurementLevel", i]] %in% list("scale", "unknown")) {
                    values[i] = as.numeric(allargs[["plotdata"]][[i]])
                } else {
                    values[i] = allargs[["plotdata"]][[i]]
                }
            }
        }, error=function(e) warns$warn(gtxt("Independent variable value for survival plot has wrong type or too many values specified"),
                                        dostop=TRUE)
        )
        pdataframe = data.frame(values)
        if (ncol(pdataframe) != length(allargs[["indep"]])) {
            warns$warn(gtxt("Number of variable values for survival plot differs from variables in the model"),
                dostop=TRUE)
        }
        names(pdataframe) = allargs[["indep"]]
        sfobj = tryCatch(survfit(res, newdata=pdataframe, censor = anytimes),
                    error=function(e) {warns$warn(e$message, dostop=TRUE)}
        )
        title = sprintf(gtxt("Survival function at\n%s"), 
            paste(allargs[["indep"]], allargs[["plotdata"]], sep="=", collapse=" "))
    } else {
        sfobj = tryCatch(survfit(res, censor=anytimes),
            error=function(e) {warns$warn(e$message, dostop=TRUE)}
        )
        title = gtxt("Survival Function at Covariate Means")
    }
    if (!is.null(allargs[["strata"]])) {
        subtitle = paste(gtxt("Strata:"), allargs[["strata"]])
    } else {
        subtitle = gtxt("Broken lines show 95% confidence interval")
    }
    plot(sfobj,
         main=title, 
         sub=subtitle,
         xlab=allargs[["time1"]], 
         ylab=gtxt("Cumulative Survival")
    )
    }
}

doproptest <- function (allargs, res) {
    # Do Cox proportionality assumption test, including plot
    # allargs is all the input arguments
    # res is the Cox result
    
    if (allargs[["proportest"]]) {
        zph = cox.zph(res, allargs[["testtransform"]])
        zlen = length(dimnames(zph$table)[[1]])
        if (zlen > 1) {
            dimnames(zph$table)[[1]][zlen] = gtxt("<GLOBAL>")
        }
        ttf = data.frame(zph$table)
        names(ttf) = c(gtxt("Rho"), gtxt("Chi-Squared"), gtxt("Sig."))
        spsspivottable.Display(ttf, 
            title=gtxt("Proportional Hazards Test"),
            templateName="COXREGTEST",
            caption=paste(gtxt("Survival times transform"), allargs[["testtransform"]], sep=": ")
        )
        plot(zph, main=gtxt("Proportional Hazards Test: Schoenfeld Residuals"))   
    }
}
# variable labels for residuals and predicted values
labels = list(res_martingale=gtxt("Residuals - Martingale"),
              res_deviance = gtxt("Residuals - Deviance"),
              res_score = gtxt("Residuals - Score"),
              res_schoenfeld = gtxt("Residuals - Schoenfeld"),
              res_ldresp = gtxt("Residuals - Ldresponse"),
              res_scaledsch = gtxt("Residuals - Scaled Schoenfeld"),
              res_dfbeta = gtxt("Dfbeta"),
              res_dfbetas = gtxt("Stdized Dfbeta"),
              pred_risk = gtxt("Predicted - Risk Score exp(linear pred.)"),
              pred_lp = gtxt("Predicted - Linear Predictor"),
              pred_expected = gtxt("Predicted - Expected Number of Events")
)

docaseresults = function(caseresults, res, resids, pred, predref, dta, idname, warns) {
    # Produce dataset of selected residuals and predicted values
    # caseresults is the dataset name
    # res is the coxph results
    # resids is the list of residuals specifications
    # pred is the predicted specification
    # predref is the prediction centering reference
    # dta is the input data (extract row.names(dta) for complete cases)
    # idname is the name of the id variable
    # warns is the error message object

    resbyevent = list("schoenfeld", "scaledsch")
    resbycase = list("martingale", "score", "dfbeta", "dfbetas", "score")
    if (any(match(resids, resbyevent), na.rm=TRUE) && any(match(resids, resbycase), na.rm=TRUE)) {
        warns$warn(gtxt("Casewise and eventwise residuals cannot be combined in the same dataset"),
            dostop=TRUE)
    }
    if (any(match(resids, resbyevent), na.rm=TRUE) && length(pred) > 0) {
        warns$warn(gtxt("Eventwise residuals cannot be combined with predicted values in the same dataset"),
            dostop=TRUE)
    }
    dict = list()
    if (any(match(resids, resbyevent), na.rm=TRUE)) {
        reqlength = res$nevent
        idname = NULL
    } else {
        reqlength = nrow(res$y)
    }
    if (is.null(idname)) {
        df = data.frame(row.names=1:reqlength)  # empty data frame with the right number of rows
        dict["ID"] = list(c("ID", "", 0, "F8.0", "nominal"))
    } else {
        df = data.frame(row.names=row.names(dta[complete.cases(dta),]))
        #df = data.frame(row.names=id)  # empty data frame with the right number of rows
        dict[["ID"]] = spssdictionary.GetDictionaryFromSPSS(idname)
        dict[["ID"]][1,] = "ID"     # ensure that we can't have a duplicate name
    }
    resultcount = 0

    # first build residual results; then prediction results
    for (rtype in resids) {
        e = tryCatch(residuals(res, type=rtype), # some types just fail with this function for some models :-)
                error = function(ea) {
                    warns$warn(gtxtf("Residual statistic %s not available for this model", rtype), dostop=FALSE)
                    return(NULL)
                }
            )
        # add successfull results to the data frame

        if (!is.null(e)) {
            if (!(rtype %in% list("dfbeta", "dfbetas", "score", "schoenfeld", "scaledsch"))) {  # dfbeta/s are multicolumn
                resultcount = resultcount + 1
                df = data.frame(df, e)
                itemname = paste("res", rtype, sep ="_")
                dict[itemname] = list(c(itemname, labels[[itemname]], 0, "F8.2", "scale"))
            } else {
                for (v in 1:ncol(e)) {
                    resultcount = resultcount + 1
                    df = data.frame(df, e[,v])
                    itemname = paste(rtype, v, sep="_")  # improve this
                    if (v < ncol(e)) {
                        itemlabel = names(res$coefficients)[[v]]
                    } else {
                        itemlabel = gtxt("Other")
                    } 
                    dict[itemname] = list(c(itemname, itemlabel, 0, "F8.2", "scale"))
                }
            }
        }
    }
    
    # predictions
    for (rtype in pred) {
        if (!(rtype %in% list("quantile", "uquantile"))) {   # always true
            e = tryCatch(predict(res, type=rtype, se.fit=TRUE, reference=predref), 
                    error = function(ea) {
                        warns$warn(gtxtf("Prediction statistic %s not available for this model", rtype), dostop=FALSE)
                        return(NULL)
                    }
                )
            if (!is.null(e)){
                resultcount = resultcount + 1
                df = data.frame(df, e)
                itemname = paste("pred", rtype, sep="_")
                dict[itemname] = list(c(itemname, labels[[itemname]], 0, "F8.2", "scale"))
                resultcount = resultcount + 1   # std error
                itemname = paste(itemname, "SE", sep="_")
                dict[itemname] = list(c(itemname, itemname, 0, "F8.2", "scale"))
            }
        } 
#         else {  # quantile and uquantile.  se.fit does not work for these (undocumented)
#             e = tryCatch(predict(res, type = rtype, se.fit=FALSE, reference=predref),
#                     error = function(ea) {
#                         warns$warn(gtxtf("Prediction statistic %s not available for this model", rtype), 
#                             dostop=FALSE)
#                         return(NULL)
#                     }
#                 )
#             if (!is.null(e)) {
#                 df = data.frame(df, e)   # as many columns as quantile values
#                 for (q in quants) {
#                     resultcount = resultcount + 1
#                     itemname = paste("pred", 
#                         ifelse(rtype == "quantile", gtxt("quantile"), gtxt("uquantile")),
#                         q,
#                         sep = "_"
#                         )
#                     dict[itemname] = list(c(itemname, itemname, 0, "F8.2", "scale"))
#                 }
#             }
#         }
    }
    
    # finally! build the dataset
    if (resultcount == 0) {
        warns$warn(gtxt("No casewise results could be computed"), dostop=TRUE)
    }

    tryCatch({
        spssdictionary.SetDictionaryToSPSS(caseresults, 
            spssdictionary.CreateSPSSDictionary(dict))
        spssdata.SetDataToSPSS(caseresults, data.frame(row.names(df), df))
        },
        error=function(e) {warns$warn(e)
          warns$warn(gtxt("Failed to create caseresults dataset."), dostop=FALSE)}
    )
    spssdictionary.EndDataStep()
}

Warn = function(procname, omsid) {
    # constructor (sort of) for message management
    lcl = list(
        procname=procname,
        omsid=omsid,
        msglist = list(),  # accumulate messages
        msgnum = 0
    )
    # This line is the key to this approach
    lcl = list2env(lcl) # makes this list into an environment

    lcl$warn = function(msg=NULL, dostop=FALSE, inproc=FALSE) {
        # Accumulate messages and, if dostop or no message, display all
        # messages and end procedure state
        # If dostop, issue a stop.

        if (!is.null(msg)) { # accumulate message
            assign("msgnum", lcl$msgnum + 1, envir=lcl)
            # There seems to be no way to update an object, only replace it
            m = lcl$msglist
            m[[lcl$msgnum]] = msg
            assign("msglist", m, envir=lcl)
        } 

        if (is.null(msg) || dostop) {
            lcl$display(inproc)  # display messages and end procedure state
            if (dostop) {
                stop(gtxt("End of procedure"), call.=FALSE)  # may result in dangling error text
            }
        }
    }
    
    lcl$display = function(inproc=FALSE) {
        # display any accumulated messages as a warnings table or as prints
        # and end procedure state, if any

        if (lcl$msgnum == 0) {   # nothing to display
            if (inproc) {
                spss.EndProcedure()
            }
        } else {
            if (!inproc) {
                procok =tryCatch({
                    StartProcedure(lcl$procname, lcl$omsid)
                    TRUE
                    },
                    error = function(e) {
                        FALSE
                    }
                )
            }
            if (procok) {  # build and display a Warnings table if we can
                table = spss.BasePivotTable("Warnings ","Warnings") # do not translate this
                rowdim = BasePivotTable.Append(table,Dimension.Place.row, 
                    gtxt("Message Number"), hideName = FALSE,hideLabels = FALSE)

                for (i in 1:lcl$msgnum) {
                    rowcategory = spss.CellText.String(as.character(i))
                    BasePivotTable.SetCategories(table,rowdim,rowcategory)
                    BasePivotTable.SetCellValue(table,rowcategory, 
                        spss.CellText.String(lcl$msglist[[i]]))
                }
                spsspkg.EndProcedure()   # implies display
            } else { # can't produce a table
                for (i in 1:lcl$msgnum) {
                    print(lcl$msglist[[i]])
                }
            }
        }
    }
    return(lcl)
}

# localization initialization
setuplocalization = function(domain) {
    # find and bind translation file names
    # domain is the root name of the extension command .R file, e.g., "SPSSINC_BREUSCH_PAGAN"
    # This would be bound to root location/SPSSINC_BREUSCH_PAGAN/lang

    fpath = Find(file.exists, file.path(.libPaths(), paste(domain, ".R", sep="")))
    bindtextdomain(domain, file.path(dirname(fpath), domain, "lang"))
} 


Run = function(args) {
    #Execute the STATS COXREGR command

    cmdname = args[[1]]
    args = args[[2]]
    oobj = spsspkg.Syntax(list(
        spsspkg.Template("CENSORTYPE", subc="", ktype="str", var="censortype",
            vallist=list("right", "left", "interval", "counting")),
        spsspkg.Template("EVENT", subc="",  ktype="existingvarlist", var="event",
            islist=FALSE),
        spsspkg.Template("TIME1", subc="", ktype="existingvarlist", var="time1",
            islist=FALSE),
        spsspkg.Template("TIME2", subc="",  ktype="existingvarlist", var="time2"),
        spsspkg.Template("INDEP", subc="", ktype="existingvarlist", var="indep", islist=TRUE),
        spsspkg.Template("STRATA", subc="", ktype="existingvarlist", var="strata", islist=TRUE),
        spsspkg.Template("ID", subc="", ktype="existingvarlist", var="id"),
        spsspkg.Template("CLUSTER", subc="", ktype="existingvarlist", var="cluster", islist=FALSE),
        
        spsspkg.Template("X1", subc="TIME", ktype="existingvarlist", var="x1", islist=TRUE),
        spsspkg.Template("T1", subc="TIME", ktype="literal", var="t1", islist=FALSE),
        spsspkg.Template("X2", subc="TIME", ktype="existingvarlist", var="x2", islist=TRUE),
        spsspkg.Template("T2", subc="TIME", ktype="literal", var="t2", islist=FALSE),
        spsspkg.Template("X3", subc="TIME", ktype="existingvarlist", var="x3", islist=TRUE),
        spsspkg.Template("T3", subc="TIME", ktype="literal", var="t3", islist=FALSE),
        spsspkg.Template("X4", subc="TIME", ktype="existingvarlist", var="x4", islist=TRUE),
        spsspkg.Template("T4", subc="TIME", ktype="literal", var="t4", islist=FALSE),
        spsspkg.Template("X5", subc="TIME", ktype="existingvarlist", var="x5", islist=TRUE),
        spsspkg.Template("T5", subc="TIME", ktype="literal", var="t5", islist=FALSE),
        
        spsspkg.Template("MAXITER", subc="OPTIONS", ktype="int", var="maxiter",
            vallist=list(1)),
        spsspkg.Template("PCONVERGE", subc="OPTIONS", ktype="float", var="pconverge",
            vallist=list(0)),
        spsspkg.Template("TIES", subc="OPTIONS", ktype="str", var="ties",
            vallist=list("efron", "breslow", "exact")),
        spsspkg.Template("MISSING", subc="OPTIONS", ktype="str", var="missingaction",
            vallist=list("listwise", "fail")),
        
        spsspkg.Template("SURVSPECS", subc="OUTPUT", ktype="bool", var="survspecs"),
        spsspkg.Template("PROPORTEST", subc="OUTPUT", ktype="bool", var="proportest"),
        spsspkg.Template("TESTTRANSFORM", subc="OUTPUT", ktype="str", var="testtransform",
            vallist=list("km", "rank", "identity")),
        spsspkg.Template("SURVIVALPLOT", subc="OUTPUT", ktype="bool", var="survivalplot"),
        spsspkg.Template("PLOTDATA", subc="OUTPUT", ktype="literal", var="plotdata", islist=TRUE), #was float
        
        spsspkg.Template("CASERESULTS", subc="SAVE", ktype="varname", var="caseresults"),
        #,"ldcase","ldresp","shape" sometimes do not work in survreg
        spsspkg.Template("RESIDUALS", subc="SAVE", ktype="str", var="resids", islist=TRUE,
            vallist = list("martingale","deviance","score", "dfbeta","dfbetas","schoenfeld","scaledsch")),  
        spsspkg.Template("PREDICT", subc="SAVE", ktype="str", var="pred", islist=TRUE,
            vallist=list("lp","risk","expected")),
        spsspkg.Template("PREDREF", subc="SAVE",ktype="str", var="predref", islist=TRUE,
            vallist=list("strata", "sample"))
    ))

    # A HELP subcommand overrides all else
    if ("HELP" %in% attr(args,"names")) {
        #writeLines(helptext)
        helper(cmdname)
    }
    else {
        res <- spsspkg.processcmd(oobj, args, "docox")
    }
}

helper = function(cmdname) {
    # find the html help file and display in the default browser
    # cmdname may have blanks that need to be converted to _ to match the file
    
    fn = gsub(" ", "_", cmdname, fixed=TRUE)
    thefile = Find(file.exists, file.path(.libPaths(), fn, "markdown.html"))
    if (is.null(thefile)) {
        print("Help file not found")
    } else {
        browseURL(paste("file://", thefile, sep=""))
    }
}
if (exists("spsspkg.helper")) {
assign("helper", spsspkg.helper)
}
