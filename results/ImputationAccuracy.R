library(dplyr)

rd <- "/project/flatiron_ucc/programs/kylie/RunMe3"
cdir = file.path (rd, 'datasets_simple/cDats')
rootdir <- file.path (rd, "final_simple")
opdir = file.path (rootdir, "density_plots")
iter = 500

dens <- function (i, mdm, pc){
    cDat = read.csv (file.path (cdir, mdm, pc, paste0 ('cData', i, '.csv')))
    mice = read.csv (file.path (rootdir, mdm, pc, paste0 ('mice_', i, '.csv')))
    rf = read.csv (file.path (rootdir, mdm, pc, paste0 ('rf_', i, '.csv')))
    daeList = lapply (1:10, function (x){
        dat = read.csv (file.path (rootdir, paste0 ('AE', mdm), pc, paste0 ('result', i, '_num_', x, '.csv')))
        .imp = rep (x, nrow (dat))
        newdat = cbind (.imp, dat)
    })
    dae = do.call ('rbind', daeList)
    missEcog = mice$b.ecogvalue[which (mice$.imp==0)]
    missInd = which (is.na (missEcog))
    miceImputed = lapply (1:10, function (x){
        dat = filter (mice, .imp==x)
        dat = select (dat[missInd,], treat, b.ecogvalue)
        dat$variable = 'MICE'
        return (dat)
    })
    rfImputed = lapply (1:10, function (x){
        dat = filter (rf, .imp==x)
        dat = select (dat[missInd,], treat, b.ecogvalue)
        dat$variable = 'RF'
        return (dat)
    })
    daeImputed = lapply (1:10, function (x){
        dat = filter (dae, .imp==x)
        dat = select (dat[missInd,], treat, b.ecogvalue)
        dat$variable = 'DAE'
        return (dat)
    })
    miceFinal = do.call ('rbind', miceImputed)
    rfFinal = do.call ('rbind', rfImputed)
    daeFinal = do.call ('rbind', daeImputed)
    true = data.frame ('treat'=cDat$treat[missInd], 'b.ecogvalue'=cDat$b.ecogvalue[missInd], 'variable'='True')
    final = rbind (true, miceFinal, rfFinal, daeFinal)
    write.csv(final, file=file.path (opdir, paste0 (mdm, '_', pc, '_', 'iter', i, '.csv')), row.names=FALSE)
}

for (j in 21:50){
    print (j)
    dens(j, 'MCAR', 50)
    dens(j, 'MAR', 50)
    dens(j, 'MNAR1', 50)
}

getDiff = function (i, mdm, pc){
    cDat = read.csv (file.path (cdir, mdm, pc, paste0 ('cData', i, '.csv')))
    mice = read.csv (file.path (rootdir, mdm, pc, paste0 ('mice_', i, '.csv')))
    rf = read.csv (file.path (rootdir, mdm, pc, paste0 ('rf_', i, '.csv')))
    daeList = lapply (1:10, function (x){
        dat = read.csv (file.path (rootdir, paste0 ('AE', mdm), pc, paste0 ('result', i, '_num_', x, '.csv')))
        dat = cbind (x, dat)
        colnames (dat)[1] = '.imp'
    })
    dae = do.call ('rbind', daeList)

    cDat = select (cDat, genderf, reth_black, reth_hisp, reth_oth, age, practypec, site_ureter, site_renal, site_urethra, smokey, dgradeh, surgery, b.ecogvalue, var1, var2)

    mice = select (mice, .imp, genderf, reth_black, reth_hisp, reth_oth, age, practypec, site_ureter, site_renal, site_urethra, smokey, dgradeh, surgery, b.ecogvalue, var1, var2)

    diffs = lapply (1:10, function (x){
        imputed = filter (mice, .imp==x)
        imputed = select (imputed, -.imp)
        diff = as.matrix (imputed) - as.matrix (cDat)
        return (diff)
    })

    mde = lapply (diffs, function (x){
        apply (x, 2, mean)
    })
    mae = lapply (diffs, function (x){
        err = abs (x)
        apply (err, 2, mean)
    })
    mse = lapply (diffs, function (x){
        err = x^2
        apply (err, 2, mean)
    })
    avgAvgDiff = Reduce ('+', mde)/length (mde)
    avgMAE = Reduce ('+', mae)/length (mae)
    avgMSE = Reduce ('+', mse)/length (mse)

    opMICE = cbind (avgAvgDiff, avgMAE, avgMSE)
    colnames (opMICE) = c ("miceAvgMD", "miceAvgMAE", "miceAvgMSE")

    rf = select (rf, .imp, genderf, reth_black, reth_hisp, reth_oth, age, practypec, site_ureter, site_renal, site_urethra, smokey, dgradeh, surgery, b.ecogvalue, var1, var2)

    diffs = lapply (1:10, function (x){
        imputed = filter (rf, .imp==x)
        imputed = select (imputed, -.imp)
        diff = as.matrix (imputed) - as.matrix (cDat)
        return (diff)
    })

    mde = lapply (diffs, function (x){
        apply (x, 2, mean)
    })
    mae = lapply (diffs, function (x){
        err = abs (x)
        apply (err, 2, mean)
    })
    mse = lapply (diffs, function (x){
        err = x^2
        apply (err, 2, mean)
    })
    avgAvgDiff = Reduce ('+', mde)/length (mde)
    avgMAE = Reduce ('+', mae)/length (mae)
    avgMSE = Reduce ('+', mse)/length (mse)

    opRF = cbind (avgAvgDiff, avgMAE, avgMSE)
    colnames (opRF) = c ("rfAvgMD", "rfAvgMAE", "rfAvgMSE")

    dae = select (dae, .imp, genderf, reth_black, reth_hisp, reth_oth, age, practypec, site_ureter, site_renal, site_urethra, smokey, dgradeh, surgery, b.ecogvalue, var1, var2)

    diffs = lapply (1:10, function (x){
        imputed = filter (dae, .imp==x)
        imputed = select (imputed, -.imp)
        diff = as.matrix (imputed) - as.matrix (cDat)
        return (diff)
    })

    mde = lapply (diffs, function (x){
        apply (x, 2, mean)
    })
    mae = lapply (diffs, function (x){
        err = abs (x)
        apply (err, 2, mean)
    })
    mse = lapply (diffs, function (x){
        err = x^2
        apply (err, 2, mean)
    })
    avgAvgDiff = Reduce ('+', mde)/length (mde)
    avgMAE = Reduce ('+', mae)/length (mae)
    avgMSE = Reduce ('+', mse)/length (mse)

    opDAE = cbind (avgAvgDiff, avgMAE, avgMSE)
    colnames (opDAE) = c ("daeAvgMD", "daeAvgMAE", "daeAvgMSE")

    X = rownames (opRF)

    op = cbind (opMICE, opRF, opDAE)
    rownames (op) = NULL
    op = cbind (X, as.data.frame (op))
    op$mech = mdm
    op$prop = pc
    return (op)
}

MCAR10 = lapply (1:iter, function (x){
    getDiff(x, 'MCAR', 10)
})
MCAR30 = lapply (1:iter, function (x){
    getDiff(x, 'MCAR', 30)
})
MCAR50 = lapply (1:iter, function (x){
    getDiff(x, 'MCAR', 50)
})

MAR10 = lapply (1:iter, function (x){
    getDiff(x, 'MAR', 10)
})
MAR30 = lapply (1:iter, function (x){
    getDiff(x, 'MAR', 30)
})
MAR50 = lapply (1:iter, function (x){
    getDiff(x, 'MAR', 50)
})

MNAR10a = lapply (1:iter, function (x){
    getDiff(x, 'MNAR1', 10)
})
MNAR30a = lapply (1:iter, function (x){
    getDiff(x, 'MNAR1', 30)
})
MNAR50a = lapply (1:iter, function (x){
    getDiff(x, 'MNAR1', 50)
})

MNAR10b = lapply (1:iter, function (x){
    getDiff(x, 'MNAR2', 10)
})
MNAR30b = lapply (1:iter, function (x){
    getDiff(x, 'MNAR2', 30)
})
MNAR50b = lapply (1:iter, function (x){
    getDiff(x, 'MNAR2', 50)
})

MCAR10long = do.call ('rbind', MCAR10)
MCAR30long = do.call ('rbind', MCAR30)
MCAR50long = do.call ('rbind', MCAR50)
MCARall = rbind (MCAR10long, MCAR30long, MCAR50long)

MAR10long = do.call ('rbind', MAR10)
MAR30long = do.call ('rbind', MAR30)
MAR50long = do.call ('rbind', MAR50)
MARall = rbind (MAR10long, MAR30long, MAR50long)

MNAR10along = do.call ('rbind', MNAR10a)
MNAR30along = do.call ('rbind', MNAR30a)
MNAR50along = do.call ('rbind', MNAR50a)
MNAR1all = rbind (MNAR10along, MNAR30along, MNAR50along)

MNAR10blong = do.call ('rbind', MNAR10b)
MNAR30blong = do.call ('rbind', MNAR30b)
MNAR50blong = do.call ('rbind', MNAR50b)
MNAR2all = rbind (MNAR10blong, MNAR30blong, MNAR50blong)

write.csv (MCARall, file.path (rootdir, 'evalImputations_MCAR.csv'))
write.csv (MARall, file.path (rootdir, 'evalImputations_MAR.csv'))
write.csv (MNAR1all, file.path (rootdir, 'evalImputations_MNAR1.csv'))
write.csv (MNAR2all, file.path (rootdir, 'evalImputations_MNAR2.csv'))
