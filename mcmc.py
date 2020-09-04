import numpy as np
import pandas as pd

from define_alleles import *
from calculate_frequencies import *
from recode_alleles import *
from switch_hidden import *
from findposteriorfrequencies import *

# MOI = multiplicity of infection
maxMOI = np.nanmax( # Return array max, ignoring NaNs
    # NOTE: Assuming genotypedata_RR is a pandas dataframe
    # Split string like so: https://cmdlinetips.com/2018/06/how-to-split-a-column-or-column-names-in-pandas-and-get-part-of-it/
    # Gets the
    pd.to_numeric(genotypedata_RR.columns.str.split("_").str[1])
)

# Get the unique Sample IDs in the dataset
ids = np.unique(genotypedata_RR[genotypedata_RR["Sample.ID"].str.contains("Day 0")]["Sample.ID"].str.replace(" Day 0", ""))
# Ditto, the unique loci for the set
locinames = np.unique(genotypedata_RR.columns[1:].str.split("_").str[0])

nids = ids.size
nloci = locinames.size

maxalleles=30

k = np.repeat(maxalleles, nloci)

alleles_definitions_RR  = define_alleles(pd.concat([genotypedata_RR,additional_neutral]), locirepeats, k)

##### calculate MOI
# TODO: What do MOI0 and MOIf stand for? Multiplicity of infection on first day/day of failure? What are they used for/mean?
MOI0 = np.repeat(0,nids)
MOIf = np.repeat(0,nids)
for i, ID in enumerate(ids):
    for lociname in locinames:
        locicolumns = genotypedata_RR.columns.str.contains(f"{lociname}_")

        nalleles0 = np.count_nonzero(~genotypedata_RR.loc[
            genotypedata_RR["Sample.ID"].str.contains(f"{ID} Day 0"),
            locicolumns]
            .isna())
        nallelesf = np.count_nonzero(~genotypedata_RR.loc[
            genotypedata_RR["Sample.ID"].str.contains(f"{ID} Day Failure"),
            locicolumns]
            .isna())

        MOI0[i] = np.max([MOI0[i],nalleles0])
        MOIf[i] = np.max([MOIf[i],nallelesf])

##### define statevector (i.e. initial MCMC state)
# TODO: What are each of these arrays used for? Should they not be initialized until closer to where they're actually used?
alleles0        = np.zeros((nids, maxMOI*nloci))
recoded0        = np.zeros((nids, maxMOI*nloci))
hidden0         = np.full_like(np.empty((nids, maxMOI*nloci)), np.nan)
recr0           = np.full_like(np.empty((nids, nloci)), np.nan)
recr_repeats0   = np.full_like(np.empty((nids, nloci)), np.nan) # number of times recrudescing allele is repeated on day 0
recr_repeatsf   = np.full_like(np.empty((nids, nloci)), np.nan) # number of times recrudescing allele is repeated on day 0
allelesf        = np.zeros((nids, maxMOI*nloci))
recodedf        = np.zeros((nids, maxMOI*nloci))
hiddenf         = np.full_like(np.empty((nids, maxMOI*nloci)), np.nan)
recrf           = np.full_like(np.empty((nids, nloci)), np.nan)
mindistance     = np.zeros((nids, nloci))
alldistance     = np.full_like(np.empty((nids, nloci, maxMOI**2)), np.nan)
allrecrf        = np.full_like(np.empty((nids, nloci, maxMOI**2)), np.nan)
classification = np.repeat(0, nids)

##### create state 0 (essentially, set NAN values to 0 and fill alleles with the appropriate initial/failure data)
for i, locus in enumerate(locinames):
    # locicolumns code is duplicated from MOI calculations
    locicolumns = genotypedata_RR.columns.str.contains(f"{locus}_")

    oldalleles = genotypedata_RR.loc[:,locicolumns].to_numpy()
    '''
    # TODO: What is this code doing? Doesn't seem necessary?
    if (len(oldalleles.shape[1]) == 0) {
        oldalleles = matrix(oldalleles,length(oldalleles),1)
    }
    '''
    newalleles = np.copy(oldalleles)
    ncolumns = oldalleles.shape[1]
    '''
    for j in range(ncolumns):
        # TODO: Can't test until I have a recodeallele implementation
        newalleles[:,j] = np.array(list(map(
            range(0, oldalleles.shape[0]),
            lambda x: recodeallele(alleles_definitions_RR[i], oldalleles[x,j]))))
    '''
    newalleles[np.isnan(newalleles)] = 0
    oldalleles[np.isnan(oldalleles)] = 0

    oldalleles[newalleles == 0] = 0

    startColumn = maxMOI*(i-1)  # TODO: Subtracted 1 for indexing reasons in Python vs R, but not for endColumn; double-check that's valid
    endColumnOldAllele = maxMOI*(i-1) + oldalleles.shape[1]
    endColumnNewAllele = maxMOI*(i-1) + newalleles.shape[1]
    alleles0[:, startColumn:endColumnOldAllele] = oldalleles[genotypedata_RR["Sample.ID"].str.contains("Day 0"),:]
    allelesf[:, startColumn:endColumnOldAllele] = oldalleles[genotypedata_RR["Sample.ID"].str.contains("Day Failure"),:]
    recoded0[:, startColumn:endColumnNewAllele] = newalleles[genotypedata_RR["Sample.ID"].str.contains("Day 0"),:]
    recodedf[:, startColumn:endColumnNewAllele] = newalleles[genotypedata_RR["Sample.ID"].str.contains("Day Failure"),:]

##### recode additional_neutral, if needed
# TODO: What does recoding do? Why is it needed?
# TODO: Seems to copy-paste much of the previous code section
recoded_additional_neutral = None
if additional_neutral.size > 0 and additional_neutral.shape[0] > 0:
    recoded_additional_neutral = np.zeros((additional_neutral.shape[0], maxMOI*nloci))
    for i, locus in enumerate(locinames):
        locicolumns = genotypedata_RR.columns.str.contains(f"{locus}_")

        oldalleles = additional_neutral[:,locicolumns] # TODO: stub
        '''
        # TODO: What is this code doing?
        if (len(oldalleles.shape[1]) == 0) {
            oldalleles = matrix(oldalleles,length(oldalleles),1)
        }
        '''
        newalleles = np.copy(oldalleles)
        ncolumns = oldalleles.shape[1]
        '''
        for j in range(ncolumns):
            # TODO: Can't test until I have a recodeallele implementation
            newalleles[:,j] = np.array(list(map(
                range(0, oldalleles.shape[0]),
                lambda x: recodeallele(alleles_definitions_RR[i], oldalleles[x,j]))))
        '''
        newalleles[np.isnan(newalleles)] = 0
        oldalleles[np.isnan(oldalleles)] = 0

        oldalleles[newalleles == 0] = 0

        startColumn = maxMOI*(i-1)  # TODO: Subtracted 1 for indexing reasons in Python vs R, but not for endColumn; double-check that's valid
        endColumnOldAllele = maxMOI*(i-1) + oldalleles.shape[1]
        recoded_additional_neutral[:, startColumn:endColumnOldAllele] = newalleles

## estimate frequencies
frequencies_RR = calculate_frequencies3(pd.concat(genotypedata_RR,additional_neutral),alleles_definitions_RR)

## assign random hidden alleles
# TODO: Figure out what this code is overall trying to do (replace non-0 elements with random values? Why is it important that the values are assigned from frequencies_RR?)
for i in range(nids):
    for j in range(nloci):
        # TODO: Code almost duplicated between top/bottom portions; refactor into single function? (needs 9 inputs: maxMOI, nids, nloci, MOIarray, alleles/recoded/hidden array, alleles_definitions_RR, frequencies_RR)
        # TODO: Start/end of what? The portion of the row w/ this locus information?
        start = maxMOI*j
        end = maxMOI*(j+1)

        nalleles0 = np.count_nonzero(alleles0[i, start:end])
        nmissing0 = MOI0[i] - nalleles0

        # TODO: Rename "nonzero_indices" and "zero_indices"?
        whichnotmissing0 = np.arange(start, end)[np.where(alleles0[i, start:start+MOI0[i]] != 0)]
        whichmissing0 = np.arange(start, end)[np.where(alleles0[i, start:start+MOI0[i]] == 0)]

        # Sample to randomly initialize the alleles/hidden variables
        if nalleles0 > 0:
            hidden0[i,whichnotmissing0] = 0
        if nmissing0 > 0:
            newhiddenalleles0 = np.random.choice(
                np.arange(0, int(frequencies_RR[0, j, 0])), # Select from first row (count of how many probabilities they are)
                size=nmissing0,
                replace=True,
                p=frequencies_RR[1, j, 0:int(frequencies_RR[0, j, 0])]
                    / frequencies_RR[1, j, 0:int(frequencies_RR[0, j, 0])].sum()) # Sum so probabilities add up to 1 (TODO: Can remove this when using real data and not just stubbing)
            recoded0[i,whichmissing0] = newhiddenalleles0
            # calculate row means
            alleles0[i,whichmissing0] = np.mean(alleles_definitions_RR[j], axis=1)[newhiddenalleles0] # hidden alleles get mean allele length
            hidden0[i,whichmissing0] = 1

        nallelesf = np.count_nonzero(allelesf[i, start:end])
        nmissingf = MOIf[i] - nallelesf

        # TODO: Rename "nonzero_indices" and "zero_indices"?
        whichnotmissingf = np.arange(start, end)[np.where(allelesf[i, start:start+MOIf[i]] != 0)]
        whichmissingf = np.arange(start, end)[np.where(allelesf[i, start:start+MOIf[i]] == 0)]

        if nallelesf > 0:
            hiddenf[i,whichnotmissingf] = 0
        if nmissingf > 0:
            newhiddenallelesf = np.random.choice(
                np.arange(0, int(frequencies_RR[0, j, 0])), # Select from first row (count of how many probabilities they are)
                size=nmissingf,
                replace=True,
                p=frequencies_RR[1, j, 0:int(frequencies_RR[0, j, 0])]
                    / frequencies_RR[1, j, 0:int(frequencies_RR[0, j, 0])].sum()) # Sum so probabilities add up to 1 (TODO: Can remove this when using real data and not just stubbing)
            recodedf[i,whichmissingf] = newhiddenallelesf
            # calculate row means
            allelesf[i,whichmissingf] = np.mean(alleles_definitions_RR[j], axis=1)[newhiddenallelesf] # hidden alleles get mean allele length
            hiddenf[i,whichmissingf] = 1

## initial estimate of q (probability of an allele being missed)
# TODO: What is qq? Is there a better name for this?
qq = np.nanmean(np.concatenate[hidden0, hiddenf])

## initial estimate of dvect (likelihood of error in analysis)
# TODO: What does dvect stand for?
dvect = np.zeros(1 + int(np.rint(
    # Get the range (max-min) of the first "nloci" rows, then the max of all those
    np.ptp(alleles_definitions_RR[0:nloci], axis=1).max()
)))
dvect[1] = 0.75
dvect[2] = 0.2
dvect[3] = 0.05

## randomly assign recrudescences/reinfections
for i in range(nids):
    z = np.random.uniform(size=1)
    if z < 0.5:
        classification[i] = 1
    for j in range(nloci): # determine which alleles are recrudescing (for beginning, choose closest pair)
        allpossiblerecrud = np.stack(np.meshgrid(
            np.arange(MOI0[i]),
            np.arange(MOIf[i]))).T.reshape(-1, 2)

        allele0_col_indices = maxMOI*j + allpossiblerecrud[:,0]
        allelef_col_indices = maxMOI*j + allpossiblerecrud[:,1]

        recrud_distances = np.abs(alleles0[i, allele0_col_indices] - allelesf[i, allelef_col_indices])
        # rename to "closest_recrud_index"?
        closestrecrud = np.argmin(recrud_distances)

        mindistance[i,j] = recrud_distances[closestrecrud]
        alldistance[i,j,:recrud_distances.size] = recrud_distances

        allrecrf[i,j,:allpossiblerecrud.shape[0]] = recodedf[i, maxMOI*j + allpossiblerecrud[:,1]]
        recr0[i,j] = maxMOI*j + allpossiblerecrud[closestrecrud,0]
        recrf[i,j] = maxMOI*j + allpossiblerecrud[closestrecrud,1]

        recr_repeats0[i,j] = np.sum(recoded0[i, maxMOI*j:maxMOI*(j+1)] == recoded0[i, int(recr0[i,j])]) # TODO: int() only needed for stub
        recr_repeatsf[i,j] = np.sum(recodedf[i, maxMOI*j:maxMOI*(j+1)] == recodedf[i, int(recrf[i,j])]) # TODO: int() only needed for stub
#===============================================================================
#   THE LINE OF SANITY
#   (code below this point has NOT been converted from R to Python)
#===============================================================================

#### correction factor (reinfection)
correction_distance_matrix = list() # for each locus, matrix of distances between each allele
for (i in 1:nloci) {
    correction_distance_matrix[[i]] = as.matrix(dist(rowMeans(alleles_definitions_RR[[i]])))
}


state_classification = matrix(NA,nids,(nruns-burnin)/record_interval)
state_alleles0 = array(NA,c(nids,maxMOI*nloci,(nruns-burnin)/record_interval))
state_allelesf = array(NA,c(nids,maxMOI*nloci,(nruns-burnin)/record_interval))
state_parameters = matrix(NA,2+2*nloci,(nruns-burnin)/record_interval)

count = 1
dposterior = 0.75
runmcmc = function() {
    # propose new classification
    # rellikelihood_reinfection = sapply(1:nids, function (x) (sum(log(frequencies_RR[[2]][cbind(1:nloci,recoded0[x,recrf[x,]])]))))
    #rellikelihood_recr = sapply(1:nids, function (x) (sum(log(dvect[round(mindistance[x,]+1)]))))
    # likelihoodratio = exp(rellikelihood_recr - rellikelihood_reinfection)
    # adjust for multiple corrections (ratio of multinomial coefficients)
    #likelihoodratio = sapply(1:nids, function (x) likelihoodratio[x]/exp(nloci*log(MOI0[x])+nloci*log(MOIf[x])-sum(log(recr_repeats0[x,]))-sum(log(recr_repeatsf[x,]))))
    #likelihoodratio = sapply(1:nids, function (x) exp(sum(log(sapply(1:nloci, function (y) mean(dvect[round(alldistance[x,y,])+1]/frequencies_RR[[2]][y,allrecrf[x,y,]],na.rm=TRUE))))))
    #likelihoodratio = sapply(1:nids, function (x) exp(sum(log(sapply(1:nloci, function (y) mean(dvect[round(alldistance[x,y,])+1]/colSums(frequencies_RR[[2]][y,1:frequencies_RR[[1]][y]]*matrix(dvect[correction_distance_matrix[[y]][,allrecrf[x,y,]]+1],frequencies_RR[[1]][y],frequencies_RR[[1]][y])),na.rm=TRUE))))))
    likelihoodratio = sapply(1:nids, function (x) exp(sum(log(sapply(1:nloci, function (y) mean(dvect[round(alldistance[x,y,])+1]/sapply(1:(maxMOI*maxMOI), function (z) sum(frequencies_RR[[2]][y,1:frequencies_RR[[1]][y]]*dvect[correction_distance_matrix[[y]][,allrecrf[x,y,z]]+1])),na.rm=TRUE))))))

    z = runif(nids)
    newclassification = classification
    newclassification[classification == 0 & z < likelihoodratio] = 1
    newclassification[classification == 1 & z < 1/likelihoodratio] = 0
    classification <<- newclassification

    # propose new hidden states
    sapply(1:nids, function (x) switch_hidden(x))

    # propose q (beta distribution is conjugate distribution for binomial process)
    q_prior_alpha = 0;
    q_prior_beta = 0;
    q_posterior_alpha = q_prior_alpha + sum(c(hidden0,hiddenf) == 1,na.rm=TRUE)
    q_posterior_beta = q_prior_beta + sum(c(hidden0,hiddenf)==0,na.rm=TRUE)
    if (q_posterior_alpha == 0) {
        q_posterior_alpha =1
    }
    qq <<- rbeta(1, q_posterior_alpha , q_posterior_beta)

    #  update dvect (approximate using geometric distribution)
    # only if there is at least 1 recrudescent infection
    if (sum(classification==1) >= 1) {
    d_prior_alpha = 0;
    d_prior_beta = 0;
    d_posterior_alpha = d_prior_alpha + length(c(mindistance[classification==1,]))
    d_posterior_beta = d_prior_beta + sum(c(round(mindistance[classification==1,])))
    if (d_posterior_beta == 0) {
        d_posterior_beta = sum(c((mindistance[classification==1,])))
    }
    if (d_posterior_beta == 0) { ## algorithm will get stuck if dposterior is allowed to go to 1
        d_posterior_beta = 1
    }


    dposterior <<- rbeta(1, d_posterior_alpha , d_posterior_beta)
    dvect = (1-dposterior) ^ (1:length(dvect)-1) * dposterior
    dvect <<- dvect / (sum(dvect))
    }

    # update frequencies
    # remove recrudescing alleles from calculations
    tempdata = recoded0
    sapply(which(classification == 1), function (x) tempdata[x,recr0[x,]] <<- 0)
    tempdata = rbind(tempdata, recodedf)
    sapply(1:nloci, function (x) findposteriorfrequencies(x,rbind(tempdata,recoded_additional_neutral)))

    # record state
    if (count > burnin & count %% record_interval == 0) {
        print(count)
        state_classification[,(count-burnin)/record_interval] <<- classification
        state_alleles0[,,(count-burnin)/record_interval] <<- alleles0
        state_allelesf[,,(count-burnin)/record_interval] <<- allelesf
        state_parameters[1,(count-burnin)/record_interval] <<- qq
        state_parameters[2,(count-burnin)/record_interval] <<- dposterior
        state_parameters[3:(3+nloci-1),(count-burnin)/record_interval] <<- apply(frequencies_RR[[2]],1,max)
        state_parameters[(3+nloci):(3+2*nloci-1),(count-burnin)/record_interval] <<- sapply(1:nloci,function (x) sum(frequencies_RR[[2]][x,]^2))

    }
    count <<- count + 1
}

replicate(nruns,runmcmc())

## make sure no NAs in result matrices
state_parameters = state_parameters[,!is.na(colSums(state_parameters))]
state_classification = state_classification[,!is.na(colSums(state_classification))]

## find mode of hidden alleles
modealleles = matrix("",2*nids,maxMOI*nloci)
for (i in 1:nids) {
    for (j in 1:nloci) {
        modealleles[2*(i-1)+1,((j-1)*maxMOI+1):(j*maxMOI)] = sapply(1:maxMOI, function (x) names(table(state_alleles0[i,(j-1)*maxMOI+x,]))[table(state_alleles0[i,(j-1)*maxMOI+x,])== max(table(state_alleles0[i,(j-1)*maxMOI+x,]))][1])
        modealleles[2*(i-1)+2,((j-1)*maxMOI+1):(j*maxMOI)] = sapply(1:maxMOI, function (x) names(table(state_allelesf[i,(j-1)*maxMOI+x,]))[table(state_allelesf[i,(j-1)*maxMOI+x,])== max(table(state_allelesf[i,(j-1)*maxMOI+x,]))][1])
    }
}

rowMeans2 = function(x){
    if (length(dim(x)) == 0) {
        ret = mean(x)
    } else {
        ret = rowMeans(x)
    }
    ret
}

temp_combined = c(sapply(1:length(ids), function (x) rep(rowMeans2(state_classification)[x],2)))
outputmatrix = cbind(temp_combined,modealleles)
colnames(outputmatrix) = c("Prob Rec",c(sapply(1:nloci, function (x) paste(locinames[x],"_",1:maxMOI,sep=""))))
write.csv(outputmatrix,paste(jobname,"_posterior",".csv",sep=""))


# summary statistics of parameters
write.csv(state_parameters,paste(jobname,"_state_parameters",".csv",sep=""))

summary_statisticsmatrix = cbind(format(rowMeans(state_parameters),digits=2),
                       apply(format(t(sapply(1:dim(state_parameters)[1], function (x) quantile(state_parameters[x,],c(0.25,0.75)))),digits=2),1, function (x) paste(x,collapse="�")))
summary_statisticsmatrix = rbind(summary_statisticsmatrix, c(format(mean(state_parameters[(3+nloci):(3+2*nloci-1),]),digits = 2),paste(format(quantile(state_parameters[(3+nloci):(3+2*nloci-1),],c(0.25,0.75)),digits=2),collapse="�")))
summary_statisticsmatrix = as.matrix(sapply(1:dim(summary_statisticsmatrix)[1], function (x) paste(summary_statisticsmatrix[x,1], " (",summary_statisticsmatrix[x,2],")",sep="")))
rownames(summary_statisticsmatrix) = c("q","d",locinames,locinames,"Mean diversity")
write.csv(summary_statisticsmatrix,paste(jobname,"_summarystatistics",".csv",sep=""))
