#Functions

#Login details for database that hide password and username
getLoginDetails <- function(){
  ## Based on code by Barry Rowlingson
  ## http://r.789695.n4.nabble.com/tkentry-that-exits-after-RETURN-tt854721.html#none
  require(tcltk)
  tt <- tktoplevel()
  tkwm.title(tt, "Get login details")
  Name <- tclVar("Login ID")
  Password <- tclVar("Password")
  entry.Name <- tkentry(tt,width="20", textvariable=Name)
  entry.Password <- tkentry(tt, width="20", show="*", 
                            textvariable=Password)
  tkgrid(tklabel(tt, text="Please enter your login details."))
  tkgrid(entry.Name)
  tkgrid(entry.Password)
  
  OnOK <- function()
  { 
    tkdestroy(tt) 
  }
  OK.but <-tkbutton(tt,text=" OK ", command=OnOK)
  tkbind(entry.Password, "<Return>", OnOK)
  tkgrid(OK.but)
  tkfocus(tt)
  tkwait.window(tt)
  
  invisible(c(loginID=tclvalue(Name), password=tclvalue(Password)))
}




#pedigree to prune pedigree
prune_pedigree=function(x="original germplasm ids", y="pedigree table or dataframe"){
      current_id=unique(x) #get all original id
      pruned_id=current_id #initialize a vector to collected germplasm related to pedigree put in original plants
      while(length(current_id) > 0){ #as long as there are unique ids to find
            iteration=y[y$germplasm_id %in% current_id,] #get parents out
            parents=c(iteration$male_parent, iteration$female_parent) #join father and mother
            parents=unique(parents) #only get unique
            parents_to_find=parents[!(parents%in%pruned_id)]
            pruned_id=c(pruned_id, parents_to_find) #update the individuals that need to be kept
            current_id=parents_to_find #update iterator
      }
      out=y[y$germplasm_id %in% pruned_id,] #only select individuals
      out=out[order(out$record_id),] #order by record id which is chronological
      return(out) #return pruned dataframe ready to make pedigree file
}

#function to fill any missing row columns to complete grid, values are missing
fill_gap=function(x, y="row max", z="range max"){
  combs=expand.grid(1:y, 1:z) #make all combinations
  colnames(combs)=c("row", "range")
  #join with x
  combs_join=merge(combs, x, by=c("row", "range"), all.x=TRUE)
  return(combs_join)
}

#Function to calculate standard error of heritabilty similar to pin function in standalone
#http://www.homepages.ed.ac.uk/iwhite/asreml/uop
pin <- function (object, transform) {
pframe <- as.list(object$gammas)
names(pframe) <- paste("V", seq(1, length(pframe)), sep = "")
tvalue <- eval(deriv(transform[[length(transform)]], names(pframe)), 
pframe)
X <- as.vector(attr(tvalue, "gradient"))
tname <- if (length(transform) == 3) 
transform[[2]]
else ""
n <- length(pframe)
i <- rep(1:n, 1:n)
j <- sequence(1:n)
k <- 1 + (i > j)
Vmat <- object$ai
se <- sqrt(sum(Vmat * X[i] * X[j] * k))
data.frame(row.names = tname, Estimate = tvalue, SE = se)
}


###The optimized "moving grid" function.### 
movingGrid2 <- function(Coords, obsPhe, radius, excludeCenter=TRUE){
  N    = nrow(Coords)
  idx  = 1:N
  dMat = sapply(idx, function(x) distHaversine(Coords, Coords[x,], r = 6378137))
  incl = 'if'(excludeCenter, sapply(idx, function(x) which(dMat[x,] <= radius  & idx!=x)),
                             sapply(idx, function(x) which(dMat[x,] <= radius)))
  
  xi   = sapply(idx, function(x) mean(obsPhe[incl[[x]]]))
  nVal = sapply(incl, length)
  
  b      = coef(lm(obsPhe ~ xi, na.action = "na.exclude"))[2]
  adjPhe = obsPhe - b*(xi - mean(xi, na.rm=T))
  ret    = cbind(Coords, 
                 adjustedPhe = adjPhe,
                 observedPhe = obsPhe,
                 movingMean  = xi, 
                 nValues     = nVal)
  return(ret)
}


###################Barplot function########################
############Credit https://github.com/mrxiaohe/R_Functions/blob/master/functions/bar#########
bar <- function(dv, factors, dataframe, percentage=FALSE, errbar=!percentage, half.errbar=TRUE, conf.level=.95, 
        xlab=NULL, ylab=NULL, main=NULL, names.arg=NULL, bar.col="black", whisker=.015,args.errbar=NULL,
        legend=TRUE, legend.text=NULL, args.legend=NULL,legend.border=FALSE, box=TRUE, args.yaxis=NULL, 
        mar=c(5,4,3,2),...){
    axes=!percentage
    dv.name<-substitute(dv)
    if(length(dv.name)>1) stop("'dv' only takes one variable")
    dv.name<-as.character(dv.name)
    dv<-dataframe[[dv.name]]
    fnames<-substitute(factors)
    if(length(fnames)==1){
        factors<-as.character(fnames)
        nf<-1
    }else{
        factors<-as.character(fnames[-1L])
        nf<-length(factors)
    }
    if(nf>2) stop("This function accepts no more than 2 factors \n",
            "\t-i.e., it only plots one-way or two-way designs.")
    if(percentage & errbar){
        warning("percentage=TRUE; error bars were not plotted")
        errbar<-FALSE
    }
    if(!percentage) xbars<-tapply(dv, dataframe[,factors], mean, na.rm=TRUE)
    else {
        xbars<-tapply(dv, list(interaction(dataframe[,factors], lex.order=TRUE)), mean, na.rm=TRUE)
        if(sum(na.omit(dv)!=0&na.omit(dv)!=1)>0) 
            stop("Data points in 'dv' need to be 0 or 1 in order to set 'percentage' to TRUE")
        xbars<-rbind(xbars, 1-xbars)*100
    }
    if(errbar){
        #se<-tapply(dv, dataframe[,factors], sd, na.rm=TRUE)/sqrt(tapply(dv, dataframe[,factors], length))
        #conf.level=1-(1-conf.level)/2
        lo.bar<-xbars-dataframe$se
        hi.bar<-xbars+dataframe$se
          #lo.bar<-xbars-se*qnorm(conf.level)
        #hi.bar<-xbars+se*qnorm(conf.level) 
    }
#   if(errbar){
#       se<-tapply(dv, dataframe[,factors], sd, na.rm=TRUE)/sqrt(tapply(dv, dataframe[,factors], length))
#       conf.level=1-(1-conf.level)/2
#       lo.bar<-xbars-se*qnorm(conf.level)
#       hi.bar<-xbars+se*qnorm(conf.level)  
#   }
    extras<-list(...)
    if(legend & !percentage){
        if(is.null(legend.text))
            legend.text<-sort(unique(dataframe[[factors[1]]]))
        args.legend.temp<-list(x="topright", bty=if(!legend.border)"n" else "o",
                               inset=c(0,0))
        if(is.list(args.legend))
            args.legend<-modifyList(args.legend.temp, args.legend)
        else 
            args.legend<-args.legend.temp
    } else if(legend & percentage){
        if(is.null(legend.text)) 
            legend.text<-c("1", "0")
        args.legend.temp<-list(x="topright", bty=if(!legend.border)"n" else "o",
                               inset=c(0,0))
        if(is.list(args.legend))
            args.legend<-modifyList(args.legend.temp, args.legend)
        else 
            args.legend<-args.legend.temp
    } else if(!legend){
        args.legend<-NULL
        legend.text<-NULL
    }
    if(errbar && legend && !percentage) ymax<-max(hi.bar)+max(hi.bar)/20
    else if(errbar && legend && percentage) ymax<-115
    else if(errbar && !legend) ymax <- max(xbars)
    else if(!errbar && legend && percentage) ymax<-110  
    else if(!errbar) ymax<-max(xbars) + max(xbars)/20
    if(!percentage){
        args.barplot<-list(beside=TRUE, height=xbars, ylim=c(0, ymax), main=main, names.arg=names.arg,
                col=hcl(h=seq(0,270, 270/(length(unique(dataframe[[factors[1]]]))))[-length(unique(dataframe[[factors[1]]]))]),
                legend.text=legend.text, args.legend=args.legend, xpd=TRUE,
                xlab=if(is.null(xlab)) factors[length(factors)] else xlab,
                ylab=if(is.null(ylab)) dv.name else ylab, axes=axes)
    }else{
        args.barplot<-list(beside=TRUE, height=xbars, ylim=c(0, ymax),  main=main, names.arg=names.arg,
                col=hcl(h=seq(0,270, 270/(length(unique(dataframe[[factors[1]]]))))[-length(unique(dataframe[[factors[1]]]))]),
                legend.text=legend.text, args.legend=args.legend, xpd=TRUE,
                xlab=if(is.null(xlab)) " "[length(factors)] else xlab,
                ylab=if(is.null(ylab)) "percentage" else ylab, axes=axes)       
    }
    args.barplot<-modifyList(args.barplot, extras)
    errbars = function(xvals, cilo, cihi, whisker, nc, args.errbar = NULL, half.errbar=TRUE) {
        if(half.errbar){
            cilo<-(cihi+cilo)/2
        }
        fixedArgs.bar = list(matlines, x=list(xvals), 
                             y=lapply(split(as.data.frame(t(do.call("rbind", 
                             list(cihi, cilo)))),1:nc),matrix, 
                             nrow=2, byrow=T))
        allArgs.bar = c(fixedArgs.bar, args.errbar)
        whisker.len = whisker*(par("usr")[2] - par("usr")[1])/2
        whiskers = rbind((xvals - whisker.len)[1,],
                         (xvals + whisker.len)[1,])
        fixedArgs.lo = list(matlines, x=list(whiskers),     
                            y=lapply(split(as.data.frame(t(do.call("rbind", 
                            list(cilo, cilo)))), 1:nc), matrix, nrow=2, byrow=T))
        allArgs.bar.lo = c(fixedArgs.lo, args.errbar)
        fixedArgs.hi = list(matlines, x=list(whiskers), 
                            y=lapply(split(as.data.frame(t(do.call("rbind", 
                            list(cihi, cihi)))), 1:nc), matrix, nrow=2, byrow=T))
        allArgs.bar.hi = c(fixedArgs.hi, args.errbar)  
        invisible(do.call(mapply, allArgs.bar))
        if(!half.errbar) invisible(do.call(mapply, allArgs.bar.lo))
        invisible(do.call(mapply, allArgs.bar.hi))
    }
    par(mar=mar)
    errloc<-as.vector(do.call(barplot, args.barplot))
    if(errbar){
        errloc<-rbind(errloc, errloc)
        lo.bar<-matrix(as.vector(lo.bar))
        hi.bar<-matrix(as.vector(hi.bar))
        args.errbar.temp<-list(col=bar.col, lty=1)
        args.errbar<-if(is.null(args.errbar)|!is.list(args.errbar)) 
                        args.errbar.temp
                     else if(is.list(args.errbar)) 
                        modifyList(args.errbar.temp, args.errbar)
        errbars(errloc, cilo=lo.bar, cihi=hi.bar, nc=1, whisker=whisker, 
                args.errbar=args.errbar, half.errbar=half.errbar)
    }
    if(box) box()
    if(percentage){
        args.yaxis.temp<-list(at=seq(0,100, 20), las=1)
        args.yaxis<-if(!is.list(args.yaxis)) args.yaxis.temp else modifyList(args.yaxis.temp, args.yaxis)
        do.call(axis, c(side=2, args.yaxis))
    }
}


#Function that counts allele depth by percentile for depth output
adepth=function(x="name of samples", y="depth dataframe"){ #function that calculates percentiles of tag depth across individuals and sites
  s_name=unique(x$FullSampleName)
  dframe=y[,colnames(y)%in%s_name]
  dpcall=apply(dframe, 1, FUN=function(x){quantile(x,c(0.05,0.1, 0.5, 0.9,.95))})
  dpcall=as.data.frame(t(dpcall))
  colnames(dpcall)=c("fifth percentile", "tenth percentile", "fiftieth percentile", "ninetieth percentile", "ninety-fifth percentile")
  dpM=melt(dpcall, id.vars=NULL)
  dpicall=apply(dframe, 2, FUN=function(x){quantile(x,c(0.05,0.1, 0.5, 0.9,.95))})
  dpicall=as.data.frame(t(dpicall))
  colnames(dpicall)=c("fifth percentile", "tenth percentile", "fiftieth percentile", "ninetieth percentile", "ninety-fifth percentile")
  dpiM=melt(dpicall)
  
  return(list(dpM, dpiM))
  }
  
  ####Function to count alleles and populations parameters###
tassel5_to_params=function(x="hap matrix", y="columns to skip", z="population number"){
      geno=x
      #recount allele A and B and het
      alleleA=rowSums(geno[,(y+1):ncol(geno)]!=substring(geno$alleles, 3, 3) & geno[,(y+1):ncol(geno)]!="N") #only counts what is not allele B and missing.  i.e. counts allele A and various calls for heterozygous
      alleleB=rowSums(geno[,(y+1):ncol(geno)]!=substring(geno$alleles, 1, 1) & geno[,(y+1):ncol(geno)]!="N")
      het=rowSums(geno[,(y+1):ncol(geno)] == "M") + rowSums( geno[,(y+1):ncol(geno)] ==   "R") + rowSums(geno[,(y+1):ncol(geno)] ==  "W") + rowSums(geno[,(y+1):ncol(geno)] ==  "K") + rowSums(geno[,(y+1):ncol(geno)] ==  "S") + rowSums(geno[,(y+1):ncol(geno)] ==  "Y")
      present=1-(rowSums(geno[,(y+1):ncol(geno)]=="N")/z)
      MAF=apply(cbind(((alleleA-het)*2+het), (alleleB-het)*2+het), 1, min)/apply(
    cbind(((alleleA-het)*2+het), ((alleleB-het)*2+het)), 1, sum) 
      percentHet=het/apply(cbind(alleleA-het, alleleB-het, het), 1, sum)
      return(cbind.data.frame(geno[,1:y], "alleleA"=alleleA, "alleleB"=alleleB, "het"=het, "present"= present, "MAF"=MAF, "percentHET"=percentHet, geno[,(y+1):ncol(geno)]))
}

##function to convert hap to 0 and 1
hap_to_G=function(x="hap matrix", y="number of columns of information"){
  ##From Prasana, pulls out first allele for a and second for b
  a = substring(x$alleles,1,1)
  #Checks the frequency of the alleles if the second allele is more frequent it is substitued
  a[x$alleleA<x$alleleB] = substring(x$alleles,3,3)[x$alleleA<x$alleleB]
  #Same thing with the second allele
  b = substring(x$alleles,3,3)
  b[x$alleleA<x$alleleB] = substring(x$alleles,1,1)[x$alleleA<x$alleleB]
  #Checks to make sure all alleles are one or the other
  #print(paste("If 0 all alleles are accounted for: ", sum(a == b), sep=""))
  
  ## Turn into letter matrix for mapping
  #makes a copy of the hap matrix
  hap01 = x
  #sets all allele values to NA
  hap01[,(y+1):ncol(hap01)]=NA
  
  ## Turn allele a and allele b into 1 and -1.  Het into 0
  #line by line if a line is a then it places 1 in hap01 for the allele
  hap01[x == a] = 1
  hap01[x == b] = -1
  hap01[x == "M"] = 0
  hap01[x == "Y"] = 0
  hap01[x == "K"] = 0
  hap01[x == "R"] = 0
  hap01[x == "W"] = 0
  hap01[x == "S"] = 0
  hap01[x== "N"]=NA
  
  return(hap01)}
  
#function to make dendrogram mapping matrix
dendrogram=function(x="genotype matrix, pass SNP and then individual columns"){
#Make original dendrogram using the binary matrix
library(ape)
library(dendextend)
library(stringr)
tree=t(x)
colnames(tree)=tree[1,]
#Tree is binary matrix can be used for all graphs
tree=tree[-1,]
hc=hclust(dist(tree))
a=dist(tree)
hc=hclust(a)
return(hc)}


#function to make sampled cross fold validation
cv_fold <- function(x = 'data frame', cv = 'cv number', ngroup = 'number per group', ss = 'starting seed' ){
  counter <- 0 #start a dummy counter
  for(i in 1:cv){ #for iteration of the cv
    n_group <- ngroup #get the number of cross fold goups
    set.seed(ss + counter) #set the seed, different each time
    ord <- sample(1:n_group, replace = TRUE, size = nrow(x)) #randomly sample the values there
    ord <- as.data.frame(ord) #made a dataframe
    colnames(ord) <- paste('CV_', n_group, '_Iter_', i, sep = '') #give the dataframe a heading
    x <- cbind.data.frame(x, ord) #join to dataframe
    counter <- counter + 1 #update the counter so the seed changes
  }
  return(x) #return the dataframe
}
