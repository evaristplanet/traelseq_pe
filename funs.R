preprocess <- function(file.fwd, file.rev, cross.out, rcross.out.negpos, data.out) {
    cat('read in\n')
    fw <- read.csv(file.fwd, sep='\t', header=F)
    re <- read.csv(file.rev, sep='\t', header=F)

    cat('merge\n')
    fw.idx <- paste0(fw$V1, ':', fw$V2)
    re.idx <- paste0(re$V1, ':', re$V2)
    all <- unique(c(fw.idx, re.idx))
    fw <- fw[match(all, fw.idx), ]
    re <- re[match(all, re.idx), ]
    chr <- unlist(lapply(strsplit(all, ':'), function(x) x[1]))
    pos <- unlist(lapply(strsplit(all, ':'), function(x) x[2]))
    pos <- as.numeric(pos) + as.numeric(pos[2])/2

    cat('put blanks to 0\n')
    fw$V4[is.na(fw$V4)] <- 0
    re$V4[is.na(re$V4)] <- 0

    cat('compute rp\n')
    rp <- (fw$V4 - re$V4) / (fw$V4 + re$V4)

    cat('export bedGraph of rp\n')
    xout <- data.frame(fw[, 1:3], rp)
    fout <- gsub('_fwd.bg', '_rp.bg', file.fwd)
    write.table(xout, fout, sep='\t', quote=F, row.names=F, col.names=F)

    cat('remove extremes\n')
    rp[abs(rp)==1] <- NA

    cat('get crossing points\n')
    tmp <- sign(rp)
    start <- which(!is.na(tmp))[1]
    bef <- tmp[start]
    id <- vector('integer', 0)
    for (i in (start +1):length(rp)) {
        if (!is.na(tmp[i]) & tmp[i]!=0) {
            if (tmp[i] != bef) {
                id <- c(id, i)
            }
            bef <- tmp[i]
        }
    }
    firstInChr <- which(!duplicated(chr))
    id <- id[!(id %in% firstInChr)]
    cp <- pos[id]
    xout <- data.frame(chr[id], as.integer(cp), as.integer(cp))
    write.table(xout, cross.out, sep='\t', col.names=F, row.names=F, quote=F)

    cat('get crossing negpos points\n')
    sel <- sign(rp[id])==1
    xout <- xout[sel, ] 
    write.table(xout, cross.out.negpos, sep='\t', col.names=F, row.names=F, quote=F)

    cat('save output\n')
    save(chr, pos, rp, file=rdata.out, compress='gzip')
}

getPerc <- function(file.fwd, file.rev) {
    cat('read in\n')
    fw <- read.csv(file.fwd, sep='\t', header=F)
    re <- read.csv(file.rev, sep='\t', header=F)

    cat('merge\n')
    fw.idx <- paste0(fw$V1, ':', fw$V2)
    re.idx <- paste0(re$V1, ':', re$V2)
    all <- unique(c(fw.idx, re.idx))
    fw <- fw[match(all, fw.idx), ]
    re <- re[match(all, re.idx), ]
    chr <- unlist(lapply(strsplit(all, ':'), function(x) x[1]))
    pos <- unlist(lapply(strsplit(all, ':'), function(x) x[2]))
    pos <- as.numeric(pos) + as.numeric(pos[2])/2

    cat('put blanks to 0\n')
    fw$V4[is.na(fw$V4)] <- 0
    re$V4[is.na(re$V4)] <- 0

    cat('compute perc\n')
    addUp <- fw$V4 + re$V4
    perc <- re$V4 / addUp

    ans <- data.frame(chr, pos, addUp, perc)
    ans
}

sexyPlot <- function(chrom, post.st, pos.en, rdata, xlab, fileout, main, width=4*7, height=2*7, rmOnes=T, pt.cex=4) {
    cat('load data\n')
    dummy <- load(rdata)

    cat('subset region\n')
    sel <- chr==chrom & pos>=pos.st  & pos<=pos.en
    pos.s <- pos[sel]
    rp.s <- rp[sel]
    o <- order(pos.s)
    pos.s <- pos.s[o]
    rp.s <- rp.s[o]
    pos.s <- pos.s / 1e+06

    write.table(xout, cross.out, sep='\t', col.names=F, row.names=F, quote=F)

    cat('remove extremes\n')
    rp.s[abs(rp.s)==1] <- NA

    cat('get crossing points\n')
    tmp <- sign(rp.s)
    while (any(is.na(tmp) | tmp==0)) {
        id <- which(is.na(tmp) | tmp==0)
        tmp[id] <- tmp[id -1]
    }
    id <- which(tmp[-length(tmp)]!=tmp[-1])
    id <- id[tmp[id +1]==1]
    cp <- pos.s[id] + (pos.s[2] - pos.s[1]) / 2

    cat('plot\n')
    ylim <- c(-.5, .5)
    ylab <- 'Read polarity'
    pdf(fileout, width=width, height=height)
    sel.nnan <- !is.nan(rp.s) #when all 0, we remove 
    col <- ifelse(rp.s[sel.nnan]>0, 'firebrick3', 'deepskyblue3')
    plot(pos.s[sel.nnan], rp.s[sel.nnan], axes=T, ylim=ylim, xlab=xlab, ylab=ylab, pch=20, col=col, main=main, cex=pt.cex)
    abline(h=0)
    abline(v=cp, col='darkgrey', lty=1)
    dev.off()

}

getNBases <- function(file) {
    tmp <- read.csv(file, header=F, sep='\t')
    ans <- sum(tmp$V3 - tmp$V2)
    ans
}

addCor <- function(x, y, legend.y=.95, legend.col=1) {
    if (length(x)>10) {
        mycor <- cor.test(x, y, method='pearson') 
        legend <- paste0('n=', length(x), ' | rho=', formatC(mycor$estimate), ', | pvalue=', formatC(mycor$p.value))
        ans <- text(.2, legend.y, legend, col=legend.col, cex=.7)
    } else {
        ans <- ''
    }
    ans
}

getVals <- function(sample, file, each=500) {
    inPos <- read.csv(file, sep='\t', header=F)
    inPos <- data.frame(chr=inPos[, 1], 
                        pos=inPos[, 2] + (inPos[, 3] - inPos[, 2]) / 2)
    ans <- lapply(as.list(1:nrow(inPos)), function(i) {
                      selchr <- inPos[i, 1]
                      selpos <- inPos[i, 2]
                      selpos <- round(selpos/each, 0) * each
                      length <- 14000 / (each *2) +1
                      selpos <- sort(unique(c(seq(selpos, by=-(each *2), length=length), 
                                              seq(selpos, by=(each *2), length=length))))
                      selpos[selpos <0] <- NA
                      sample.chr <- sample[sample$chr==selchr, ]
                      idx <- match(selpos, sample.chr$pos)
                      ans <- sample.chr$perc[idx]
                      names(ans) <- seq(-14000, 14000, by=each *2)
                      ans
                        })
    ans <- do.call(rbind, ans)
    ans
}

getDist <- function(bed1, bed2, middle=F, ssamplebed1, ssamplebed2) {
    tmp <- read.csv(bed1, sep='\t', header=F)
    if (middle) {
        middle <- tmp[, 2] + (tmp[, 3] - tmp[, 2]) /2
        tmp <- data.frame(V1=tmp[, 1], V2=as.integer(middle), V3=as.integer(middle +1)) 
    }
    if (!missing(ssamplebed1)) {
        tmp <- tmp[sort(sample(nrow(tmp))[1:ssamplebed1]), ]
    }
    if (!missing(ssamplebed2)) {
        cmd <- paste0('shuf ', bed2, ' | head -n ', ssamplebed2, 
                      ' | sort -k1,1 -k2,2n > /tmp/shufBed.bed')
        system(cmd)
        bed2 <- '/tmp/shufBed.bed'
    }
    ans <- addDist2bed(tmp, bed2, sortBed=T)$dist2feature
    ans
}

getPvalsDist <- function(s1.up, s2.up, s1.do, s2.do, s1.ot, s2.ot) {
    t.txt <- paste0('K9.up. S1 vc S2 pvalue=', formatC(wilcox.test(s1.up, s2.up)$p.value), '\n',
                    'K9.down. S1 vc S2 pvalue=', formatC(wilcox.test(s1.do, s2.do)$p.value), '\n',
                    'K9.other. S1 vc S2 pvalue=', formatC(wilcox.test(s1.ot, s2.ot)$p.value))
    t.txt.ss <- paste0('K9.up. S1 vc S2 pvalue=', formatC(t.test.ss(s1.up, s2.up, sssize=500)), '\n',
                       'K9.down. S1 vc S2 pvalue=', formatC(t.test.ss(s1.do, s2.do, sssize=500)), '\n',
                       'K9.other. S1 vc S2 pvalue=', formatC(t.test.ss(s1.ot, s2.ot, sssize=500)))
    wilcox.txt <- paste0('K9.up. S1 vc S2 pvalue=', formatC(wilcox.test(s1.up, s2.up)$p.value), '\n',
                         'K9.down. S1 vc S2 pvalue=', formatC(wilcox.test(s1.do, s2.do)$p.value), '\n',
                         'K9.other. S1 vc S2 pvalue=', formatC(wilcox.test(s1.ot, s2.ot)$p.value))
    wilcox.txt.ss <- paste0('K9.up. S1 vc S2 pvalue=', formatC(wilcox.test.ss(s1.up, s2.up, sssize=500)), '\n',
                            'K9.down. S1 vc S2 pvalue=', formatC(wilcox.test.ss(s1.do, s2.do, sssize=500)), '\n',
                            'K9.other. S1 vc S2 pvalue=', formatC(wilcox.test.ss(s1.ot, s2.ot, sssize=500)))
    txt <- paste0('t.test\n', t.txt, '\n\nt.test with subsampling (500)\n', t.txt.ss, 
                  '\n\nwilcox.test\n', wilcox.txt, '\n\nwilcox.test with subsampling (500)\n', wilcox.txt.ss)
    txt
}
