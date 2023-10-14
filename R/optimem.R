ss_kb <- function(M, y, eta, n.k, n.b, n.r, min.nonDA, ss.br.mn.kb,
                  ss.br.ul.kb, ss.br.ll.kb, ss.ul.obj, early.exit, kmin.prop){
        n.sample <- nrow(M); n.taxa <- ncol(M)
        gid <- unique(y); n.gid <- length(gid)
        g.mtx <- matrix(0, nrow=n.sample, ncol=n.gid)
        for(i in 1:n.gid) g.mtx[,i] <- ifelse(y==gid[i], 1, 0)
        j.G <- rep(1, n.gid); n.g <- colSums(g.mtx)
        if(is.null(n.k)) n.k <- ceiling(log(min.nonDA+2/n.taxa)/log(1-eta))
        n.cores <- detectCores()
        optimem.cl <- makeCluster(n.cores-1)
        registerDoParallel(cl=optimem.cl)
        taxa.in.M <- colnames(M); f.rslt <- NULL
        n.sel <- n.taxa; cut.off.n.taxa <- 0.2*n.taxa
        cut.off.n.taxa2 <- min(log(n.taxa, 2)*0.1/eta, cut.off.n.taxa)
        tmp.b.hat <- tmp.b.hat.exit <- NULL
        k.min.taxa.cutoff <- round(kmin.prop*n.taxa)
        if(eta*n.taxa < 1){
                for(k in 1:n.k){
                        rslt.br <- NULL
                        rslt.br <- foreach(b=1:n.sel, .combine="rbind") %dopar% {
                                sel.taxa <- taxa.in.M[-b]
                                M.tld <- M[,sel.taxa]
                                n.sel.r <- floor(n.sel/2)
                                x.g <- matrix(0, nrow=n.r, ncol=n.gid)
                                n.excl.r <- 0
                                for(r in 1:n.r){
                                        num.taxa.r <- sample(sel.taxa, n.sel.r)
                                        den.taxa.r <- setdiff(sel.taxa, num.taxa.r)
                                        x.r <- log(rowSums(M.tld[,num.taxa.r,drop=FALSE])/
                                                           rowSums(M.tld[,den.taxa.r,drop=FALSE]))
                                        esi <- which(x.r %in% c(-Inf, Inf, NaN))
                                        if(length(esi)>0){
                                                n.g.xr <- colSums(g.mtx[-esi,])
                                                if(min(n.g.xr) > 0){
                                                        x.g[r,] <- crossprod(x.r[-esi], g.mtx[-esi,])/n.g.xr
                                                } else{
                                                        x.g[r,] <- 0
                                                        n.excl.r <- n.excl.r + 1
                                                }
                                        } else{
                                                x.g[r,] <- crossprod(x.r, g.mtx)/n.g
                                        }
                                }
                                H.kb <- crossprod(x.g)
                                D.kb <- (tcrossprod(diag(H.kb), j.G) + tcrossprod(j.G, diag(H.kb)) - 2*H.kb)/(n.r-n.excl.r)
                                rslt.b.ss <- sum(D.kb)
                                rslt.taxa <- sel.taxa
                                return(c(rslt.b.ss, rslt.taxa))
                        }
                        kb.ss <- as.numeric(rslt.br[,1])
                        min.b.ss.indx <- which.min(kb.ss)
                        taxa.in.M <- rslt.br[min.b.ss.indx,-1]
                        f.rslt[[k]] <- list(b.hat=min(kb.ss), b.hat.taxa=taxa.in.M, is.all.taxa.nonDA=FALSE)
                        if(early.exit==TRUE){
                                tmp.b.hat <- max(ss.ul.obj, tmp.b.hat, f.rslt[[k]]$b.hat)
                                tmp.b.hat.exit <- min(tmp.b.hat.exit, f.rslt[[k]]$b.hat)
                                if(k > cut.off.n.taxa & tmp.b.hat < ss.br.ul.kb) {
                                        f.rslt[[k]]$is.all.taxa.nonDA <- TRUE
                                        break
                                }
                                if(k > cut.off.n.taxa & f.rslt[[k]]$b.hat < ss.br.ll.kb*0.9) break
                                if(tmp.b.hat.exit < ss.br.mn.kb & f.rslt[[k]]$b.hat > ss.br.ul.kb) break
                        }
                        n.taxa <- length(taxa.in.M)
                        n.sel <- n.sel-1
                        if(n.sel < k.min.taxa.cutoff) break
                }
        } else{
                for(k in 1:n.k){
                        rslt.br <- NULL
                        rslt.br <- foreach(b=1:n.b, .combine="rbind") %dopar% {
                                sel.taxa <- sample(taxa.in.M, n.sel)
                                M.tld <- M[,sel.taxa]
                                n.sel.r <- floor(n.sel/2)
                                x.g <- matrix(0, nrow=n.r, ncol=n.gid)
                                n.excl.r <- 0
                                for(r in 1:n.r){
                                        num.taxa.r <- sample(sel.taxa, n.sel.r)
                                        den.taxa.r <- setdiff(sel.taxa, num.taxa.r)
                                        x.r <- log(rowSums(M.tld[,num.taxa.r,drop=FALSE])/
                                                           rowSums(M.tld[,den.taxa.r,drop=FALSE]))
                                        esi <- which(x.r %in% c(-Inf, Inf, NaN))
                                        if(length(esi)>0){
                                                n.g.xr <- colSums(g.mtx[-esi,])
                                                if(min(n.g.xr) > 0){
                                                        x.g[r,] <- crossprod(x.r[-esi], g.mtx[-esi,])/n.g.xr
                                                } else{
                                                        x.g[r,] <- 0
                                                        n.excl.r <- n.excl.r + 1
                                                }
                                        } else{
                                                x.g[r,] <- crossprod(x.r, g.mtx)/n.g
                                        }
                                }
                                H.kb <- crossprod(x.g)
                                D.kb <- (tcrossprod(diag(H.kb), j.G) + tcrossprod(j.G, diag(H.kb)) - 2*H.kb)/(n.r-n.excl.r)
                                rslt.b.ss <- sum(D.kb)
                                rslt.taxa <- sel.taxa
                                return(c(rslt.b.ss, rslt.taxa))
                        }
                        kb.ss <- as.numeric(rslt.br[,1])
                        min.b.ss.indx <- which.min(kb.ss)
                        taxa.in.M <- rslt.br[min.b.ss.indx,-1]
                        f.rslt[[k]] <- list(b.hat=min(kb.ss), b.hat.taxa=taxa.in.M, is.all.taxa.nonDA=FALSE)
                        if(early.exit==TRUE){
                                tmp.b.hat <- max(ss.ul.obj, tmp.b.hat, f.rslt[[k]]$b.hat)
                                tmp.b.hat.exit <- min(tmp.b.hat.exit, f.rslt[[k]]$b.hat)
                                if(k > cut.off.n.taxa2 & tmp.b.hat < ss.br.ul.kb) {
                                        f.rslt[[k]]$is.all.taxa.nonDA <- TRUE
                                        break
                                }
                                if(k > cut.off.n.taxa2 & f.rslt[[k]]$b.hat < ss.br.ll.kb*0.95) break
                                if(tmp.b.hat.exit < ss.br.mn.kb & f.rslt[[k]]$b.hat > ss.br.ul.kb) break
                        }
                        n.taxa <- length(taxa.in.M)
                        n.sel <- floor((1-eta)*n.taxa)
                        if(n.sel < k.min.taxa.cutoff) break
                }
        }
        stopCluster(cl=optimem.cl)
        return(f.rslt)
}

ss_br <- function(M, y, n.perm, n.b, n.r){
        n.sample <- nrow(M); n.taxa <- ncol(M)
        gid <- unique(y); n.gid <- length(gid)
        g.mtx <- matrix(0, nrow=n.sample, ncol=n.gid)
        for(i in 1:n.gid) g.mtx[,i] <- ifelse(y==gid[i], 1, 0)
        j.G <- rep(1, n.gid); n.g <- colSums(g.mtx)
        n.cores <- detectCores()
        optimem.cl <- makeCluster(n.cores-1)
        registerDoParallel(cl=optimem.cl)
        taxa.in.M <- colnames(M); f.rslt0 <- NULL
        f.rslt.obj <- foreach(b=1:n.b, .combine="c") %dopar% {
                M.tld <- M
                n.taxa.r <- floor(n.taxa/2)
                x.r <- matrix(0, nrow=n.sample, ncol=n.r)
                for(r in 1:n.r){
                        num.taxa.r <- sample(taxa.in.M, n.taxa.r)
                        den.taxa.r <- setdiff(taxa.in.M, num.taxa.r)
                        x.r[,r] <- log(rowSums(M.tld[,num.taxa.r])/rowSums(M.tld[,den.taxa.r]))
                        n.non.num.indx <- which(x.r[,r] %in% c(-Inf, Inf, NaN))
                        x.r[n.non.num.indx, r] <- 0
                        n.non.num.r <- length(n.non.num.indx)
                }
                x.g <- sweep(crossprod(x.r, g.mtx), 2, n.g, "/")
                H.kb <- crossprod(x.g)
                D.kb <- (tcrossprod(diag(H.kb), j.G) + tcrossprod(j.G, diag(H.kb)) - 2*H.kb)/n.r
                return(sum(D.kb))
        }
        for(k in 1:n.perm){
                g.mtx <- matrix(0, nrow=n.sample, ncol=n.gid)
                perm.y <- sample(y)
                for(i in 1:n.gid) g.mtx[,i] <- ifelse(perm.y==gid[i], 1, 0)
                rslt.br <- NULL
                rslt.br <- foreach(b=1:n.b, .combine="c") %dopar% {
                        M.tld <- M
                        n.taxa.r <- floor(n.taxa/2)
                        x.r <- matrix(0, nrow=n.sample, ncol=n.r)
                        for(r in 1:n.r){
                                num.taxa.r <- sample(taxa.in.M, n.taxa.r)
                                den.taxa.r <- setdiff(taxa.in.M, num.taxa.r)
                                x.r[,r] <- log(rowSums(M.tld[,num.taxa.r])/rowSums(M.tld[,den.taxa.r]))
                                n.non.num.indx <- which(x.r[,r] %in% c(-Inf, Inf, NaN))
                                x.r[n.non.num.indx, r] <- 0
                                n.non.num.r <- length(n.non.num.indx)
                        }
                        x.g <- sweep(crossprod(x.r, g.mtx), 2, n.g, "/")
                        H.kb <- crossprod(x.g)
                        D.kb <- (tcrossprod(diag(H.kb), j.G) + tcrossprod(j.G, diag(H.kb)) - 2*H.kb)/n.r
                        return(sum(D.kb))
                }
                f.rslt0 <- c(f.rslt0, rslt.br)
        }
        stopCluster(cl=optimem.cl)
        return(list(f.rslt.obj=f.rslt.obj, f.rslt.null=f.rslt0))
}

optimem <- function(M, y, eta=0.1, alpha=0.05, n.k=NULL, n.b=200, n.r=300,
                    min.nonDA=0.1, n.perm=20, k.sel.plot=TRUE, early.exit=TRUE, k.min.prop=0.1){
        grp.id <- unique(y); n.grp <- length(grp.id)
        n.obs.grp <- table(y); excl.taxa <- NULL
        for(i in 1:n.grp){
                excl.taxa <- which(apply(M[y==grp.id[i],], 2, function(x) sum(x>0)) < n.obs.grp[i]*0.1)
        }
        excl.taxa <- sort(unique(excl.taxa))
        if(length(excl.taxa)>0){
                M.ex <- M[,-excl.taxa]
        } else{
                M.ex <- M
        }
        ref.rslt <- ss_br(M.ex, y, n.perm=n.perm, n.b=200, n.r=100)
        mean.kb <- mean(ref.rslt$f.rslt.null)
        ul.kb <- quantile(ref.rslt$f.rslt.null, 1-alpha/2)
        ll.kb <- quantile(ref.rslt$f.rslt.null, alpha/2)
        ss.obj <- quantile(ref.rslt$f.rslt.obj, 0.6)
        tmp.rslt <- ss_kb(M.ex, y, eta=eta, n.k=n.k, n.b=n.b, n.r=n.r, min.nonDA=min.nonDA,
                          ss.br.mn.kb=mean.kb, ss.br.ul.kb=ul.kb, ss.br.ll.kb=ll.kb, ss.ul.obj=ss.obj,
                          early.exit=early.exit, kmin.prop=k.min.prop)
        b.hat <- sapply(tmp.rslt,"[[",1)
        opt.min.kb <- which(b.hat == min(b.hat))
        n.rs <- length(tmp.rslt)
        s <- mss <- NULL
        k.sel.dat <- data.frame(s=0:n.rs, mss=c(ss.obj, b.hat))
        k.sel.plt <- ggplot(k.sel.dat, aes(x=s, y=mss)) +
                geom_point(color="gray40") +
                theme_bw() +
                theme(panel.grid = element_blank()) +
                labs(x="Removal Step (k)", y="MSS") +
                geom_abline(slope=0, intercept=mean.kb, color="gray60", linetype="twodash") +
                geom_abline(slope=0, intercept=ul.kb, color="gray70", linetype="dotted") +
                geom_abline(slope=0, intercept=ll.kb, color="gray70", linetype="dotted") +
                geom_point(aes(x=opt.min.kb, y=b.hat[opt.min.kb]), color="coral4")
        if(k.sel.plot==TRUE) print(k.sel.plt)
        if(tmp.rslt[[n.rs]]$is.all.taxa.nonDA==TRUE){
                print("All taxa are non-DA or noise levels are larger than signals!!!")
                return(NULL)
        } else{
                nonDA.min <- matrix(gtools::mixedsort(tmp.rslt[[opt.min.kb]]$b.hat.taxa), ncol=1)
                colnames(nonDA.min) <- "Remaining Taxa at min(MSS)"; rownames(nonDA.min) <- 1:length(nonDA.min)
                optimem.rslt <- list(nonDAtaxa_min=nonDA.min, mss_taxa=tmp.rslt, min_MSS=opt.min.kb,
                                     mean_null=mean.kb, upper_limit_null=ul.kb, lower_limit_null=ll.kb)
                class(optimem.rslt) <- "optimem"
                return(optimem.rslt)
        }
}

print.optimem <- function(x, ...){
        cat("Remaining Taxa at min(MSS): ", x$nonDAtaxa_min)
}

plot.optimem <- function(x, ...){
        b.hat <- sapply(x$mss_taxa, "[[", 1)
        opt.min.kb <- which(b.hat == min(b.hat))
        mean.kb <- x$mean_null
        ul.kb <- x$upper_limit_null
        ll.kb <- x$lower_limit_null
        n.rs <- length(b.hat)
        s <- mss <- NULL
        k.sel.dat <- data.frame(s=1:n.rs, mss=b.hat)
        ggplot(k.sel.dat, aes(x=s, y=mss)) +
                geom_point(color="gray40") +
                theme_bw() +
                theme(panel.grid = element_blank()) +
                labs(x="Removal Step (k)", y="MSS") +
                geom_abline(slope=0, intercept=mean.kb, color="gray60", linetype="twodash") +
                geom_abline(slope=0, intercept=ul.kb, color="gray70", linetype="dotted") +
                geom_abline(slope=0, intercept=ll.kb, color="gray70", linetype="dotted") +
                geom_point(aes(x=opt.min.kb, y=b.hat[opt.min.kb]), color="coral4")
}
