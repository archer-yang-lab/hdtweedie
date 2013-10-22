! --------------------------------------------------------------------------
! tweediegrpnet.f90: the IRLS-BMD algorithm for Tweedie model with grouped elastic net penalized learning.
! --------------------------------------------------------------------------
! 
! USAGE:
! 
! SUBROUTINE tweediegrpnet(alpha,rho,vt,bn,bs,ix,iy,gam,nobs,nvars,x,y,pf,dfmax,pmax,nlam,flmin,ulam,eps,maxit,nalam,b0,beta,idx,nbeta,alam,npass,jerr)
! 
! 
! INPUT ARGUMENTS:
!    alpha = group-elastic net quadratic mixing parameter, 0 < alpha <= 1
!    rho = power for variance-mean relationship
!    vt(nobs) = observation weight vector
!    bn = number of groups
!    bs(bn) = size of each group
!    ix(bn) = first index for each group
!    iy(bn) = last index for each group
!    gam(bn) = second derivatives
!    nobs = number of observations
!    nvars = number of predictor variables
!    x(nobs, nvars) = matrix of predictors, of dimension N * p; each row is an observation vector.
!    y(nobs) = response variable. This argument should be in {0, inf} 
!    pf(bn) = relative penalties for each group
!                pf(k) = 0 => kth group unpenalized
!    dfmax = limit the maximum number of variables in the model.
!            (one of the stopping criterion)
!    pmax = limit the maximum number of variables ever to be nonzero. 
!           For example once beta enters the model, no matter how many 
!           times it exits or re-enters model through the path, it will 
!           be counted only once. 
!    nlam = the number of lambda values
!    flmin = user control of lambda values (>=0)
!            flmin < 1.0 => minimum lambda = flmin*(largest lambda value)
!            flmin >= 1.0 => use supplied lambda values (see below)
!    ulam(nlam) = user supplied lambda values (ignored if flmin < 1.0)
!    eps = convergence threshold for coordinate majorization descent. 
!          Each inner coordinate majorization descent loop continues 
!          until the relative change in any coefficient is less than eps.
!    isd = standarization flag:
!          isd = 0 => regression on original predictor variables
!          isd = 1 => regression on standardized predictor variables
!          Note: output solutions always reference original
!                variables locations and scales.
!    maxit = maximum number of outer-loop iterations allowed at fixed lambda value. 
!            (suggested values, maxit = 100000)
! 
! OUTPUT:
! 
!    nalam = actual number of lambda values (solutions)
!    b0(nvars) = intercept values for each solution
!    beta(nvars, nlam) = compressed coefficient values for each solution
!    idx(pmax) = pointers to compressed coefficients
!    nbeta(nlam) = number of compressed coefficients for each solution
!    alam(nlam) = lambda values corresponding to each solution
!    npass = actual number of passes over the data for all lambda values
!    jerr = error flag:
!           jerr  = 0 => no error
!           jerr > 0 => fatal error - no output returned
!                    jerr < 7777 => memory allocation error
!                    jerr = 10000 => maxval(pf) <= 0.0
!                    jerr =15000 => minval(vt) < 0
!                    jerr = 20000+ => predictor matrix standardization fails
!                    jerr = 30000+ => SVD fails
!           jerr < 0 => non fatal error - partial output:
!                    Solutions for larger lambdas (1:(k-1)) returned.
!                    jerr = -k => convergence for kth lambda value not reached
!                           after maxit (see above) iterations.
!                    jerr = -10000-k => number of non zero coefficients along path
!                           exceeds pmax (see above) at kth lambda value.
! 
! LICENSE: GNU GPL (version 2 or later)
! 
! AUTHORS: Wei Qian, Yi Yang and Hui Zou
!     
!    School of Statistics, University of Minnesota.
! 
! REFERENCES:
!    Qian, W., Yang, Y., Yang, Y. and Zou, H. (2013). 
!    Tweedie's Compound Poisson Model With Grouped Elastic Net
!    submitted to Journal of Computational and Graphical Statistics.


! --------------------------------------------------
SUBROUTINE tweediegrpnet(alpha,rho,vt,bn,bs,ix,iy,gam,nobs,nvars,x,y,pf,&
dfmax,pmax,nlam,flmin,ulam,eps,isd,maxit,nalam,b0,beta,&
idx,nbeta,alam,npass,jerr)
! --------------------------------------------------
    IMPLICIT NONE
    ! - - - arg types - - -
    DOUBLE PRECISION, PARAMETER :: big=9.9E30 
    DOUBLE PRECISION, PARAMETER :: small=1.0D-16 
    DOUBLE PRECISION, PARAMETER :: mfl = 1.0E-6 
    INTEGER, PARAMETER :: mnlam = 6 
    INTEGER:: mnl
    INTEGER:: bn
    INTEGER::bs(bn)
    INTEGER::ix(bn)
    INTEGER::iy(bn)
    INTEGER:: nobs
    INTEGER::nvars
    INTEGER::dfmax
    INTEGER::pmax
    INTEGER::nlam
    INTEGER::isd
    INTEGER::maxit
    INTEGER::nalam
    INTEGER::npass
    INTEGER::jerr
    INTEGER:: idx(pmax)  
    INTEGER::nbeta(nlam)
    DOUBLE PRECISION :: alpha
    DOUBLE PRECISION :: rho
    DOUBLE PRECISION :: vt(nobs)
    DOUBLE PRECISION::gam(bn) 
    DOUBLE PRECISION:: x(nobs,nvars)
    DOUBLE PRECISION::y(nobs)
    DOUBLE PRECISION::pf(bn)
    DOUBLE PRECISION:: flmin
    DOUBLE PRECISION::ulam(nlam)
    DOUBLE PRECISION::eps
    DOUBLE PRECISION:: b0(nlam)
    DOUBLE PRECISION::beta(nvars,nlam)
    DOUBLE PRECISION::alam(nlam)
    ! - - - local declarations - - -
    DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: xmean
    DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: xnorm
    INTEGER :: ju=0 
    DOUBLE PRECISION:: vtt(nobs) 
    DOUBLE PRECISION:: yt(nobs) 
    DOUBLE PRECISION:: max_gam 
    CHARACTER:: jobu = 'N'
    CHARACTER:: jobvt = 'N'
    DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: hj 
    INTEGER:: lda 
    DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: uj 
    INTEGER:: ldu
    DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: vj 
    INTEGER:: ldvt
    DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: work 
    INTEGER:: lwork 
    INTEGER:: info  
    DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: eig
    DOUBLE PRECISION::d 
    DOUBLE PRECISION::t 
    DOUBLE PRECISION::dif  
    DOUBLE PRECISION::unorm 
    DOUBLE PRECISION::al 
    DOUBLE PRECISION::alf 
    DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: b
    DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::oldbeta
    DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::r 
    DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: r1 
    DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: r2 
    DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::oldb
    DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::u 
    DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::v 
    DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::dd 
    INTEGER, DIMENSION (:), ALLOCATABLE :: oidx 
    INTEGER:: g
    INTEGER::j
    INTEGER::l
    INTEGER::ierr
    INTEGER::ni 
    INTEGER::me 
    INTEGER::start
    INTEGER::end
    INTEGER::ii 
    INTEGER::jj  
    ! - - - local declarations - - -
    DOUBLE PRECISION:: tlam 
    INTEGER:: jx 
    INTEGER:: jxx(bn) 
    DOUBLE PRECISION:: ga(bn) 
    DOUBLE PRECISION:: vl(nvars) 
    DOUBLE PRECISION:: al0 

! - - - begin - - -
! - - - allocate variables - - -
    ALLOCATE(b(0:nvars),STAT=jerr)
    ALLOCATE(oldbeta(0:nvars),STAT=ierr)
    jerr=jerr+ierr
    ALLOCATE(r(1:nobs),STAT=ierr)
    jerr=jerr+ierr
    ALLOCATE(r1(1:nobs),STAT=ierr)
    jerr=jerr+ierr
    ALLOCATE(r2(1:nobs),STAT=ierr)
    jerr=jerr+ierr
    ALLOCATE(oidx(1:bn),STAT=ierr)
    jerr=jerr+ierr
    ALLOCATE(xmean(1:nvars), STAT=ierr)
    jerr=jerr+ierr
    ALLOCATE(xnorm(1:nvars),STAT=ierr)
    jerr=jerr+ierr
    IF(jerr/=0) RETURN
! - - - checking pf - - -
    IF(maxval(pf) <= 0.0D0) THEN
        jerr=10000
        RETURN
    ENDIF
    pf=max(0.0D0,pf)

    IF (minval(vt) < 0) THEN
        jerr = 15000
        RETURN
    ENDIF
    vt = vt / sum(vt) 
    call standard(nobs, nvars, x, isd, xmean, xnorm, ju)
    IF(ju > 0) THEN
        jerr = 20000+ju
        RETURN
    ENDIF

! - - - some initial setup - - -
    jxx = 0
    al = 0.0D0
    mnl = Min (mnlam, nlam)
    r = 0.0D0
    b=0.0D0
    oldbeta=0.0D0
    idx=0
    oidx=0
    npass=0
    ni=npass
! --------- lambda loop ----------------------------
    IF(flmin < 1.0D0) THEN
        flmin = Max (mfl, flmin)
        alf=flmin**(1.0D0/(nlam-1.0D0))
    ENDIF
    r1 = vt*y*exp(-(rho-1.0D0)*r)
    r2 = vt*exp((2.0D0-rho)*r)
    vl = matmul(r1-r2, x)
    DO g = 1,bn
            ALLOCATE(u(bs(g)),STAT=ierr)
            IF(ierr/=0) RETURN
            u = vl(ix(g):iy(g))
            ga(g) = sqrt(dot_product(u,u))
            DEALLOCATE(u)
    END DO
    DO l=1, nlam
        al0 = al 
        IF(flmin>=1.0D0) THEN
            al=ulam(l)
        ELSE
            IF(l > 2) THEN
                al=al*alf
            ELSE IF(l==1) THEN
                al=big
            ELSE IF(l==2) THEN 
                al0 = 0.0D0
                r1 = vt*y*exp(-(rho-1.0D0)*r)
                r2 = vt*exp((2.0D0-rho)*r)
                vl = matmul(r1-r2, x)
                DO g = 1,bn
                    ALLOCATE(u(bs(g)),STAT=ierr)
                    IF(ierr/=0) RETURN
                    u = vl(ix(g):iy(g))
                    ga(g) = sqrt(dot_product(u,u))
                    IF(pf(g)>1.0D-16) THEN
                        al0 = max(al0, ga(g)/pf(g)/alpha)
                    ENDIF
                    DEALLOCATE(u)
                END DO
                al = al0 * alf 
            ENDIF
        ENDIF
        tlam = (2.0*al-al0)
        DO g = 1, bn
            IF(jxx(g) == 1) CYCLE
            IF(ga(g) > pf(g) * tlam) jxx(g) = 1
        ENDDO
! --------- outer loop ----------------------------
        DO
            oldbeta(0)=b(0)
            IF(ni>0) THEN
                DO j=1,ni
                    g=idx(j)
                    oldbeta(ix(g):iy(g))=b(ix(g):iy(g))
                ENDDO
            ENDIF
            r1 = vt*y*exp(-(rho-1.0D0)*r)
            r2 = vt*exp((2.0D0-rho)*r)
            vtt = (rho-1.0D0)*r1 + (2.0D0-rho)*r2
            yt = r + (r1-r2)/vtt
            vtt = vtt
            DO g=1,bn
                IF(jxx(g) == 0) CYCLE
                start=ix(g)
                end=iy(g)
                ALLOCATE(hj(1:bs(g),1:bs(g)),STAT=ierr)
                jerr=jerr+ierr
                lda = bs(g)
                ALLOCATE(uj(1:bs(g),1:bs(g)),STAT=ierr)
                jerr=jerr+ierr
                ldu = bs(g)
                ALLOCATE(vj(1:bs(g),1:bs(g)),STAT=ierr)
                jerr=jerr+ierr
                ldvt = bs(g)
                lwork = 5*bs(g)
                ALLOCATE(work(lwork),STAT=ierr)
                jerr=jerr+ierr
                ALLOCATE(eig(bs(g)),STAT=ierr)
                jerr=jerr+ierr
                IF(jerr/=0) RETURN
                hj = 0.0D0
                DO ii=1,bs(g)
                    DO jj=1,ii
                        hj(ii,jj) = sum(vtt*x(:,start+ii-1)*x(:,start+jj-1))
                        hj(jj,ii) = hj(ii,jj)
                    ENDDO
                ENDDO
                call dgesvd(jobu,jobvt,bs(g),bs(g),hj,lda,eig,uj,ldu,vj,ldvt,work,lwork,info)
                IF (info /=0) THEN
                    jerr = 30000+info
                    RETURN
                ENDIF
                gam(g) = eig(1)
                DEALLOCATE(hj,uj,vj,work,eig)
            ENDDO

! --middle loop-------------------------------------
            DO
                npass=npass+1
                dif=0.0D0
                DO g=1,bn
                    IF(jxx(g) == 0) CYCLE
                    start=ix(g)
                    end=iy(g)
                    ALLOCATE(u(bs(g)),STAT=ierr)
                    jerr=jerr+ierr
                    ALLOCATE(dd(bs(g)),STAT=ierr)
                    jerr=jerr+ierr
                    ALLOCATE(oldb(bs(g)),STAT=ierr)
                    jerr=jerr+ierr
                    IF(jerr/=0) RETURN
                    oldb=b(start:end)
                    u=matmul(vtt*(yt-r),x(:,start:end))
                    u=gam(g)*b(start:end)+u
                    unorm=sqrt(dot_product(u,u))
                    t=unorm-pf(g)*al*alpha
                    IF(t>0.0D0) THEN
                        b(start:end)=u*t/((gam(g)+al*(1.0D0-alpha))*unorm)
                    ELSE
                        b(start:end)=0.0D0
                    ENDIF
                    dd=b(start:end)-oldb
                    IF(any(abs(dd)>small)) THEN
                        dif=max(dif,dot_product(dd,dd))
                        r=r+matmul(x(:,start:end),dd)
                        IF(oidx(g)==0) THEN
                            ni=ni+1
                            IF(ni>pmax) EXIT
                            oidx(g)=ni
                            idx(ni)=g
                        ENDIF
                    ENDIF
                    DEALLOCATE(u,dd,oldb)
                ENDDO
                d = sum(vtt*(yt-r))/sum(vtt)
                IF(abs(d) > small) THEN
                    b(0)=b(0)+d
                    r=r+d
                    dif=max(dif,d**2)
                ENDIF
                IF (ni > pmax) EXIT
                IF (dif < eps) EXIT 
                IF(npass > maxit) THEN
                    jerr=-l
                    RETURN
               ENDIF
                
! --inner loop----------------------
                DO
                    npass=npass+1
                    dif=0.0D0
                    DO j=1,ni
                        g=idx(j)
                        start=ix(g)
                        end=iy(g)
                        ALLOCATE(u(bs(g)),STAT=ierr)
                        jerr=jerr+ierr
                        ALLOCATE(dd(bs(g)),STAT=ierr)
                        jerr=jerr+ierr
                        ALLOCATE(oldb(bs(g)),STAT=ierr)
                        jerr=jerr+ierr
                        IF(jerr/=0) RETURN
                        oldb=b(start:end)
                        u=matmul(vtt*(yt-r),x(:,start:end))
                        u=gam(g)*b(start:end)+u
                        unorm=sqrt(dot_product(u,u))
                        t=unorm-pf(g)*al*alpha
                        IF(t>0.0D0) THEN
                            b(start:end)=u*t/((gam(g)+al*(1.0D0-alpha))*unorm)
                        ELSE
                            b(start:end)=0.0D0
                        ENDIF
                        dd=b(start:end)-oldb
                        IF(any(abs(dd)>small)) THEN
                            dif=max(dif,dot_product(dd,dd))
                            r=r+matmul(x(:,start:end),dd)
                        ENDIF
                        DEALLOCATE(u,dd,oldb)
                    ENDDO
                    d = sum(vtt*(yt-r))/sum(vtt)
                    IF(abs(d) > small) THEN
                        b(0)=b(0)+d
                        r=r+d
                        dif=max(dif, d**2)
                    ENDIF
                    IF(dif<eps) EXIT
                    IF(npass > maxit) THEN
                        jerr=-l
                        RETURN
                    ENDIF
                ENDDO
            ENDDO
            IF(ni>pmax) EXIT
!--- final check ------------------------
            jx = 0
            max_gam = maxval(gam)
            IF(any((b-oldbeta)**2 >= eps)) jx = 1
            IF (jx /= 0) CYCLE
            r1 = vt*y*exp(-(rho-1.0D0)*r)
            r2 = vt*exp((2.0D0-rho)*r)
            vl = matmul(r1-r2, x)
            DO g = 1, bn
                IF(jxx(g) == 1) CYCLE
                ALLOCATE(u(bs(g)),STAT=ierr)
                IF(ierr/=0) RETURN
                u = vl(ix(g):iy(g))
                ga(g) = sqrt(dot_product(u,u))
                IF(ga(g) > al*pf(g)*alpha) THEN 
                    jxx(g) = 1
                    jx = 1
                ENDIF
                DEALLOCATE(u)
            ENDDO
            IF(jx == 1) CYCLE
            EXIT
        ENDDO
!---------- final update variable and save results------------
        IF(ni>pmax) THEN
            jerr=-10000-l
            EXIT
        ENDIF
        IF(ni>0) THEN
            DO j=1,ni
                g=idx(j)
                beta(ix(g):iy(g),l)=b(ix(g):iy(g))
            ENDDO
        ENDIF
        nbeta(l)=ni
        b0(l)=b(0)
        alam(l)=al
        nalam=l
        IF (l < mnl) CYCLE 
        me=0
        DO j=1,ni
            g=idx(j)
            IF(any(abs(beta(ix(g):iy(g),l))>small)) me=me+1
        ENDDO
        IF(me>dfmax) EXIT
    ENDDO
    
    DO l=1,nalam
        IF (isd == 1) THEN
            DO j=1,nvars
                beta(j,l) = beta(j,l) / xnorm(j)
            ENDDO
            b0(l) = b0(l) - dot_product(beta(:,l),xmean)
        ENDIF
    ENDDO
    DEALLOCATE(b,oldbeta,r,r1,r2,oidx,xmean,xnorm)
    RETURN
END SUBROUTINE tweediegrpnet


! --------------------------------------------------
SUBROUTINE standard(nobs, nvars, x, isd, xmean, xnorm, ju)
! --------------------------------------------------
    IMPLICIT NONE
    ! - - - arg types - - -
    INTEGER :: nobs
    INTEGER :: nvars
    INTEGER :: isd
    DOUBLE PRECISION :: x(nobs,nvars)
    DOUBLE PRECISION :: xmean(nvars)
    DOUBLE PRECISION :: xnorm(nvars)
    INTEGER :: ju
    INTEGER :: j
    DOUBLE PRECISION :: maj

    DO j=1,nvars
        IF(isd==1) THEN
            xmean(j) = sum(x(:,j))/nobs
            x(:,j) = x(:,j)-xmean(j)
            maj = dot_product(x(:,j),x(:,j))/nobs
            IF(maj < 1.0D-16) THEN 
                ju = j
                RETURN
            ENDIF
            xnorm(j) = sqrt(maj)
            x(:,j) = x(:,j)/xnorm(j)
        ENDIF
    ENDDO
END SUBROUTINE standard


