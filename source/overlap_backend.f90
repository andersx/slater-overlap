module overlap_backend

    implicit none

    double precision, parameter :: zero = 0.0d0
    double precision, parameter :: half = 0.5d0
    double precision, parameter :: one = 1.0d0
    double precision, parameter :: two = 2.0d0

    ! Use this constant instead of global, because of compatibility with MNDO etc
    double precision, parameter :: AU_TO_EV = 27.21d0

    ! contains the factorials of i-1
    double precision, parameter :: fc(1:25) =&
        (/ 1.0d0,1.0d0, 2.0d0, 6.0d0, 24.0d0, &
        120.0d0, 720.0d0, 5040.0d0, 40320.0d0, 362880.0d0, &
        3628800.0d0, 39916800.0d0, 4.790016d+08, 6.2270208d+09, 8.71782912d+10, &
        1.307674368d+12, 2.092278989d+13, 3.55687428096d+14, 6.402373705728d+15, 1.21645100408832d+17, &
        2.43290200817664d+18, 5.109094217170944d+19, 1.12400072777760768d+21, 2.585201673888497664d+22, &
        6.2044840173323943936d+23 /)

    double precision, parameter :: logfc(1:25) = (/ 0.0d0, 0.0d0, 0.6931471805599d0,       &
        &  1.7917594692281d0,  3.1780538303479d0,  4.7874917427820d0, &
        &  6.5792512120101d0,  8.5251613610654d0, 10.6046029027453d0, &
        & 12.8018274800815d0, 15.1044125730755d0, 17.5023078458739d0, &
        & 19.9872144956619d0, 22.5521638531234d0, 25.1912211827387d0, &
        & 27.8992713838409d0, 30.6718601061763d0, 33.5050734501369d0, &
        & 36.3954452080331d0, 39.3398841871995d0, 42.3356164607535d0, &
        & 45.3801388984769d0, 48.4711813518352d0, 51.6066755677644d0, &
        & 54.7847293981123d0 /)

    !     define c coefficients for associate legendre polynomials.
    double precision, parameter::cc(1:21,1:3) = reshape ( (/   &
        8.0d0,   8.0d0,   4.0d0,  -4.0d0,  4.0d0,   &
        4.0d0, -12.0d0,  -6.0d0,  20.0d0,  5.0d0,   &
        3.0d0, -30.0d0, -10.0d0,  35.0d0,  7.0d0,   &
        15.0d0,   7.5d0, -70.0d0, -17.5d0, 63.0d0,  &
        10.5d0,                                     &
        0.0d0,   0.0d0,   0.0d0,  12.0d0,  0.0d0,   &
        0.0d0,  20.0d0,  30.0d0,   0.0d0,  0.0d0,   &
        -30.0d0,  70.0d0,  70.0d0,   0.0d0,  0.0d0, &
        -70.0d0, -105.d0, 210.0d0, 157.5d0,  0.0d0, &
        0.0d0,                                      &
        0.0d0,   0.0d0,   0.0d0,   0.0d0,  0.0d0,   &
        0.0d0,   0.0d0,   0.0d0,   0.0d0,  0.0d0,   &
        35.0d0,   0.0d0,   0.0d0,   0.0d0,  0.0d0,  &
        63.0d0, 157.5d0,   0.0d0,   0.0d0,  0.0d0,  &
        0.0d0/), (/ 21, 3 /) )

    double precision, save::A(15),B(15)

    ! Constants
    !     THE ARRAY B0(I) CONTAINS THE B INTEGRALS FOR ZERO ARGUMENT.
    double precision, parameter :: B0(1:15)= &
        &      (/ 2.0D0,0.0D0,0.666666666666667D0,0.0D0,0.4D0,0.0D0,       &
        &         0.285714285714286D0,0.0D0,0.222222222222222D0,0.0D0,     &
        &         0.181818181818182D0,0.0D0,0.153846153846154D0,0.0D0,     &
        &         0.133333333333333D0 /)


    !     define addresses for index pairs (00,10,20,30,40,50,60,70).
    integer,parameter:: iad(1:8)= (/ 1,2,4,7,11,16,22,29 /)

    !     define binomial coefficients (00,10,11,20,...,77).
    integer,parameter:: ibinom(1:36)= (/&
        &            1,1,1,1,2,1,1,3,3,1,1,4,6,4,1,1,5,10,10,5,1,          &
        &            1,6,15,20,15,6,1,1,7,21,35,35,21,7,1 /)
contains

function binomialcoefficient(m, n) result (biocoeff)

    integer, intent(in)::m,n
    double precision::biocoeff

    integer, parameter::size=30
    integer, save::bc(size,size)=0
    logical, save::initialized=.false.;

    integer::i,j,k

    if (.not.initialized) then
        do i=1,size
        bc(i,1)=one
        bc(i,2:size)=zero
        end do
        do i=2, size
        do j=2, i
        bc(i,j)=bc(i-1,j-1)+bc(i-1,j)
        end do
        end do

        initialized=.true.
    end if

    biocoeff=one*bc(m,n)

end function binomialcoefficient

function getslateroverlap(na, la, nb, lb, mm, zeta_a, zeta_b, rab) result (overlap)

    implicit none
    integer, intent(in)::na, la, nb, lb, mm
    double precision, intent(in)::zeta_a,zeta_b, rab
    double precision::overlap

    ! local variables
    integer, save::ntotal=-1
    double precision, save:: zeta_a_old=-1.0d99, zeta_b_old=-1.0d99
    double precision, save:: rab_old=-1.0d99
    double precision, parameter:: tolerance=1.d0-16

    logical::resetup

    if ((la.ge.na) .or. (lb.ge.nb)) then
        overlap=0.d0
        return
    end if
    resetup=.true.
    if (ntotal.eq. (na+nb)) then
        if (abs(rab_old/rab-1.d0)< tolerance) then
            if (abs(zeta_a_old/zeta_a-1.d0) < tolerance) then
                if (abs(zeta_b_old/zeta_b-1.d0) < tolerance) then
                    resetup=.false.
                end if
            end if
        end if
    end if

    if (resetup) then
        ntotal=(na+nb)
        rab_old=rab
        zeta_a_old=zeta_a
        zeta_b_old=zeta_b
        call setupslaterauxiliary(ntotal, zeta_a, zeta_b, rab)
    end if

    overlap=calculateoverlap (na,la,mm,nb,lb,zeta_a*rab,zeta_b*rab)

end function getslateroverlap


subroutine setupslaterauxiliary(n,sa,sb,rab)
    !     *
    !     calculation of auxiliary integrals for sto overlaps.
    !     *

    implicit double precision (a-h,o-z)
    implicit integer (i-n)

    !explicit type to satisfy pgi v8 compiler
    double precision betpow(17)


    ! *** initialization.
    alpha  = 0.5d0*rab*(sa+sb)
    beta   = 0.5d0*rab*(sa-sb)
    ! *** auxiliary a integrals for calculation of overlaps.
    c      = exp(-alpha)
    ralpha = 1.0d0/alpha
    a(1)   = c*ralpha
    do i=1,n
        a(i+1) = (a(i)*i+c)*ralpha
    end do
    ! *** auxiliary b integrals for calculation of overlaps.
    !     the code is valid only for n.le.14, i.e. for overlaps
    !     involving orbitals with main quantum numbers up to 7.
    !     branching depending on absolute value of the argument.
    absx   = abs(beta)
    !     zero argument.
    if(absx.lt.1.0d-06) then
        do i=1,n+1
            b(i) = b0(i)
        enddo
        return
    endif
    !     large argument.
    if((absx.gt.0.5d0 .and. n.le.5) .or.                              &
        &   (absx.gt.1.0d0 .and. n.le.7) .or.                              &
        &   (absx.gt.2.0d0 .and. n.le.10).or.                              &
        &    absx.gt.3.0d0) then
        expx   = exp(beta)
        expmx  = 1.0d0/expx
        rx     = 1.0d0/beta
        b(1)   = (expx-expmx)*rx
        do i=1,n
            expx   = -expx
            b(i+1)= (i*b(i)+expx-expmx)*rx
        enddo
        return
    endif
    !     small argument.
    if(absx.le.0.5d0) then
        last = 6
    else if(absx.le.1.0d0) then
        last = 7
    else if(absx.le.2.0d0) then
        last = 12
    else
        last = 15
    endif
    betpow(1) = 1.0d0
    do m=1,last
        betpow(m+1) = -beta*betpow(m)
    enddo
    do i=1,n+1
        y      = 0.0d0
        ma     = 1-mod(i,2)
        do m=ma,last,2
            y      = y+betpow(m+1)/(fc(m+1)*(m+i))
        enddo
        b(i)   = y*2.0d0
    enddo

    return
end subroutine setupslaterauxiliary

function calculateoverlap (na,la,mm,nb,lb,alpha,beta)
    !     *
    !     overlap integrals between slater type orbitals.
    !     *
    !     quantum numbers (na,la,mm) and (nb,lb,mm).
    !     na and nb must be positive and less than or equal to 7.
    !     la, lb and abs(mm) must be less than or equal to 5.
    !     further restrictions are la.le.na, lb.le.nb,
    !     mm.le.la, and mm.le.lb.
    !     *

    implicit double precision (a-h,o-z)
    implicit integer (i-n)


    ! *** initialization.
    m      = abs(mm)
    nab    = na+nb+1
    x      = 0.0d0
    ! *** find a and b integrals.
    !     p      = (alpha + beta)*0.5
    !     pt     = (alpha - beta)*0.5
    !     call aintgs(a,p,na+nb)
    !     call bintgs(b,pt,na+nb)
    ! *** section used for overlap integrals involving s functions.
    if((la.gt.0).or.(lb.gt.0)) go to 20
    iada   = iad(na+1)
    iadb   = iad(nb+1)
    do 10 i=0,na
    iba    = ibinom(iada+i)
    do 10 j=0,nb
    ibb    = iba*ibinom(iadb+j)
    if(mod(j,2).eq.1) ibb=-ibb
    ij     = i+j
    x      = x+ibb*a(nab-ij)*b(ij+1)
    10 continue
    ss     = x  * 0.5d0
    ss     = ss * sqrt( alpha**(2*na+1)*beta**(2*nb+1)/               &
        &                    (fc(2*na+1)*fc(2*nb+1)) )
    calculateoverlap = ss
    return
    ! *** section used for overlap integrals involving p functions.
    ! *** special case m=0, s-p(sigma), p(sigma)-s, p(sigma)-p(sigma).
    20 if(la.gt.1 .or. lb.gt.1) go to 320
    if(m.gt.0) go to 220
    iu     = mod(la,2)
    iv     = mod(lb,2)
    namu   = na-iu
    nbmv   = nb-iv
    iadna  = iad(namu+1)
    iadnb  = iad(nbmv+1)
    do 130 kc=0,iu
    ic     = nab-iu-iv+kc
    jc     = 1+kc
    do 130 kd=0,iv
    id     = ic+kd
    jd     = jc+kd
    do 130 ke=0,namu
    ibe    = ibinom(iadna+ke)
    ie     = id-ke
    je     = jd+ke
    do 130 kf=0,nbmv
    ibf    = ibe*ibinom(iadnb+kf)
    if(mod(kd+kf,2).eq.1) ibf=-ibf
    x      = x+ibf*a(ie-kf)*b(je+kf)
    130 continue
    ss     = x  * sqrt( (2*la+1)*(2*lb+1)*0.25d0 )
    !     compute overlap integral from reduced overlap integral.
    ss     = ss * sqrt( alpha**(2*na+1)*beta**(2*nb+1)/               &
        &                    (fc(2*na+1)*fc(2*nb+1)) )
    if(mod(lb,2).eq.1) ss=-ss
    calculateoverlap = ss
    return
    ! *** section used for overlap integrals involving p functions.
    ! *** special case la=lb=m=1, p(pi)-p(pi).
    220 iadna  = iad(na)
    iadnb  = iad(nb)
    do 230 ke=0,na-1
    ibe    = ibinom(iadna+ke)
    ie     = nab-ke
    je     = ke+1
    do 230 kf=0,nb-1
    ibf    = ibe*ibinom(iadnb+kf)
    if(mod(kf,2).eq.1) ibf=-ibf
    i      = ie-kf
    j      = je+kf
    x      = x+ibf*(a(i)*b(j)-a(i)*b(j+2)-a(i-2)*b(j)+a(i-2)*b(j+2))
    230 continue
    ss     = x  * 0.75d0
    !     compute overlap integral from reduced overlap integral.
    ss     = ss * sqrt( alpha**(2*na+1)*beta**(2*nb+1)/               &
        &                    (fc(2*na+1)*fc(2*nb+1)) )
    if(mod(lb+mm,2).eq.1) ss=-ss
    calculateoverlap = ss
    return
    ! *** section used for overlap integrals involving non-s functions.
    ! *** general case la.gt.1 or lb.gt.1, m.ge.0.
    320 lam    = la-m
    lbm    = lb-m
    iada   = iad(la+1)+m
    iadb   = iad(lb+1)+m
    iadm   = iad(m+1)
    iu1    = mod(lam,2)
    iv1    = mod(lbm,2)
    iuc    = 0
    do 340 iu=iu1,lam,2
    iuc    = iuc+1
    cu     = cc(iada,iuc)
    namu   = na-m-iu
    iadna  = iad(namu+1)
    iadu   = iad(iu+1)
    ivc    = 0
    do 340 iv=iv1,lbm,2
    ivc    = ivc+1
    nbmv   = nb-m-iv
    iadnb  = iad(nbmv+1)
    iadv   = iad(iv+1)
    sum    = 0.0d0
    do 330 kc=0,iu
    ibc    = ibinom(iadu+kc)
    ic     = nab-iu-iv+kc
    jc     = 1+kc
    do 330 kd=0,iv
    ibd    = ibc*ibinom(iadv+kd)
    id     = ic+kd
    jd     = jc+kd
    do 330 ke=0,namu
    ibe    = ibd*ibinom(iadna+ke)
    ie     = id-ke
    je     = jd+ke
    do 330 kf=0,nbmv
    ibf    = ibe*ibinom(iadnb+kf)
    iff    = ie-kf
    jff    = je+kf
    do 330 ka=0,m
    iba    = ibf*ibinom(iadm+ka)
    i      = iff-2*ka
    do 330 kb=0,m
    ibb    = iba*ibinom(iadm+kb)
    if(mod(ka+kb+kd+kf,2).eq.1) ibb=-ibb
    j      = jff+2*kb
    sum    = sum+ibb*a(i)*b(j)
    330 continue
    x      = x+sum*cu*cc(iadb,ivc)
    340 continue
    ss     = x*(fc(m+2)/8.0d0)**2* sqrt( (2*la+1)*fc(la-m+1)*         &
        &         (2*lb+1)*fc(lb-m+1)/(4.0d0*fc(la+m+1)*fc(lb+m+1)))
    !     compute overlap integral from reduced overlap integral.
    ss     = ss * sqrt( alpha**(2*na+1)*beta**(2*nb+1)/               &
        &                    (fc(2*na+1)*fc(2*nb+1)) )

    if(mod(lb+mm,2).eq.1) ss=-ss
    calculateoverlap = ss
    return
end function calculateoverlap

! Calculate the radial part of one-center two-electron integrals (slater-condon parameter).
function getslatercondonparameter(k,na,ea,nb,eb,nc,ec,nd,ed) result(slatercondon)


    ! Type of integral, can be equal to 0,1,2,3,4 in spd-basis
    integer, intent(in) :: k

    ! Principle quantum number of ao, electron 1
    integer, intent(in) :: na
    integer, intent(in) :: nb

    ! Principle quantum number of ao, electron 2
    integer, intent(in) :: nc
    integer, intent(in) :: nd

    ! Exponents of ao, electron 1
    double precision, intent(in) :: ea
    double precision, intent(in) :: eb

    ! Exponents of ao, electron 2
    double precision, intent(in) :: ec
    double precision, intent(in) :: ed

    ! The resulting Slater-Condon parameter
    double precision :: slatercondon

    ! Local variables
    integer :: nab, ncd, n, i, m, m1, m2
    double precision :: eab, ecd, e, c, s0, s1, s2, s3
    double precision :: aea, aeb, aec, aed, ae, a2, acd, aab

    aea    = log(ea)
    aeb    = log(eb)
    aec    = log(ec)
    aed    = log(ed)
    nab    = na+nb
    ncd    = nc+nd
    ecd    = ec+ed
    eab    = ea+eb
    e      = ecd+eab
    n      = nab+ncd
    ae     = log(e)
    a2     = log(two)
    acd    = log(ecd)
    aab    = log(eab)
    c      = exp(logfc(n)+na*aea+nb*aeb+nc*aec+nd*aed &
        +half*(aea+aeb+aec+aed)+a2*(n+2) &
        -half*(logfc(2*na+1)+logfc(2*nb+1)   &
        +logfc(2*nc+1)+logfc(2*nd+1))-ae*n)
    c      = c*au_to_ev
    s0     = one/e
    s1     = zero
    s2     = zero
    m      = ncd-k

    do i=1,m
    s0     = s0*e/ecd
    s1     = s1+s0*(binomialcoefficient(ncd-k,i)-binomialcoefficient(ncd+k+1,i))/binomialcoefficient(n,i)
    enddo

    m1     = m+1
    m2     = ncd+k+1

    do i=m1,m2
    s0     = s0*e/ecd
    s2     = s2+s0*binomialcoefficient(m2,i)/binomialcoefficient(n,i)
    enddo

    s3     = exp(ae*n-acd*m2-aab*(nab-k))/binomialcoefficient(n,m2)
    slatercondon = c*(s1-s2+s3)
    return
end function getslatercondonparameter

end module overlap_backend
