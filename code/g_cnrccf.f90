subroutine cnrccf(zdat,kzdt,mzdt,nzdt,flow,fhgh,finc,nset,nhgh, &
    ndsh,smooth,lines,icolfg)
    !
    !  This is a version of the old conrec which does color-filled contours.  this
    !  is a modified cpcnrc which was written by dave kennison.  it requires
    !  the conpack library.   the first ten arguments are identical to the
    !  old conrec call.
    !
    !  smooth - zero for no smoothing; small values (e.g. .001) yield
    !     approximately cubic splines; large values (e.g. 50.) yield nearly
    !     polygonal curves; suggested starting value is 2.5; if negative,
    !     smoothing is done before coordinate mapping
    !
    !  lines - if nonzero, contour lines are drawn in the foreground color
    !
    !  icolfg - if zero, a default color table is used; if nonzero, user
    !     must set color indices 2 to ncl+2, where ncl = number of contour
    !     levels.  the number of colors required is always one more than the
    !     number of contour levels.  color index numbering for associated
    !     contour fill colors must begin with 2 since color indices 0 and 1
    !     are reserved for background and foreground, respectively.  these
    !     may be set by the user regardless of the value of icolfg.
    !     choosing the default color table will reset indices 2 to maxcol.
    !     a subset of these will actually be used.  see routine dfclrs for
    !     the value of maxcol.
    !
    dimension zdat(kzdt,*)
    !
    ! define some needed dimensions.
    !
    parameter (lama=100000,lrwk=2000,liwk=1000,ncra=2000,ngrps=10, &
    locv=10)
    !
    ! define required workspace arrays.
    !
    dimension rwrk(lrwk),iwrk(liwk),iama(lama),xcra(ncra),ycra(ncra), &
    iaia(ngrps),igia(ngrps),rec(4)
    !
    ! logical to fix spaghetti code - mat
    !
    logical lexit
    !
    ! define a character variable to use for point-value labelling.
    !
    character*(locv) croz
    character*100 ilts
    common/mapawm/klat,icmp,icon
    !
    ! declare the contour-line drawing routine.
    !
    external colram
        lexit = .false.
    !
    !  ..set gks internal parameters
    !
    call gqclip(ier,iclip,rec)
    call gqfais(ier,ifais)
    call gsclip(0)
    call gsfais(1)
    !
    !  ..set the tension on the two dimensional smoother
    !
    call cpsetr('t2d',smooth)
    !
    ! arrange for the selection of contour levels as desired by the user.
    !
    if (finc.lt.0.) then
        call cpseti('cls - contour level selector',max(1,int(-finc)))
        call cpsetr('cis - contour interval specifier',0.)
        elseif (finc.eq.0.) then
        call cpseti('cls - contour level selector',16)
        call cpsetr('cis - contour interval specifier',0.)
    else
        call cpseti('cls - contour level selector',1)
        call cpsetr('cis - contour interval specifier',finc)
        if (flow.lt.fhgh) then
            call cpsetr('cmn - contour minimum',flow)
            call cpsetr('cmx - contour maximum',fhgh)
        endif
    endif
    !
    ! set up the desired mapping of output onto the plotter frame.
    !
    if (nset.lt.0) then
        call cpseti('set - do-set-call flag',1)
        call getset(xvpl,xvpr,yvpb,yvpt,xwdl,xwdr,ywdb,ywdt,lnlg)
        call cpsetr('vpl - viewport left edge',xvpl)
        call cpsetr('vpr - viewport right edge',xvpr)
        call cpsetr('vpb - viewport bottom edge',yvpb)
        call cpsetr('vpt - viewport top edge',yvpt)
        call cpseti('vps - viewport shape',0)
        elseif (nset.eq.0) then
        call cpseti('set - do-set-call flag',1)
        call cpsetr('vpl - viewport left edge',.05)
        call cpsetr('vpr - viewport right edge',.95)
        call cpsetr('vpb - viewport bottom edge',.05)
        call cpsetr('vpt - viewport top edge',.95)
        call cpseti('vps - viewport shape',4)
    else
        call cpseti('set - do-set-call flag',0)
    endif
    !
    ! decide what dash pattern to use.
    !
    if (lines.ne.0) then
        idsh=abs(ndsh)
        if (idsh.eq.0.or.idsh.eq.1.or.idsh.eq.1023) then
            idsh=ior(ishift(32767,1),1)
        else
            idsh=ior(ishift(idsh,6),iand(ishift(idsh,-4),63))
        endif
    endif
    !
    ! decide whether to label highs and lows or not.
    !
    if (nhgh.eq.0) then
        call cpsetc('hlt - high/low label text', &
        'h:b:$zdv$:e:''l:b:$zdv$:e:')
    else
        call cpsetc('hlt - high/low label text',' ')
    endif
    !
    ! initialize conpack and give it all array dimensions.
    !
    call cprect(zdat,kzdt,mzdt,nzdt,rwrk,lrwk,iwrk,liwk)
    !
    ! pick contour levels.
    !
    call cppkcl(zdat,rwrk,iwrk)
    !
    ! retrieve the contour levels selected, one at a time.  discard levels
    ! which are outside the range, if any, specified by the user-supplied
    ! values of flow and fhgh, and move the parameters for all remaining
    ! levels to the beginning of the parameter arrays.  set dash patterns
    ! for all levels.  the value of 'ciu' must be saved for possible
    ! restoration if it gets clobbered as a side effect of setting contour
    ! level 1.
    !
    call cpgetr('ciu - contour interval used',cinu)
    call cpgeti('ncl - number of contour levels',nclo)
    ncln=0
    do iclo=1,nclo
        call cpseti('pai - parameter array index',iclo)
        call cpgetr('clv - contour level',clev)
        if (flow.ge.fhgh.or.(clev.ge.flow.and.clev.le.fhgh)) then
            ncln=ncln+1
            if (ncln.ne.iclo) then
                call cpgeti('clu - contour level use flag',iclu)
                call cpseti('pai - parameter array index',ncln)
                call cpsetr('clv - contour level',clev)
                call cpseti('clu - contour level use flag',iclu)
                call cpseti('aia - area identifier above level',ncln+1)
                call cpseti('aib - area identifier below level',ncln)
                call cpseti('clc - contour line color index',-1)
                call cpsetc('cld - contour line dash pattern', &
                '$$$$$$$$$$$$$$$$')
                call cpseti('cll - contour line line width',-1)
                call cpseti('llc - line label color index',-1)
                call cpsetc('llt - line label text',' ')
            endif
        endif
        if (ndsh.gt.0.or.(ndsh.lt.0..and.clev.lt.0.)) &
        call cpseti('cld - contour line dash pattern',idsh)
    enddo
    !
    ! log an error if no contour levels were within the user's bounds.
    !
    if (ncln.eq.0) then
        call seter('cnrccf - no contour levels in specified range', &
        1,2)
        return
    endif
    !
    ! if the number of contour levels decreased, reset parameters affected.
    !
    if (ncln.lt.nclo) then
        call cpseti('ncl - number of contour levels',ncln)
        call cpsetr('ciu - contour interval used',cinu)
    endif
    !
    !  ..default color table
    !
    if (icolfg.eq.0) then
        ncol=ncln+1
        call dfclrs(ncol)
    endif
    !
    !  ..color filled contour
    !
    call arinam(iama,lama)
    call cpclam(zdat,rwrk,iwrk,iama)
    call arscam(iama,xcra,ycra,ncra,iaia,igia,ngrps,colram)
    !
    !  ..map
    !
    !     call getusv('ii',icolor)
    !     call setusv('ii',icmp)
    !     call mapgrd
    !     call maplbl
    !     call setusv('lw',2000)
    !     call maplot
    !     call setusv('lw',1000)
    !     call maplbm
    !     call setusv('ii',icolor)
    !
    ! if requested, put out a simple background.
    !
    if (nset.eq.0) call cpback(zdat,rwrk,iwrk)
    !
    ! see how the user has chosen to position contour levels.
    !
    call cpgeti('llp - line label positioning flag',llpf)
    !
    ! draw the contour lines, masking them if necessary.
    !
    if (lines.ne.0) then
        if (llpf.le.1) then
            call gqplci(ier,icolor)
            call gsplci(0)
            call cpcldr(zdat,rwrk,iwrk)
            call gsplci(icolor)
        else
            call arinam(iama,lama)
            call cplbam(zdat,rwrk,iwrk,iama)
            call cpcldm(zdat,rwrk,iwrk,iama,cpdrpl)
        endif
    endif
    !
    ! plot labels.
    !
    call cpgetc('ilt - informational label text',ilts)
    !     print*,'informational label is ',ilts
    call cplbdr(zdat,rwrk,iwrk)
    !
    ! if requested, label every point on the grid.
    !
    if (nhgh.gt.0) then
        call getset(xvpl,xvpr,yvpb,yvpt,xwdl,xwdr,ywdb,ywdt,lnlg)
        call cpgetr('cwm - character width multiplier',chwm)
        call cpgetr('hla - high/low label angle',angd)
        call cpgetr('hls - high/low label size',size)
        call cpgeti('map - mapping flag',imap)
        call cpgetr('orv - out-of-range value',orva)
        call cpgetr('spv - special value',spva)
        call cpgetr('xc1 - x coordinate at i = 1',xca1)
        call cpgetr('xcm - x coordinate at i = m',xcam)
        call cpgetr('yc1 - y coordinate at j = 1',yca1)
        call cpgetr('ycn - y coordinate at j = n',ycan)
        size=(xvpr-xvpl)*chwm*size
        if (xca1.eq.xcam) then
            xca1=1.
            xcam=real(mzdt)
        endif
        if (yca1.eq.ycan) then
            yca1=1.
            ycan=real(nzdt)
        endif
        do j=1,nzdt
            ypos=yca1+real(j-1)*(ycan-yca1)/real(nzdt-1)
            do i=1,mzdt
                xpos=xca1+real(i-1)*(xcam-xca1)/real(mzdt-1)
                if (spva.eq.0..or.zdat(i,j).ne.spva) then
                    call cpsetr('zdv - z data value',zdat(i,j))
                    call cpgetc('zdv - z data value',croz)
                    do k=locv,2,-1
                        if (croz(k:k).ne.' ') then
                            lcrz=k
                            lexit = .true.
                            exit
                        endif
                    enddo
                    if (.not.lexit) lcrz=1
                    if (imap.eq.0) then
                        call plchhq(xpos,ypos,croz(1:lcrz),size,angd,0.)
                    else
                        call cpmpxy(imap,xpos,ypos,xmpd,ympd)
                        if (orva.eq.0..or.xmpd.ne.orva) &
                        call plchhq(xmpd,ympd,croz(1:lcrz),size,angd,0.)
                    endif
                endif
            enddo
        enddo
    endif
    !
    !  ..reset gks internal parameters
    !
    call gsclip(iclip)
    call gsfais(ifais)
    return
end
