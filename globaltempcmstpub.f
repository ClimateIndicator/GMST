      program globalmeancmst
c     program to calculate spatial means from China-MST gridded data. Note
c     that only land values are used in paper (and in IPCC) due to use of
c     climatological values over sea ice. Full documentation in 
c     globaltemphadcrupub.f - this code is equivalent except for data ingest
      integer yearvals(200)
      real annobs(200,72,36)
      real obs(200,12,72,36)
      real obsline(200)
      real trendvals(72,36,3)
      integer sigmat(72,36)
      real spatmean(200,12,3,3)
      real annspatmean(200,3,3)
      real gridmean(72,36)
      real landmask(360,180)
      real avlandmask(72,36)
      integer gridcount(200,12,3,3)
      integer startyr
      integer endyr
      integer trendlim(2)
      integer trendstart
      integer trendend
      real regvals(2)
      integer length
      integer trendlen
      real residuals(200)
      real corrval
      real efflen
      real temp
      real errval
      integer year
      integer compflag
      do 10 i=1,200
       yearvals(i)=0
       obsline(i)=-999.0
10    continue
20    format(i4)
22    format('Enter first year of time series')
23    format('Enter last year of time series')
24    format('Enter first year for trend calculation')
25    format('Enter last year for trend calculation')
      write(6,22)
      read(5,20) startyr
      write(6,23)
      read(5,20) endyr
      write(6,24)
      read(5,20) trendstart
      write(6,25)
      read(5,20) trendend
      length=endyr-startyr+1
      trendlen=trendend-trendstart+1
      trendlim(1)=trendstart-startyr+1
      trendlim(2)=trendend-startyr+1
      do 28 i=1,length
       yearvals(i)=i+startyr-1
28    continue
29    format(f6.2)
      call obsingest(obs)
      call landmaskingest(landmask)
      call landmaskav(landmask,avlandmask)
      call anngridcalc(obs,annobs)
      call spatmeancalc(obs,length,spatmean,gridcount,avlandmask)
      call annspatcalc(spatmean,annspatmean,avlandmask)
      call gridmeancalc(annobs,gridmean)
      do 100 j=1,72
       do 90 k=1,36
        do 40 l=1,length
         obsline(l)=annobs(l,j,k)
40      continue
        call compcheck(obsline,trendlim,trendlen,compflag)
      if (compflag.eq.1) then
       call regcalc(yearvals,obsline,trendlim,regvals)
       call rescalc(yearvals,obsline,regvals,residuals,trendlim)
       call corrcalc(residuals,trendlim,corrval)
       efflen=trendlen*(1-corrval)/(1+corrval)
       call errcalc(residuals,trendlim,efflen,yearvals,errval)
       trendvals(j,k,1)=regvals(2)*10.0
       trendvals(j,k,2)=(regvals(2)-errval*1.645)*10.0
       trendvals(j,k,3)=(regvals(2)+errval*1.645)*10.0
        if ((trendvals(j,k,2)*trendvals(j,k,3)).ge.0) then
         sigmat(j,k)=1
        else
         sigmat(j,k)=0
        endif 
      else
       do 50 m=1,3
        trendvals(j,k,m)=-999.0
50     continue
       sigmat(j,k)=-99
55    format(2i4,f6.2)
      endif
90    continue
100   continue
160    format(72(1x,f9.4,1x))
161   format(72i4)
      open(unit=1,file='trendvals')
      open(unit=2,file='trendvalsext')
      open(unit=3,file='sigmat')
      do 170 k=1,3
       do 165 j=1,36
        write(2,160) (trendvals(i,j,k),i=1,72)
        if (k.eq.1) then
         write(1,160) (trendvals(i,j,k),i=1,72)
         write(3,161) (sigmat(i,j),i=1,72)
        endif
165     continue
170    continue
      close(3)
      close(2)
      close(1)
180   format(2i5,3(f6.2),3i6)
181   format(i5,3(f6.2))
182   format(72(f8.2))
183   format(2i5,6(f6.2),6i6)
184   format(i5,6(f6.2))
      open(unit=1,file='spatmean')
      do 200 i=1,length
       do 190 j=1,12
        write(1,180) (i+startyr-1),j,(spatmean(i,j,k,1),k=1,3),
     *(gridcount(i,j,k,1),k=1,3)
190    continue
200    continue
      close(1)
      open(unit=1,file='annspatmean')
      do 220 i=1,length
       write(1,181) (i+startyr-1),(annspatmean(i,k,1),k=1,3)
220   continue
      close(1)
      open(unit=1,file='gridmean')
      write(1,182) ((gridmean(j,k),j=1,72),k=1,36)
      close(1)
      open(unit=1,file='spatmeanlo')
      do 230 i=1,length
       do 225 j=1,12
        write(1,183) (i+startyr-1),j,((spatmean(i,j,k,l),k=1,3),l=2,3),
     *((gridcount(i,j,k,l),k=1,3),l=2,3)
225    continue
230   continue
      close(1)
      open(unit=1,file='annspatmeanlo')
      do 240 i=1,length
       write(1,184) (i+startyr-1),((annspatmean(i,k,l),k=1,3),l=2,3)
240   continue
      close(1)
      stop
      end
      subroutine regcalc(yearvals,obs,trendlim,regvals)
      integer yearvals(200)
      real obs(200)
      integer trendlim(2)
      real regvals(2)
      real xysum
      integer count
      real xsum
      real x2sum
      real ysum
      real means(2)
995   format(f5.2)
      count=0
      xsum=0.0
      ysum=0.0
      xysum=0.0
      x2sum=0.0
      do 1000 i=trendlim(1),trendlim(2)
       if (obs(i).gt.-999) then
        xsum=xsum+yearvals(i)+0.0
        ysum=ysum+obs(i)
        xysum=xysum+yearvals(i)*obs(i)
        x2sum=x2sum+(yearvals(i)**2)+0.0
        count=count+1
       endif
1000  continue
      means(1)=xsum/(count+0.0)
      means(2)=ysum/(count+0.0) 
      regvals(2)=(xysum-(xsum*ysum)/(count+0.0))/(x2sum-(xsum**2)/
     *(count+0.0))
      regvals(1)=means(2)-regvals(2)*means(1)
1010  format(i6)
      return
      end
      subroutine rescalc(yearvals,obs,regvals,residuals,trendlim)
      integer yearvals(200)
      real obs(200)
      real regvals(2)
      real residuals(200)
      integer trendlim(2)
      do 1100 i=1,200
       residuals(i)=-999.0
1100  continue
      do 1110 i=trendlim(1),trendlim(2)
       if (obs(i).gt.-999) then
        residuals(i)=obs(i)-(regvals(1)+regvals(2)*(yearvals(i)+0.0))
       else
        residuals(i)=-999.0
       endif
1110  continue
      return
      end
      subroutine corrcalc(residuals,trendlim,corrval)
      real residuals(200)
      integer trendlim(2)
      real corrval
      real corrmat(200,2)
      integer corrlen
      real xsum
      real ysum
      real x2sum
      real y2sum
      real xysum
      real corrmeans(2)
      corrlen=0
      do 1200 i=trendlim(1),(trendlim(2)-1)
      if ((residuals(i).gt.-999).and.(residuals((i+1)).gt.-999)) then
       corrlen=corrlen+1
       corrmat(corrlen,1)=residuals(i)
       corrmat(corrlen,2)=residuals((i+1))
      endif
1200  continue
      xsum=0.0
      ysum=0.0
      do 1210 i=1,corrlen
       xsum=xsum+corrmat(i,1)
       ysum=ysum+corrmat(i,2)
1210  continue
      corrmeans(1)=xsum/(corrlen+0.0)
      corrmeans(2)=ysum/(corrlen+0.0)
      xsum=0.0
      ysum=0.0
      x2sum=0.0
      y2sum=0.0
      xysum=0.0
      do 1220 i=1,corrlen
       xsum=xsum+(corrmat(i,1)-corrmeans(1))
       ysum=ysum+(corrmat(i,2)-corrmeans(2))
       x2sum=x2sum+(corrmat(i,1)-corrmeans(1))**2
       y2sum=y2sum+(corrmat(i,2)-corrmeans(2))**2
       xysum=xysum+(corrmat(i,1)-corrmeans(1))*(corrmat
     *(i,2)-corrmeans(2))
1220  continue
1229  format(2(f7.4))
1230  format(i6,5(f7.2))
      corrval=xysum/(sqrt((x2sum*y2sum)))
      return
      end
      subroutine errcalc(residuals,trendlim,efflen,yearvals,errval)
      real residuals(200)
      integer trendlim(2)
      real efflen
      integer count
      integer yearvals(200)
      real errval
      real sqsum
      real xsum
      real x2sum
      real xmean
      real s2eval
      sqsum=0.0
      xsum=0.0
      x2sum=0.0
      count=0
      do 1300 i=trendlim(1),trendlim(2)
       if (residuals(i).gt.-999) then
        sqsum=sqsum+residuals(i)**2
        xsum=xsum+yearvals(i)+0.0
        count=count+1
       endif
1300  continue
      xmean=xsum/(count+0.0)
      do 1310 i=trendlim(1),trendlim(2)
       if (residuals(i).gt.-999) then
        x2sum=x2sum+(yearvals(i)-xmean)**2
       endif
1310  continue
      s2eval=sqsum/(efflen-2.0)
      errval=sqrt((s2eval/x2sum))
      return
      end
      subroutine compcheck(obsline,trendlim,trendlen,compflag)
      real obsline(200)
      integer trendlim(2)
      integer trendlen
      integer compflag
      integer count(3)
      do 1500 i=1,3
       count(i)=0
1500  continue
      do 1600 i=trendlim(1),trendlim(2)
       if (obsline(i).gt.-999) then
        count(1)=count(1)+1
        if (i.lt.(trendlim(1)+0.1*(trendlen-1))) then
         count(2)=count(2)+1
        elseif (i.gt.(trendlim(1)+0.9*(trendlen-1))) then
         count(3)=count(3)+1
        endif
       endif
1600  continue
1610  format(4i6)
      if (count(1).ge.(0.7*trendlen)) then
        compflag=1
        if (count(2).lt.(0.02*trendlen)) then
         compflag=0
        endif
        if (count(3).lt.(0.02*trendlen)) then
         compflag=0
        endif
      else
        compflag=0
      endif
      return
      end
      subroutine obsingest(obs)
      real obs(200,12,72,36)
      do 1703 i=1,200
       do 1702 j=1,12
        do 1701 k=1,72
         do 1700 l=1,36
          obs(i,j,k,l)=-999.0
1700  continue
1701  continue
1702  continue
1703  continue
1705  format(i3)
c     replace ###DIRECTORY### with working directory
      open(unit=1,file=
     *'###DIRECTORY/cmsttest2022')
      do 1730 i=1,173
      write(6,1705) i
       do 1720 j=1,12
        do 1710 l=1,36
          read(1,*,end=1740) (obs(i,j,k,l),k=1,72)
1710    continue
1720   continue
1730  continue
1740  continue
1750  format(72(f8.2))
1751  format(2i4,f8.2)
      close(1)
      return
      end
      subroutine anngridcalc(obs,annobs)
      real obs(200,12,72,36)
      real annobs(200,72,36)
      real obssum
      integer obscount
      do 1902 i=1,200
       do 1901 j=1,72
        do 1900 k=1,36
         annobs(i,j,k)=-999.0
1900  continue
1901  continue
1902  continue
      do 1950 i=1,200
       do 1940 k=1,72
        do 1930 l=1,36
         obssum=0.0
         obscount=0
         do 1920 j=1,12
          if (obs(i,j,k,l).gt.-999) then
           obssum=obssum+obs(i,j,k,l)
           obscount=obscount+1
          endif
1920     continue
         if (obscount.eq.12) then
          annobs(i,k,l)=obssum/12.0
         else
          annobs(i,k,l)=-999.0
         endif
1930  continue
1940  continue
1950  continue
      return
      end
      subroutine spatmeancalc(obs,length,spatmean,gridcount,avlandmask)
      real obs(200,12,72,36)
      real avlandmask(72,36)
      integer length
      real spatmean(200,12,3,3)
      integer gridcount(200,12,3,3)
      real pi
      real spatsum(2,3)
      integer count(2,3)
      real weight
      real wsum(2,3)
      real lat(2)
      pi=3.14159
      do 2050 i=1,length
       do 2040 j=1,12
        do 2001 n=1,3
        do 2000 m=1,2
         spatsum(m,n)=0.0
         count(m,n)=0
         wsum(m,n)=0.0
2000    continue
2001    continue
2005  format(2(f6.2))
        do 2030 l=1,36
         lat(1)=((5.0*l)-95.0)*pi/180.0
         lat(2)=((5.0*l)-90.0)*pi/180.0
         weight=abs((sin(lat(1))-sin(lat(2))))
         do 2020 k=1,72
          if (obs(i,j,k,l).gt.-999) then
           if (lat(1).lt.0) then
            spatsum(1,1)=spatsum(1,1)+weight*obs(i,j,k,l)
            count(1,1)=count(1,1)+1
            wsum(1,1)=wsum(1,1)+weight
            if (avlandmask(k,l).ge.0.5) then
             spatsum(1,2)=spatsum(1,2)+weight*obs(i,j,k,l)
             count(1,2)=count(1,2)+1
             wsum(1,2)=wsum(1,2)+weight
            else
             spatsum(1,3)=spatsum(1,3)+weight*obs(i,j,k,l)
             count(1,3)=count(1,3)+1
             wsum(1,3)=wsum(1,3)+weight
            endif
           else
            spatsum(2,1)=spatsum(2,1)+weight*obs(i,j,k,l)
            count(2,1)=count(2,1)+1
            wsum(2,1)=wsum(2,1)+weight
            if (avlandmask(k,l).ge.0.5) then
             spatsum(2,2)=spatsum(2,2)+weight*obs(i,j,k,l)
             count(2,2)=count(2,2)+1
             wsum(2,2)=wsum(2,2)+weight
            else
             spatsum(2,3)=spatsum(2,3)+weight*obs(i,j,k,l)
             count(2,3)=count(2,3)+1
             wsum(2,3)=wsum(2,3)+weight
            endif
           endif
          endif
2020     continue
2030    continue
        do 2036 n=1,3
        do 2035 m=1,2
         if (wsum(m,n).gt.0) then
           spatmean(i,j,m,n)=spatsum(m,n)/wsum(m,n)
         else
           spatmean(i,j,m,n)=-999.0
         endif
         gridcount(i,j,m,n)=count(m,n)
2035    continue
        gridcount(i,j,3,n)=gridcount(i,j,1,n)+gridcount(i,j,2,n)
        if ((wsum(1,n).gt.0).and.(wsum(2,n).gt.0)) then
         spatmean(i,j,3,n)=(spatmean(i,j,1,n)+spatmean(i,j,2,n))/2.0
        else
         spatmean(i,j,3,n)=-999.0
        endif
2036  continue
2040  continue
2050  continue
      return
      end   
      subroutine annspatcalc(spatmean,annspatmean)
      real spatmean(200,12,3,3)
      real annspatmean(200,3,3)
      real sum
      integer count
      do 2110 i=1,200
       do 2100 k=1,3
        do 2095 n=1,3
        sum=0.0
        count=0
        do 2090 j=1,12
         if (spatmean(i,j,k,n).gt.-999) then
          sum=sum+spatmean(i,j,k,n)
          count=count+1
         endif
2090    continue
        if (count.eq.12) then
         annspatmean(i,k,n)=sum/12.0
        else
         annspatmean(i,k,n)=-999.0
        endif
2095  continue
2100   continue
2110  continue
      return
      end
      subroutine gridmeancalc(annobs,gridmean)
      real annobs(200,72,36)
      real gridmean(72,36)
      real sum
      integer count
      do 2220 j=1,72
       do 2210 k=1,36
        sum=0.0
        count=0
        do 2200 i=92,121
         if (annobs(i,j,k).gt.-999) then
          sum=sum+annobs(i,j,k)
          count=count+1
         endif
2200    continue
        if (count.ge.24) then
          gridmean(j,k)=sum/(count+0.0)
        else
          gridmean(j,k)=-999.0
        endif
2210  continue
2220  continue
      return
      end      
      subroutine landmaskingest(landmask)
      real landmask(360,180)
      do 2301 i=1,360
       do 2300 j=1,180
        landmask(i,j)=-999.0
2300   continue
2301  continue
      open(unit=1,file='berklandmask')
      do 2320 i=1,180
       read(1,*,end=2330) (landmask(j,i),j=1,360)
2320  continue
2330  continue
      close(1)
      return
      end
      subroutine landmaskav(landmask,avlandmask)
      real landmask(360,180)
      real avlandmask(72,36)
      real sum
      integer pos
      do 2430 i=1,72
       if (i.le.36) then
        pos=i+36
       else
        pos=i-36
       endif
       do 2420 j=1,36
        sum=0.0
        do 2410 k=1,5
         do 2400 l=1,5
          sum=sum+landmask((i*5+k-5),(j*5+l-5))
2400     continue
2410    continue
        avlandmask(pos,j)=sum/25.0
2420   continue
2430  continue
      return
      end
