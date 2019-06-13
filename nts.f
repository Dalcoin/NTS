        
c      Version 2.5.0
        
       program nts_generator
       implicit real*8 (a-h,o-z)
       real*8 ndens,last_density 
       parameter (number_central_densities=30)                
       parameter (neos=97)                
       dimension cen_dens(number_central_densities)
       dimension radius(10000)
       dimension smass(number_central_densities)
       dimension ndens(neos)
       dimension eps(neos)
       dimension pr(neos)
       dimension y(2), ynew(2), ders(2)
       dimension gamout(100)            
       
       integer iii
       logical format_val
       logical ntsprint
             
       dimension gam(100), alf(100), c(100), alf2(100), c2(100)
       dimension eng(100), prs(100), den(100)
       dimension rho(100), drucke(100), dichte(100)
       dimension drucke1(100), dichte1(100)
       
       dimension cnd14(100), rad14(100), smm14(100)
       dimension dehors14(100), facteur14(4,100)
       dimension dehors14red(100), facteur14red(4,100)
       dimension dehors14prs(100), facteur14prs(4,100)
       dimension dehors14eng(100), facteur14eng(4,100)
       
       dimension cnd(100), rad(100), smm(100)
       dimension cnds(100), rads(100), smms(100)
       dimension ced(100), ceds(100)
        
       dimension dehors(100), facteur(4,100)
       dimension dehors1(100), facteur1(4,100)
       dimension dehors2(100), facteur2(4,100)
       dimension dehors3(100), facteur3(4,100)
       
       dimension outptd1(100), outptd2(100)
       dimension coeffd1(100), coeffd2(100)
       dimension dpde(100)
       
       dimension set1d(100)
       dimension set1r(100),set2r(100),set3r(100),set4r(100)
       dimension set5r(100),set6r(100),set7r(100),set0r(100)
       
       dimension set1m(100),set2m(100),set3m(100),set4m(100)
       dimension set5m(100),set6m(100),set7m(100),set0m(100)
       
       dimension den66(100), prs66(100), eng66(100)
        
       dimension outpt4(100), coeff4(4,100), css(100), css_comp(100)
       dimension outpt5(100), coeff5(4,100)
       dimension denex(100), engex(100), prsex(100)
       
       dimension vsl(100), denvsl(100), outpt6(100), coeff6(4,100)
       
       common/eosval1/outpt1(neos), coeff1(4,neos)
       common/eosval2/outpt2(neos), coeff2(4,neos)
       common/eosval3/outpt3(neos), coeff3(4,neos)
       common/block1/g,fourpi,n2,sm_solar
       common/block2/number_differential_eqs,max_rk4_steps,diff_eq_step
       common/block3/first_density,last_density, xkg, density_step
       common/block4/central_energy_density,const_1,const_2
       
       ncd=number_central_densities
         
       
       open(300,file='eosin.don')
       open(223,file='in.don')
       open(323,file='parameos.don')
       open(424,file='paramnts.don')
       open(525,file='gamval.don')
       open(626,file='maxmass.srt')
       open(727,file='mas14.srt')
       open(000,file='dump.don')
       open(900,file='eosall.don')
       open(666,file='ntsall.don')
       open(888,file='causnts.don')    

       g=6.67259D-45*197.327
       pi = 3.14597d0
       fourpi=12.56637061     
       n2= neos-1    

       read(323,*) ngval, m, n, rho_i, rho_f, rho_fin, match
       read(424,*) nde, mrk4s, des, fd, dd, nstep
       read(223,*) ncontrol, format_val, ntsprint, noption14, ncaus

       if(nstep .EQ. 1) then
          open(101,file='nset_1.txt')
          open(102,file='nset_2.txt')
          open(103,file='nset_3.txt')
          open(104,file='nset_4.txt')
          open(105,file='nset_5.txt')
          open(106,file='nset_6.txt')
          open(107,file='nset_7.txt')
       end if   
         
       
c      reference values
c      ngval=7 
c      m=22
c      n=97
c      rho_i = 0.36d0
c      rho_f = 0.51d0
c      rho_fin = 1.8d0
c      match = 1
c      ncontrol = 2


       if (match .EQ. 1) then
          valle = rho_i
       else if (match .NE. 1) then
          valle = rho_f
       end if

       do i=1,ngval
       read(525,*) gam(i)
       end do
        
       rho_inc = (rho_f - rho_i)/4.d0  
       rho_pas = (rho_fin - rho_f)/m
       rho_com = rho_f
        
       inac = 0
       do i=1,n
          read(300,*) eng(i), prs(i), den(i)
          if (den(i) .LT. valle) then
             rhobc = den(i)
             pbc = prs(i)
             ebc = eng(i) 
             inac = inac +1
c          write(000,*) ebc, pbc, rhobc
          end if     
       end do

c      Generates a list of Gammas and hybrid EOS all in the same file


c---------------------------------------------------------------------------------------------------

       if (ncontrol .NE.  1 .AND. ncontrol .NE. 0) then
          go to 111   
       end if

c---------------------------------------------------------------------------------------------------        
       
       do i=1,ngval
          alf(i) = pbc/(rhobc**gam(i))
          c(i) = (1/rhobc)*(ebc-((alf(i)/(gam(i)-1))*(rhobc**gam(i))))
c          write(000,*) gam(i), alf(i), c(i) 
       end do
       
c       jj=1
       do jj=1,ngval
          do j=1,ngval
             write(900,*) gam(jj), gam(j) 
             write(900,*)"         "
             do k=1,inac
                write(900,*) eng(k), prs(k), den(k)
             end do
             if (match .EQ. 1) then
                do i=1,4  
                   rho(i)=rho_i+rho_inc*(i-1.d0)
                   drucke1(i)=druck(rho(i),gam(jj),alf(jj))
                   dichte1(i)=avg(rho(i),gam(jj),alf(jj),c(jj))
                   write(900,*) dichte1(i), drucke1(i), rho(i)
                end do             
                alf2(j) = drucke1(4)/(rho(4)**gam(j))
             c2(j) = (1.d0/rho(4))*(dichte1(4)-((alf2(j)/(gam(j)-1.d0))*
     1             (rho(4)**gam(j))))             
c               write(000,*) rho(4), drucke1(4), dichte1(4)
c               write(000,*) gam(j), alf2(j), c2(j)             
                do i=1,22
                   rho(i)=rho_com+rho_pas*(i-1)
                   drucke(i)=druck(rho(i),gam(j),alf2(j))
                   dichte(i)=avg(rho(i),gam(j),alf2(j),c2(j))
                   write(900,*) dichte(i), drucke(i), rho(i)
                end do
             else if (match .NE. 1) then
                do i=1,26
                   rho(i)=rho_i+rho_inc*(i-1)
                   drucke1(i)=druck(rho(i),gam(jj),alf(jj))
                   dichte1(i)=avg(rho(i),gam(jj),alf(jj),c(jj))
                   write(900,*) dichte1(i), drucke1(i), rho(i)
                end do
             end if 
             write(900,*)"---------------------------------------------"
          end do
       end do 
        
111    continue 

c---------------------------------------------------------------------------------------------------

       if (ncontrol .NE. 2 .AND. ncontrol .NE. 0) then
          go to 1445
       end if

c---------------------------------------------------------------------------------------------------
              
       do i=1,ngval
          alf(i) = pbc/(rhobc**gam(i))
      c(i) = (1.d0/rhobc)*(ebc-((alf(i)/(gam(i)-1.d0))*(rhobc**gam(i))))
c          write(000,*) gam(i), alf(i), c(i)
       end do


       do jj=1,ngval
          do j=1,ngval
             iii = 900 + j + jj*10
             do k=1,inac
                write(iii,*) eng(k), prs(k), den(k)
             end do
             if (match .EQ. 1) then
                do i=1,4  
                   rho(i)=rho_i+rho_inc*(i-1.d0)
                   drucke1(i)=druck(rho(i),gam(jj),alf(jj))
                   dichte1(i)=avg(rho(i),gam(jj),alf(jj),c(jj))
                   write(iii,*) dichte1(i), drucke1(i), rho(i)
                end do             
                alf2(j) = drucke1(4)/(rho(4)**gam(j))
             c2(j) = (1.d0/rho(4))*(dichte1(4)-((alf2(j)/(gam(j)-1.d0))*
     1             (rho(4)**gam(j))))                     
                do i=1,22
                   rho(i)=rho_com+rho_pas*(i-1)
                   drucke(i)=druck(rho(i),gam(j),alf2(j))
                   dichte(i)=avg(rho(i),gam(j),alf2(j),c2(j))
                   write(iii,*) dichte(i), drucke(i), rho(i)
                end do
             else if (match .NE. 1) then
                do i=1,26
                   rho(i)=rho_i+rho_inc*(i-1)
                   drucke1(i)=druck(rho(i),gam(jj),alf(jj))
                   dichte1(i)=avg(rho(i),gam(jj),alf(jj),c(jj))
                   write(iii,*) dichte1(i), drucke1(i), rho(i)
                end do
             end if 
          end do
       end do 
             
1445   continue

c---------------------------------------------------------------------------------------------------
       
       if (ncontrol .NE. 1 .AND. ncontrol .NE. 0) then
          go to 666
       end if

c---------------------------------------------------------------------------------------------------


       do i=1,ngval
          alf(i) = pbc/(rhobc**gam(i))
          c(i) = (1/rhobc)*(ebc-((alf(i)/(gam(i)-1))*(rhobc**gam(i))))
       end do

       do jj=1,ngval
          do j=1,ngval
             write(666,*) gam(jj), gam(j) 
             write(666,*)"         "
             ijk = 800 + j + jj*10
             if (match .EQ. 1) then
                do i=1,4
                   den(i+inac)=rho_i+rho_inc*(i-1.d0)
                   prs(i+inac)=druck(den(i+inac),gam(jj),alf(jj))
                   eng(i+inac)=avg(den(i+inac),gam(jj),alf(jj),c(jj))
                end do
                alf2(j) = prs(inac+4)/(den(inac+4)**gam(j))
                c2(j) = (1/den(inac+4))*(eng(inac+4)-((alf2(j)/
     1                  (gam(j)-1))*(den(inac+4)**gam(j))))
                do i=1,22
                   den(i+4+inac)=rho_com+rho_pas*(i-1)
                   prs(i+4+inac)=druck(den(i+4+inac),gam(j),alf2(j))
                   eng(i+4+inac)=avg(den(i+4+inac),gam(j),alf2(j),c2(j))
                end do
             else if (match .NE. 1) then
                do i=1,26
                   den(i+inac)=rho_i+rho_inc*(i-1.d0)
                   prs(i+inac)=druck(den(i+inac),gam(jj),alf(jj))
                   eng(i+inac)=avg(den(i+inac),gam(jj),alf(jj),c(jj))
                end do
             end if

c            troubleshoot sequence

c             do i = 1,inac+26
c                write(000,*) den(i), prs(i), eng(i)
c             end do

c---------------------------------------------------------------------------------------------------

c             Cluckly Clumper Group

c---------------------------------------------------------------------------------------------------

c            
c              number_differential_eqs=2   
c number of Runge-Kutta iterations 
c              max_rk4_steps=100000        
c dimensionless step in RK4 procedure 
c              diff_eq_step=0.001     
c chosen first and last number densities 
c              first_density=0.01      
c              last_density=1.8


             number_differential_eqs=nde
c number of Runge-Kutta iterations
             max_rk4_steps=mrk4s
c dimensionless step in RK4 procedure
             diff_eq_step=des
c chosen first and last number densities
             first_density=fd
             last_density=dd

c conversion factor
             xkg=1.E+30/1.78266270
c solar mass in MeV/c^2
             sm_solar=1.989E+30*xkg
             density_step=(last_density - first_density)/
     1             number_central_densities

c---------------------------------------------------------------------------------------------------

c             Splunky Spliner Group

c---------------------------------------------------------------------------------------------------

             call dcsakm(neos,prs,eng,outpt1,coeff1)
             call dcsakm(neos,den,prs,outpt2,coeff2)
             call dcsakm(neos,den,eng,outpt3,coeff3)

c---------------------------------------------------------------------------------------------------


              
c---------------------------------------------------------------------------------------------------
c                                                                                                  |
c             Just Cause Group: Causality Limit                                                    |
c                                                                                                  |
c---------------------------------------------------------------------------------------------------

             if (ncaus .EQ. 1) then
                bvr=1.d0

c                 call dcsakm(neos,eng,prs,outptd1,coeffd1)

                n73=24 !24
                n273=n73-1
                n37=73 !73

                do i=1,n73
                   den66(i)=den(i+n37)
                   eng66(i)=eng(i+n37)
                   prs66(i)=prs(i+n37)
                end do

                nst=25 !24
                nst2=nst-1
                ninv=60 !73

                do i=1,nst
                   denex(i)=den(i+ninv)
                   engex(i)=eng(i+ninv)
                   prsex(i)=prs(i+ninv)
                end do
                 
                 
                call dcsakm(nst,engex,prsex,outpt4,coeff4)
                call dcsakm(nst,denex,engex,outpt5,coeff5)

                call dcsakm(n73,eng66,prs66,outptd1,coeffd1)

      
                do i=1,n73
                   dpde(i)=dcsder(1,eng66(i),n273,outptd1,coeffd1)
                   write(000,*) den66(i), eng66(i), dpde(i)
                end do
                 
                call dcsakm(n73,dpde,den66,outptd2,coeffd2)
                denlim=dcsval(bvr,n273,outptd2,coeffd2)
                write(000,*) denlim
             end if
              
c---------------------------------------------------------------------------------------------------

c             Nunts Crumper Group: Neutron Star Mass, Radius, CenDen

c---------------------------------------------------------------------------------------------------           

             DO 200 i=1, number_central_densities
                ync=first_density+density_step*i
c                ync=last_density
c    modification
                central_energy_density=dcsval(ync,n2,outpt3,coeff3) 
c                central_energy_density=dcsval(0.3944d0,n2,outpt3,coeff3)
                ced = central_energy_density
c    end modification                    
                pressure=dcsval(ync,n2,outpt2,coeff2) 
                const_1=1.E-19/sqrt(fourpi*g*central_energy_density)
                const_2=fourpi*central_energy_density/   
     1             (sqrt(fourpi*g*central_energy_density)**3)
                y(1)=pressure/central_energy_density
                star_radius=diff_eq_step/10.
                y(2)=((star_radius)**3)/3.
c    Start Runge-Kutta solution for fixed central density 
                k=0
                do while((pressure > 0.).and.(k<= max_rk4_steps))
                   k=k+1 
                   star_radius=star_radius+diff_eq_step       
                   call derivatives(star_radius,y,ders)
                   call runge_kutta_4(star_radius,y,ynew,ders)
                   y(1)=ynew(1)
                   y(2)=ynew(2)
                   pressure=y(1)*central_energy_density
c                   radius(i) = star_radius*const_1*10.
c                   ebdb = ced*ders(2)/(star_radius * star_radius)
                end do 
                tp = 10.d0
                cen_dens(i)=ync
                radius(i)=star_radius*const_1
                smass(i)=y(2)*const_2/sm_solar

                cnd(i)=cen_dens(i)
                rad(i)=radius(i)*tp  
                smm(i)=smass(i)  

                neos1=neos-1
                css(i)= dsqrt(dcsder(1,dcsval(ync,nst2,outpt5,coeff5)
     1                        ,nst2,outpt4,coeff4))         
                 
c                Printing subgroup
                 
c                 write(000,*) ntsprint
                if (ntsprint .EQ. .TRUE.) then
                   WRITE(ijk,*) cen_dens(i),radius(i)*tp, smass(i)
                end if
  
                if (format_val .EQ. .TRUE.) then
                   WRITE(666,'(E12.3,6E12.4)') cen_dens(i),
     1                                          radius(i)*tp, smass(i)
                else if (format_val .NE. .TRUE.) then
                WRITE(666,500) cen_dens(i), radius(i)*tp, 
     1                          smass(i),css(i)
                end if
                 
c                49ers Subgroup: Prints for easy graphing
                 
                if (nstep .EQ. 1) then
                   if(j .EQ. 1) then
                      set1d(i)=cen_dens(i) 
                      set1r(i)=radius(i)*tp
                      set1m(i)=smass(i)
                   end if
                   if(j .EQ. 2) then
                      set2r(i)=radius(i)*tp
                      set2m(i)=smass(i)
                   end if
                   if(j .EQ. 3) then
                      set3r(i)=radius(i)*tp
                      set3m(i)=smass(i)
                   end if
                   if(j .EQ. 4) then
                      set4r(i)=radius(i)*tp
                      set4m(i)=smass(i)
                   end if
                   if(j .EQ. 5) then
                      set5r(i)=radius(i)*tp
                      set5m(i)=smass(i)
                   end if
                   if(j .EQ. 6) then
                      set6r(i)=radius(i)*tp
                      set6m(i)=smass(i)
                   end if
                   if(j .EQ. 7) then
                      set7r(i)=radius(i)*tp
                      set7m(i)=smass(i)
                   end if
                end if
                    

c                Just Cause subgroup: 

                if (ncaus .EQ. 1) then
                   if(cen_dens(i) .LT. denlim) then
c                      write(888,500) cen_dens(i), radius(i)*tp,smass(i)
                       causden=cen_dens(i) 
                       causrad=radius(i)*tp
                       causmass=smass(i)
                   end if
                end if              
200          continue 
500          format(2x,E10.4,2x,E10.4,2x,F6.4,2x,F10.4) 
c      stop      

c                Just Cause subgroup

             WRITE(666,*)'-----------------------------------'
             if(ncaus .EQ. 1) then
                write(888,*) gam(jj),gam(j),causmass,causrad,causden
             end if

             ivsl = 2
             nedy = ncd-2
             nedy1 = nedy-1
             do i=1,nedy
                denvsl(i) = cnd(i+ivsl)
                vsl(i) = css(i+ivsl)
             end do
             call dcsakm(nedy,denvsl,vsl,outpt6,coeff6) 
c---------------------------------------------------------------------------------------------------
             
c             Splunky Spliners Mad Max and 14 Group
             
c---------------------------------------------------------------------------------------------------
             
             if (noption14 .EQ. 1) then                          
                smmmax=3.0d0
                smmmin=0.9d0

                n2=n-1
                ncd2=ncd-1

                smmtemp=0.d0
                mm_c=0
                do i=1,ncd2
                   mm_c=mm_c+1
                   if(smm(i) .LE. smmmax .AND. smm(i) .GE. smmmin) then
                      if(smm(i) .GT. smmtemp) then
                         xmaxim=smm(i)
                         nmax=mm_c 
c                         write(000,*) gam(jj), gam(j), nmax
                      end if
                   end if
                   smmtemp=smm(i) 
c                       write(000,*) gam(jj), gam(j), smmtemp
                end do   

                h_mms=smm(nmax) 
c                write(000,*) nmax, h_mms
                  
                call dcsakm(ncd,rad,cnd,dehors1,facteur1)
                cen_num_den=dcsval(rad(nmax),ncd2,dehors1,facteur1)
                                  
                call dcsakm(n,den,eng,dehors2,facteur2)
                cen_eng_den=dcsval(cen_num_den,n2,dehors2,facteur2)
                   
                call dcsakm(n,den,prs,dehors3,facteur3)
                cen_prs=dcsval(cen_num_den,n2,dehors3,facteur3)

                cen_vsl = dcsval(cen_num_den,nedy1,outpt6,coeff6)

                write(626,8934) gam(jj),gam(j),h_mms,rad(nmax),
     1                          cen_num_den,cen_eng_den,cen_prs,cen_vsl

8934            FORMAT(2x,F3.1,2x,F3.1,2x,F9.4,2x,F9.4,2x,
     $                 F9.4,2x,F9.4,2x,F9.4,2x,F9.4)           

c               Mas-14 Subgroup
            
                nmax1=nmax-1 
                nmax2=nmax1-1                   
                     
                if(smm(nmax) .LT. 1.4d0) then
                   cen_den_14=0.d0
                   rad_14=0.d0
                   go to 223
                end if
                                           
                rad_cutoff=50.d0

                do i=1,nmax1
                   cnd14(i)=cnd(i)
                   rad14(i)=rad(i)
                   smm14(i)=smm(i)
                   if(rad14(i) .GT. rad_cutoff) then
                      smm14(i)=0.05d0
                      rad14(i)=30.d0
                   end if
                end do
                   
                
                call dcsakm(nmax1,smm14,cnd14,dehors14,facteur14)
                cen_den_14 = dcsval(1.4d0,nmax2,dehors14,facteur14)
       
                call dcsakm(nmax1,cnd14,rad14,
     1                      dehors14red,facteur14red)
                rad_14 = dcsval(cen_den_14,nmax2,
     1                          dehors14red,facteur14red)
                       
                eng_14 = dcsval(cen_den_14,n2,dehors2,facteur2) 
                   
                cen_vsl_14=dcsval(cen_den_14,nedy1,outpt6,coeff6)
              
223             continue 
                    
                write(727,8933) gam(jj), gam(j),
     1          rad_14, cen_den_14, eng_14, cen_vsl_14
                   
8933            FORMAT(2x,F3.1,2x,F3.1,2x,F9.4,2x,F9.4,2x,F9.4,2x,F9.4)

                else 
                continue
             end if
          end do

          if (noption14 .EQ. 1) then
             write(626,*)"                                "
             write(727,*)"                                "
          else
             continue
          end if
          
          if (nstep .EQ. 1) then
             do i=1,ncd
                write(100+jj,669) set1d(i),set1r(i),set1m(i),
     1                            set2r(i),set2m(i),set3r(i),
     1                            set3m(i),set4r(i),set4m(i),
     1                            set5r(i),set5m(i),set6r(i),
     1                            set6m(i),set7r(i),set7m(i)
             end do
          end if

669       format(E10.4,2x,E10.4,2x,F6.4,2x,E10.4,2x,F6.4,2x,E10.4,2x,
     $           F6.4,2x,E10.4,2x,F6.4,2x,E10.4,2x,F6.4,
     $           2x,E10.4,2x,F6.4,2x,E10.4,2x,F6.4)
   

          if(ncaus .EQ. 1) then
             write(888,*)"   "
          end if
       end do
666    continue

       end

c      This is Program's End. Here after are the Subroutines and Functions. 

c---------------------------------------------------------------------------------------------------

c      Function Set 1

c---------------------------------------------------------------------------------------------------


       function avg(r,gam,alfa,cc)
          implicit real*8(a-h,o-z)
          avg=(alfa/(gam-1.d0))*r**gam+(cc*r) 
          return
       end

       function druck(r,gam,alfa)
          implicit real*8(a-h,o-z) 
          druck=(alfa*(r**gam))
          return 
       end       

c---------------------------------------------------------------------------------------------------

c      Function\Subroutine Set 2

c---------------------------------------------------------------------------------------------------


       SUBROUTINE runge_kutta_4(x,y,yout,dydx)
          implicit real*8(a-h,o-z) 
          real*8 last_density 
          dimension yt(2), dyt(2), dym(2) 
          dimension y(2), dydx(2), yout(2) 
          common/block1/g,fourpi,n2,sm_solar 
          common/block2/number_differential_eqs,max_rk4_steps,
     1                  diff_eq_step 
          common/block3/first_density,last_density,xkg,density_step 
          common/block4/central_energy_density,const_1,const_2  

          hh=diff_eq_step*0.5                                  
          h6=diff_eq_step/6.                 
          xh=x+hh 
          yt(1)=y(1)+hh*dydx(1)
          yt(2)=y(2)+hh*dydx(2) 
          call derivatives(xh,yt,dyt)
          yt(1)=y(1)+hh*dyt(1) 
          yt(2)=y(2)+hh*dyt(2) 
          call derivatives(xh,yt,dym)
          yt(1)=y(1)+diff_eq_step*dym(1)              
          yt(2)=y(2)+diff_eq_step*dym(2)              
          dym(1)=dyt(1)+dym(1) 
          dym(2)=dyt(2)+dym(2) 
          call derivatives(x+diff_eq_step,yt,dyt)
          yout(1)=y(1) + h6*(dydx(1)+dyt(1)+2.*dym(1))
          yout(2)=y(2) + h6*(dydx(2)+dyt(2)+2.*dym(2))
       return 
       end                              


       subroutine derivatives(r,y,ders)                 
          implicit real*8(a-h,o-z) 
          real*8 last_density 
          parameter (neos=97)                
          dimension y(2), ders(2)                      
          common/eosval1/outpt1(neos), coeff1(4,neos)                  
          common/eosval2/outpt2(neos), coeff2(4,neos)                  
          common/eosval3/outpt3(neos), coeff3(4,neos)                  
          common/block1/g,fourpi,n2,sm_solar 
          common/block2/number_differential_eqs,max_rk4_steps,
     1                  diff_eq_step 
          common/block3/first_density,last_density, xkg,density_step
          common/block4/central_energy_density,const_1,const_2   

          if(y(1).gt.0.) then 
             p1=y(1)*central_energy_density
             e_rho=dcsval(p1,n2,outpt1,coeff1)
             e_rho=e_rho/central_energy_density
c
             ders(1)=tov(r,e_rho,y)
             ders(2)=(r**2)*e_rho
          ENDIF
         return                      
       end 



       DOUBLE PRECISION FUNCTION tov(r,e_rho,y)
          implicit real*8(a-h,o-z)
          dimension y(2)   
          tov=-(e_rho+y(1))*                 
     1       (y(2) + (r**3)*y(1))/(r*r-2*r*y(2))
          return            
       end     
       function gamval(l)
          implicit real*8(a-h,o-z)
          gamval = 3.0 + 0.5*(l-1) 
          return 
       end
