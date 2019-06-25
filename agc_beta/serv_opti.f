
       PROGRAM server
       implicit real*8 (a-h, o-z)
       common/block1/n,m,pr2o,r2o,dentemp(100),
     1               den(100),pres(100),eng(100),prssec(100)

c SERVER INTIALIZATION

       dimension :: rho(100), rho0(100), eps(100)
       dimension :: druke(100), epfrac(100)
       dimension :: break(100), coeff(4,100), depdp(100)

c       open(000,file='dump.don')
       open(100,file='betaeos.don')
       open(400,file='par.don')

       read(400,*) n
       n2=n-1
       m=5

       do i=1,n
          read(100,*) eng(i), pres(i), den(i)
       end do
 
       call dcsakm(n,den,pres,break,coeff)       
       r2o = 0.32d0
       pr2o  = dcsval(r2o,n2,break,coeff)
 
       do i=1,m
          dentemp(i) = (0.32d0)+(0.48d0-0.32d0)*((i-1)/4.d0)
          prssec(i)=dcsval(dentemp(i),n2,break,coeff)
c          write(000,*) dentemp(i), prssec(i)
       end do
       
c SERVER MAIN LOOP

c read the command and parameters
c if command is 0 (or rather "not 1"), evaluate energy and repeat
c if command is 1, stop the server
       icmd = 0
       do 100 while (icmd .ne. 1)
          
c         read the command and parameters of rho() from STDIN
c         Begin: Originalcode: Server Main Loop: read(*,*)
          read(*,*) icmd, gam 
c         End: Originalcode: Server Main Loop: read(*,*)
c         Begin: Modifiedcode: Server Main Loop: read(*,*)
c         read(*,*) icmd, rp, cp, wp, rn, cn, wn
c         End: Modifiedcode: Server Main Loop: read(*,*)

          if (icmd .ne. 1) then
c            calculate energy
             dev = sigma(gam)
c            write(000,*) alp,gam,dev
             write(*, *) dev
          endif 

 100   continue           
            
       stop
       END PROGRAM
   
       
       FUNCTION sigma(gam) 
       implicit real*8 (a-h,o-z)
       common/block1/n,m,pr2o,r2o,dentemp(100),
     1               den(100),pres(100),eng(100),prssec(100)
       sum=0.d0
       do 20 i=1,m  
          functfit= (pr2o/(r2o**gam))*(dentemp(i))**gam
          err=dabs((prssec(i)-functfit)/prssec(i)) 
   20     sum=sum+err
       sum=sum/5.d0   
       sigma=sum               
       return
       end 

