    
 
       program rho_0_prs
       implicit real*8(a-h,o-z)

       dimension :: eng(100), prs(100), den(100)
       dimension dehorsprs(100), facteurprs(4,100)
       dimension dehorseng(100), facteureng(4,100)


       open(555,file='eos.don')
       open(858,file='prs.don')
       open(999,file='eng_prs.don')

       rho0 = 0.16

       n=45
       n2=n-1
       do i=1,n
          read(555,*) eng(i), prs(i), den(i)
       end do
 
       call dcsakm(n,den,prs,dehorsprs,facteurprs) 
       call dcsakm(n,den,eng,dehorseng,facteureng) 

       prs0 = dcsval(rho0,n2,dehorsprs,facteurprs)
       eng0 = dcsval(rho0,n2,dehorseng,facteureng)

       write(858,5000) prs0
       write(999,6000) prs0, eng0

5000  format(2x,F7.3)
6000  format(2x,F7.3,2x,F7.3)

       end
