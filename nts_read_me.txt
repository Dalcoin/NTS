
      Documentation 2.0.0

      Previous versions:
 
      0.5.0: NTS data generated into one file.
 
      0.5.5: Options 1) and 2), generating NTS data into multiple seperate files or one file.

      1.0: Generating EOS, NTS data.
 
      1.0.5: Generating EOS, NTS data into different files or all togeather with formatting options. 

      1.5.0: Generating EOS, NTS data and extracting maximum mass and corrosponding radius, central density and central energy density.

      1.5.5: Adding 1.4 SM values to 1.5.0

      1.5.5.m: Enhanced 1.4 SM and MM routines with ECC.

      Update:
      
      2.0.0: Adding tabled values for various polytropic expansions, complete with 1.4 SM and MM SM data, and causality limit. [Current Version]

      Upcomming: 

      2.0.5: Cleaning option for single gamma match.  

------------------------------------------------------------------------------------------------------------------------------------------------------|
                                                                                                                                                      |
      Simple Directions:                                                                                                                              |  
                                                                                                                                                      |
      1.) Inport EOS that is exactly 97 values long with Energy density, Pressure, and Density in that order. Put in file titled 'eosin.don'          |
                                                                                                                                                      | 
      2.) Set values in file titled 'in.don' Default values should be '3 false false 1'                                                               |
                                                                                                                                                      |
      3.) Compile and run program the bash script "./nts.sh"                                                                                          |
                                                                                                                                                      |
      4.) NTS values for maximum mass and mass at 1.4 solar masses printed in 'maxmass.srt' and 'mas14.srt' respectively                              |
                                                                                                                                                      |
          'ntsall.don' contains the output of NTS values for each polytropic parameter change.                                                        | 
                                                                                                                                                      |
          'eosall.don' contains the output of EOS values for each polytropic parameter change.                                                        |   
                                                                                                                                                      |
          'causnts.don' containts the density for which causality is violated for each polytropic parameter change.                                   |
                                                                                                                                                      |
          'nset_x.txt' x is an integer from 1 to ngval. These are formated NTS values for easy graphing.                                              |
                                                                                                                                                      | 
------------------------------------------------------------------------------------------------------------------------------------------------------|       
      Input Documentation:

      Parameters in 'paramnlo.don'. -> 'paramnlo.don' has been renamed 'parameos.don' in version 1.5.5.

      Parameter defaults are below.

      Gamma values in 'gamval.don'.

      Default values are set to double matching seven by seven
      polytropic expansion starting with a adiabatic parameter
      of 1.5 and running to 4.5 with increments of 0.5.

      Number density matches at 0.33 and again at 0.46 using
      from the data found from asymmetric isospin EOS derived using
      chiral nuclear force models. 

      The 'eosin.don' file contains the EOS. N2LO and N3LO have been 
      verified with a lambda parameter of 500 MeV. 


      Program Directions:

      1.) Set values in 'parameos.don'. (These are the EOS parameters)

      Parameter list reads as below: (Note: in 1.5.5 ncontrol is no longer an EOS parameter)

      ngval, m, n, rho_i, rho_f, rho_fin, match

      ngval: Number of gamma values to be read, starting from first value: select 1 or higher integer
      m: number of values to be matched after second matching: select between 1 and 26
      n: number of values in EOS file: select between 1 and 97
      rho_i: first matching point: select between 0.01 and 1.5
      rho_f: second matching point: select between 0.01 and 1.5
      match: Determines whether or not to match once or twice: select 1 to match twice, 0 to match once.
      *ncontrol: Determines whether to output EOS to single file, individual files, out put NTS data or all of the above.
        

      Reference values:
      ngval=7           
      m=22
      n=97
      rho_i = 0.36d0
      rho_f = 0.51d0
      rho_fin = 1.8d0
      match = 1
      ncontrol = 2 (3 since version 1.0)

      *Warning: It is not recommended that you change anything except ngval, match and ncontrol.
                The system hasn't been tested extensively with any other changes.

      2.) Set gamma values in 'gamval.don'.

      e.g.):
      
      0.5
      1.0
      1.5
      .
      .
      .
      4.5

      Be sure that you do not have fewer gamma values than the integer you set as ngval!

      3.) Set values in 'paramnts.don'. (These are the nts parameters) 
          
          (Note: 'pnts.don' has been renamed 'paramnts.don' in version 1.5.5, format_val and ntsprint have been removed)

          Parameter list reads as below:

          nde, mrk4s, des, fd, dd, format_val, ntsprint

          nde: Number of differential equations: 2 (Don't change this values unless serious modifications are made to the TOV equations!)        

          mrk4s: Number of Runge-Kutta interations: 100000 (This shouldn't need to be changed either)

          des: Star central density increment value: 0.001 (This shouldn't need to be changed either)

          fd: First central number density value: 0.01 (This value may be changed as needed)
         
          dd: Last (dernier) number density value: 1.8 (This values may be changed as needed)     

          *format_val: Boolean, Set to True for expoential formatting of nts values, else False for long form. 
 
          *ntsprint: Boolean, Set to True if each nts values are to be printed in their own seperate files (per gamma-gamma value) 
                              as well as in the complete file.

          nstep: Set to 1 to print out seven files each with seven NTS values. Makes for easy graphing, 
                 only use if there are seven gamma values.

          *Warning: It is not recommended that you change anything except fd, dd and format_val and ntsprint. 
                    The system hasn't been tested extensively with any other changes.   

      4.) Set values in 'in.don' (New in version 1.5.5)

          Set values of ncontrol, format_val, ntsprint, noption14

          ncontrol: Determines whether to output EOS to single file, individual files, or not output EOS at all. [1, 2, 3 respectively]     

          format_val: Boolean, Set to True for expoential formatting of nts values, else False for long form.

          ntsprint: Boolean, Set to True if each nts values are to be printed in their own seperate files (per gamma-gamma value) 
                             as well as in the complete file.
          
          noption14: If 1 then the central density and radius at 1.4 solar masses will be printed to 'mas14.srt' 
                     and maximum mass values to 'maxmass.srt', else nothing happens. 

          ncaus: If 1 then the causality limits will be computed and printed to the file 'causnts.don'

      3.) Compile the program: 'f90 $F90FLAGS -o run -s -w nts.f  $LINK_FNL'

      3.5) Run the program: './run'

      3.m) Run bash script which compiles and runs: './nts.sh'

      4.) Output: 

      If ncontrol = 3: Polytropic/NXLO hybrid EOS will appear altogether in 'nloall.don', 
      individually in files labled fort.xxx, and the NTS data will appear in 'ntsall.don'
         
      If format_val was set to true then the formated NTS values will be printed, else the unformated will be printed in a single file.
 
      If ntsprint was set to true then each nts set will be printed seperately in it's own file. 
         
      If ncontrol = 1,2 Only the first or second above will be initalized, 1 corosponding to the first choice, 2 to the second.

      If ncontrol = 0 all choices will be initalized.   
      
      If noption = 1 Then the values corrposonding to maximum mass and 1.4 solar mass will be printed to 'maxmass.srt' and 'mas14.srt' respectively.


-----------------------------------------------------------------------------------------------------------------------|
                                                                                                                       |
                                                                                                                       |
                                                                                                                       |
  ______        _    _           _                 ____          ____                                                  |
 |  ____|      | |  | |         | |               |_   _|       |  __ \                                                |  
 | |__  __  __ | |  | |_ __ ___ | |__  _ __ __ _    | |  _ __   | |__) | __ __  _ ____   __ _ ___                      |  
 |  __| \ \/ / | |  | | '_ ` _ \| '_ \| '__/ _` |   | | | '_ \  |  ___/ '__| | | | '_ \ / _` / __|                     |
 | |____ >  <  | |__| | | | | | | |_) | | | (_| |  _| |_| | | | | |   | |  | |_| | | | | (_| \__ \                     |
 |______/_/\_\  \____/|_| |_| |_|_.__/|_|  \__,_| |_____|_| |_| |_|   |_|   \__,_|_| |_|\__,_|___/                     |
                                                                                                                       |
                                                                                                                       |
                                                                                                                       |
-----------------------------------------------------------------------------------------------------------------------|
