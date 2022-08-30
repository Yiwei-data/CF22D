# UF22D
One may need an electronic structure software to run the UF22D method and there is a quick way to implement the UF22D functional in Gaussian 16. 

  First, make a new copy of the source codes (name the new folder as revG16), then modify the parameters of MN15 in Gaussian 16 source codes and recompile the codes as revG16. (To find out where the MN15 parameters are, one may simply grep the parameters of MN15 in the original Gaussian 16 source code.) When one utilizes the MN15 functional in revG16, it will instead use the UF22D functional; in the input file, one may used the keywords “MN15 em=GD3 iop(3/174=1000000,3/175=-1,3/176=1530000)” to get UF22D.

  The parameters of UF22D are available in the Supplementary Table S2 of this paper. And the corresponding codes can be found in ‘Code’ directory. In addition, the percentage of the HF exchange energy is 0.462806 for UF22D, and this needs to be changed accordingly in the following lines of the original source code utilam.F:
     else if(IMOpt.eq.-31) then
C       MN15
C       IMHF = 440000
C       UF22D
        IMHF = 462806
        
  Remember to run the test file in the ‘Example’ directory to make sure you get the same single-point energies as the benchmark results of UF22D for CH3 and H2O; it is important to test both CH3 and H2O. Also remember to use the original (unmodified) version if you want to use MN15.
