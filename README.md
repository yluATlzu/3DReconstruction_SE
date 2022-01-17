# How to do 3D Reconstruction with Spherical Embedding

----------------------------------------------------
Download Aspire 0.14 
----------------------------------------------------
Download Aspire 0.14 Matlab code from http://spr.math.princeton.edu/

Assuming that the Aspire package has been extracted to a folder named $ASPIRE.

Assuming that the files in 3DReconstruction_SE have been copied to a folder named $SE.

Start Matlab and do the following:


----------------------------------------------------
Installation
----------------------------------------------------
1. Change into the directory $ASPIRE
2. Run 'initpath'
3. Run 'install' to install ASPIRE (This only needs to be run once)


----------------------------------------------------
Initialization
----------------------------------------------------
1. Change into the directory $ASPIRE
2. Run 'initpath' (This needs to be run each time you start Matlab session)
3. Change into the directory $SE
4. Run 'initSEPath'  (This needs to be run each time you start Matlab session)


----------------------------------------------------
To run the experiments with simulation data
----------------------------------------------------
1. Change into the directory $SE/SimulatedData
2. Run 'produceSimulatedProjections(NumP, SNR)', where the values of NumP and SNR can be choosen from the following list
2. Run 'testSimulatedData(NumP, SNR)', where the values of NumP and SNR can be choosen from the following list

NumP 	  SNR

100     0.2

500     0.2

1000	  0.2

2000	  0.2

1000	  0.1

1000	  0.3

1000	  0.4

1000	  0.5

1000	  1

An Example: to test reconstruction from 100 images with SNR=0.2

 1. run  'produceSimulatedProjections(100,0.2)' 
 2. run 'testSimulatedData(100, 0.2)'

To computie RSME in degrees from MSE in radians, the following formula can be used:

RMSE = sqrt(MSE)*180/pi


----------------------------------------------------
To run the experiments with the real data
----------------------------------------------------
1. Change into the directory $SE/RealData/EMPIAR-10028
2. Run 'testRealData('EMPIAR-10028')
3. Change into the directory $SE/RealData/EMPIAR-10328
4. Run 'testRealData('EMPIAR-10328')



