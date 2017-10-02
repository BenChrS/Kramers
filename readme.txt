PROGRAM INFORMATION

program version 5
18.02.2015

Run the program "colnoise" in the build folder to run the program with standard settings.

xml-files can be used to load predefined settings. A xml-file "simulationOptions.xml" with the used options will be created at each run and stored in the folder of the current run. In order to load settings at the beginning of the program the xml-file with all predefined variables and the file "simulationOptionsXML.xsd" need to be in the same folder the program is called from. When this file can not be found, the program will use the standard hard coded variables instead. For the exact same results, the program can also load the seed of the random generator from the file called "randGen.xml", also stored each time in the folder of the current run.

Here is the list of all available potentials and their corresponding numbers (potNr).

potentials
0) U(x) = 0
1) U(x) = 1/2*k^2*x
2) mexican hat
  potDepth: well height
  potWidth: distance between minimum and maximum of the potential
3) U(x) = a*y^4 + b*y^3 + c*y^2 -- assymetric mexican hat
4)not working
5) time dependent mexican hat
  startPotDepth: well height before startTime
  endPotDepth: well height after endTime
  potWidth: distance between minimum and maximum of the potential
  startTime: potential well starts to move
  endTime: potential well stops to move
6) potential of double well potential with no potential to the right
   potDepth is the well height and 
   potWidth is the distance between minimum and maximum of the potential
  
correlation Functions
0) I(t) = 0.5*D/tau * exp(-t/tau)
1) I(t) = 1/(a*sqrt(Pi) * e^(-t/a)^2

noise Functions
0) no noise
1) white Noise (only works properly with sdeSolver 0)
2) colored Noise

sdeSolver
0) adamsBWhiteNoise
1) eulerColNoise
2) eulerCromerColNoise
3) adamsBColNoise (recommended)
4) eulerWhiteNoise


kinEnerAv:
parameters 
	Corr0:
		temperature=1
		mass=1
		D=10
		T=30
		nSteps=30*2**7
		np=8000
		x0=0
		v0=0
		duration: 50 minutes

	Corr3:
		(1)
		temperature=1
		mass=0.1
		T=5
		alpha=2*sqrt(m)
		gammMin -0.01
		MeijerGMin -0.06
		
		(2)
		temperature=1
		mass=0.1
		T=1
		alpha=10*sqrt(m)
		gammMin -0.3
