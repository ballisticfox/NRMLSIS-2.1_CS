# NRMLSIS-2.1_CS
Transferring NRMLSIS to C# and refactoring it for use in games like KSP and KSA.

Translated by ballisticfox.

Based on the original Fortran 90 Code hosted here: https://map.nrl.navy.mil/map/pub/nrl/NRLMSIS/NRLMSIS2.1/

---

Branches:

- main - stable state of the project, may not be a 1:1 code mapping from F90 to C#, but retains the same outputs.
- archive - an archive of the 1:1 F90 to C# mapping created at the beginning of the project.

Future branches:
- refactor - refactoring the code to be more human readable, adding new calc methods to sample only temperature, only gas density of a specific, ambient pressure, etc.
- modularization - moving the raw data out of splines into an XML or KSP-alike format.
- interface - adding an interface so external applications can call and access data from inside NRMLSIS.
- burst - second refactor to convert the code into burst compiled binaries.
