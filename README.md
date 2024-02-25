# Huckel method code
# Description
This FORTRAN 90 code is developed for the computation of the eigenvalues and eigenvectors of polyenes with Huckel method.
The Huckel method is a semi-empirical approach that provides a simplified yet insightful understanding of the electronic structure and properties of conjugated organic molecules. It neglects the σ-electrons and considers only the π-orbitals that are involved in the delocalized electron cloud. By doing so, the computational complexity is significantly reduced. The Hückel method employs the secular determinant approach to solve the Schrödinger equation for the molecular orbitals in the π-system. 

This program aims to build the Hamiltonian matrix (Huckel Hamiltonian) for the system specified in the input file and diagonalize it, in order to obtain the aigenvalues and the eigenvectors.

# System Requirements
- Fortran 90 Compiler
- LAPACK library

# Folder Structure
1. Input file `huckel.in`
2. Main program `huckel.f90`
3. Executable file `huckel`
4. Output file `huckel.dat` generated after the execution of the program
5. `eigenvectors.dat` and `eigenvalues.dat` files containing the eigenvectors and the eigenvalues in a format suitable for plotting with gnuplot

# Installation
1. Clone the repository or downlnad the source code
2. Install LAPACK library, following the instruction given [here](https://www.netlib.org/lapack/)
3. Compile the program using the Fortran 90 compiler
   ```
   gfortran huckel.f90 -o huckel -llapack

4. Execute the program 
   ```
   ./huckel

  # Usage
  1. Insert the parameters of the simulation in the input file, following the example given in `huckel.in` file
  2. Compile and execute the program, as explained in Installation section
  3. The program will initiate the calculations, providing the output file (huckel.out) containing the Huckel matrix, the eigenvectors and the eigenvalues for the chosen system and the files (eigenvectors.dat and eigenvalues.dat) for the plots.

# Licence
This program is released under the GPL-3.0 Licence.

# Contact
For inquiries, bug reports, or suggestions, please reach out to martina.colucci@studenti.unipg.it
