--------------------------------------------------------------------------------------------------
					Exam PPNM
---------------------------------------------------------------------------------------------------
Name: Anders Bak-Nyhus
studentnumber: 201707103
Assignment: 3. Bilinear interpolation

---------------------------------------------------------------------------------------------------
					Method

There are several ways to do bilinear interpolation. I chose to interpolate along x first, and then along y with the functionvalues attained from the first two interpolations along x.

---------------------------------------------------------------------------------------------------
					Input files

dataPts.txt:
Is the file that has the vectors of points in x and y, where we need to interpolate between them. x is the first collum and y is the second.


dataVals.txt:
Is the file containing the input values for the tabulated function values in the points mentioned in dataPts.txt. At the moment the layout of the values is first take all the values with x = -1 and then take y from -1 to 1, and then change x to -0.5 and do the same again. as is have taken values of x and y between -1 and 1 in intervals of 1/2, there will be 25 function values. The function chosen for this exercise is f(x,y)=x^2-y^2.
I chose to have a second collum of zeros, that are not used in dataVals.txt so as to be able to use the same functions.

---------------------------------------------------------------------------------------------------
					basic code structure
main.c creates a vector and a matrix with data points and functionvalues respectively. This is then sent to a double sum in steps of 0.1, where each x and y will get a functionvalue that has been estimated by the interpolation by first doing a binsearch, to see what vector elements are the closest.
