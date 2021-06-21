Here i will describe the basics of the code as well as the way to change input files.

Input files:

x_i_y_i.txt is the file that has the four points x1, x2, y1 and y2 as well as the point (x,y), whichs needs to be interpolated.

layout:
x1 y1
x2 y2
x y


F_ijVals.txt is the file containing the input values for the tabulated function values in the points
(x1,y1), (x1,y2), (x2,y1) and (x2,y2)

layout:
F11 F12
F21 F22

code:
the code interpolates along x at first, where it finds the function values at (x,y1) and (x,y2), then it does
a final interpolation along y, to find the function value at (x,y).
