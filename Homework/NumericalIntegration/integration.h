#ifndef NUMERICAL_INTEGRATION_INTEGRATION_H
#define NUMERICAL_INTEGRATION_INTEGRATION_H

double adapt24 (  double func(double), double leftEndPt, double rightEndPt, double absAcc, double relAcc, double second_funcVal, double third_funcVal, int numOfRecursions, double* integrationError  );
double adapt (    double func(double), double leftEndPt, double rightEndPt, double absAcc, double relAcc, double* integrationError );
double open_quad( double func(double), double leftEndPt, double rightEndPt, double absAcc, double relAcc, double* integrationError );
double integrate( double func(double), double leftEndPt, double rightEndPt, double absAcc, double relAcc, double* integrationError );

#endif 
