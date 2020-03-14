#include <cmath>
using namespace std;

void RGBtoXYZ(double R, double G, double B, double &X, double &Y, double &Z)
{
	double sR;
	double sG;
	double sB;

	if (R > 0.04045)
		sR = pow((R + 0.055) / 1.055, 2.4);
	else
		sR = R / 12.92;

	if (G > 0.04045)
		sG = pow((G + 0.055) / 1.055, 2.4);
	else
		sG = G / 12.92;

	if (B > 0.04045)
		sB = pow((B + 0.055) / 1.055, 2.4);
	else
		sB = B / 12.92;

	sR *= 100;
	sG *= 100;
	sB *= 100;

	X = 0.412424 * sR + 0.357570 * sG + 0.180464 * sB;
	Y = 0.212656 * sR + 0.715158 * sG + 0.072186 * sB;
	Z = 0.019332 * sR + 0.119193 * sG + 0.950444 * sB;
}

void XYZtoRGB(double X, double Y, double Z, double &R, double &G, double &B)
{
	double sR;
	double sG;
	double sB;

	sR =  3.24067947165339 * X + -1.53720386341109 * Y + -0.49857053794902 * Z;
	sG = -0.96924887667406 * X +  1.87597888576149 * Y +  0.04155456150049 * Z;
	sB =  0.05563619566875 * X + -0.20399638719236 * Y +  1.05707393711408 * Z;

	sR /= 100;
	sG /= 100;
	sB /= 100;

	if (sR > 0.0031308049535603715170278637770898)
		R = 1.055*pow(sR, 1/2.4) - 0.055;
	else
		R = 12.92 * sR;

	if (sG > 0.0031308049535603715170278637770898)
		G = 1.055*pow(sG, 1/2.4) - 0.055;
	else
		G = 12.92 * sG;

	if (sB > 0.0031308049535603715170278637770898)
		B = 1.055*pow(sB, 1/2.4) - 0.055;
	else
		B = 12.92 * sB;
}

void XYZtoLUV(double X, double Y, double Z, double &L, double &U, double &V)
{
	double XN =  95.0458;
	double ZN = 108.8969;
	double YN = 100.0000;

	double uprimen = 4*XN/(XN + 15*YN + 3*ZN);
	double vprimen = 9*YN/(XN + 15*YN + 3*ZN);

	double uprime = 4*X/(X + 15*Y + 3*Z);
	double vprime = 9*Y/(X + 15*Y + 3*Z);

	Y /= YN;

	if (Y > 0.008856)
		L = 116 * pow(Y, 1/3.0) - 16;
	else
		L = 903.3 * Y;

	if (isfinite(uprime))
		U = 13*L*(uprime - uprimen);
	else
		U = 0.0;

	if (isfinite(vprime))
		V = 13*L*(vprime - vprimen);
	else
		V = 0.0;
}

void LUVtoXYZ(double L, double U, double V, double &X, double &Y, double &Z)
{
	double XN =  95.0458;
	double ZN = 108.8969;
	double YN = 100.0000;

	double uprimen = 4*XN/(XN + 15*YN + 3*ZN);
	double vprimen = 9*YN/(XN + 15*YN + 3*ZN);

	double uprime = U/(13*L) + uprimen;
	double vprime = V/(13*L) + vprimen;

	if (L > 7.9996248)
		Y = pow((L + 16)/116, 3.0);
	else
		Y = L / 903.3;

	Y *= YN;

	if (vprime == 0.0)
		X = 0.0;
	else
		X = 9*Y*uprime/(4*vprime);

	if (uprime == 0.0)
		Z = 0.0;
	else
		Z = (4 - uprime)*X/(3*uprime) - 5*Y;
}

void RGBtoLUV(double R, double G, double B, double &L, double &U, double &V)
{
	double X;
	double Y;
	double Z;

	RGBtoXYZ(R, G, B, X, Y, Z);
	XYZtoLUV(X, Y, Z, L, U, V);
}

void LUVtoRGB(double L, double U, double V, double &R, double &G, double &B)
{
	double X;
	double Y;
	double Z;

	LUVtoXYZ(L, U, V, X, Y, Z);
	XYZtoRGB(X, Y, Z, R, G, B);
}