#include "interp1.h"


void Interp1Wrapper::interp1(std::vector<double>* Y, std::vector<int>* X, std::vector<double>* V, std::vector<int>* Xq, std::string method) {

	auto X_length = X->size();
	auto Xq_length = Xq->size();
	std::vector<double> x_i(X->begin(), X->end());
	std::vector<double> y_i(V->begin(), V->end());
	std::vector<int>* x_o;
	std::vector<double>* y_o;
	std::vector<int> ParaSet(Xq_length);
	std::vector<int> ParaSetIndex(Xq_length);

	x_o = Xq;
	y_o = Y;

	int initIndex = 1;
	for (int i = 0; i < Xq_length; i++) {
		float CurParaSetCur = x_i[initIndex];
		if (i >= CurParaSetCur) {
			initIndex = initIndex + 1;
		}
		if (initIndex >= X_length) {
			initIndex = X_length - 1;
		}

		ParaSet[i] = x_i[initIndex - 1] + 1;
		ParaSetIndex[i] = initIndex - 1;
	}

	SplineInterpolation(X_length, Xq_length, &x_i, &y_i, x_o, y_o, &ParaSet, &ParaSetIndex, 0);

}

// the FindPositionID is the function to find the index i where xn >= x_i.at(i)
// the UserPara the para to send to FindPositionID
// the input of xi must be sorted
void Interp1Wrapper::SplineInterpolation(int InputDataLen, int OutputDataLen, std::vector<double>* x_i, std::vector<double>* y_i, std::vector<int>* x_o, std::vector<double>* y_o, std::vector<int>* ParaSet, std::vector<int>* ParaSetIndex, double UserPara)
{
	std::vector<double> An(InputDataLen);
	std::vector<double> Bn(InputDataLen);
	std::vector<double> Cn(InputDataLen);
	std::vector<double> Dn(InputDataLen);

	// initilize
	for (int i = 0; i < InputDataLen; i++)
	{
		An.at(i) = y_i->at(i);
	}

	// sove linear equations
	std::vector<std::vector<double>> A;
	for (int i = 0; i < InputDataLen; i++) {
		A.push_back(std::vector<double>(InputDataLen));
	}
	std::vector<double> b(InputDataLen);

	ConstractSplineLinearEquations(InputDataLen, &A, &b, x_i, y_i);

	// solve Cn
	SolveSplineEquations(InputDataLen, &A, &b, &Cn);

	// solve Bn and Dn
	for (int i = 0; i < InputDataLen - 1; i++)
	{
		double hi = x_i->at(i + 1) - x_i->at(i);
		Bn.at(i) = 1 / hi * (An.at(i + 1) - An.at(i)) - hi / 3 * (2 * Cn.at(i) + Cn.at(i + 1));
		Dn.at(i) = (Cn.at(i + 1) - Cn.at(i)) / (3 * hi);
	}

	// yn(i) = An(ParaSet) + Bn(ParaSet)*(curx-xi) + Cn(ParaSet)*(curx-xi)^2 + Dn(ParaSet)*(curx-xi)^3;
	InterplotData(InputDataLen, OutputDataLen, x_i, x_o, y_o, &An, &Bn, &Cn, &Dn, ParaSet, ParaSetIndex, UserPara);

}



void Interp1Wrapper::InterplotData(int InputDataLen, int OutputDataLen, std::vector<double>* x_i, std::vector<int>* x_o, std::vector<double>* y_o, std::vector<double>* An, std::vector<double>* Bn, std::vector<double>* Cn, std::vector<double>* Dn, std::vector<int>* ParaSet, std::vector<int>* ParaSetIndex, double UserPara)
{
	int CurParaSet = 0;

	for (int i = 0; i < OutputDataLen; i++)
	{
		double xi = ParaSet->at(i);
		int CurParaSet = ParaSetIndex->at(i);
		int xn = i + 1;
		y_o->at(i) = An->at(CurParaSet) + Bn->at(CurParaSet) * (xn - xi) + Cn->at(CurParaSet) * (xn - xi) * (xn - xi) + Dn->at(CurParaSet) * (xn - xi) * (xn - xi) * (xn - xi);

	}
}

void Interp1Wrapper::ConstractSplineLinearEquations(int InputDataLen, std::vector<std::vector<double>>* A, std::vector<double>* b, std::vector<double>* x_i, std::vector<double>* y_i)
{
	std::vector<std::vector<double>>* pA = A; // for parameter array

														  // 
	for (int r = 0; r < InputDataLen; r++)
	{
		for (int c = 0; c < InputDataLen; c++)
		{
			pA->at(r).at(c) = 0;
		}
	}
	// construct A

	pA->at(0).at(0) = 1;
	pA->at(InputDataLen - 1).at(InputDataLen - 1) = 1;


	for (int i = 1; i < InputDataLen - 1; i++)
	{
		double h0 = x_i->at(i) - x_i->at(i - 1);
		double h1 = x_i->at(i + 1) - x_i->at(i);
		pA->at(i).at(i - 1) = h0;
		pA->at(i).at(i) = 2 * (h0 + h1);
		pA->at(i).at(i + 1) = h1;
	}

	// construct b
	b->at(0) = 0;
	b->at(InputDataLen - 1) = 0;

	for (int i = 1; i < InputDataLen - 1; i++)
	{
		double h0 = x_i->at(i) - x_i->at(i - 1);
		double h1 = x_i->at(i + 1) - x_i->at(i);

		double a0 = y_i->at(i - 1);
		double a1 = y_i->at(i);
		double a2 = y_i->at(i + 1);

		b->at(i) = 3 / h1 * (a2 - a1) - 3 / h0 * (a1 - a0);
	}
}



void Interp1Wrapper::SolveSplineEquations(int DataLen, std::vector<std::vector<double>>* A, std::vector<double>* b, std::vector<double>* x)
{
	std::vector<double> y(DataLen);
	std::vector<double> Beta(DataLen);

	std::vector<std::vector<double>>* pA = A; // for parameter array


	for (int i = 0; i < DataLen - 1; i++)
	{
		double ci = pA->at(i).at(i + 1);
		double bi = pA->at(i).at(i);

		if (i == 0)
		{
			Beta.at(i) = ci / bi;
		}
		else
		{
			double ai = pA->at(i).at(i - 1);
			Beta.at(i) = ci / (bi - ai * Beta.at(i - 1));
		}
	}

	for (int i = 0; i < DataLen; i++)
	{
		double bi = pA->at(i).at(i);

		if (i == 0)
		{
			y.at(i) = b->at(i) / bi;
		}
		else
		{
			double ai = pA->at(i).at(i - 1);
			y.at(i) = (b->at(i) - ai * y.at(i - 1)) / (bi - ai * Beta.at(i - 1));
		}
	}

	for (int i = DataLen - 1; i >= 0; i--)
	{
		if (i == DataLen - 1)
		{
			x->at(i) = y.at(i);
		}
		else
		{
			x->at(i) = y.at(i) - Beta.at(i) * x->at(i + 1);
		}
	}
}

// find the index i where xn >= x_i->at(i)
int Interp1Wrapper::FindPositionID_0(std::vector<double>*x_i, double xn, int InputDataLen, double InterplotGap)
{
#define X0_i_FirstPixel		1.0f

	int i = xn - X0_i_FirstPixel;

	return i;
}



// find the index i where xn >= x_i->at(i)
int Interp1Wrapper::FindPositionID_search(std::vector<double>*x_i, double xn, int InputDataLen, double UserPara)
{
	int i;
	for (i = 0; i < InputDataLen; i++)
	{
		if (xn >= x_i->at(i))
		{
			break;
		}
	}

	return i;
}


