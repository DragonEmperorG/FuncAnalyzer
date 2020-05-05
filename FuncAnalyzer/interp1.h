#pragma once
#include <vector>
#include <string>


class Interp1Wrapper {
public:

	Interp1Wrapper() {

	}

	~Interp1Wrapper() {

	}

	void interp1(std::vector<double>* Y, std::vector<int>* X, std::vector<double>* V, std::vector<int>* Xq, std::string method);

private:

	void SplineInterpolation(int InputDataLen, int OutputDataLen, std::vector<double>* x_i, std::vector<double>* y_i, std::vector<int>* x_o, std::vector<double>* y_o, std::vector<int>* ParaSet, std::vector<int>* ParaSetIndex, double UserPara);

	// find the index i where xn >= x_i.at(i)
	int FindPositionID_search(std::vector<double>* x_i, double xn, int InputDataLen, double UserPara);

	int FindPositionID_0(std::vector<double>* x_i, double xn, int InputDataLen, double InterplotGap);

	void InterplotData(int InputDataLen, int OutputDataLen, std::vector<double>* x_i, std::vector<int>* x_o, std::vector<double>* y_o, std::vector<double>* An, std::vector<double>* Bn, std::vector<double>* Cn, std::vector<double>* Dn, std::vector<int>* ParaSet, std::vector<int>* ParaSetIndex, double UserPara);

	void ConstractSplineLinearEquations(int InputDataLen, std::vector<std::vector<double>>* A, std::vector<double>* b, std::vector<double>* x_i, std::vector<double>* y_i);

	void SolveSplineEquations(int DataLen, std::vector<std::vector<double>>* A, std::vector<double>* b, std::vector<double>* x);

};
