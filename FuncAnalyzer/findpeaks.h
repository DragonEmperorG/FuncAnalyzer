#pragma once
#include <vector>
#include <string>
class FindpeaksWrapper {
public:

	FindpeaksWrapper();

	FindpeaksWrapper(double min_peak_distance) {
		this->min_peak_distance_ = min_peak_distance;
	}

	~FindpeaksWrapper();

	std::vector<int> findpeaks(const std::vector<double>* x);

private:

	double min_peak_distance_ = 0.0;

	bool validateInputParameters(const std::vector<double>* x);

	void getAllPeaks(std::vector<int>* iPk, std::vector<double>* iInf, std::vector<int>* iInflect, const std::vector<double>* y);

	void findLocalMaxima(std::vector<int>* iPk, std::vector<int>* iInflect, std::vector<double>* yTemp);
	
	std::vector<int> removePeaksBelowMinPeakHeight(const std::vector<double>* Y, const std::vector<int>* iPk, const double Ph, const std::string widthRef);
	
	std::vector<int> removePeaksBelowThreshold(const std::vector<double>* Y, const std::vector<int>* iPk, const double Th);

	void combinePeaks(std::vector<int>* iPkOut, const std::vector<int>* iPk, const std::vector<double>* iInf);

	std::vector<int> findPeaksSeparatedByMoreThanMinPeakDistance(const std::vector<double>* y, const std::vector<int>* x, const std::vector<int>* iPk, const double Pd);

	std::vector<int> orderPeaks(const std::vector<double>* Y, const std::vector<int>* iPk, const std::vector<int>* idx, const std::string Str);

	std::vector<int> keepAtMostNpPeaks(const std::vector<int>* idx, const int Np);

	std::vector<int> assignOutputs(const std::vector<double>* y, const std::vector<int>* x, const std::vector<int>* iPk, const bool yIsRow, const bool xIsRow);
};

