#include "envelope.h"
#include "findpeaks.h"
#include "interp1.h"
#include <numeric>
#include <iostream>


std::vector<double> EnvelopeWrapper::envelope(const std::vector<double> x, const int np, const std::string method) {
	if (!validateInputParameters(&x, np, method)) {
		return std::vector<double>();
	}

	auto x_size = x.size();

	std::vector<double> y_upper(x_size);
	if (method == "peak") {
		envPeak(&y_upper, &x, np);
	}

	return y_upper;

}

bool EnvelopeWrapper::validateInputParameters(const std::vector<double> *x, const int n, const std::string method) {
	if ((n <= 0) || (method != "peak")) {
		return false;
	}
	return true;
}

void EnvelopeWrapper::envPeak(std::vector<double>* y_upper, const std::vector<double>* x, const int n) {

	const auto nx = x->size();

	// handle default case where not enough input is given
	if (nx < 2) {
		y_upper->push_back(x->at(0));
		return;
	}

    // compute upper envelope
	std::vector<int> iPk;
	if (nx > (n + 1)) {
		FindpeaksWrapper findpeaks_wrapper(n);
		iPk = findpeaks_wrapper.findpeaks(x);
	} 

	std::vector<int> iLocs = std::vector<int>(iPk.begin(), iPk.end());
	if (iPk.size() < 2) {
		iLocs.push_back(nx);
		iLocs.insert(iLocs.begin(), 1);
	}

	auto iLoc_length = iLocs.size();
	std::vector<double> x_intermediate(iLoc_length);
	for (int i = 0; i < iLoc_length; i++) {
		x_intermediate.at(i) = x->at(iLocs.at(i) - 1);
	}

	std::vector<int> nx_index(nx);
	std::iota(nx_index.begin(), nx_index.end(), 1);

	std::cout << "iLocs: " << std::endl;
	for (int i = 0; i < iLocs.size(); i++) {
		std::cout << iLocs.at(i) << ", " << x_intermediate.at(i) << std::endl;
	}
	
	Interp1Wrapper interp1_wrapper;
	interp1_wrapper.interp1(y_upper, &iLocs, &x_intermediate, &nx_index, "spline");

	y_upper->erase(y_upper->begin(), y_upper->begin() + 1);
	y_upper->push_back(0);

	std::cout << "y_upper: " << std::endl;
	for (int i = 0; i < y_upper->size(); i++) {
		std::cout << y_upper->at(i) << std::endl;
	}

}
