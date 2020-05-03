#include "envelope.h"
#include "findpeaks.h"


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
	if (nx > (n + 1)) {
		FindpeaksWrapper findpeaks_wrapper(n);
		std::vector<int> iPk = findpeaks_wrapper.findpeaks(x);
	} 

	iLocs
	if (true) {

	}
	else {

	}
}
