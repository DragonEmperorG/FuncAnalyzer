#pragma once

#include <string>
#include <vector>


class EnvelopeWrapper {
public:

	EnvelopeWrapper() {

	}

	~EnvelopeWrapper() {

	}	

	std::vector<double> envelope(const std::vector<double> x, const int np, const std::string method);

private:

	bool validateInputParameters(const std::vector<double>* x, const  int n, const std::string method);

	void envPeak(std::vector<double>* y_upper, const std::vector<double>* x, const  int n);
};

