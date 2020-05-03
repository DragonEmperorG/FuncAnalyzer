#include "findpeaks.h"
#include <algorithm>
#include <limits>
#include <numeric>


FindpeaksWrapper::FindpeaksWrapper() {
}

FindpeaksWrapper::~FindpeaksWrapper() {
}

std::vector<int> FindpeaksWrapper::findpeaks(const std::vector<double>* Yin) {

	std::vector<double> y(Yin->begin(), Yin->end());
	std::vector<int> x;
	for (auto i = 1; i <= Yin->size(); i++) {
		x.push_back(i);
	}

	// find indices of all finiteand infinite peaksand the inflection points
	std::vector<int> iFinite;
	std::vector<double> iInfinite;
	std::vector<int> iInflect;
	getAllPeaks(&iFinite, &iInfinite, &iInflect, &y);

	// keep only the indices of finite peaks that meet the required
	// minimum heightand threshold
	double minH = -std::numeric_limits<double>::infinity();
	std::string	refW = "halfprom";
	std::vector<int> iPk;
	std::vector<int> iPkOut;
	iPk = removePeaksBelowMinPeakHeight(&y, &iFinite, minH, refW);

	double minT = 0;
	iPkOut = removePeaksBelowThreshold(&y, &iPk, minT);
	iPk = iPkOut;

	// combine finiteand infinite peaks into one list
	iPkOut.clear();
	combinePeaks(&iPkOut, &iPk, &iInfinite);
	iPk = iPkOut;

	// find the indices of the largest peaks within the specified distance
	std::vector<int> idx;
	std::vector<int> idxOut;
	idx = findPeaksSeparatedByMoreThanMinPeakDistance(&y, &x, &iPk, min_peak_distance_);

	// re - order and bound the number of peaks based upon the index vector
	idxOut = orderPeaks(&y, &iPk, &idx, "none");
	idx = idxOut;
	idxOut = keepAtMostNpPeaks(&idx, Yin->size());
	idx = idxOut;

	// use the index vector to fetch the correct peaks.
	auto idx_length = idx.size();
	iPkOut = std::vector<int>(idx_length, 0);
	for (int i = 0; i < idx_length; i++) {
		iPkOut.at(i) = iPk.at(idx.at(i) - 1);
	}
	iPk = iPkOut;

	return assignOutputs(&y, &x, &iPk, false, false);
}

bool FindpeaksWrapper::validateInputParameters(const std::vector<double>* x) {
	return false;
}

void FindpeaksWrapper::getAllPeaks(std::vector<int>* iPk, std::vector<double>* iInf, std::vector<int>* iInflect, const std::vector<double>* y) {

	// fetch indices all infinite peaks
	

	// temporarily remove all + Inf values
	std::vector<double> yTemp(y->begin(), y->end());

	// determine the peaks and inflection points of the signal
	findLocalMaxima(iPk, iInflect, &yTemp);
}

void FindpeaksWrapper::findLocalMaxima(std::vector<int>* iPk, std::vector<int>* iInflect, std::vector<double>* yTemp) {
	// bookend Y by NaNand make index vector
	yTemp->push_back(NAN);
	yTemp->insert(yTemp->begin(), NAN);

	auto yTemp_length = yTemp->size();
	std::vector<int> iTemp(yTemp_length, 0);
	for (int i = 0; i < yTemp_length; i++)	{
		iTemp.at(i) = i + 1;
	}
	
	// keep only the first of any adjacent pairs of equal values(including NaN).
	std::vector<int> yFinite(yTemp_length, 0);
	for (int i = 0; i < yTemp_length; i++) {
		if (isnan(yTemp->at(i))) {
			yFinite.at(i) = 0;
		} else {
			yFinite.at(i) = 1;
		}
	}

	std::vector<int> yTemp_intermediate;
	std::vector<int> yFinite_intermediate;
	std::vector<int> iNeq;
	for (int i = 0; i < (yTemp_length - 1); i++) {
		if (yTemp->at(i) == yTemp->at(i + 1)) {
			yTemp_intermediate.push_back(0);
		} else {
			yTemp_intermediate.push_back(1);
		}
		if ((yFinite.at(i) == 0) && (yFinite.at(i + 1) == 0)) {
			yFinite_intermediate.push_back(0);
		}
		else {
			yFinite_intermediate.push_back(1);
		}
		if ((yTemp_intermediate.at(i) == 1) && (yFinite_intermediate.at(i) == 1)) {
			iNeq.push_back(i + 2);
		}
	}
	iNeq.insert(iNeq.begin(), 1);

	std::vector<int> iTemp_output;
	for (int i = 0; i < iNeq.size(); i++) {
		iTemp_output.push_back(iTemp.at(iNeq.at(i) - 1));
	}
	iTemp = iTemp_output;

	// take the sign of the first sample derivative
	auto iTemp_length_minus_1 = iTemp.size() - 1;
	std::vector<double> s(iTemp_length_minus_1, 0);
	for (int i = 0; i < iTemp_length_minus_1; i++) {
		double yTemp_n = yTemp->at(iTemp.at(i + 1) - 1);
		double yTemp_n_minus_1 = yTemp->at(iTemp.at(i) - 1);
		if (isnan(yTemp_n) || isnan(yTemp_n_minus_1)) {
			s.at(i) = NAN;
		}
		else {
			double diff = yTemp_n - yTemp_n_minus_1;
			if (diff > 0) {
				s.at(i) = 1;
			}
			else if (diff == 0) {
				s.at(i) = 0;
			} else {
				s.at(i) = -1;
			}
		}
	}

	// find local maxima
	std::vector<int> iMax;
	for (int i = 1; i < s.size(); i++) {
		double s_n = s.at(i);
		double s_n_minus_1 = s.at(i - 1);
		double diff = s_n - s_n_minus_1;
		if (diff < 0) {
			iMax.push_back(i + 1);
		}
	}

	// find all transitions from rising to falling or to NaN
	std::vector<int> iAny;
	for (int i = 1; i < s.size(); i++) {
		double s_n = s.at(i);
		double s_n_minus_1 = s.at(i - 1);
		if (isnan(s_n) || isnan(s_n_minus_1)) {
			iAny.push_back(i + 1);
		}
		else {
			if (s_n != s_n_minus_1) {
				iAny.push_back(i + 1);
			}
		}
	}

	// index into the original index vector without the NaN bookend.
	for (int i = 0; i < iAny.size(); i++) {
		iInflect->push_back(iTemp.at(iAny.at(i) - 1) - 1);
	}
	for (int i = 0; i < iMax.size(); i++) {
		iPk->push_back(iTemp.at(iMax.at(i) - 1) - 1);
	}

}

std::vector<int> FindpeaksWrapper::removePeaksBelowMinPeakHeight(const std::vector<double>* Y, const std::vector<int>* iPk, const double Ph, const std::string widthRef) {
	
	if (!iPk->empty()) {
		auto iPk_length = iPk->size();
		std::vector<int> iPK_out(iPk_length, 0);
		for (int i = 0; i < iPk_length; i++) {
			if (Y->at(iPk->at(i)) > Ph) {
				iPK_out.at(i) = iPk->at(i);
			}
		}
		return iPK_out;
	}
	
	return std::vector<int>(iPk->begin(), iPk->end());
}

std::vector<int> FindpeaksWrapper::removePeaksBelowThreshold(const std::vector<double>* Y, const std::vector<int>* iPk, const double Th) {
	
	std::vector<double> base;
	auto iPk_length = iPk->size();
	for (int i = 0; i < iPk_length; i++) {
		base.push_back(std::max(Y->at(iPk->at(i) - 1 - 1), Y->at(iPk->at(i) - 1 + 1)));
	}

	std::vector<int> iPk_out;
	for (int i = 0; i < iPk_length; i++) {
		if ((Y->at(iPk->at(i) - 1) - base.at(i)) >= Th) {
			iPk_out.push_back(iPk->at(i));
		}
	}

	return iPk_out;
}

void FindpeaksWrapper::combinePeaks(std::vector<int>* iPkOut, const std::vector<int>* iPk, const std::vector<double>* iInf) {
	for (auto i = 0; i < iPk->size(); i++) {
		iPkOut->push_back(iPk->at(i));
	}
}

std::vector<int> FindpeaksWrapper::findPeaksSeparatedByMoreThanMinPeakDistance(const std::vector<double>* y, const std::vector<int>* x, const std::vector<int>* iPk, const double Pd) {
	
	// Start with the larger peaks to make sure we don't accidentally keep a
	// small peakand remove a large peak in its neighborhood.

	// copy peak valuesand locations to a temporary place
	auto iPk_length = iPk->size();
	std::vector<double> pks(iPk_length);
	std::vector<int> locs(iPk_length);
	for (auto i = 0; i < iPk_length; i++) {
		pks.at(i) = y->at(iPk->at(i) - 1);
		locs.at(i) = x->at(iPk->at(i) - 1);
	}
    
	// Order peaks from large to small
	std::vector<int> sortIdx(iPk_length);
	std::iota(sortIdx.begin(), sortIdx.end(), 1);
	std::sort(sortIdx.begin(), sortIdx.end(), [pks](int i1, int i2) {return pks[i1 - 1] > pks[i2 - 1]; });

	std::vector<int> locs_temp(iPk_length);
	for (auto i = 0; i < sortIdx.size(); i++) {
		locs_temp.at(i) = locs.at(sortIdx.at(i) - 1);
	}

	std::vector<int> idelete(iPk_length);
	for (auto i = 0; i < iPk_length; i++) {
		if (idelete.at(i) == 0) {
			// If the peak is not in the neighborhood of a larger peak, find
			// secondary peaks to eliminate.
			int locs_temp_i = locs_temp.at(i);
			std::vector<int> locs_temp_intermediate(iPk_length, 0);
			for (int j = 0; j < iPk_length; j++) {
				locs_temp_intermediate.at(j) = (locs_temp.at(j) >= locs_temp_i - Pd) & (locs_temp.at(j) <= locs_temp_i + Pd);
				idelete.at(j) = idelete.at(j) | locs_temp_intermediate.at(j);
			}
			idelete.at(i) = 0;    // Keep current peak
		}
	}

	// report back indices in consecutive order
	std::vector<int> idx;
	for (auto i = 0; i < idelete.size(); i++) {
		if (idelete.at(i) == 0) {
			idx.push_back(sortIdx.at(i));
		}
	}
	std::sort(idx.begin(), idx.end());
	
	return idx;
}

std::vector<int> FindpeaksWrapper::orderPeaks(const std::vector<double>* Y, const std::vector<int>* iPk, const std::vector<int>* idx, const std::string Str) {
	if (idx->empty() || Str == "none") {
		return std::vector<int>(idx->begin(), idx->end());
	}
	
	return std::vector<int>(idx->begin(), idx->end());
}

std::vector<int> FindpeaksWrapper::keepAtMostNpPeaks(const std::vector<int>* idx, const int Np) {

	if (idx->size() > Np) {
		return std::vector<int>(idx->begin(), idx->begin() + Np);
	}

	return std::vector<int>(idx->begin(), idx->end());
}

std::vector<int> FindpeaksWrapper::assignOutputs(const std::vector<double>* y, const std::vector<int>* x, const std::vector<int>* iPk, const bool yIsRow, const bool xIsRow) {
	
	auto iPk_length = iPk->size();
	std::vector<int> XpkOut(iPk_length);
	for (auto i = 0; i < iPk_length; i++) {
			XpkOut.at(i) = iPk->at(i);
	}
	
	return XpkOut;
}




