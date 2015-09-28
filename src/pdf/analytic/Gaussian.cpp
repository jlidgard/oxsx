#include <Gaussian.h>
#include <math.h>
#include <PdfExceptions.h>

// Set/Get Parameters are the protected interface from base class Integrable Pdf

Gaussian::Gaussian(size_t nDims_){
    fNDims = nDims_;
}

Gaussian::Gaussian(double mean_, double stdDev_){
    if(stdDev_ <= 0)
        throw ParameterError("Gaussian standard deviation must be greater than 0!");
    fNDims = 1;

    std::vector<double> params;
    params.push_back(mean_);
    params.push_back(stdDev_);
    SetParameters(params);
    
}

Gaussian::Gaussian(const std::vector<double>& means_, const std::vector<double>& stdDevs_){
    if(means_.size() != stdDevs_.size())
        throw DimensionError("Gaussian: Unequal number of means and std devs");

    fNDims = means_.size();
    std::vector<double> params(fNDims * 2, 0);
    for(size_t i = 0; i < means_.size(); i++){
        params.push_back(means_.at(i));
        params.push_back(stdDevs_.at(i));
    }
    SetParameters(params);
}

Gaussian::Gaussian(const Gaussian& other_) : IntegrablePdf(other_){
    fNDims = other_.fNDims;
    SetParameters(other_.GetParameters());
}

double 
Gaussian::operator() (const std::vector<double>& vals_) const{
    if (vals_.size() != fNDims)
        throw DimensionError("Gaussian dimensionality does not match the observable vector passed!");

    if(GetParameters().size() != fNDims)
        throw DimensionError("Gaussian: Attempted to set wrong number of params");

    double exponent = 0;
    for(size_t i = 0; i < fNDims; i++){
        double mean  = GetMean(i);
        double stDev = GetStDev(i);

        double nDevs = (vals_[i] - mean)/stDev;
        exponent += nDevs * nDevs;
    }
    double norm = sqrt(2*M_PI);
    for(size_t i = 0; i < fNDims; i++)
        norm *= GetStDev(i);
    return exp(- 0.5 * exponent) / norm; 
}

double 
Gaussian::Integral(const std::vector<double>& mins_, const std::vector<double>& maxs_) const{
    if(mins_.size() != fNDims || maxs_.size() != fNDims)
        throw DimensionError("Gaussian, tried to integrate over interval of wrong dimensionality");

    double integral = 1;
    for(size_t i = 0; i < mins_.size(); i++)
        integral *= ( Cdf(i, maxs_[i]) - Cdf(i, mins_[i]));
  
    return integral;  
}

double 
Gaussian::Integral(double min_, double max_) const{
    return Integral(std::vector<double>(1, min_), std::vector<double>(1,max_) );
}

Pdf* 
Gaussian::Clone() const{
    return static_cast<Pdf*> (new Gaussian(*this));
}

// double 
// Gaussian::Cdf(size_t dim_, double val_) const{
//     return gsl_cdf_gaussian_P(val_ - GetMean(dim_), GetStDev(dim_));
// }

double 
Gaussian::GetMean(size_t dimension_) const{
    return GetParameter(dimension_ * 2);
}

double 
Gaussian::GetStDev(size_t dimension_) const{
    return GetParameter(dimension_ * 2 + 1);
}
