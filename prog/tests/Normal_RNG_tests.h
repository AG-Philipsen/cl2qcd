#include <cmath>
#include <vector>
#include <numeric>

/**
 * This function calculates the mean of a vector entries
 */
hmc_float mean(std::vector<hmc_float> a){
  return std::accumulate(a.begin(),a.end(),0.0)/a.size();
}

/**
 * This function calculates the variance of a vector entries
 */
hmc_float square_add (hmc_float x, hmc_float y) {return x+y*y;}
hmc_float variance(std::vector<hmc_float> a){
  hmc_float out=std::accumulate(a.begin(),a.end(),0.0,square_add)/a.size();
  out-=mean(a)*mean(a);
  return out;
}

/**
 * This function calculates the standard deviation of a vector entries
 */
hmc_float std_deviation(std::vector<hmc_float> a){
  return sqrt(variance(a));
}

/**
 * This function makes the mean test on a single set of samples. The average of a
 * sequence of independent normal values is itself a normal random value with the same mean.
 * So if we average a million samples from a normal with mean 4 and standard deviation 3,
 * the average will have a normal distribution with mean 4.
 * But the standard deviation of the mean is smaller than the standard
 * deviation of the individual samples by a factor of 1/√n, where n is the number of samples.
 * So in this example the average of our samples will have standard deviation 3/√10^6 = 0.003. An
 * important rule of thumb about normal distributions is that samples will lie within 2 standard
 * deviations of the mean about 95% of the time. This function calculates the mean of
 * the set of samples, then calculate the theoretical standard deviation of the mean
 * (mean_std_dev). This function returns true if the mean of the data is between mu-num_sigma*mean_std_dev
 * and mu+num_sigma*mean_std_dev
 * @param a The set of data that are supposed to be normally distributed
 * @param mu The mean of the theoretical distribution according which the data have been drawn
 * @param num_sigma The number of standard deviation to take into account in the test
 */
bool mean_test_single_set(vector<hmc_float> a, hmc_float num_sigma, hmc_float mu=0, hmc_float sigma=1){
  hmc_float data_mean=mean(a);
  hmc_float mean_std_dev=sigma/sqrt(a.size());
  if(data_mean>=(mu-num_sigma*mean_std_dev) && data_mean<=(mu+num_sigma*mean_std_dev))
    return true;
  else{
    if(logger.beDebug()){
      logger.trace() << "mean_test_single_set failed since " << mu-num_sigma*mean_std_dev << "<=" << data_mean << "<=" << mu+num_sigma*mean_std_dev << " is not true!";
    }
    return false;
  }
}

/**
 * This function repeats the mean test on an ensemble of set of samples. We know that the 
 * mean_test_single_set should success with a probability equal to the integral between
 * mu-num_sigma*sigma and mu+num_sigma*sigma of the normal distribution used to draw data.
 * This integral is equal to Erf[num_sigma/sqrt(2)]. This function print to the shell
 * how often the mean_test_single_set has passed together with the theoretical frequency.
 */
void mean_test_multiple_set(vector<vector<hmc_float>> samples, hmc_float num_sigma, hmc_float mu=0, hmc_float sigma=1){
  hmc_float success_exp;
  hmc_float success_teo=erf(num_sigma/sqrt(2))*100;
  for(uint i=0, k=0; i<samples.size(); i++){
    if(mean_test_single_set(samples[i],num_sigma,mu,sigma))
      k++;
    else{
      if(logger.beDebug())
	logger.trace() << "The " << i << "th mean_test_single_set failed!";
    }
    if(i==samples.size()-1)
      success_exp=((hmc_float)k)/samples.size()*100;
  }
  logger.info() << "The mean_test_single_set up to " << num_sigma << "sigma has passed the " <<setprecision(8) <<  success_exp << "% of the times.";
  logger.info() << "It should pass the " << setprecision(8) << success_teo << "% of the times.";
  //The following tests are meaningfull with a big number of sets of samples (bigger than 2000)
  if(num_sigma<1.9)
    BOOST_CHECK(fabs(success_exp-success_teo)<2);
  else if(num_sigma>=1.9 && num_sigma<2.9)
    BOOST_CHECK(fabs(success_exp-success_teo)<1);
  else if(num_sigma>=2.9 && num_sigma<3.9)
    BOOST_CHECK(fabs(success_exp-success_teo)<0.25);
  else 
    BOOST_CHECK(fabs(success_exp-success_teo)<0.1);
}



 /**
 * This function makes the variance test on a single set of samples in a similar way in which
 * we tested the mean. Suppose an RNG produces values from a normal distribution with variance
 * σ^2. Let S^2 be the sample variance based on n values from the RNG. If n is very large,
 * then S2 approximately has a normal distribution with mean σ^2 and variance 2σ^4/(n−1).
 * As before, we apply the idea that anything with a normal distribution will lie within two
 * standard deviations of its mean 95% of the time (and in general the probability is equal to
 * Erf[num_sigma/sqrt(2)]). This function calculates the variance of the set of samples,
 * then calculate the theoretical standard deviation of the variance (variance_std_dev).
 * This function returns true if the variance of the data is between sigma^2-num_sigma*variance_std_dev
 * and mu+num_sigma*variance_std_dev
 * @param a The set of data that are supposed to be normally distributed
 * @param mu The mean of the theoretical distribution according which the data have been drawn
 * @param num_sigma The number of standard deviation to take into account in the test
 */
bool variance_test_single_set(vector<hmc_float> a, hmc_float num_sigma, hmc_float sigma=1){
  hmc_float data_variance=variance(a);
  hmc_float variance_std_dev=sqrt(2./(a.size()-1))*sigma*sigma;
  if(data_variance>=(sigma*sigma-num_sigma*variance_std_dev) && 
     data_variance<=(sigma*sigma+num_sigma*variance_std_dev))
    return true;
  else{
    if(logger.beDebug()){
      logger.trace() << "variance_test_single_set failed since " << sigma*sigma-num_sigma*variance_std_dev << "<=" << data_variance << "<=" << sigma*sigma+num_sigma*variance_std_dev << " is not true!";
    }
    return false;
  }
}
  
/**
 * This function repeats the variance test on an ensemble of set of samples. We know that the 
 * variance_test_single_set should success with a probability equal to the integral between
 * mu-num_sigma*sigma and mu+num_sigma*sigma of the normal distribution used to draw data.
 * This integral is equal to Erf[num_sigma/sqrt(2)]. This function print to the shell
 * how often the variance_test_single_set has passed together with the theoretical frequency.
 */
void variance_test_multiple_set(vector<vector<hmc_float>> samples, hmc_float num_sigma, hmc_float sigma=1){
  hmc_float success_exp;
  hmc_float success_teo=erf(num_sigma/sqrt(2))*100;
  for(uint i=0, k=0; i<samples.size(); i++){
    if(variance_test_single_set(samples[i],num_sigma,sigma))
      k++;
    else{
      if(logger.beDebug())
	logger.trace() << "The " << i << "th variance_test_single_set failed!";
    }
    if(i==samples.size()-1)
      success_exp=((hmc_float)k)/samples.size()*100;
  }
  logger.info() << "The variance_test_single_set up to " << num_sigma << "sigma has passed the " <<setprecision(8) <<  success_exp << "% of the times.";
  logger.info() << "It should pass the " << setprecision(8) << success_teo << "% of the times.";
  //The following tests are meaningfull with a big number of sets of samples (bigger than 2000)
  if(num_sigma<1.9)
    BOOST_CHECK(fabs(success_exp-success_teo)<2);
  else if(num_sigma>=1.9 && num_sigma<2.9)
    BOOST_CHECK(fabs(success_exp-success_teo)<1);
  else if(num_sigma>=2.9 && num_sigma<3.9)
    BOOST_CHECK(fabs(success_exp-success_teo)<0.5);
  else 
    BOOST_CHECK(fabs(success_exp-success_teo)<0.25);
  //Note that "Typically, sample variances will be more variable than the normal
  //approximation predicts. This means our tests will fail more often than predicted. But
  //as before, we may not need to be too careful. Coding errors are likely to cause the tests to fail
  //every time, not just a little more often than expected. Also, the tests in the following sections
  //do a more careful job of testing the distribution of the samples and have a solid theoretical
  //justification."
}
  