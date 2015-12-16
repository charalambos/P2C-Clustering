////////////////////////////////////////////////////////////////////////////////////
// Copyright © Charalambos "Charis" Poullis, charalambos@poullis.org    	  //
// This work can only be used under an exclusive license of the author.           //
////////////////////////////////////////////////////////////////////////////////////
#ifndef __GAUSSIAN_DISTRIBUTION_H__
#define __GAUSSIAN_DISTRIBUTION_H__

#include <vector>

#include "../Vector.h"
#include "../Matrix.h"

using namespace std;

typedef enum DISTANCE_METRIC {
	KULLBACK_LEIBLER_DISTANCE,
	X2_DISTANCE,
   	DELTA_DISTANCE,
    	HELLINGER_MATSUSITA_BHATTACHARYA_DISTANCE,
     	BHATTACHARYA_DISTANCE
} DISTANCE_METRIC;

template<typename T, int D>
class GaussianDistribution	{
	public:
		///Constructor
		GaussianDistribution<T,D>()	{
			means.Zero();
			covariance_matrix.Zero();
			inv_covariance_matrix.Zero();
			unnorm_covariance_matrix.Zero();
			samples.clear();
			sample_diffs.clear();
			determinant = 0.0f;
		}


  		///Constructor
		GaussianDistribution<T,D>(Vector<T,D> const &_means, Matrix<T,D,D> const &_covariance_matrix)	{
			means = _means;
			covariance_matrix = _covariance_matrix;
			for (int i=0;i<D;i++)	{
				if (covariance_matrix(i,i) < EPSILON)	{
					covariance_matrix(i,i) = EPSILON;
				}
			}
			inv_covariance_matrix = covariance_matrix;
			unnorm_covariance_matrix = covariance_matrix;
			if (inv_covariance_matrix.invert() < 0)	{
				std::cout << "Inverse covariance matrix failed" << std::endl;
				std::cout << covariance_matrix << std::endl;
			}
			determinant = (float) max(EPSILON,covariance_matrix.determinant());
		}

		///Destructor
		~GaussianDistribution()	{
		}

		///Given a set of samples of type T and dimensions D
		///this function fits a gaussian of the same dimensions
		void fit(std::vector<Vector<T,D> > const &_samples)	{

			if (_samples.size() == 0) return;

			samples.clear();
			if (_samples.size() == 1)	samples.push_back(_samples[0]);
			else	samples.insert(samples.begin(), _samples.begin(), _samples.end());

			sample_diffs.clear();

			///First, compute the mean of the samples
			means.Zero();
			for (int i=0;i<samples.size();i++)	{
				means += samples[i];
			}
			sum_of_samples = means;
			means /= float(samples.size());

			///Second, use the mean to compute the covariance matrix
			covariance_matrix.identity();
			for (int i=0;i<samples.size();i++)	{
				///Compute the difference of the sample from the mean
				Vector<T,D> diffs = samples[i] - means;
				sample_diffs.push_back(diffs);

				///Compute the values for each element of the matrix
				///Do it symmetrically
				for (unsigned int j=0;j<D;j++)	{
					for (unsigned int k=j;k<D;k++)	{
						float diff = diffs(j)*diffs(k);
						covariance_matrix(j,k) += diff;
						covariance_matrix(k,j) += diff;
					}
				}
			}
			//covariance_matrix /= float(samples.size());
			///Normalize the matrix. This is to avoid cases where a distribution with more samples
			///appears to be different from an identical distribution of less samples because of
			///the above division by the samples.size()
			unnorm_covariance_matrix = covariance_matrix;
			covariance_matrix.normalize();

			///Compute the inverse covariance matrix
			inv_covariance_matrix = covariance_matrix;
			if (inv_covariance_matrix.invert() < 0)	{
				std::cout << "Inverse covariance matrix failed" << std::endl;
				std::cout << covariance_matrix << std::endl;
				std::cout << means << std::endl;
			}

			///Compute the determinant of the covariance matrix
			determinant = (float) max(EPSILON,covariance_matrix.determinant());

			///Set the recompute flag to true
			recompute = true;

			return;
		}

		///This is an update function for processing speed-up. All necessary quantities are kept in memory.
		///So when a new sample is added there is no need to re-fit the distribution on ALL samples
		void update(Vector<T,D> const &sample)	{
			///if it's the first time then initialize
			if (samples.size() == 0)	{
				unnorm_covariance_matrix.identity();
			}
			///Add the new sample to it
			samples.push_back(sample);
			///Add to the sum of the samples
			sum_of_samples += sample;
			means = sum_of_samples/float(samples.size());
			///Compute the difference from the means
			Vector<T,D> diffs = sample - means;
			sample_diffs.push_back(diffs);

			///Compute the values for each element of the matrix
			///Do it symmetrically
			for (unsigned int j=0;j<D;j++)	{
				for (unsigned int k=j;k<D;k++)	{
					float diff = diffs(j)*diffs(k);
					unnorm_covariance_matrix(j,k) += diff;
					unnorm_covariance_matrix(k,j) += diff;
				}
			}

			///Copy it to the covariance matrix
			covariance_matrix = unnorm_covariance_matrix;
			///Normalize
			covariance_matrix.normalize();
			///Compute the inverse covariance matrix
			inv_covariance_matrix = covariance_matrix;
			if (inv_covariance_matrix.invert() < 0)	{
				std::cout << "Inverse covariance matrix failed" << std::endl;
				std::cout << covariance_matrix << std::endl;
				std::cout << means << std::endl;
			}

			///Compute the determinant of the covariance matrix
			determinant = (float) max(EPSILON,covariance_matrix.determinant());

			///Set the recompute flag to true
			recompute = true;

			return;
		}

		///Same as above but for multiple samples
		///This is an update function for processing speed-up. All necessary quantities are kept in memory.
		///So when a new sample is added there is no need to re-fit the distribution on ALL samples
		void update(std::vector<Vector<T,D> >const &_samples)	{
			///if it's the first time then initialize
			if (samples.size() == 0)	{
				unnorm_covariance_matrix.identity();
			}
			for (int i=0;i<_samples.size();i++)	{
				Vector<T,D> sample = _samples[i];
				///Add the new sample to it
				samples.push_back(sample);
				///Add to the sum of the samples
				sum_of_samples += sample;
				means = sum_of_samples/float(samples.size());
				///Compute the difference from the means
				Vector<T,D> diffs = sample - means;
				sample_diffs.push_back(diffs);

				///Compute the values for each element of the matrix
				///Do it symmetrically
				for (unsigned int j=0;j<D;j++)	{
					for (unsigned int k=j;k<D;k++)	{
						float diff = diffs(j)*diffs(k);
						unnorm_covariance_matrix(j,k) += diff;
						unnorm_covariance_matrix(k,j) += diff;
					}
				}
			}

			///Copy it to the covariance matrix
			covariance_matrix = unnorm_covariance_matrix;
			///Normalize
			covariance_matrix.normalize();
			///Compute the inverse covariance matrix
			inv_covariance_matrix = covariance_matrix;
			if (inv_covariance_matrix.invert() < 0)	{
				std::cout << "Inverse covariance matrix failed" << std::endl;
				std::cout << covariance_matrix << std::endl;
				std::cout << means << std::endl;
			}

			///Compute the determinant of the covariance matrix
			determinant = (float) max(EPSILON,covariance_matrix.determinant());

			///Set the recompute flag to true
			recompute = true;

			return;
		}

		///Evaluates the gaussian function at the given sample
		float evaluate(Vector<T,D> const &sample)	{
			///If it's the first time then compute the constant
			///otherwise use the computed value
			if (recompute)	{
 				constant = 1.0f/(pow(2.0f*float(M_PI),float(D)/2.0f) * std::sqrt(determinant));
				recompute = false;
			}

			///Compute the difference from the mean and its transpose
			Vector<T,D> diffs = (sample-means);
			Matrix<T,1,D> diffsT;
			for (unsigned int i=0;i<D;i++)	{
				diffsT(0,i) = diffs(i);
			}

			Vector<float,1> arg_vec = -0.5f*diffsT*(inv_covariance_matrix*diffs);
			float arg = arg_vec(0);

			float val = (float) max(EPSILON,double(constant * exp(arg)));

			/*if (val < float(EPSILON))	{
				std::cout << " The evaluation is less than EPSILON" << std::endl;
			}*/

			return val;
		}

		///Returns the probability of the given sample in this distribution
		///This is the negative log probability value
		float prob(Vector<T,D> const &sample)	{
			///Evaluate the gaussian function at this point
			float eval = evaluate(sample);
			///Compute the negative log probability
			float neg_log_prob = -log(eval);
			return neg_log_prob;
		}

		///Compute the bhattacharya distance between two gaussian distributions
		///If two distributions are identical then the similarity is 0
		///Otherwise is a positive real number
		friend float bhattacharya_distance( GaussianDistribution<T,D> *gd1,
					 	    GaussianDistribution<T,D> *gd2)	{
			///Get the means for the two distributions
			Vector<T,D> means1 = gd1->getMeans();
			Vector<T,D> means2 = gd2->getMeans();
			Vector<T,D> diffs = means2-means1;
			Matrix<T,1,D> diffsT;
			for (unsigned int i=0;i<D;i++)	{
				diffsT(0,i) = diffs(i);
			}

			///Get the covariance matrices for the two distributions
			Matrix<T,D,D> cov_matrix1 = gd1->getCovarianceMatrix();
			Matrix<T,D,D> cov_matrix2 = gd2->getCovarianceMatrix();
			Matrix<T,D,D> mat = 0.5f*(cov_matrix1 + cov_matrix2);
			Matrix<T,D,D> inv_mat = mat;
			if (inv_mat.invert()<0)	{
				std::cout << "Matrix inversion failed." << std::endl;
				std::cout << mat << std::endl;
			}

			///Get the determinants for the covariance matrices
			float det1 = cov_matrix1.determinant();
			float det2 = cov_matrix2.determinant();
			float det_mat = mat.determinant();
//std::cout << "det_mat: " << det_mat << std::endl;
			///Compute the constants
			Vector<float,1> first_part = 0.125f*diffsT*inv_mat*diffs;
			float denom = max(float(EPSILON), std::sqrt(det1*det2));

			float similarity = first_part(0) + 0.5f*log(det_mat/denom);
/*
			///Go through the points and compute the probability in both distributions
			for (unsigned int i=0;i<samples.size();i++)	{
				float eval1 = gd1->evaluate(samples[i]);
				float eval2 = gd2->evaluate(samples[i]);

				switch(distance_type)	{
					case KULLBACK_LEIBLER_DISTANCE: similarity += eval1 * log(eval1/eval2);
									break;
					case X2_DISTANCE: similarity += (eval1-eval2)*(eval1-eval2)/eval2;
							  break;
					case DELTA_DISTANCE: similarity += (eval1-eval2)*(eval1-eval2)/(eval1+eval2);
							     break;
					case HELLINGER_MATSUSITA_BHATTACHARYA_DISTANCE: value = sqrtf(eval1) - sqrtf(eval2);
											similarity += value*value;
											break;
					case BHATTACHARYA_DISTANCE: similarity += sqrtf(eval1*eval2);
								    break;
				}
			}

			switch(distance_type)	{
				case BHATTACHARYA_DISTANCE:	similarity = -log(similarity);
								break;
			}
			if (distance_type == HELLINGER_MATSUSITA_BHATTACHARYA_DISTANCE)	{
				similarity = sqrtf(similarity);
			}*/

			return similarity;
		}

		///Compute the mahalanobis distance between a multivariate vector and a gaussian distribution
		friend float mahalanobis_distance( GaussianDistribution<T,D> *gd,
					 	    Vector<T,D> const &mv)	{
			Matrix<T,D,D> cov_mat =gd->getCovarianceMatrix();
			Matrix<T,D,D> inv_cov_mat = cov_mat;
			if (inv_cov_mat.invert()<0)	{
				std::cout << "Matrix inversion failed." << std::endl;
			}

			Vector<T,D> diffs = mv - gd->getMeans();
			Matrix<T,1,D> diffsT;
			for (unsigned int i=0;i<D;i++)	{
				diffsT(0,i) = diffs(i);
			}

			Vector<float,1> val = diffsT*inv_cov_mat*diffs;
			float distance = std::sqrt(val(0));

			return distance;
		}

		Vector<T,D> getMeans()	{
			return means;
		}

		///HACKING function. DO NOT USE unless π.χ.
		void setMeans(Vector<T,D> const &_means)	{
			means = _means;
		}

		Matrix<T,D,D> getCovarianceMatrix()	{
			return covariance_matrix;
		}

		std::vector<Vector<T,D> > getSamples()	{
			return samples;
		}

	private:
		///The mean values
		Vector<T,D>	means;
		///The covariance matrix
		Matrix<T,D,D>	covariance_matrix;
		///The inverse covariance matrix (kept for speed optimization)
		Matrix<T,D,D> inv_covariance_matrix;
		///Flag that determines whether the constant should be recomputed
		bool recompute;
		///The determinant of the covariance matrix
		float determinant;
		///The constant
		float constant;

		///Private information used for speeding up the process
		///The samples
		std::vector<Vector<T,D> > samples;
		///The diffs of the samples from the means
		std::vector<Vector<T,D> > sample_diffs;
		///The unnormalized covariance matrix
		Matrix<T,D,D> unnorm_covariance_matrix;
		///The sum of the samples
		Vector<T,D> sum_of_samples;
};

typedef	GaussianDistribution<float,1>	GaussianDistribution1f;
typedef	GaussianDistribution<float,2>	GaussianDistribution2f;
typedef	GaussianDistribution<float,3>	GaussianDistribution3f;
typedef	GaussianDistribution<float,4>	GaussianDistribution4f;
typedef	GaussianDistribution<float,5>	GaussianDistribution5f;
typedef	GaussianDistribution<float,6>	GaussianDistribution6f;
typedef	GaussianDistribution<float,7>	GaussianDistribution7f;

#endif
