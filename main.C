#include <iostream>
#include <fstream>
#include <cstdlib>
#include <math.h>
#include "random.h"

using namespace std;

void calc_mean_sigma(int N, double& mean, double& sigma);

//This program is a simple model that predicts what percentage of the population
//of a set of cities will be college graduates, given a parameter X related
//to how much college graduates prefer to live in cities with other college
//graduates.  X<=1 means they don't care at all if there are other college
//graduates (just about the size of the total population), and for X>1, the 
//larger the value of X, the more they are attracted to cities with other 
//college graduates.
//
//Specifically we start with data for total city populations and percentage
//college graduates.
//We remove the entire college graduate population from all of the cities and
//redistribute them one-by-one back into the cities randomly with probably
//proportional to the "effective population" of the city.  The effective
//population is just X*(number of college grads) + (number of non-college
//grads).
//
//Thus X=1 just corresponds to the case where college grads only
//decide based on the current population of the city independent of number
//of college grads already there, X=0 is the case where they only decide
//based on the initial population, and X>1 is the case where college grads
//are weighted more strongly than non-college grads in their decision.
int main(int argc,char *argv[])
{
  //Read parameters from the command line.
  if (argc!=3) {
    cout << "Usage: cities X Nsamp\n";
    return 1;
  }
  const double X=atof(argv[1]);
  const int Nsamp=atoi(argv[2]);

  //Number of cities
  const int ncity=100;
  //City populations from 2010 census from
  //http://en.wikipedia.org/wiki/List_of_U.S._cities_by_population
  //plus some I filled in myself from looking up the city's wiki page.
  //See data/populations_college_order_number_only.txt
  int citypop[ncity]={601723, 945942, 144229, 805235, 233209, 617594, 403892, 790390, 600158, 382578, 608660, 8175133, 112488, 620961, 124775, 420003, 416427, 2695598, 1307402, 97856, 1526006, 408958, 583776, 210565, 181045, 787033, 459787, 731424, 203433, 120083, 337256, 129779, 594833, 204214, 1197816, 3792621, 32736, 197899, 820445, 82825, 520116, 319294, 129272, 601222, 466488, 296943, 545852, 145170, 153060, 305704, 186440, 178874, 49528, 437994, 199110, 178042, 2099451, 205671, 173514, 261310, 399457, 238300, 396815, 579999, 713777, 1445632, 382368, 821784, 58409, 55759, 343829, 229493, 118032, 212237, 188040, 335709, 193524, 741096, 269666, 1327407, 646889, 103190, 391906, 195844, 141527, 287208, 156633, 154305, 167674, 76089, 583756, 494665, 649121, 303871, 66982, 94406, 291707, 201165, 129877, 347483};
  //Fraction of college graduates for year 2010 from 
  //http://www.nytimes.com/interactive/2012/05/31/us/education-in-metro-areas.html
  //See data/fractions_2010.txt
  double frac_college[ncity]={0.468, 0.453, 0.44, 0.434, 0.433, 0.43, 0.41, 0.394, 0.382, 0.379, 0.37, 0.36, 0.352, 0.351, 0.346, 0.341, 0.341, 0.34, 0.337, 0.332, 0.331, 0.33, 0.33, 0.33, 0.327, 0.325, 0.325, 0.322, 0.32, 0.319, 0.319, 0.318, 0.317, 0.317, 0.311, 0.31, 0.309, 0.308, 0.307, 0.301, 0.3, 0.299, 0.298, 0.297, 0.294, 0.293, 0.293, 0.292, 0.291, 0.291, 0.29, 0.288, 0.287, 0.285, 0.285, 0.285, 0.284, 0.283, 0.283, 0.283, 0.281, 0.281, 0.277, 0.276, 0.273, 0.272, 0.271, 0.269, 0.269, 0.269, 0.268, 0.266, 0.266, 0.263, 0.262, 0.262, 0.262, 0.258, 0.256, 0.254, 0.251, 0.251, 0.248, 0.245, 0.244, 0.243, 0.241, 0.234, 0.232, 0.222, 0.216, 0.201, 0.196, 0.195, 0.193, 0.179, 0.177, 0.16, 0.158, 0.15};
  //Use the above to calculate number of college grads in each city.
  int college[ncity];
  int non_college[ncity];
  int ntot_college=0;
  for (int i=0; i<ncity; i++) {
    college[i]=int(frac_college[i]*citypop[i]+0.5);
    ntot_college+=college[i];
    non_college[i]=citypop[i]-college[i];
  }

  //Find the percent of college grads in each city and calculate the mean
  //and sigma percent.
  double percent, real_percent_mean=0.0, real_percent_sigma=0.0;
  for (int i=0; i<ncity; i++) {
    percent=double(college[i])/double(citypop[i]);
    real_percent_mean+=percent;
    real_percent_sigma+=percent*percent;
  }
  calc_mean_sigma(ncity,real_percent_mean,real_percent_sigma);
    
  
  //Main loop.  Remove college grads and redistribute them back amongst the
  //cities according to the rule.  Loop over this Nsamp times to get
  //statistics.
  int college_samp[ncity], i;
  double weight[ncity], sum_weight, x;
  double samp_percent_mean, samp_percent_sigma;
  double mean_percent_mean=0.0, sigma_percent_mean=0.0, mean_percent_sigma=0.0, sigma_percent_sigma=0.0;
  int n_above=0; //Number of samples with percent sigma above real_percent_sigma
  for (int count=0; count<Nsamp; count++) {
    
    //The array non_college gives the number of non-college graduates in 
    //each city.  This is fixed.  We will redistribute the ntot_college
    //college graduates into an array college_samp, which we initialize 
    //to 0.
    for (i=0; i<ncity; i++)
      college_samp[i]=0;
    
    //The array weight gives the probability weight for the distribution
    //of grads into each city.  It is given by
    //weight[i]=non_college[i]+X*college_samp[i]
    //Since we haven't distributed anyone into college_samp yet, weight[i]
    //is initialized to non_college[i].
    sum_weight=0.0;
    for (i=0; i<ncity; i++) {
      weight[i]=double(non_college[i]);
      sum_weight+=weight[i];
    }

    //Loop over the ntot_college grads, distributing them into college_samp.
    for (int grad=0; grad<ntot_college; grad++) {
      
      //Random number 0<=x<sum_weight
      x=ran3()*sum_weight;
      
      //See which city x corresponds to.
      double tmp_sum=0.0;
      for (i=0; i<ncity; i++) {
	tmp_sum+=weight[i];
	if ( x<tmp_sum ) break;
      }
      if (i==ncity) i--;
      
      //i is now the city number the current college grad moves to.
      //Update weight[i], sum_weight, and college_samp[i] appropriately.
      weight[i]+=X;
      sum_weight+=X;
      college_samp[i]++;

    }

    //College grads have been distributed.  Collect statistics.
    samp_percent_mean=0.0;
    samp_percent_sigma=0.0;
    for (i=0; i<ncity; i++) {
      percent=double(college_samp[i])/double(college_samp[i]+non_college[i]);
      samp_percent_mean+=percent;
      samp_percent_sigma+=percent*percent;
    }
    calc_mean_sigma(ncity,samp_percent_mean,samp_percent_sigma);

    mean_percent_mean+=samp_percent_mean;
    mean_percent_sigma+=samp_percent_sigma;
    sigma_percent_mean+=samp_percent_mean*samp_percent_mean;
    sigma_percent_sigma+=samp_percent_sigma*samp_percent_sigma;
    if (samp_percent_sigma>real_percent_sigma) n_above++;
    
  }

  //Get means and sigmas from sums and sums of squares.
  calc_mean_sigma(Nsamp,mean_percent_mean,sigma_percent_mean);
  calc_mean_sigma(Nsamp,mean_percent_sigma,sigma_percent_sigma);

  //Output results
  cout << "The mean percentage of college grads across all cities in the real data is " << real_percent_mean << "\n";
  cout << "The average and standard deviation of this quantity from the simulation is " << mean_percent_mean << " " << sigma_percent_mean << "\n";
  cout << "\n";
  cout << "The standard deviation of the percentage of college grads across all cities in the real data is " << real_percent_sigma << "\n";
  cout << "The average and standard deviation of this quantity from the simulation is " << mean_percent_sigma << " " << sigma_percent_sigma << "\n";
  cout << "\n";
  cout << "The number of simulation samples with greater standard deviation of the percentage is " << n_above << "\n";
  cout << "The total number of simulation samples is " << Nsamp << "\n";
  cout << "Number with greater standard devation / total number of samples = " << double(n_above)/double(Nsamp) << "\n";
  
  return 0;
}

//Calculates the mean and the standard deviation.  "mean" is input
//as sum(x_i,i=1..N) and "sigma" is input as sum(x_i^2,i=1..N), and
//they will be output as the appropriate quantities.
void calc_mean_sigma(int N, double& mean, double& sigma)
{
  mean/=double(N);
  sigma=sqrt(1.0/double(N-1)*(sigma-double(N)*mean*mean));
}
