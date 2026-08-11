#include<iostream>
#include<cmath>
#include<random>
#include<sys/time.h>

std::mt19937 mt_rand(0);
std::uniform_real_distribution<double> d(0.0,1.0);

double rand_u(){
  return d(mt_rand);
}


class SIR{
 private:
  double a;
  double oneOverN;
  int N;
  double S, I, R;
  double averageS, averageI;
  
 public:
  SIR(double a, int N) : a(a), N(N){
    S=.5; I=.5; R=0;
    oneOverN = 1./(double)N;
  }
  void simulate(int T, bool print=false){
    double t=0;
    averageS=0;
    averageI=0;
    if(print) std::cout << t << " " << S << " " << I << "\n";
    while(t<T) {
      double rateSI = S*(1+10*I/(a+S));
      double rateIR = 5*I;
      double rateRI = (10*S+0.1)*R;
      double totalRate = rateSI+rateIR+rateRI; 
      double u = rand_u()*totalRate;
      double delta_t = -log(rand_u()) / totalRate / N; 
      t += delta_t;
      if (t >= 100){
	averageS += S*delta_t / (T-100);
	averageI += I*delta_t / (T-100); 
      }
      
      if (u < rateSI){
	S-=oneOverN; 
	I+=oneOverN; 
      }
      else if (u<rateSI+rateIR){
	I-=oneOverN; 
	R+=oneOverN; 
      }
      else {
	R-=oneOverN; 
	S+=oneOverN; 
      }
      if (print) std::cout << t << " " << S << " " << I << "\n";
    } ;
  }
  void steady_state(){
    simulate(1000,false);
    std::cout << averageS << " " << averageI << "\n";
  }
};

int main(int argc, char**argv){
  struct timeval tv;
  struct timezone tz;
  gettimeofday(&tv, &tz);
  mt_rand.seed(tv.tv_usec);
  
  int N=1000;
  double a = 0.1;
  bool print_only_one_trajectory = false;
  for(int i=1;i<argc;i++){
    switch(* (argv[i]) ){
    case 'N': N=atoi(argv[i]+1); break;
    case 'a': a=atof(argv[i]+1); break;
    case 't': print_only_one_trajectory=true; break;
    default: std::cerr << "unkown option : "<<argv[i]<<"\n";exit(1);
    }
  }
  if (print_only_one_trajectory){
    SIR sir(a,N);
    sir.simulate(10,true);
  }
  else {
    for (int i=0;i<100;i++){
      SIR sir(a,N);
      sir.steady_state();
    }
  }
}
