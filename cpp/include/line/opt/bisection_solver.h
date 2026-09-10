#ifndef LINE_OPT_BISECTION_SOLVER_H
#define LINE_OPT_BISECTION_SOLVER_H
#include <chrono>
#include <cmath>
#include <map>
#include <string>
#include <utility>
#include <vector>
#include "line/opt/problem.h"
#include "line/util/error.h"
namespace line { namespace opt {
class BisectionSolver{
 struct Probe{bool feasible;EvaluationResult base;std::map<std::string,double>violations;};
public:BisectionSolver(const OptimizationProblem&p,std::string direction="min_feasible",LineEvaluator::SolveFunction solve={}):problem_(p),direction_(std::move(direction)),solve_(std::move(solve)){
  if(!problem_.validate().empty())throw InputError("BisectionSolver: invalid optimization problem");if(problem_.variables().size()!=1||problem_.variables()[0]->dimension()!=1)throw InputError("BisectionSolver requires exactly one dimension-1 decision variable");var_=problem_.variables()[0];
  const double lo=scalar_value(var_->decode({0})),hi=scalar_value(var_->decode({1}));if(!std::isfinite(lo)||!std::isfinite(hi)||std::fmod(lo,1.0)!=0.0||std::fmod(hi,1.0)!=0.0)throw InputError("BisectionSolver requires an integer-valued variable");lo_=long(std::llround(lo));hi_=long(std::llround(hi));
  evaluators_.emplace_back(problem_.model_variant(),problem_.variables(),problem_.fixed_variables(),solve_);for(const auto&s:problem_.scenarios())evaluators_.emplace_back(s.model,problem_.variables(),problem_.fixed_variables(),solve_);
 }
 OptimizationResult solve(){auto start=std::chrono::steady_clock::now();long lo=lo_,hi=hi_;std::size_t it=0;if(direction_=="min_feasible")while(lo<hi){long m=(lo+hi)/2;if(probe(m).feasible)hi=m;else lo=m+1;++it;}else if(direction_=="max_feasible")while(lo<hi){long m=(lo+hi+1)/2;if(probe(m).feasible)lo=m;else hi=m-1;++it;}else throw InputError("BisectionSolver: unknown direction '"+direction_+"'");auto&p=probe(lo);OptimizationResult out;out.variable_values[var_->name()]={double(lo)};out.feasible=p.feasible;out.constraint_violations=p.violations;out.iterations=it;for(const auto&e:evaluators_)out.model_evaluations+=e.evaluation_count();out.terminated_by="bisection";VariableValues all=fixed_values();all[var_->name()]={double(lo)};out.objective_value=problem_.objective()->evaluate(p.base,all);out.solve_time=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();return out;}
private:std::vector<ConstraintPtr>all_constraints()const{auto c=problem_.objective()->constraints;c.insert(c.end(),problem_.constraints().begin(),problem_.constraints().end());return c;}VariableValues fixed_values()const{VariableValues v;for(const auto&f:problem_.fixed_variables())v[f.first->name()]=f.second;return v;}
 Probe&probe(long value){auto found=cache_.find(value);if(found!=cache_.end())return found->second;VariableValues v;v[var_->name()]={double(value)};VariableValues all=fixed_values();all.insert(v.begin(),v.end());Probe p{true,EvaluationResult(),{}};bool first=true;for(auto&e:evaluators_){auto r=e.evaluate(v);if(first){p.base=r;first=false;}if(!r.feasible){p.feasible=false;continue;}for(const auto&c:all_constraints()){double x=c->evaluate(r,all);if(x>0){p.feasible=false;p.violations[c->name()]=std::max(p.violations[c->name()],x);}}}return cache_.emplace(value,std::move(p)).first->second;}
 const OptimizationProblem&problem_;std::string direction_;LineEvaluator::SolveFunction solve_;VariablePtr var_;long lo_=0,hi_=0;std::vector<LineEvaluator>evaluators_;std::map<long,Probe>cache_;
};
} }
#endif
