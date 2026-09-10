#ifndef LINE_OPT_RESULTS_H
#define LINE_OPT_RESULTS_H

#include <cmath>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <string>
#include <vector>

namespace line { namespace opt {

using Value = std::vector<double>;
using VariableValues = std::map<std::string, Value>;

inline double scalar_value(const Value& v) { return v.empty() ? 0.0 : v[0]; }
inline std::string metric_key(const std::string& station, const std::string& jobclass) {
    return station + "||" + jobclass;
}

class SensitivityData {
public:
    using ParameterMap = std::map<std::string, double>;
    using MetricMap = std::map<std::string, ParameterMap>;
    using Data = std::map<std::string, MetricMap>;

    void add(const std::string& kind, const std::string& metric,
             const std::string& parameter, double value) {
        data_[kind][metric][parameter] += value;
    }
    const MetricMap* for_kind(const std::string& kind) const {
        const auto found = data_.find(kind);
        return found == data_.end() ? nullptr : &found->second;
    }
    bool empty() const { return data_.empty(); }
    const Data& data() const { return data_; }
    static std::string metric_key(const std::string& station,
                                  const std::string& jobclass) {
        return station + "||" + jobclass;
    }
    static std::string parameter_key(const std::string& station,
                                     const std::string& jobclass) {
        return "rate||" + station + "||" + jobclass;
    }
private:
    Data data_;
};

struct EvaluationResult {
    bool feasible = true;
    std::map<std::string,double> response_times, throughputs, queue_lengths, utilizations;
    std::map<std::string,double> system_response_times, system_throughputs;
    double solve_time = 0.0;
    std::string solver_used;
    std::shared_ptr<SensitivityData> sensitivities;

    void set_response_time(const std::string& s,const std::string& c,double v){response_times[metric_key(s,c)]=v;}
    void set_throughput(const std::string& s,const std::string& c,double v){throughputs[metric_key(s,c)]=v;}
    void set_queue_length(const std::string& s,const std::string& c,double v){queue_lengths[metric_key(s,c)]=v;}
    static double exact(const std::map<std::string,double>& m,const std::string& k,double d){auto i=m.find(k);return i==m.end()?d:i->second;}
    static double aggregate(const std::map<std::string,double>& m,const std::string& s,bool mean,double d){
        const std::string p=s+"||"; double total=0; std::size_t n=0;
        for(const auto& kv:m) if(kv.first.compare(0,p.size(),p)==0){total+=kv.second;++n;}
        return n==0?d:(mean?total/n:total);
    }
    double response_time(const std::string&s,const std::string&c="")const{return c.empty()?aggregate(response_times,s,true,std::numeric_limits<double>::infinity()):exact(response_times,metric_key(s,c),std::numeric_limits<double>::infinity());}
    double throughput(const std::string&s,const std::string&c="")const{return c.empty()?aggregate(throughputs,s,false,0):exact(throughputs,metric_key(s,c),0);}
    double queue_length(const std::string&s,const std::string&c="")const{return c.empty()?aggregate(queue_lengths,s,false,0):exact(queue_lengths,metric_key(s,c),0);}
    double utilization(const std::string&s)const{return exact(utilizations,s,0);}
    double system_response_time(const std::string& c="")const{
        if(!c.empty()) return exact(system_response_times,c,std::numeric_limits<double>::infinity());
        if(system_response_times.empty()) return std::numeric_limits<double>::infinity();
        double w=0,x=0; for(const auto& kv:system_response_times){double t=exact(system_throughputs,kv.first,0);w+=kv.second*t;x+=t;}
        if(x>0)return w/x; double s=0;for(const auto&kv:system_response_times)s+=kv.second;return s/system_response_times.size();
    }
    double system_throughput(const std::string& c="")const{if(!c.empty())return exact(system_throughputs,c,0);double x=0;for(const auto&kv:system_throughputs)x+=kv.second;return x;}
};

struct OptimizationResult {
    double objective_value=std::numeric_limits<double>::infinity();
    VariableValues variable_values;
    std::map<std::string,double> constraint_violations;
    bool feasible=false; std::size_t iterations=0,model_evaluations=0; double solve_time=0;
    std::vector<double> convergence_history; std::string terminated_by;
    double total_violation()const{double x=0;for(const auto&kv:constraint_violations)x+=kv.second;return x;}
};

} }
#endif
