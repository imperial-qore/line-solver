#ifndef LINE_OPT_VARIABLES_H
#define LINE_OPT_VARIABLES_H

#include <algorithm>
#include <cmath>
#include <memory>
#include <optional>
#include <string>
#include <vector>
#include "line/lang/lqn/lqn_reader.h"
#include "line/lang/qn/network_builder.h"
#include "line/opt/results.h"
#include "line/util/error.h"

namespace line { namespace opt {

class DecisionVariable {
public:
    explicit DecisionVariable(std::string n,std::size_t d=1):name_(std::move(n)),dimension_(d){}
    virtual ~DecisionVariable()=default;
    const std::string& name()const{return name_;} std::size_t dimension()const{return dimension_;}
    virtual Value decode(const std::vector<double>& x)const=0;
    virtual void apply(qn::Network<double>&,const Value&)const {
        throw InputError("line-opt: variable '" + name_ + "' cannot be applied to a flat Network");
    }
    virtual void apply(lqn::LqnModel<double>&,const Value&)const {
        throw InputError("line-opt: variable '" + name_ + "' cannot be applied to a LayeredNetwork");
    }
    virtual std::string type()const=0;
    virtual std::vector<std::string> layers(const lqn::LqnModel<double>&) const { return {}; }
    virtual std::optional<Value> current_value(const lqn::LqnModel<double>&) const {
        return std::nullopt;
    }
    virtual bool supports_sensitivity() const { return false; }
    virtual std::string sensitivity_key(const lqn::LqnModel<double>&) const { return {}; }
    virtual std::map<std::string,std::string> sensitivity_metric_targets(
        const lqn::LqnModel<double>&) const { return {}; }
    virtual double rate_jacobian(const Value&) const { return 0.0; }
    virtual double decode_jacobian(double) const { return 0.0; }
protected:
    static std::size_t node(const qn::NetworkStruct<double>&s,const std::string&n){for(std::size_t i=0;i<s.nodes.size();++i)if(s.nodes[i].name==n)return i+1;throw InputError("line-opt: unknown node '"+n+"'");}
    static std::size_t cls(const qn::NetworkStruct<double>&s,const std::string&n){for(std::size_t i=0;i<s.classes.size();++i)if(s.classes[i].name==n)return i+1;throw InputError("line-opt: unknown class '"+n+"'");}
    std::string name_; std::size_t dimension_;
};

class ServerAllocation final:public DecisionVariable{
public: ServerAllocation(std::string station,int lo,int hi,std::string n=""):DecisionVariable(n.empty()?station+"_servers":n),station_(std::move(station)),lo_(lo),hi_(hi){if(lo>hi)throw InputError("ServerAllocation: invalid bounds");}
 Value decode(const std::vector<double>&x)const override{return {double(std::clamp<int>(int(std::llround(lo_+x.at(0)*(hi_-lo_))),lo_,hi_))};}
 void apply(qn::Network<double>&m,const Value&v)const override{m.set_number_of_servers(node(m.raw_struct(),station_),scalar_value(v));}
 std::string type()const override{return "server_allocation";} private:std::string station_;int lo_,hi_;
};
class StationReplicas final:public DecisionVariable{
public: StationReplicas(std::string station,int lo,int hi,std::string n=""):DecisionVariable(n.empty()?station+"_replicas":n),station_(std::move(station)),lo_(lo),hi_(hi){}
 Value decode(const std::vector<double>&x)const override{return {double(std::clamp<int>(int(std::llround(lo_+x.at(0)*(hi_-lo_))),lo_,hi_))};}
 void apply(qn::Network<double>&m,const Value&v)const override{auto&sn=m.raw_struct();auto nd=node(sn,station_);auto st=sn.nodes.at(nd-1).station;double base=sn.stations.at(st-1).nservers;if(!std::isfinite(base)||base<1)base=1;m.set_number_of_servers(nd,base*scalar_value(v));}
 std::string type()const override{return "station_replicas";} private:std::string station_;int lo_,hi_;
};
class ServiceRate final:public DecisionVariable{
public: ServiceRate(std::string station,std::string jobclass,double lo,double hi,std::string n=""):DecisionVariable(n.empty()?station+"_"+jobclass+"_rate":n),station_(std::move(station)),class_(std::move(jobclass)),lo_(lo),hi_(hi){}
 Value decode(const std::vector<double>&x)const override{return {lo_+x.at(0)*(hi_-lo_)};}
 void apply(qn::Network<double>&m,const Value&v)const override{auto&s=m.raw_struct();m.set_service(node(s,station_),cls(s,class_),lang::Distrib<double>::exp_rate(scalar_value(v)));}
 std::string type()const override{return "service_rate";} private:std::string station_,class_;double lo_,hi_;
};
class JobPopulation final:public DecisionVariable{
public: JobPopulation(std::string c,int lo,int hi,std::string n=""):DecisionVariable(n.empty()?c+"_population":n),class_(std::move(c)),lo_(lo),hi_(hi){}
 Value decode(const std::vector<double>&x)const override{return {double(std::clamp<int>(int(std::llround(lo_+x.at(0)*(hi_-lo_))),lo_,hi_))};}
 void apply(qn::Network<double>&m,const Value&v)const override{auto&s=m.raw_struct();s.classes.at(cls(s,class_)-1).population=scalar_value(v);}
 std::string type()const override{return "job_population";} private:std::string class_;int lo_,hi_;
};
class RoutingProbabilities final:public DecisionVariable{
public: RoutingProbabilities(std::string c,std::string source,std::vector<std::string>targets,std::string n=""):DecisionVariable(n.empty()?c+"_routing_from_"+source:n,std::max<std::size_t>(1,targets.size()-1)),class_(std::move(c)),source_(std::move(source)),targets_(std::move(targets)){if(targets_.empty())throw InputError("RoutingProbabilities: no targets");}
 Value decode(const std::vector<double>&x)const override{Value p(targets_.size());if(p.size()==1){p[0]=1;return p;}double rem=1;for(std::size_t i=0;i+1<p.size();++i){p[i]=rem*x.at(i);rem-=p[i];}p.back()=rem;return p;}
 void apply(qn::Network<double>&m,const Value&v)const override{auto&s=m.raw_struct();auto c=cls(s,class_),src=node(s,source_);for(std::size_t i=0;i<targets_.size();++i)s.set_route(c,c,src,node(s,targets_[i]),v.at(i));}
 std::string type()const override{return "routing";} private:std::string class_,source_;std::vector<std::string>targets_;
};
class ClassPriority final:public DecisionVariable{
public: ClassPriority(std::vector<std::string>c,std::string mode="levels",int lo=1,int hi=10,std::string n="class_priorities"):DecisionVariable(std::move(n),mode=="levels"?c.size():std::max<std::size_t>(1,c.size()-1)),classes_(std::move(c)),mode_(std::move(mode)),lo_(lo),hi_(hi){}
 Value decode(const std::vector<double>&x)const override{Value v;if(mode_=="levels"){for(double z:x)v.push_back(std::round(lo_+z*(hi_-lo_)));return v;}std::vector<std::pair<double,std::size_t>>k;for(std::size_t i=0;i<classes_.size();++i)k.push_back({-(i<x.size()?x[i]:0),i});std::stable_sort(k.begin(),k.end());for(auto&p:k)v.push_back(double(p.second));return v;}
 void apply(qn::Network<double>&m,const Value&v)const override{auto&s=m.raw_struct();if(mode_=="levels")for(std::size_t i=0;i<classes_.size();++i)s.classes.at(cls(s,classes_[i])-1).prio=int(v.at(i));else for(std::size_t rank=0;rank<v.size();++rank)s.classes.at(cls(s,classes_.at(std::size_t(v[rank])))-1).prio=int(v.size()-rank);}
 std::string type()const override{return "class_priority";} private:std::vector<std::string>classes_;std::string mode_;int lo_,hi_;
};

class ClassServiceMapping final : public DecisionVariable {
public:
    ClassServiceMapping(std::string jobclass, std::vector<std::string> stations,
                        std::string name = "")
        : DecisionVariable(name.empty() ? jobclass + "_mapping" : std::move(name)),
          class_(std::move(jobclass)), stations_(std::move(stations)) {}

    Value decode(const std::vector<double>& x) const override {
        if (stations_.empty()) return {-1.0};
        const std::size_t index = std::min(
            static_cast<std::size_t>(std::floor(std::clamp(x.at(0), 0.0, 1.0) * stations_.size())),
            stations_.size() - 1);
        return {static_cast<double>(index)};
    }

    void apply(qn::Network<double>& model, const Value& value) const override {
        if (stations_.empty()) return;
        auto& sn = model.raw_struct();
        const std::size_t mapped_class = cls(sn, class_);
        const std::size_t selected = std::min(
            static_cast<std::size_t>(std::max(0.0, std::floor(scalar_value(value)))),
            stations_.size() - 1);
        std::vector<bool> blocked(sn.nodes.size() + 1, false);
        for (std::size_t k = 0; k < stations_.size(); ++k) {
            const std::size_t station = node(sn, stations_[k]);
            if (k != selected) blocked[station] = true;
        }

        std::vector<std::vector<bool>> connected(
            sn.nodes.size() + 1, std::vector<bool>(sn.nodes.size() + 1, false));
        for (const auto& block : sn.P)
            for (std::size_t i = 1; i <= sn.nodes.size(); ++i)
                for (std::size_t j = 1; j <= sn.nodes.size(); ++j)
                    if (block.second(i - 1, j - 1) > 0.0) connected[i][j] = true;

        qn::RoutingMatrix<double> routing;
        for (std::size_t r = 1; r <= sn.classes.size(); ++r) {
            for (std::size_t i = 1; i <= sn.nodes.size(); ++i) {
                std::size_t degree = 0;
                for (std::size_t j = 1; j <= sn.nodes.size(); ++j)
                    if (connected[i][j] && (r != mapped_class || !blocked[j])) ++degree;
                if (degree == 0) continue;
                const double probability = 1.0 / static_cast<double>(degree);
                for (std::size_t j = 1; j <= sn.nodes.size(); ++j)
                    if (connected[i][j] && (r != mapped_class || !blocked[j]))
                        routing.set(r, r, i, j, probability);
            }
        }
        sn.P.clear();
        sn.Peff.clear();
        model.link(routing);
    }

    std::string type() const override { return "class_mapping"; }

private:
    std::string class_;
    std::vector<std::string> stations_;
};

using VariablePtr=std::shared_ptr<DecisionVariable>;
} }
#endif
