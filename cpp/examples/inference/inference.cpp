/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * `python/examples/inference/est_*.py`: the seventeen demand-estimation examples.
 *
 * WHAT THE C++ PORT CARRIES. `cpp/include/line/api/infer/` holds the
 * estimator leaves -- `infer_qmle`, `infer_rps`, `infer_gibbs`, the
 * arrival-instant queue-length replay (`infer_compute_ql_at_arrival`,
 * `infer_get_qlen_arrival`) and the LQN identification set (`infer_lqn_ekf`,
 * `infer_lqn_jacobian`, `infer_lqn_getobs`, `infer_lqn_findbyname`).
 * `line/inference/sampled_metric.h` carries the sampled-data contract and
 * `line/inference/param_estimator.h` carries collection, lookup, interpolation,
 * automatic selection, typed dispatch, model updates, and all eleven MATLAB
 * estimator methods.
 *
 * Every example below now exercises `ParamEstimator`, including the QMLE,
 * Gibbs, MLPS and FMLPS adapters around their existing API leaves.
 *
 * THE SYNTHETIC DATA. Every reference script draws its samples from numpy,
 * whose Mersenne Twister has no C++ counterpart, and eight of the seventeen do
 * not even seed it, so their own numbers change from run to run. Where an
 * estimator needs individual observations the realized draws are carried as
 * literals, with the call and seed named. Otherwise the examples carry the
 * deterministic mean structure the reference perturbs, since inventing a
 * second random stream would make the port look like it reproduced numbers it
 * did not.
 *
 * `Exp(nan)`, the placeholder the reference parks on an unestimated station,
 * has no C++ spelling: `Distrib::exp_rate` resolves a non-positive or NaN rate
 * to the Immediate rate rather than propagating the NaN. It is written here as
 * `unestimated()` and is overwritten by the estimate before a solver sees it.
 */

#include <cstdio>
#include <limits>
#include <string>
#include <vector>

#include "examples_common.h"
#include "line/inference/param_estimator.h"
#include "line/solvers/mva/solver_mva_runner.h"

namespace line {
namespace examples {

namespace {

// ---------------------------------------------------------------------------
// Sample series
// ---------------------------------------------------------------------------

std::vector<double> filled(std::size_t n, double v) { return std::vector<double>(n, v); }

double mean_of(const std::vector<double>& v) {
    double s = 0.0;
    for (std::size_t i = 0; i < v.size(); ++i) s += v[i];
    return v.empty() ? 0.0 : s / static_cast<double>(v.size());
}

/** The model the reference estimates on: its stations, servers and classes. */
void model_line(Net& m) {
    const Sn& sn = m.get_struct();
    std::printf("MODEL: %s,", sn.name.c_str());
    for (std::size_t i = 0; i < sn.nstations; ++i)
        std::printf(" %s(servers=%g)", sn.stations[i].name.c_str(), sn.stations[i].nservers);
    for (std::size_t r = 0; r < sn.nclasses; ++r)
        std::printf(" %s(N=%g)", sn.classes[r].name.c_str(), sn.classes[r].population);
    std::printf("\n");
}

/** One `SampledMetric` handed to the estimator, as a refusal can report it. */
void sampled(const char* metric, const char* node, const char* jobclass,
             const std::vector<double>& v) {
    std::printf("  data %-6s %-8s %-8s n=%zu mean=%.6g\n", metric, node, jobclass, v.size(),
                mean_of(v));
}

// ---------------------------------------------------------------------------
// The two model shapes the directory uses
// ---------------------------------------------------------------------------

/** `Exp(nan)`: the placeholder on a station whose demand is to be estimated. */
D unestimated() { return Exp(std::numeric_limits<double>::quiet_NaN()); }

std::vector<double> timestamps(std::size_t n, std::size_t first = 0) {
    std::vector<double> t(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) t[i] = static_cast<double>(first + i);
    return t;
}

Matrix<double> run_ubr(Net& model, std::size_t queue, const std::vector<double>& t,
                       const std::vector<std::vector<double>>& arrivals,
                       const std::vector<double>& aggregate_util,
                       const std::vector<std::vector<double>>& class_util = {}) {
    infer::ParamEstimator estimator(model);
    for (std::size_t r = 0; r < arrivals.size(); ++r) {
        estimator.add_samples(
            infer::SampledMetric(lang::MetricType::ArvR, t, arrivals[r], queue, r + 1));
        if (r < class_util.size() && !class_util[r].empty())
            estimator.add_samples(
                infer::SampledMetric(lang::MetricType::Util, t, class_util[r], queue, r + 1));
    }
    estimator.add_samples(
        infer::SampledMetric(lang::MetricType::Util, t, aggregate_util, queue));
    return estimator.estimate_at({queue});
}

Matrix<double> run_observation_estimator(
    Net& model, const std::string& method, std::size_t queue, const std::vector<double>& t,
    const std::vector<std::vector<double>>& arrivals,
    const std::vector<std::vector<double>>& response, const std::vector<double>& aggregate,
    const std::vector<double>& x0 = {}) {
    infer::EstimatorOptions options;
    options.method = method;
    options.x0 = x0;
    infer::ParamEstimator estimator(model, options);
    for (std::size_t r = 0; r < arrivals.size(); ++r) {
        estimator.add_samples(
            infer::SampledMetric(lang::MetricType::ArvR, t, arrivals[r], queue, r + 1));
        estimator.add_samples(
            infer::SampledMetric(lang::MetricType::RespT, t, response[r], queue, r + 1));
    }
    estimator.add_samples(
        infer::SampledMetric(lang::MetricType::Util, t, aggregate, queue));
    return estimator.estimate_at({queue});
}

Matrix<double> run_erps(Net& model, std::size_t queue, const std::vector<double>& t,
                        const std::vector<std::vector<double>>& response,
                        const std::vector<std::vector<double>>& conditional_qlen) {
    infer::EstimatorOptions options;
    options.method = "erps";
    infer::ParamEstimator estimator(model, options);
    for (std::size_t r = 0; r < response.size(); ++r) {
        estimator.add_samples(
            infer::SampledMetric(lang::MetricType::RespT, t, response[r], queue, r + 1));
        infer::SampledMetric q(lang::MetricType::QLen, t, conditional_qlen[r], queue);
        q.set_conditional(infer::ConditionEvent(queue, r + 1, lang::EventType::ARV));
        estimator.add_samples(q);
    }
    return estimator.estimate_at({queue});
}

void solve_after_estimate(Net& model) {
    mva::MvaOptions options;
    Matrix<double> init;
    const mva::AvgResult<double> result = mva::solver_mva_run_analyzer(model.get_struct(), options, init);
    section("MVA");
    print_avg(model.get_struct(), result);
}

/**
 * The four-node open model: Source -> Delay -> Queue1 -> Sink.
 *
 * Node order is the reference's -- Delay, Queue1, Source, Sink -- because the
 * scripts index `node[1]` for the queue.
 */
Net open4(SchedStrategy queue_sched, const std::vector<D>& delay_svc,
          const std::vector<D>& arrivals, std::size_t& delay, std::size_t& queue) {
    Net m("model");
    delay = m.add_delay("Delay");
    queue = m.add_queue("Queue1", queue_sched);
    Source source(m, "Source");
    Sink sink(m, "Sink");
    Routing P;
    for (std::size_t r = 0; r < delay_svc.size(); ++r) {
        OpenClass c(m, "Class" + std::to_string(r + 1), 0);
        m.set_service(delay, c, delay_svc[r]);
        m.set_service(queue, c, unestimated());
        source.set_arrival(c, arrivals[r]);
        P.set(c, c, delay, queue, 1.0);
        P.set(c, c, queue, sink, 1.0);
        P.set(c, c, source, delay, 1.0);
    }
    m.link(P);
    return m;
}

/** The two-node closed model: Delay -> Queue1 -> Delay, one class per population. */
Net closed2(const std::vector<double>& populations, std::size_t& delay, std::size_t& queue) {
    Net m("model");
    delay = m.add_delay("Delay");
    queue = m.add_queue("Queue1", SchedStrategy::PS);
    Routing P;
    for (std::size_t r = 0; r < populations.size(); ++r) {
        ClosedClass c(m, "Class" + std::to_string(r + 1), populations[r], delay, 0);
        m.set_service(delay, c, D::exp_mean(1.0));
        m.set_service(queue, c, unestimated());
        P.set(c, c, delay, queue, 1.0);
        P.set(c, c, queue, delay, 1.0);
    }
    m.link(P);
    return m;
}

/** `HyperExp(p, lambda1, lambda2)` as the reference spells it. */
D hyper(double p, double r1, double r2) { return HyperExp(p, r1, r2); }

// ---------------------------------------------------------------------------
// est_trace_closed: the realized trace, so infer_gibbs sees the reference data
// ---------------------------------------------------------------------------

/**
 * `np.sort(np.random.rand(200) * 100)` after `np.random.seed(1)`, in seconds.
 * The Gibbs estimator consumes these sample by sample, so they are the realized
 * draws and not a restatement of the distribution they came from.
 */
const double kArrivalTimes[200] = {
    0.011437481734488664, 0.28703270311589701, 1.2555980159115854, 1.3951572975597015,
    1.5821242846556283, 1.8288277344191806, 1.8576202177409518, 1.8647289372943021,
    1.9366957870297075, 1.9880133839795588, 2.6210986877719278, 2.7387593197926163,
    2.8306488020794607, 3.9054783232882362, 4.455187854476172, 4.9953458946087164,
    5.336254511708038, 5.9917689512211663, 6.5961090684023782, 6.6000172722062489,
    6.6334834428441569, 6.653648135411494, 7.0022143719222329, 7.1974279689486771,
    8.5044211369777916, 9.233859476879779, 9.8346833833050091, 10.233442882782583,
    10.322600657764202, 10.749412910609291, 11.474597295337519, 12.134345574073734,
    12.417331511991115, 12.42709619721647, 13.002857211827767, 13.645522566068502,
    13.713574962887776, 13.747470414623752, 13.927634725075855, 14.038693859523377,
    14.672857490581015, 14.675589081711305, 15.679139464608427, 16.53541971169328,
    16.983041956456891, 17.234050834532855, 18.62602113776709, 19.343428262332772,
    19.542948110931878, 19.810148908487879, 20.329323466099048, 20.445224973151742,
    21.01740099148396, 21.162811600005902, 22.57093386078547, 23.297427384102043,
    23.43620861214205, 23.702698024302769, 23.984775914758615, 24.621106760304588,
    25.232574457032342, 26.031509857854097, 26.329677048711098, 26.491955766280938,
    26.554665937222623, 26.992789176502608, 27.91836790111395, 28.044399206440517,
    28.777533858634875, 29.361414837367949, 30.233257263183976, 31.342417815924286,
    31.551563100606295, 31.736240932216077, 32.664490177209615, 34.556072704304775,
    34.776585974550656, 34.889834197784253, 35.726976000249977, 37.008419791410631,
    38.014117262355043, 38.786064406417175, 39.676747423066992, 39.767683698553355,
    40.813680276128117, 41.405598781956833, 41.417926952690266, 41.702200470257402,
    41.730480236712694, 41.919451440329482, 42.110762500505217, 42.809118987129494,
    44.789352617590517, 44.991213347994055, 48.634511093703182, 49.157315928033832,
    49.376971426872998, 51.48891120583086, 52.467030912373367, 52.705810225760928,
    53.316528497301704, 53.589640591551159, 53.881673400335693, 53.883106434165285,
    55.09482191178968, 55.282197868576588, 55.624023399041889, 55.868982844575164,
    55.971698205414242, 56.103021925570992, 56.810046191994211, 56.885143708648137,
    57.367948667228589, 57.411760549201304, 57.838961438713177, 57.97452192457969,
    58.135892727325775, 58.575927145828786, 58.655504050199291, 58.930553690328423,
    60.632946165333038, 61.677835700165758, 61.714491362072387, 61.995571838137977,
    62.169572020912177, 62.336011579180273, 62.367220705560889, 62.971750702156449,
    63.946088087994013, 66.344149781844806, 66.379464521978875, 66.923289345318466,
    67.046751017840222, 67.883553293989095, 68.521950039675943, 68.650092768158373,
    69.089691751692399, 69.18771139504733, 69.232261566931413, 69.440015772774515,
    69.681816148990023, 69.975836002093118, 71.152475862847169, 71.298898038267666,
    72.032449344215806, 72.599798535045153, 73.506596328866948, 74.382585407509296,
    74.533443090650209, 74.712164273718457, 74.81656543798394, 75.014431494496748,
    75.081210313615557, 75.094243402733724, 75.275555373881389, 75.387618846124639,
    75.546305260246641, 77.217802954324682, 78.927932845148845, 80.063267268061637,
    80.074456867553664, 80.475456374334541, 80.710519561877916, 80.739128870952385,
    82.898089955017866, 83.462567189737285, 84.203089235960576, 84.682880149003523,
    86.002794868288802, 86.354185455942869, 87.638915229603825, 87.811743639094544,
    87.814250342941307, 88.330609120580988, 88.594209931077444, 89.460666350384727,
    89.588621819606686, 90.337952056225376, 90.340191528788353, 90.781585250352407,
    90.853515091979915, 90.859550309309554, 92.302453554648338, 92.480797039935069,
    92.750858039603386, 92.943723374376134, 93.197206919683723, 93.259546303716363,
    94.45947559908133, 94.901632068761643, 94.948925870707129, 95.788953015050197,
    96.484004714838562, 96.727633000027197, 96.826157571939746, 96.959574831967458,
    97.00199890883124, 98.861615441244894, 98.886108890649467, 99.732285045148046,
};

/** `0.5 + np.random.rand(200) * 0.3`, the response times of the same draw. */
const double kResponseTimes[200] = {
    0.78505283577412388, 0.66699595645854692, 0.77468190492988231, 0.69246986268390121,
    0.61700231424237384, 0.64579720012907293, 0.68129314487599202, 0.66486437645256879,
    0.77785442801193616, 0.77562003069008179, 0.61846268387706649, 0.78897875853221355,
    0.55218670000413927, 0.53789885583189123, 0.54052374741401943, 0.65169864970306901,
    0.50645744158225936, 0.78439106336530429, 0.74813464135121976, 0.50450569422263669,
    0.55285887667251654, 0.59961907231005107, 0.53929905344327511, 0.74284720763797463,
    0.60342099580498809, 0.78203224470001009, 0.67460425398412405, 0.76364959532355314,
    0.75342033361766658, 0.77161769561259475, 0.6379640797450421, 0.66390404480611198,
    0.73958107734561174, 0.58571565552024429, 0.64707605678597835, 0.67973309229376344,
    0.50465998266525069, 0.67804442245899121, 0.6301029046968375, 0.74220815866545242,
    0.59457344092861186, 0.7678666125575454, 0.67335716458536066, 0.55520306048823953,
    0.7363787701476513, 0.68360935311348014, 0.51617278162212588, 0.62605810400033968,
    0.70372065096962888, 0.77558053339325772, 0.50012060746740727, 0.79302774470931792,
    0.6129740944237323, 0.79213506150749668, 0.68141483029221595, 0.7486537423942079,
    0.6724134514124307, 0.68842285949220505, 0.58567288450870814, 0.67605000219682487,
    0.72500652911079799, 0.75749415092871275, 0.72652465654030407, 0.70941717453419084,
    0.75934382901637987, 0.59680429905102395, 0.70123663723627616, 0.63526218092400477,
    0.61463082560945514, 0.62324340497665565, 0.62044387504086218, 0.59521518378748306,
    0.68657581037609039, 0.62907418124637926, 0.79214062337817559, 0.70334026743029332,
    0.55957096652813321, 0.62801030280440984, 0.60300387193232696, 0.73929164118756996,
    0.7639994865690295, 0.77115258674791165, 0.69881594371257871, 0.58106247860892735,
    0.57571001045137693, 0.75646938280922071, 0.65831439389262403, 0.74064832520137935,
    0.67174655515748194, 0.71994275758625337, 0.65570348823921676, 0.73126517315056638,
    0.67065739721141471, 0.63971296357759422, 0.60280667238598451, 0.52046280452501126,
    0.61337725379842989, 0.5238878233094757, 0.79484513411913349, 0.55448385539922918,
    0.74355760931616188, 0.7624884934867695, 0.70652397571578296, 0.6708483238236127,
    0.54829143104468259, 0.64006400682899189, 0.60355161534656454, 0.56751198734422537,
    0.67775356062973902, 0.59368095131057508, 0.77489166604050519, 0.77289065748546704,
    0.57713548813465887, 0.53326739022320879, 0.55788881960573389, 0.64987525120366585,
    0.71857570039237884, 0.56245833152263747, 0.57441006751316925, 0.75550156248091005,
    0.62475461548025824, 0.68500552014657079, 0.57009984177177497, 0.53059017782773921,
    0.65475710509055896, 0.64314229611493468, 0.54580149322794902, 0.68654186952212459,
    0.66320303564418137, 0.69624120409122325, 0.54336366203739794, 0.72545834514057317,
    0.56661474193997685, 0.65580554730980989, 0.73558880846648567, 0.50669912839754183,
    0.59730873791785599, 0.76187671292003323, 0.75341288228062087, 0.66153217777836315,
    0.75998248224995391, 0.78494179740608794, 0.74792209928880238, 0.75623463315096162,
    0.52962302054610455, 0.69539129970229729, 0.7110550964457959, 0.6830722437939889,
    0.73988457852080836, 0.51037136596149024, 0.73107162036649398, 0.71951858022185811,
    0.57790951798950596, 0.57712078964664404, 0.68969099522903832, 0.60358923847732449,
    0.73897660340218629, 0.63384386960103711, 0.73482482443525154, 0.79714153508715224,
    0.59007450186017008, 0.54290174847742978, 0.7703925309047821, 0.66246781366621199,
    0.79242211125652462, 0.69098132000056922, 0.79817390738314364, 0.66382124124288544,
    0.65792778017165632, 0.54062837092165095, 0.6067115512951643, 0.50786557018891199,
    0.5481185538558947, 0.7236911578122488, 0.50911990697863663, 0.60996292917312334,
    0.7587038758514133, 0.7078033152522063, 0.70728264261985052, 0.55659104028659323,
    0.63257128422452702, 0.67447322220328731, 0.79692551229911157, 0.56117186756993653,
    0.57431987052870648, 0.57865192513187358, 0.72505172398939788, 0.63709259823070785,
    0.51707883152332845, 0.65255487218307184, 0.56358804939310281, 0.73958127342806301,
    0.58919941445218171, 0.50828180358635111, 0.67802973483938633, 0.75315212867935666,
    0.61430483721866425, 0.72495749321784297, 0.65334244348913839, 0.66228554148865659,
};

}  // namespace

// ---------------------------------------------------------------------------
// EKF
// ---------------------------------------------------------------------------

/** `est_ekf_changepoint.py`: a sliding EKF window over a demand change point. */
void est_ekf_changepoint() {
    std::size_t delay = 0, queue = 0;
    Net m = open4(SchedStrategy::FCFS, {D::exp_mean(1.0)}, {Exp(1.0)}, delay, queue);
    model_line(m);

    const std::size_t n = 200, changeT = 100, W = 30;
    const double D1 = 0.3, D2 = 0.6, lam = 1.0;
    // Deterministic part only: the reference adds randn*0.02 from an UNSEEDED stream.
    std::vector<double> arvr = filled(n, lam), util(n, 0.0), respt(n, 0.0);
    for (std::size_t t = 0; t < n; ++t) {
        const double d = t < changeT ? D1 : D2;
        util[t] = lam * d;
        respt[t] = d / (1.0 - util[t]);
    }
    sampled("ArvR", "Queue1", "Class1", arvr);
    sampled("RespT", "Queue1", "Class1", respt);
    sampled("Util", "Queue1", "(aggr)", util);
    std::printf("  window W=%zu, warm started from the previous window's estimate\n", W);

    std::vector<double> estimates(n - W + 1, 0.0), warm{0.3};
    for (std::size_t end = W; end <= n; ++end) {
        const std::vector<double> tw = timestamps(W, end - W);
        const std::vector<double> aw(arvr.begin() + static_cast<std::ptrdiff_t>(end - W),
                                     arvr.begin() + static_cast<std::ptrdiff_t>(end));
        const std::vector<double> rw(respt.begin() + static_cast<std::ptrdiff_t>(end - W),
                                     respt.begin() + static_cast<std::ptrdiff_t>(end));
        const std::vector<double> uw(util.begin() + static_cast<std::ptrdiff_t>(end - W),
                                     util.begin() + static_cast<std::ptrdiff_t>(end));
        const Matrix<double> estimate =
            run_observation_estimator(m, "ekf", queue, tw, {aw}, {rw}, uw, warm);
        warm[0] = estimate(0, 0);
        estimates[end - W] = warm[0];
    }
    std::printf("\n=== EKF Change-Point Detection (warm-start) ===\n");
    std::printf("True demand: D=%g (t<=%zu), D=%g (t>%zu)\n", D1, changeT, D2, changeT);
    std::printf("Window size: %zu\n\n", W);
    const std::size_t cp[10] = {W, 50, 80, 100, 110, 120, 130, 150, 180, 200};
    std::printf("%6s  %10s  %10s\n", "t", "Estimated", "True");
    for (std::size_t i = 0; i < 10; ++i)
        std::printf("%6zu  %10.4f  %10.4f\n", cp[i], estimates[cp[i] - W],
                    cp[i] <= changeT ? D1 : D2);
}

/** `est_ekf_closed.py`: the EKF and the autoMethod selection on a closed model. */
void est_ekf_closed() {
    std::size_t delay = 0, queue = 0;
    Net m = closed2({2.0}, delay, queue);
    model_line(m);

    const std::size_t n = 100;
    // np.random.seed(1): arvr = 1.5 - rand(n)*0.1, so the mean structure is 1.45.
    const std::vector<double> arvr = filled(n, 1.45);
    std::vector<double> util(n, 0.0), respt(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        util[i] = 0.4 * arvr[i];
        respt[i] = 0.4 / (1.0 - util[i]);
    }
    sampled("ArvR", "Queue1", "Class1", arvr);
    sampled("RespT", "Queue1", "Class1", respt);
    sampled("Util", "Queue1", "(aggr)", util);

    std::printf("\n=== EKF Estimator ===\n");
    const Matrix<double> ekf = run_observation_estimator(
        m, "ekf", queue, timestamps(n), {arvr}, {respt}, util, {0.2});
    std::printf("Estimated demand: [[%.8f]]\n", ekf(0, 0));
    std::printf("\n=== autoMethod selection ===\n");
    infer::ParamEstimator automatic(m);
    automatic.add_samples(infer::SampledMetric(lang::MetricType::ArvR, timestamps(n), arvr,
                                                queue, 1));
    automatic.add_samples(infer::SampledMetric(lang::MetricType::RespT, timestamps(n), respt,
                                                queue, 1));
    automatic.add_samples(
        infer::SampledMetric(lang::MetricType::Util, timestamps(n), util, queue));
    const std::string selected = automatic.auto_method();
    std::printf("Selected: %s (%s)\n", selected.c_str(),
                infer::ParamEstimator::get_required_metrics(selected).c_str());
    automatic.estimate_at({queue});
    solve_after_estimate(m);
}

/** `est_ekf_open.py`: the same filter on the four-node open model. */
void est_ekf_open() {
    std::size_t delay = 0, queue = 0;
    Net m = open4(SchedStrategy::FCFS, {D::exp_mean(1.0)}, {Exp(1.0)}, delay, queue);
    model_line(m);

    const std::size_t n = 100;
    // Deterministic part only: the reference adds randn noise from an UNSEEDED stream.
    const std::vector<double> arvr = filled(n, 1.0);
    const std::vector<double> util = filled(n, 0.4);
    std::vector<double> respt(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) respt[i] = 0.4 / (1.0 - util[i]);
    sampled("ArvR", "Queue1", "Class1", arvr);
    sampled("RespT", "Queue1", "Class1", respt);
    sampled("Util", "Queue1", "(aggr)", util);

    std::printf("\n=== EKF Estimator (Open Network) ===\n");
    const Matrix<double> estimate = run_observation_estimator(
        m, "ekf", queue, timestamps(n), {arvr}, {respt}, util, {0.2});
    std::printf("Estimated demand: [[%.8f]]\n", estimate(0, 0));
    solve_after_estimate(m);
}

// ---------------------------------------------------------------------------
// ERPS
// ---------------------------------------------------------------------------

namespace {

/** The two-class ERPS/UBO data set, which four scripts share. */
void erps_data(std::size_t n, std::vector<double>& u, std::vector<double>& r1,
               std::vector<double>& r2, std::vector<double>& q) {
    // np.random.seed(1): arvr1 = 1 - rand*0.15 and arvr2 = 2 - rand*0.15, so the
    // mean structure is 0.925 and 1.925.
    u.assign(n, 0.0);
    r1.assign(n, 0.0);
    r2.assign(n, 0.0);
    q.assign(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        u[i] = 0.1 * 0.925 + 0.3 * 1.925;
        r1[i] = 0.1 / (1.0 - u[i]);
        r2[i] = 0.3 / (1.0 - u[i]);
        q[i] = 1.0 + u[i] / (1.0 - u[i]);
    }
}

}  // namespace

/** `est_erps_closed.py`: ERPS on a two-class closed PS station. */
void est_erps_closed() {
    std::size_t delay = 0, queue = 0;
    Net m = closed2({1.0, 3.0}, delay, queue);
    model_line(m);

    const std::size_t n = 1000;
    std::vector<double> u, r1, r2, q;
    erps_data(n, u, r1, r2, q);
    sampled("QLen", "Queue1", "arv Cls1", q);
    sampled("QLen", "Queue1", "arv Cls2", q);
    sampled("RespT", "Queue1", "Class1", r1);
    sampled("RespT", "Queue1", "Class2", r2);

    const Matrix<double> estimate = run_erps(m, queue, timestamps(n), {r1, r2}, {q, q});
    std::printf("Estimated demands: [[%.8f %.8f]]\n", estimate(0, 0), estimate(0, 1));
    solve_after_estimate(m);
}

/** `est_erps_open.py`: the same estimator on the open model. */
void est_erps_open() {
    std::size_t delay = 0, queue = 0;
    Net m = open4(SchedStrategy::PS, {hyper(0.5, 3.0, 10.0), hyper(0.5, 2.0, 8.0)},
                  {Exp(0.1), Exp(0.05)}, delay, queue);
    model_line(m);

    const std::size_t n = 1000;
    std::vector<double> u, r1, r2, q;
    erps_data(n, u, r1, r2, q);
    sampled("QLen", "Queue1", "arv Cls1", q);
    sampled("QLen", "Queue1", "arv Cls2", q);
    sampled("RespT", "Queue1", "Class1", r1);
    sampled("RespT", "Queue1", "Class2", r2);

    const Matrix<double> estimate = run_erps(m, queue, timestamps(n), {r1, r2}, {q, q});
    std::printf("Estimated demands: [[%.8f %.8f]]\n", estimate(0, 0), estimate(0, 1));
    solve_after_estimate(m);
}

// ---------------------------------------------------------------------------
// MCMC, MLE
// ---------------------------------------------------------------------------

/** `est_mcmc_closed.py`: grid-posterior Gibbs sampling from aggregate queue lengths. */
void est_mcmc_closed() {
    std::size_t delay = 0, queue = 0;
    Net m = closed2({2.0, 3.0}, delay, queue);
    model_line(m);

    const std::size_t n = 500;
    // np.random.seed(1): qlen = 1.5 + rand(n)*0.5, mean structure 1.75.
    const std::vector<double> qlen = filled(n, 1.75);
    sampled("QLen", "Queue1", "(aggr)", qlen);

    std::printf("\n=== MCMC Estimator ===\n");
    infer::EstimatorOptions options;
    options.method = "mcmc";
    options.random_seed = 1;
    infer::ParamEstimator estimator(m, options);
    estimator.add_samples(
        infer::SampledMetric(lang::MetricType::QLen, timestamps(n), qlen, queue));
    const Matrix<double> estimate = estimator.estimate_at({queue});
    std::printf("Estimated demands: [[%.8f %.8f]]\n", estimate(0, 0), estimate(0, 1));
    solve_after_estimate(m);
}

/** `est_mle_open.py`: maximum likelihood over an MVA prediction, two open classes. */
void est_mle_open() {
    std::size_t delay = 0, queue = 0;
    Net m = open4(SchedStrategy::FCFS, {D::exp_mean(1.0), D::exp_mean(1.0)},
                  {Exp(1.0), Exp(0.5)}, delay, queue);
    model_line(m);

    const std::size_t n = 100;
    // Deterministic part only: the reference adds randn noise from an UNSEEDED stream.
    const std::vector<double> a1 = filled(n, 1.0), a2 = filled(n, 0.5);
    const std::vector<double> util = filled(n, 0.25);
    std::vector<double> r1(n, 0.0), r2(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        r1[i] = 0.1 / (1.0 - util[i]);
        r2[i] = 0.3 / (1.0 - util[i]);
    }
    sampled("ArvR", "Queue1", "Class1", a1);
    sampled("ArvR", "Queue1", "Class2", a2);
    sampled("RespT", "Queue1", "Class1", r1);
    sampled("RespT", "Queue1", "Class2", r2);
    sampled("Util", "Queue1", "(aggr)", util);

    std::printf("\n=== MLE Estimator (Open Network) ===\n");
    const Matrix<double> estimate = run_observation_estimator(
        m, "mle", queue, timestamps(n), {a1, a2}, {r1, r2}, util, {0.1, 0.3});
    std::printf("Estimated demands: [[%.8f %.8f]]\n", estimate(0, 0), estimate(0, 1));
    solve_after_estimate(m);
}

// ---------------------------------------------------------------------------
// QMLE: the one estimator this directory runs end to end
// ---------------------------------------------------------------------------

/**
 * `est_qmle_closed.py`, ported in full through `infer_qmle`.
 *
 * The dispatcher's `qmle` arm is model bookkeeping around one call: it reads
 * the populations and the delay-station think times off the model, averages
 * each queue-length series, and hands the three to `infer_qmle`. That
 * bookkeeping is written out here; the estimator itself is the ported one.
 *
 * The two series are `0.3 + rand(1000)*0.1` and `0.8 + rand(1000)*0.1` after
 * `np.random.seed(1)`. `infer_qmle` consumes only their MEANS, so the realized
 * means of that draw are carried rather than two thousand realized samples;
 * they reproduce the reference's estimate to every digit it prints.
 */
void est_qmle_closed() {
    std::size_t delay = 0, queue = 0;
    Net m = closed2({2.0, 3.0}, delay, queue);

    const std::size_t n = 1000;
    const double qlen1_mean = 0.35006045994559054;  // mean of 0.3 + rand(1000)*0.1, seed 1
    const double qlen2_mean = 0.8518112390463691;   // mean of 0.8 + rand(1000)*0.1, seed 1
    sampled("QLen", "Queue1", "Class1", filled(n, qlen1_mean));
    sampled("QLen", "Queue1", "Class2", filled(n, qlen2_mean));

    infer::EstimatorOptions options;
    options.method = "qmle";
    infer::ParamEstimator estimator(m, options);
    estimator.add_samples(infer::SampledMetric(lang::MetricType::QLen, timestamps(n),
                                                filled(n, qlen1_mean), queue, 1));
    estimator.add_samples(infer::SampledMetric(lang::MetricType::QLen, timestamps(n),
                                                filled(n, qlen2_mean), queue, 2));
    const Matrix<double> est = estimator.estimate_at({queue});
    std::printf("Estimated demands: [[%.8f %.8f]]\n", est(0, 0), est(0, 1));

    mva::MvaOptions mopt;
    Matrix<double> init;
    const mva::AvgResult<double> res = mva::solver_mva_run_analyzer(m.get_struct(), mopt, init);
    section("MVA");
    print_avg(m.get_struct(), res);
}

// ---------------------------------------------------------------------------
// VI: variational inference over transition counts
// ---------------------------------------------------------------------------

/**
 * `est_vi_closed.m` / `est_vi_closed.py`, run end to end through the `vi`
 * estimator.
 *
 * The method (Perez-Casale, Adv. Appl. Prob. 53(3), 2021) infers service rates
 * from NOISY QUEUE-LENGTH READINGS taken over time: each reading is exact with
 * probability 1-epsilon and uniform over the remaining feasible values
 * otherwise. It is the only estimator here that returns a conjugate Gamma
 * POSTERIOR rather than a point, and the only one that reads the queue lengths
 * of EVERY station rather than of the estimated one alone -- the transition
 * counts it is written in are pinned by the whole picture.
 *
 * It carries no random-number stream, so the numbers below are the MATLAB and
 * Python ones digit for digit.
 */
void est_vi_closed() {
    const double N = 10.0;
    Net m("model");
    Delay delay(m, "Delay");
    Queue queue(m, "Queue1", SchedStrategy::FCFS);
    ClosedClass cl(m, "Class1", N, delay, 0);
    m.set_service(delay, cl, D::exp_mean(2.0));
    m.set_service(queue, cl, D::exp_mean(0.5));  // starting point of the estimate
    Routing P;
    P.set(cl, cl, delay, queue, 1.0);
    P.set(cl, cl, queue, delay, 1.0);
    m.link(P);

    // queue-length readings, one per unit time, 10% of them faulty
    const std::vector<double> ts = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    const std::vector<double> qlen = {1, 2, 3, 2, 4, 3, 5, 4, 3, 4};
    std::vector<double> dlen(qlen.size());
    for (std::size_t k = 0; k < qlen.size(); ++k) dlen[k] = N - qlen[k];
    sampled("QLen", "Delay", "Class1", dlen);
    sampled("QLen", "Queue1", "Class1", qlen);

    infer::EstimatorOptions options;
    options.method = "vi";
    options.epsilon = 0.1;       // probability that a reading is faulty
    options.prior_shape = 2.0;   // Gamma prior shape; the rate comes from the model
    options.variational.ngrid = 51;
    options.variational.nsamples = 32;
    options.variational.ymax = 60;
    options.variational.iter_max = 5;
    infer::ParamEstimator estimator(m, options);
    estimator.add_samples(
        infer::SampledMetric(lang::MetricType::QLen, ts, dlen, delay, 1));
    estimator.add_samples(
        infer::SampledMetric(lang::MetricType::QLen, ts, qlen, queue, 1));
    const Matrix<double> est = estimator.estimate_at({queue});
    std::printf("Estimated demand: %.8f\n", est(0, 0));
    std::printf("posterior service rate ~ Gamma(%.4f, %.4f), mean %.4f\n",
                estimator.options.posterior_alpha[0], estimator.options.posterior_beta[0],
                estimator.options.posterior_alpha[0] / estimator.options.posterior_beta[0]);
    std::printf("evidence lower bound over the iterations:");
    for (std::size_t i = 0; i < estimator.options.bound.size(); ++i)
        std::printf(" %.3f", estimator.options.bound[i]);
    std::printf("\n");

    mva::MvaOptions mopt;
    Matrix<double> init;
    const mva::AvgResult<double> res = mva::solver_mva_run_analyzer(m.get_struct(), mopt, init);
    section("MVA");
    print_avg(m.get_struct(), res);
}

// ---------------------------------------------------------------------------
// Trace-driven estimators: MLPS, FMLPS, Gibbs
// ---------------------------------------------------------------------------

/**
 * `est_trace_closed.py`: the three trace-driven estimators on one PS station.
 *
 * THE GIBBS NUMBER WILL NOT MATCH THE REFERENCE DIGIT FOR DIGIT and is not
 * meant to. The estimator draws its test set and each coordinate update from a
 * uniform stream; the port takes an explicit `McRng` instead of reproducing
 * numpy's Mersenne Twister, so the two agree in distribution and on every
 * deterministic intermediate rather than sample by sample. Measured: 0.1819
 * here against 0.1841 in Python, on identical input data.
 */
void est_trace_closed() {
    const double N = 5.0;
    std::size_t delay = 0, queue = 0;
    Net m = closed2({N}, delay, queue);

    std::vector<double> at(kArrivalTimes, kArrivalTimes + 200);
    std::vector<double> rt(kResponseTimes, kResponseTimes + 200);
    const std::vector<double> tput = filled(200, N / (1.0 + 0.5));
    sampled("ArvR", "Queue1", "Class1", at);
    sampled("RespT", "Queue1", "Class1", rt);
    sampled("Tput", "Queue1", "Class1", tput);

    const auto run_trace = [&](const std::string& method) {
        infer::EstimatorOptions options;
        options.method = method;
        options.random_seed = 1;
        infer::ParamEstimator estimator(m, options);
        infer::SampledMetric arrivals(lang::MetricType::ArvR, timestamps(at.size()), at, queue,
                                      1);
        arrivals.set_trace();
        infer::SampledMetric response(lang::MetricType::RespT, timestamps(rt.size()), rt, queue,
                                      1);
        response.set_trace();
        estimator.add_samples(arrivals);
        estimator.add_samples(response);
        if (method == "gibbs")
            estimator.add_samples(infer::SampledMetric(lang::MetricType::Tput,
                                                        timestamps(tput.size()), tput, queue, 1));
        return estimator.estimate_at({queue});
    };

    std::printf("\n=== MLPS Estimator ===\n");
    const Matrix<double> mlps = run_trace("mlps");
    std::printf("MLPS demand: Class1=%.4f\n", mlps(0, 0));
    std::printf("\n=== FMLPS Estimator ===\n");
    const Matrix<double> fmlps = run_trace("fmlps");
    std::printf("FMLPS demand: Class1=%.4f\n", fmlps(0, 0));

    std::printf("\n=== Gibbs Estimator ===\n");
    const Matrix<double> gibbs = run_trace("gibbs");
    std::printf("Gibbs demand: Class1=%.4f\n", gibbs(0, 0));

    std::printf("\n=== Comparison (true demand ~ 0.5) ===\n");
    std::printf("MLPS:  %.4f\n", mlps(0, 0));
    std::printf("FMLPS: %.4f\n", fmlps(0, 0));
    std::printf("Gibbs: %.4f\n", gibbs(0, 0));
    mva::MvaOptions mopt;
    Matrix<double> init;
    const mva::AvgResult<double> res = mva::solver_mva_run_analyzer(m.get_struct(), mopt, init);
    section("MVA");
    print_avg(m.get_struct(), res);
}

// ---------------------------------------------------------------------------
// UBO
// ---------------------------------------------------------------------------

/** `est_ubo_changepoint.py`: a sliding UBO window over a two-class change point. */
void est_ubo_changepoint() {
    std::size_t delay = 0, queue = 0;
    Net m = open4(SchedStrategy::FCFS, {D::exp_mean(1.0), D::exp_mean(1.0)},
                  {Exp(1.0), Exp(0.5)}, delay, queue);
    model_line(m);

    const std::size_t n = 200, changeT = 100, W = 30;
    const double D1_before = 0.1, D1_after = 0.3, D2_val = 0.2;
    const double lambda1 = 1.0, lambda2 = 0.5;
    // Deterministic part only: the reference adds randn*0.01 from an UNSEEDED stream.
    const std::vector<double> a1 = filled(n, lambda1), a2 = filled(n, lambda2);
    std::vector<double> util(n, 0.0), r1(n, 0.0), r2(n, 0.0);
    for (std::size_t t = 0; t < n; ++t) {
        const double d1 = t < changeT ? D1_before : D1_after;
        util[t] = lambda1 * d1 + lambda2 * D2_val;
        r1[t] = d1 / (1.0 - util[t]);
        r2[t] = D2_val / (1.0 - util[t]);
    }
    sampled("ArvR", "Queue1", "Class1", a1);
    sampled("ArvR", "Queue1", "Class2", a2);
    sampled("RespT", "Queue1", "Class1", r1);
    sampled("RespT", "Queue1", "Class2", r2);
    sampled("Util", "Queue1", "(aggr)", util);

    std::vector<double> estimate1(n - W + 1, 0.0), estimate2(n - W + 1, 0.0);
    for (std::size_t end = W; end <= n; ++end) {
        const std::vector<double> tw = timestamps(W, end - W);
        const std::vector<double> a1w(a1.begin() + static_cast<std::ptrdiff_t>(end - W),
                                      a1.begin() + static_cast<std::ptrdiff_t>(end));
        const std::vector<double> a2w(a2.begin() + static_cast<std::ptrdiff_t>(end - W),
                                      a2.begin() + static_cast<std::ptrdiff_t>(end));
        const std::vector<double> r1w(r1.begin() + static_cast<std::ptrdiff_t>(end - W),
                                      r1.begin() + static_cast<std::ptrdiff_t>(end));
        const std::vector<double> r2w(r2.begin() + static_cast<std::ptrdiff_t>(end - W),
                                      r2.begin() + static_cast<std::ptrdiff_t>(end));
        const std::vector<double> uw(util.begin() + static_cast<std::ptrdiff_t>(end - W),
                                     util.begin() + static_cast<std::ptrdiff_t>(end));
        const Matrix<double> estimate =
            run_observation_estimator(m, "ubo", queue, tw, {a1w, a2w}, {r1w, r2w}, uw);
        estimate1[end - W] = estimate(0, 0);
        estimate2[end - W] = estimate(0, 1);
    }
    std::printf("\n=== UBO Change-Point Detection (2-class) ===\n");
    std::printf("Class 1: D=%g (t<=%zu), D=%g (t>%zu)\n", D1_before, changeT, D1_after, changeT);
    std::printf("Class 2: D=%g (constant)\n", D2_val);
    std::printf("Window size: %zu\n\n", W);
    const std::size_t cp[10] = {W, 50, 80, 100, 110, 120, 130, 150, 180, 200};
    std::printf("%6s  %8s  %8s  %8s  %8s\n", "t", "Est D1", "True D1", "Est D2", "True D2");
    for (std::size_t i = 0; i < 10; ++i)
        std::printf("%6zu  %8.4f  %8.4f  %8.4f  %8.4f\n", cp[i],
                    estimate1[cp[i] - W], cp[i] <= changeT ? D1_before : D1_after,
                    estimate2[cp[i] - W], D2_val);
}

/** `est_ubo_closed.py`: UBO on a two-class closed PS station. */
void est_ubo_closed() {
    std::size_t delay = 0, queue = 0;
    Net m = closed2({1.0, 2.0}, delay, queue);
    model_line(m);

    const std::size_t n = 30;
    // np.random.seed(1): arvr1 = 2 - rand*0.15 and arvr2 = 1 - rand*0.15.
    const std::vector<double> a1 = filled(n, 1.925), a2 = filled(n, 0.925);
    std::vector<double> util(n, 0.0), r1(n, 0.0), r2(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        util[i] = 0.1 * a1[i] + 0.3 * a2[i];
        r1[i] = 0.1 / (1.0 - util[i]);
        r2[i] = 0.3 / (1.0 - util[i]);
    }
    sampled("ArvR", "Queue1", "Class1", a1);
    sampled("ArvR", "Queue1", "Class2", a2);
    sampled("RespT", "Queue1", "Class1", r1);
    sampled("RespT", "Queue1", "Class2", r2);
    sampled("Util", "Queue1", "(aggr)", util);

    const Matrix<double> estimate = run_observation_estimator(
        m, "ubo", queue, timestamps(n), {a1, a2}, {r1, r2}, util);
    std::printf("Estimated demands: [[%.8f %.8f]]\n", estimate(0, 0), estimate(0, 1));
    solve_after_estimate(m);
}

/** `est_ubo_open.py`: UBO on the open model with hyperexponential think times. */
void est_ubo_open() {
    std::size_t delay = 0, queue = 0;
    Net m = open4(SchedStrategy::FCFS, {hyper(0.5, 3.0, 10.0), hyper(0.5, 2.0, 8.0)},
                  {Exp(0.1), Exp(0.05)}, delay, queue);
    model_line(m);

    const std::size_t n = 1000;
    // Deterministic part only: arvr1 = 2 - rand*0.15, arvr2 = 1 - rand*0.10, UNSEEDED.
    const std::vector<double> a1 = filled(n, 1.925), a2 = filled(n, 0.95);
    std::vector<double> util(n, 0.0), r1(n, 0.0), r2(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        util[i] = 0.1 * a1[i] + 0.3 * a2[i];
        r1[i] = 0.1 / (1.0 - util[i]);
        r2[i] = 0.3 / (1.0 - util[i]);
    }
    sampled("ArvR", "Queue1", "Class1", a1);
    sampled("ArvR", "Queue1", "Class2", a2);
    sampled("RespT", "Queue1", "Class1", r1);
    sampled("RespT", "Queue1", "Class2", r2);
    sampled("Util", "Queue1", "(aggr)", util);

    const Matrix<double> estimate = run_observation_estimator(
        m, "ubo", queue, timestamps(n), {a1, a2}, {r1, r2}, util);
    std::printf("Estimated demands: [[%.8f %.8f]]\n", estimate(0, 0), estimate(0, 1));
    solve_after_estimate(m);
}

/** `est_ubo_variants_closed.py`: UBO against MLE on the same closed model. */
void est_ubo_variants_closed() {
    std::size_t delay = 0, queue = 0;
    Net m = closed2({1.0, 3.0}, delay, queue);
    model_line(m);

    const std::size_t n = 1000;
    // np.random.seed(1): arvr1 = 2 - rand*0.15 and arvr2 = 1 - rand*0.15.
    const std::vector<double> a1 = filled(n, 1.925), a2 = filled(n, 0.925);
    std::vector<double> util(n, 0.0), r1(n, 0.0), r2(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        util[i] = 0.1 * a1[i] + 0.3 * a2[i];
        r1[i] = 0.1 / (1.0 - util[i]);
        r2[i] = 0.3 / (1.0 - util[i]);
    }
    sampled("ArvR", "Queue1", "Class1", a1);
    sampled("ArvR", "Queue1", "Class2", a2);
    sampled("RespT", "Queue1", "Class1", r1);
    sampled("RespT", "Queue1", "Class2", r2);
    sampled("Util", "Queue1", "(aggr)", util);

    std::printf("\n=== UBO Estimator ===\n");
    const Matrix<double> ubo = run_observation_estimator(
        m, "ubo", queue, timestamps(n), {a1, a2}, {r1, r2}, util);
    std::printf("\n=== MLE Estimator ===\n");
    const Matrix<double> mle = run_observation_estimator(
        m, "mle", queue, timestamps(n), {a1, a2}, {r1, r2}, util, {0.1, 0.3});
    std::printf("\n=== Comparison ===\n");
    std::printf("True demands: Class1=0.1000, Class2=0.3000\n");
    std::printf("UBO:  Class1=%.4f, Class2=%.4f\n", ubo(0, 0), ubo(0, 1));
    std::printf("MLE:  Class1=%.4f, Class2=%.4f\n", mle(0, 0), mle(0, 1));
    solve_after_estimate(m);
}

// ---------------------------------------------------------------------------
// UBR
// ---------------------------------------------------------------------------

/** `est_ubr_changepoint.py`: a sliding UBR window over a demand change point. */
void est_ubr_changepoint() {
    std::size_t delay = 0, queue = 0;
    Net m = open4(SchedStrategy::FCFS, {D::exp_mean(1.0)}, {Exp(1.0)}, delay, queue);
    model_line(m);

    const std::size_t n = 200, changeT = 100, W = 30;
    const double D1 = 0.3, D2 = 0.6, lam = 1.0;
    // Deterministic part only: the reference adds randn*0.02 from an UNSEEDED stream.
    const std::vector<double> arvr = filled(n, lam);
    std::vector<double> util(n, 0.0);
    for (std::size_t t = 0; t < n; ++t) util[t] = lam * (t < changeT ? D1 : D2);
    sampled("ArvR", "Queue1", "Class1", arvr);
    sampled("Util", "Queue1", "(aggr)", util);

    std::vector<double> estimates(n - W + 1, 0.0);
    for (std::size_t end = W; end <= n; ++end) {
        const std::vector<double> tw = timestamps(W, end - W);
        const std::vector<double> aw(arvr.begin() + static_cast<std::ptrdiff_t>(end - W),
                                     arvr.begin() + static_cast<std::ptrdiff_t>(end));
        const std::vector<double> uw(util.begin() + static_cast<std::ptrdiff_t>(end - W),
                                     util.begin() + static_cast<std::ptrdiff_t>(end));
        estimates[end - W] = run_ubr(m, queue, tw, {aw}, uw)(0, 0);
    }
    std::printf("\n=== UBR Change-Point Detection ===\n");
    std::printf("True demand: D=%g (t<=%zu), D=%g (t>%zu)\n", D1, changeT, D2, changeT);
    std::printf("Window size: %zu\n\n", W);
    const std::size_t cp[10] = {W, 50, 80, 100, 110, 120, 130, 150, 180, 200};
    std::printf("%6s  %10s  %10s\n", "t", "Estimated", "True");
    std::printf("%6s  %10s  %10s\n", "------", "----------", "----------");
    for (std::size_t i = 0; i < 10; ++i)
        std::printf("%6zu  %10.4f  %10.4f\n", cp[i], estimates[cp[i] - W],
                    cp[i] <= changeT ? D1 : D2);
}

/** `est_ubr_closed.py`: UBR from per-class and aggregate utilization. */
void est_ubr_closed() {
    std::size_t delay = 0, queue = 0;
    Net m = closed2({1.0, 2.0}, delay, queue);
    model_line(m);

    const std::size_t n = 1000;
    // np.random.seed(1): 2 - rand*0.15, 3 - rand*0.25, 1 - rand*0.05, 0.4*(2 - rand*0.15).
    const std::vector<double> a1 = filled(n, 1.925), a2 = filled(n, 2.875);
    const std::vector<double> util = filled(n, 0.975), util1 = filled(n, 0.4 * 1.925);
    sampled("ArvR", "Queue1", "Class1", a1);
    sampled("ArvR", "Queue1", "Class2", a2);
    sampled("Util", "Queue1", "(aggr)", util);
    sampled("Util", "Queue1", "Class1", util1);

    const Matrix<double> estimate = run_ubr(m, queue, timestamps(n), {a1, a2}, util,
                                             {util1, std::vector<double>()});
    std::printf("Estimated demands: [[%.8f %.8f]]\n", estimate(0, 0), estimate(0, 1));
    solve_after_estimate(m);
}

/** `est_ubr_open.py`: UBR on the open model with a hyperexponential think time. */
void est_ubr_open() {
    std::size_t delay = 0, queue = 0;
    Net m = open4(SchedStrategy::FCFS, {hyper(0.5, 3.0, 10.0)}, {Exp(0.1)}, delay, queue);
    model_line(m);

    const std::size_t n = 1000;
    // Deterministic part only: 2 - rand*0.15 and 1 - rand*0.05, UNSEEDED.
    const std::vector<double> arvr = filled(n, 1.925);
    const std::vector<double> util = filled(n, 0.975);
    std::vector<double> respt(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) respt[i] = 0.7 / (1.0 - util[i]);
    sampled("ArvR", "Queue1", "Class1", arvr);
    sampled("RespT", "Queue1", "Class1", respt);
    sampled("Util", "Queue1", "(aggr)", util);

    const Matrix<double> estimate = run_ubr(m, queue, timestamps(n), {arvr}, util);
    std::printf("Estimated demand: [[%.8f]]\n", estimate(0, 0));
    solve_after_estimate(m);
}

LINE_EXAMPLE("inference", est_ekf_changepoint);
LINE_EXAMPLE("inference", est_ekf_closed);
LINE_EXAMPLE("inference", est_ekf_open);
LINE_EXAMPLE("inference", est_erps_closed);
LINE_EXAMPLE("inference", est_erps_open);
LINE_EXAMPLE("inference", est_mcmc_closed);
LINE_EXAMPLE("inference", est_mle_open);
LINE_EXAMPLE("inference", est_qmle_closed);
LINE_EXAMPLE("inference", est_trace_closed);
LINE_EXAMPLE("inference", est_ubo_changepoint);
LINE_EXAMPLE("inference", est_ubo_closed);
LINE_EXAMPLE("inference", est_ubo_open);
LINE_EXAMPLE("inference", est_ubo_variants_closed);
LINE_EXAMPLE("inference", est_ubr_changepoint);
LINE_EXAMPLE("inference", est_ubr_closed);
LINE_EXAMPLE("inference", est_ubr_open);
LINE_EXAMPLE("inference", est_vi_closed);

}  // namespace examples
}  // namespace line
