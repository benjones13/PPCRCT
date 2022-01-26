// Generated by rstantools.  Do not edit by hand.

/*
    PPCRCT is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    PPCRCT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with PPCRCT.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef MODELS_HPP
#define MODELS_HPP
#define STAN__SERVICES__COMMAND_HPP
#include <rstan/rstaninc.hpp>
// Code generated by Stan version 2.21.0
#include <stan/model/model_header.hpp>
namespace model_Hier_PP_HistoricOnly_namespace {
using std::istream;
using std::string;
using std::stringstream;
using std::vector;
using stan::io::dump;
using stan::math::lgamma;
using stan::model::prob_grad;
using namespace stan::math;
static int current_statement_begin__;
stan::io::program_reader prog_reader__() {
    stan::io::program_reader reader;
    reader.add_event(0, 0, "start", "model_Hier_PP_HistoricOnly");
    reader.add_event(35, 33, "end", "model_Hier_PP_HistoricOnly");
    return reader;
}
#include <stan_meta_header.hpp>
class model_Hier_PP_HistoricOnly
  : public stan::model::model_base_crtp<model_Hier_PP_HistoricOnly> {
private:
        int N0;
        int J0;
        int P;
        std::vector<double> y0;
        std::vector<int> SchoolCode0;
        matrix_d X0;
        double a0;
public:
    model_Hier_PP_HistoricOnly(stan::io::var_context& context__,
        std::ostream* pstream__ = 0)
        : model_base_crtp(0) {
        ctor_body(context__, 0, pstream__);
    }
    model_Hier_PP_HistoricOnly(stan::io::var_context& context__,
        unsigned int random_seed__,
        std::ostream* pstream__ = 0)
        : model_base_crtp(0) {
        ctor_body(context__, random_seed__, pstream__);
    }
    void ctor_body(stan::io::var_context& context__,
                   unsigned int random_seed__,
                   std::ostream* pstream__) {
        typedef double local_scalar_t__;
        boost::ecuyer1988 base_rng__ =
          stan::services::util::create_rng(random_seed__, 0);
        (void) base_rng__;  // suppress unused var warning
        current_statement_begin__ = -1;
        static const char* function__ = "model_Hier_PP_HistoricOnly_namespace::model_Hier_PP_HistoricOnly";
        (void) function__;  // dummy to suppress unused var warning
        size_t pos__;
        (void) pos__;  // dummy to suppress unused var warning
        std::vector<int> vals_i__;
        std::vector<double> vals_r__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        try {
            // initialize data block variables from context__
            current_statement_begin__ = 2;
            context__.validate_dims("data initialization", "N0", "int", context__.to_vec());
            N0 = int(0);
            vals_i__ = context__.vals_i("N0");
            pos__ = 0;
            N0 = vals_i__[pos__++];
            check_greater_or_equal(function__, "N0", N0, 1);
            current_statement_begin__ = 3;
            context__.validate_dims("data initialization", "J0", "int", context__.to_vec());
            J0 = int(0);
            vals_i__ = context__.vals_i("J0");
            pos__ = 0;
            J0 = vals_i__[pos__++];
            check_greater_or_equal(function__, "J0", J0, 1);
            current_statement_begin__ = 4;
            context__.validate_dims("data initialization", "P", "int", context__.to_vec());
            P = int(0);
            vals_i__ = context__.vals_i("P");
            pos__ = 0;
            P = vals_i__[pos__++];
            check_greater_or_equal(function__, "P", P, 1);
            current_statement_begin__ = 5;
            validate_non_negative_index("y0", "N0", N0);
            context__.validate_dims("data initialization", "y0", "double", context__.to_vec(N0));
            y0 = std::vector<double>(N0, double(0));
            vals_r__ = context__.vals_r("y0");
            pos__ = 0;
            size_t y0_k_0_max__ = N0;
            for (size_t k_0__ = 0; k_0__ < y0_k_0_max__; ++k_0__) {
                y0[k_0__] = vals_r__[pos__++];
            }
            current_statement_begin__ = 6;
            validate_non_negative_index("SchoolCode0", "N0", N0);
            context__.validate_dims("data initialization", "SchoolCode0", "int", context__.to_vec(N0));
            SchoolCode0 = std::vector<int>(N0, int(0));
            vals_i__ = context__.vals_i("SchoolCode0");
            pos__ = 0;
            size_t SchoolCode0_k_0_max__ = N0;
            for (size_t k_0__ = 0; k_0__ < SchoolCode0_k_0_max__; ++k_0__) {
                SchoolCode0[k_0__] = vals_i__[pos__++];
            }
            size_t SchoolCode0_i_0_max__ = N0;
            for (size_t i_0__ = 0; i_0__ < SchoolCode0_i_0_max__; ++i_0__) {
                check_greater_or_equal(function__, "SchoolCode0[i_0__]", SchoolCode0[i_0__], 0);
                check_less_or_equal(function__, "SchoolCode0[i_0__]", SchoolCode0[i_0__], J0);
            }
            current_statement_begin__ = 7;
            validate_non_negative_index("X0", "N0", N0);
            validate_non_negative_index("X0", "P", P);
            context__.validate_dims("data initialization", "X0", "matrix_d", context__.to_vec(N0,P));
            X0 = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(N0, P);
            vals_r__ = context__.vals_r("X0");
            pos__ = 0;
            size_t X0_j_2_max__ = P;
            size_t X0_j_1_max__ = N0;
            for (size_t j_2__ = 0; j_2__ < X0_j_2_max__; ++j_2__) {
                for (size_t j_1__ = 0; j_1__ < X0_j_1_max__; ++j_1__) {
                    X0(j_1__, j_2__) = vals_r__[pos__++];
                }
            }
            current_statement_begin__ = 8;
            context__.validate_dims("data initialization", "a0", "double", context__.to_vec());
            a0 = double(0);
            vals_r__ = context__.vals_r("a0");
            pos__ = 0;
            a0 = vals_r__[pos__++];
            check_greater_or_equal(function__, "a0", a0, 0);
            check_less_or_equal(function__, "a0", a0, 1);
            // initialize transformed data variables
            // execute transformed data statements
            // validate transformed data
            // validate, set parameter ranges
            num_params_r__ = 0U;
            param_ranges_i__.clear();
            current_statement_begin__ = 13;
            num_params_r__ += 1;
            current_statement_begin__ = 14;
            num_params_r__ += 1;
            current_statement_begin__ = 15;
            num_params_r__ += 1;
            current_statement_begin__ = 16;
            validate_non_negative_index("beta", "P", P);
            num_params_r__ += P;
            current_statement_begin__ = 17;
            validate_non_negative_index("eta0_raw", "J0", J0);
            num_params_r__ += J0;
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    ~model_Hier_PP_HistoricOnly() { }
    void transform_inits(const stan::io::var_context& context__,
                         std::vector<int>& params_i__,
                         std::vector<double>& params_r__,
                         std::ostream* pstream__) const {
        typedef double local_scalar_t__;
        stan::io::writer<double> writer__(params_r__, params_i__);
        size_t pos__;
        (void) pos__; // dummy call to supress warning
        std::vector<double> vals_r__;
        std::vector<int> vals_i__;
        current_statement_begin__ = 13;
        if (!(context__.contains_r("sigma")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable sigma missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("sigma");
        pos__ = 0U;
        context__.validate_dims("parameter initialization", "sigma", "double", context__.to_vec());
        double sigma(0);
        sigma = vals_r__[pos__++];
        try {
            writer__.scalar_lb_unconstrain(0, sigma);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable sigma: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 14;
        if (!(context__.contains_r("alpha")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable alpha missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("alpha");
        pos__ = 0U;
        context__.validate_dims("parameter initialization", "alpha", "double", context__.to_vec());
        double alpha(0);
        alpha = vals_r__[pos__++];
        try {
            writer__.scalar_unconstrain(alpha);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable alpha: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 15;
        if (!(context__.contains_r("sigma_eta")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable sigma_eta missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("sigma_eta");
        pos__ = 0U;
        context__.validate_dims("parameter initialization", "sigma_eta", "double", context__.to_vec());
        double sigma_eta(0);
        sigma_eta = vals_r__[pos__++];
        try {
            writer__.scalar_lb_unconstrain(0, sigma_eta);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable sigma_eta: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 16;
        if (!(context__.contains_r("beta")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable beta missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("beta");
        pos__ = 0U;
        validate_non_negative_index("beta", "P", P);
        context__.validate_dims("parameter initialization", "beta", "vector_d", context__.to_vec(P));
        Eigen::Matrix<double, Eigen::Dynamic, 1> beta(P);
        size_t beta_j_1_max__ = P;
        for (size_t j_1__ = 0; j_1__ < beta_j_1_max__; ++j_1__) {
            beta(j_1__) = vals_r__[pos__++];
        }
        try {
            writer__.vector_unconstrain(beta);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable beta: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 17;
        if (!(context__.contains_r("eta0_raw")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable eta0_raw missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("eta0_raw");
        pos__ = 0U;
        validate_non_negative_index("eta0_raw", "J0", J0);
        context__.validate_dims("parameter initialization", "eta0_raw", "vector_d", context__.to_vec(J0));
        Eigen::Matrix<double, Eigen::Dynamic, 1> eta0_raw(J0);
        size_t eta0_raw_j_1_max__ = J0;
        for (size_t j_1__ = 0; j_1__ < eta0_raw_j_1_max__; ++j_1__) {
            eta0_raw(j_1__) = vals_r__[pos__++];
        }
        try {
            writer__.vector_unconstrain(eta0_raw);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable eta0_raw: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        params_r__ = writer__.data_r();
        params_i__ = writer__.data_i();
    }
    void transform_inits(const stan::io::var_context& context,
                         Eigen::Matrix<double, Eigen::Dynamic, 1>& params_r,
                         std::ostream* pstream__) const {
      std::vector<double> params_r_vec;
      std::vector<int> params_i_vec;
      transform_inits(context, params_i_vec, params_r_vec, pstream__);
      params_r.resize(params_r_vec.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r(i) = params_r_vec[i];
    }
    template <bool propto__, bool jacobian__, typename T__>
    T__ log_prob(std::vector<T__>& params_r__,
                 std::vector<int>& params_i__,
                 std::ostream* pstream__ = 0) const {
        typedef T__ local_scalar_t__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // dummy to suppress unused var warning
        T__ lp__(0.0);
        stan::math::accumulator<T__> lp_accum__;
        try {
            stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
            // model parameters
            current_statement_begin__ = 13;
            local_scalar_t__ sigma;
            (void) sigma;  // dummy to suppress unused var warning
            if (jacobian__)
                sigma = in__.scalar_lb_constrain(0, lp__);
            else
                sigma = in__.scalar_lb_constrain(0);
            current_statement_begin__ = 14;
            local_scalar_t__ alpha;
            (void) alpha;  // dummy to suppress unused var warning
            if (jacobian__)
                alpha = in__.scalar_constrain(lp__);
            else
                alpha = in__.scalar_constrain();
            current_statement_begin__ = 15;
            local_scalar_t__ sigma_eta;
            (void) sigma_eta;  // dummy to suppress unused var warning
            if (jacobian__)
                sigma_eta = in__.scalar_lb_constrain(0, lp__);
            else
                sigma_eta = in__.scalar_lb_constrain(0);
            current_statement_begin__ = 16;
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> beta;
            (void) beta;  // dummy to suppress unused var warning
            if (jacobian__)
                beta = in__.vector_constrain(P, lp__);
            else
                beta = in__.vector_constrain(P);
            current_statement_begin__ = 17;
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> eta0_raw;
            (void) eta0_raw;  // dummy to suppress unused var warning
            if (jacobian__)
                eta0_raw = in__.vector_constrain(J0, lp__);
            else
                eta0_raw = in__.vector_constrain(J0);
            // transformed parameters
            current_statement_begin__ = 21;
            validate_non_negative_index("eta0", "J0", J0);
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> eta0(J0);
            stan::math::initialize(eta0, DUMMY_VAR__);
            stan::math::fill(eta0, DUMMY_VAR__);
            stan::math::assign(eta0,multiply(sigma_eta, eta0_raw));
            // validate transformed parameters
            const char* function__ = "validate transformed params";
            (void) function__;  // dummy to suppress unused var warning
            current_statement_begin__ = 21;
            size_t eta0_j_1_max__ = J0;
            for (size_t j_1__ = 0; j_1__ < eta0_j_1_max__; ++j_1__) {
                if (stan::math::is_uninitialized(eta0(j_1__))) {
                    std::stringstream msg__;
                    msg__ << "Undefined transformed parameter: eta0" << "(" << j_1__ << ")";
                    stan::lang::rethrow_located(std::runtime_error(std::string("Error initializing variable eta0: ") + msg__.str()), current_statement_begin__, prog_reader__());
                }
            }
            // model body
            {
            current_statement_begin__ = 25;
            validate_non_negative_index("mu0", "N0", N0);
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> mu0(N0);
            stan::math::initialize(mu0, DUMMY_VAR__);
            stan::math::fill(mu0, DUMMY_VAR__);
            stan::math::assign(mu0,add(add(stan::model::rvalue(eta0, stan::model::cons_list(stan::model::index_multi(SchoolCode0), stan::model::nil_index_list()), "eta0"), alpha), multiply(X0, beta)));
            current_statement_begin__ = 27;
            lp_accum__.add((a0 * normal_log(y0, mu0, sigma)));
            current_statement_begin__ = 28;
            lp_accum__.add(normal_log(eta0_raw, 0, 1));
            current_statement_begin__ = 29;
            lp_accum__.add(normal_log(alpha, 0, 5));
            current_statement_begin__ = 30;
            lp_accum__.add(normal_log(beta, 0, 5));
            current_statement_begin__ = 31;
            lp_accum__.add(exponential_log(sigma, 1));
            current_statement_begin__ = 32;
            lp_accum__.add((cauchy_log(sigma_eta, 0, .3) - cauchy_ccdf_log(0, 0, .3)));
            }
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
        lp_accum__.add(lp__);
        return lp_accum__.sum();
    } // log_prob()
    template <bool propto, bool jacobian, typename T_>
    T_ log_prob(Eigen::Matrix<T_,Eigen::Dynamic,1>& params_r,
               std::ostream* pstream = 0) const {
      std::vector<T_> vec_params_r;
      vec_params_r.reserve(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        vec_params_r.push_back(params_r(i));
      std::vector<int> vec_params_i;
      return log_prob<propto,jacobian,T_>(vec_params_r, vec_params_i, pstream);
    }
    void get_param_names(std::vector<std::string>& names__) const {
        names__.resize(0);
        names__.push_back("sigma");
        names__.push_back("alpha");
        names__.push_back("sigma_eta");
        names__.push_back("beta");
        names__.push_back("eta0_raw");
        names__.push_back("eta0");
    }
    void get_dims(std::vector<std::vector<size_t> >& dimss__) const {
        dimss__.resize(0);
        std::vector<size_t> dims__;
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(P);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(J0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(J0);
        dimss__.push_back(dims__);
    }
    template <typename RNG>
    void write_array(RNG& base_rng__,
                     std::vector<double>& params_r__,
                     std::vector<int>& params_i__,
                     std::vector<double>& vars__,
                     bool include_tparams__ = true,
                     bool include_gqs__ = true,
                     std::ostream* pstream__ = 0) const {
        typedef double local_scalar_t__;
        vars__.resize(0);
        stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
        static const char* function__ = "model_Hier_PP_HistoricOnly_namespace::write_array";
        (void) function__;  // dummy to suppress unused var warning
        // read-transform, write parameters
        double sigma = in__.scalar_lb_constrain(0);
        vars__.push_back(sigma);
        double alpha = in__.scalar_constrain();
        vars__.push_back(alpha);
        double sigma_eta = in__.scalar_lb_constrain(0);
        vars__.push_back(sigma_eta);
        Eigen::Matrix<double, Eigen::Dynamic, 1> beta = in__.vector_constrain(P);
        size_t beta_j_1_max__ = P;
        for (size_t j_1__ = 0; j_1__ < beta_j_1_max__; ++j_1__) {
            vars__.push_back(beta(j_1__));
        }
        Eigen::Matrix<double, Eigen::Dynamic, 1> eta0_raw = in__.vector_constrain(J0);
        size_t eta0_raw_j_1_max__ = J0;
        for (size_t j_1__ = 0; j_1__ < eta0_raw_j_1_max__; ++j_1__) {
            vars__.push_back(eta0_raw(j_1__));
        }
        double lp__ = 0.0;
        (void) lp__;  // dummy to suppress unused var warning
        stan::math::accumulator<double> lp_accum__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        if (!include_tparams__ && !include_gqs__) return;
        try {
            // declare and define transformed parameters
            current_statement_begin__ = 21;
            validate_non_negative_index("eta0", "J0", J0);
            Eigen::Matrix<double, Eigen::Dynamic, 1> eta0(J0);
            stan::math::initialize(eta0, DUMMY_VAR__);
            stan::math::fill(eta0, DUMMY_VAR__);
            stan::math::assign(eta0,multiply(sigma_eta, eta0_raw));
            if (!include_gqs__ && !include_tparams__) return;
            // validate transformed parameters
            const char* function__ = "validate transformed params";
            (void) function__;  // dummy to suppress unused var warning
            // write transformed parameters
            if (include_tparams__) {
                size_t eta0_j_1_max__ = J0;
                for (size_t j_1__ = 0; j_1__ < eta0_j_1_max__; ++j_1__) {
                    vars__.push_back(eta0(j_1__));
                }
            }
            if (!include_gqs__) return;
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    template <typename RNG>
    void write_array(RNG& base_rng,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& vars,
                     bool include_tparams = true,
                     bool include_gqs = true,
                     std::ostream* pstream = 0) const {
      std::vector<double> params_r_vec(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r_vec[i] = params_r(i);
      std::vector<double> vars_vec;
      std::vector<int> params_i_vec;
      write_array(base_rng, params_r_vec, params_i_vec, vars_vec, include_tparams, include_gqs, pstream);
      vars.resize(vars_vec.size());
      for (int i = 0; i < vars.size(); ++i)
        vars(i) = vars_vec[i];
    }
    std::string model_name() const {
        return "model_Hier_PP_HistoricOnly";
    }
    void constrained_param_names(std::vector<std::string>& param_names__,
                                 bool include_tparams__ = true,
                                 bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        param_name_stream__.str(std::string());
        param_name_stream__ << "sigma";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "alpha";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "sigma_eta";
        param_names__.push_back(param_name_stream__.str());
        size_t beta_j_1_max__ = P;
        for (size_t j_1__ = 0; j_1__ < beta_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "beta" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t eta0_raw_j_1_max__ = J0;
        for (size_t j_1__ = 0; j_1__ < eta0_raw_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "eta0_raw" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
            size_t eta0_j_1_max__ = J0;
            for (size_t j_1__ = 0; j_1__ < eta0_j_1_max__; ++j_1__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "eta0" << '.' << j_1__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
        }
        if (!include_gqs__) return;
    }
    void unconstrained_param_names(std::vector<std::string>& param_names__,
                                   bool include_tparams__ = true,
                                   bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        param_name_stream__.str(std::string());
        param_name_stream__ << "sigma";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "alpha";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "sigma_eta";
        param_names__.push_back(param_name_stream__.str());
        size_t beta_j_1_max__ = P;
        for (size_t j_1__ = 0; j_1__ < beta_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "beta" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t eta0_raw_j_1_max__ = J0;
        for (size_t j_1__ = 0; j_1__ < eta0_raw_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "eta0_raw" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
            size_t eta0_j_1_max__ = J0;
            for (size_t j_1__ = 0; j_1__ < eta0_j_1_max__; ++j_1__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "eta0" << '.' << j_1__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
        }
        if (!include_gqs__) return;
    }
}; // model
}  // namespace
typedef model_Hier_PP_HistoricOnly_namespace::model_Hier_PP_HistoricOnly stan_model;
#ifndef USING_R
stan::model::model_base& new_model(
        stan::io::var_context& data_context,
        unsigned int seed,
        std::ostream* msg_stream) {
  stan_model* m = new stan_model(data_context, seed, msg_stream);
  return *m;
}
#endif
#endif
