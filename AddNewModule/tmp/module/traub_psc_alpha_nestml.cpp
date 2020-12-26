/*
*  traub_psc_alpha_nestml.cpp
*
*  This file is part of NEST.
*
*  Copyright (C) 2004 The NEST Initiative
*
*  NEST is free software: you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation, either version 2 of the License, or
*  (at your option) any later version.
*
*  NEST is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with NEST.  If not, see <http://www.gnu.org/licenses/>.
*
*  2020-03-26 06:36:32.839665
*/

// C++ includes:
#include <limits>

// Includes from libnestutil:
#include "numerics.h"

// Includes from nestkernel:
#include "exceptions.h"
#include "kernel_manager.h"
#include "universal_data_logger_impl.h"

// Includes from sli:
#include "dict.h"
#include "dictutils.h"
#include "doubledatum.h"
#include "integerdatum.h"
#include "lockptrdatum.h"

#include "traub_psc_alpha_nestml.h"


/* ----------------------------------------------------------------
* Recordables map
* ---------------------------------------------------------------- */
nest::RecordablesMap<traub_psc_alpha_nestml> traub_psc_alpha_nestml::recordablesMap_;

namespace nest
{
  // Override the create() method with one call to RecordablesMap::insert_()
  // for each quantity to be recorded.
  template <> void RecordablesMap<traub_psc_alpha_nestml>::create(){
  // use standard names whereever you can for consistency!  

  insert_("alpha_n_init", &traub_psc_alpha_nestml::get_alpha_n_init);

  insert_("beta_n_init", &traub_psc_alpha_nestml::get_beta_n_init);

  insert_("alpha_m_init", &traub_psc_alpha_nestml::get_alpha_m_init);

  insert_("beta_m_init", &traub_psc_alpha_nestml::get_beta_m_init);

  insert_("alpha_h_init", &traub_psc_alpha_nestml::get_alpha_h_init);

  insert_("beta_h_init", &traub_psc_alpha_nestml::get_beta_h_init);

  insert_("V_m", &traub_psc_alpha_nestml::get_V_m);

  insert_("Act_n", &traub_psc_alpha_nestml::get_Act_n);

  insert_("Act_m", &traub_psc_alpha_nestml::get_Act_m);

  insert_("Inact_h", &traub_psc_alpha_nestml::get_Inact_h);

  insert_("I_syn_in", &traub_psc_alpha_nestml::get_I_syn_in);

  insert_("I_syn_in__d", &traub_psc_alpha_nestml::get_I_syn_in__d);

  insert_("I_syn_ex", &traub_psc_alpha_nestml::get_I_syn_ex);

  insert_("I_syn_ex__d", &traub_psc_alpha_nestml::get_I_syn_ex__d);
  }
}

/* ----------------------------------------------------------------
 * Default constructors defining default parameters and state
 * Note: the implementation is empty. The initialization is of variables
 * is a part of the traub_psc_alpha_nestml's constructor.
 * ---------------------------------------------------------------- */
traub_psc_alpha_nestml::Parameters_::Parameters_(){}

traub_psc_alpha_nestml::State_::State_(){}

/* ----------------------------------------------------------------
* Parameter and state extractions and manipulation functions
* ---------------------------------------------------------------- */

traub_psc_alpha_nestml::Buffers_::Buffers_(traub_psc_alpha_nestml &n):
  logger_(n), __s( 0 ), __c( 0 ), __e( 0 ){
  // Initialization of the remaining members is deferred to
  // init_buffers_().
}

traub_psc_alpha_nestml::Buffers_::Buffers_(const Buffers_ &, traub_psc_alpha_nestml &n):
  logger_(n), __s( 0 ), __c( 0 ), __e( 0 ){
  // Initialization of the remaining members is deferred to
  // init_buffers_().
}

/* ----------------------------------------------------------------
 * Default and copy constructor for node, and destructor
 * ---------------------------------------------------------------- */
traub_psc_alpha_nestml::traub_psc_alpha_nestml():Archiving_Node(), P_(), S_(), B_(*this)
{
  recordablesMap_.create();
  // use a default `good` enough value for the absolute error.
  // it cab be adjusted via `SetStatus`
  P_.__gsl_error_tol = 1e-3;
  
  P_.t_ref = 2.0*1.0; // as ms
  
  P_.g_Na = 10000.0*1.0; // as nS
  
  P_.g_K = 8000.0*1.0; // as nS
  
  P_.g_L = 10*1.0; // as nS
  
  P_.C_m = 100.0*1.0; // as pF
  
  P_.E_Na = 50.0*1.0; // as mV
  
  P_.E_K = (-100.0*1.0); // as mV
  
  P_.E_L = (-67.0*1.0); // as mV
  
  P_.V_Tr = (-20.0*1.0); // as mV
  
  P_.tau_syn_ex = 0.2*1.0; // as ms
  
  P_.tau_syn_in = 2.0*1.0; // as ms
  
  P_.I_e = 0*1.0; // as pA
  
  S_.r = 0; // as integer
  
  S_.alpha_n_init = 0.032 * (S_.ode_state[State_::V_m] / 1.0 + 52.0) / (1.0 - std::exp((-(S_.ode_state[State_::V_m] / 1.0 + 52.0)) / 5.0)); // as real
  
  S_.beta_n_init = 0.5 * std::exp((-(S_.ode_state[State_::V_m] / 1.0 + 57.0)) / 40.0); // as real
  
  S_.alpha_m_init = 0.32 * (S_.ode_state[State_::V_m] / 1.0 + 54.0) / (1.0 - std::exp((-(S_.ode_state[State_::V_m] / 1.0 + 54.0)) / 4.0)); // as real
  
  S_.beta_m_init = 0.28 * (S_.ode_state[State_::V_m] / 1.0 + 27.0) / (std::exp((S_.ode_state[State_::V_m] / 1.0 + 27.0) / 5.0) - 1.0); // as real
  
  S_.alpha_h_init = 0.128 * std::exp((-(S_.ode_state[State_::V_m] / 1.0 + 50.0)) / 18.0); // as real
  
  S_.beta_h_init = 4.0 / (1.0 + std::exp((-(S_.ode_state[State_::V_m] / 1.0 + 27.0)) / 5.0)); // as real
  
  S_.ode_state[State_::V_m] = (-70.0*1.0); // as mV
  
  S_.ode_state[State_::Act_n] = get_alpha_n_init() / (get_alpha_n_init() + get_beta_n_init()); // as real
  
  S_.ode_state[State_::Act_m] = get_alpha_m_init() / (get_alpha_m_init() + get_beta_m_init()); // as real
  
  S_.ode_state[State_::Inact_h] = get_alpha_h_init() / (get_alpha_h_init() + get_beta_h_init()); // as real
  
  S_.ode_state[State_::I_syn_in] = 0; // as real
  
  S_.ode_state[State_::I_syn_in__d] = 0; // as real
  
  S_.ode_state[State_::I_syn_ex] = 0; // as real
  
  S_.ode_state[State_::I_syn_ex__d] = 0; // as real
}

traub_psc_alpha_nestml::traub_psc_alpha_nestml(const traub_psc_alpha_nestml& __n):
  Archiving_Node(), P_(__n.P_), S_(__n.S_), B_(__n.B_, *this){
  P_.t_ref = __n.P_.t_ref;
  P_.g_Na = __n.P_.g_Na;
  P_.g_K = __n.P_.g_K;
  P_.g_L = __n.P_.g_L;
  P_.C_m = __n.P_.C_m;
  P_.E_Na = __n.P_.E_Na;
  P_.E_K = __n.P_.E_K;
  P_.E_L = __n.P_.E_L;
  P_.V_Tr = __n.P_.V_Tr;
  P_.tau_syn_ex = __n.P_.tau_syn_ex;
  P_.tau_syn_in = __n.P_.tau_syn_in;
  P_.I_e = __n.P_.I_e;
  
  S_.r = __n.S_.r;
  
  S_.ode_state[State_::V_m] = __n.S_.ode_state[State_::V_m];
  S_.ode_state[State_::Act_n] = __n.S_.ode_state[State_::Act_n];
  S_.ode_state[State_::Act_m] = __n.S_.ode_state[State_::Act_m];
  S_.ode_state[State_::Inact_h] = __n.S_.ode_state[State_::Inact_h];
  S_.ode_state[State_::I_syn_in] = __n.S_.ode_state[State_::I_syn_in];
  S_.ode_state[State_::I_syn_in__d] = __n.S_.ode_state[State_::I_syn_in__d];
  S_.ode_state[State_::I_syn_ex] = __n.S_.ode_state[State_::I_syn_ex];
  S_.ode_state[State_::I_syn_ex__d] = __n.S_.ode_state[State_::I_syn_ex__d];
  
  V_.RefractoryCounts = __n.V_.RefractoryCounts;
  
}

traub_psc_alpha_nestml::~traub_psc_alpha_nestml(){ 
  // GSL structs may not have been allocated, so we need to protect destruction
  if (B_.__s)
    gsl_odeiv_step_free( B_.__s );
  if (B_.__c)
    gsl_odeiv_control_free( B_.__c );
  if (B_.__e)
    gsl_odeiv_evolve_free( B_.__e );
}

/* ----------------------------------------------------------------
* Node initialization functions
* ---------------------------------------------------------------- */

void traub_psc_alpha_nestml::init_state_(const Node& proto){
  const traub_psc_alpha_nestml& pr = downcast<traub_psc_alpha_nestml>(proto);
  S_ = pr.S_;
}



extern "C" inline int traub_psc_alpha_nestml_dynamics(double, const double ode_state[], double f[], void* pnode){
  typedef traub_psc_alpha_nestml::State_ State_;
  // get access to node so we can almost work as in a member function
  assert( pnode );
  const traub_psc_alpha_nestml& node = *( reinterpret_cast< traub_psc_alpha_nestml* >( pnode ) );

  // ode_state[] here is---and must be---the state vector supplied by the integrator,
  // not the state vector in the node, node.S_.ode_state[].
  
  f[State_::V_m] = (-(((node.get_g_Na() * ode_state[State_::Act_m] * ode_state[State_::Act_m] * ode_state[State_::Act_m] * ode_state[State_::Inact_h] * (ode_state[State_::V_m] - node.get_E_Na())) + (node.get_g_K() * ode_state[State_::Act_n] * ode_state[State_::Act_n] * ode_state[State_::Act_n] * ode_state[State_::Act_n] * (ode_state[State_::V_m] - node.get_E_K())) + (node.get_g_L() * (ode_state[State_::V_m] - node.get_E_L())))) + node.get_I_e() + node.B_.I_stim_grid_sum_ + (ode_state[State_::I_syn_in]) + (ode_state[State_::I_syn_ex])) / node.get_C_m();
  f[State_::Act_n] = ((0.032 * (ode_state[State_::V_m] / 1.0 + 52.0) / (1.0 - std::exp(-((ode_state[State_::V_m] / 1.0 + 52.0)) / 5.0))) * (1 - ode_state[State_::Act_n]) - (0.5 * std::exp(-((ode_state[State_::V_m] / 1.0 + 57.0)) / 40.0)) * ode_state[State_::Act_n]) / 1.0;
  f[State_::Act_m] = ((0.32 * (ode_state[State_::V_m] / 1.0 + 54.0) / (1.0 - std::exp(-((ode_state[State_::V_m] / 1.0 + 54.0)) / 4.0))) * (1 - ode_state[State_::Act_m]) - (0.28 * (ode_state[State_::V_m] / 1.0 + 27.0) / (std::exp((ode_state[State_::V_m] / 1.0 + 27.0) / 5.0) - 1.0)) * ode_state[State_::Act_m]) / 1.0;
  f[State_::Inact_h] = ((0.128 * std::exp(-((ode_state[State_::V_m] / 1.0 + 50.0)) / 18.0)) * (1 - ode_state[State_::Inact_h]) - (4.0 / (1.0 + std::exp(-((ode_state[State_::V_m] / 1.0 + 27.0)) / 5.0))) * ode_state[State_::Inact_h]) / 1.0;
  f[State_::I_syn_in] = ode_state[State_::I_syn_in__d];
  f[State_::I_syn_in__d] = -(1) / pow(node.get_tau_syn_in(), 2) * ode_state[State_::I_syn_in] + -(2) / node.get_tau_syn_in() * ode_state[State_::I_syn_in__d];
  f[State_::I_syn_ex] = ode_state[State_::I_syn_ex__d];
  f[State_::I_syn_ex__d] = -(1) / pow(node.get_tau_syn_ex(), 2) * ode_state[State_::I_syn_ex] + -(2) / node.get_tau_syn_ex() * ode_state[State_::I_syn_ex__d];
  return GSL_SUCCESS;
}



void traub_psc_alpha_nestml::init_buffers_(){
  get_spikeInh().clear(); //includes resize
  get_spikeExc().clear(); //includes resize
  get_I_stim().clear(); //includes resize
  
  B_.logger_.reset(); // includes resize
  Archiving_Node::clear_history();
  
  if ( B_.__s == 0 ){
    B_.__s = gsl_odeiv_step_alloc( gsl_odeiv_step_rkf45, 8 );
  } else {
    gsl_odeiv_step_reset( B_.__s );
  }

  if ( B_.__c == 0 ){
    B_.__c = gsl_odeiv_control_y_new( P_.__gsl_error_tol, 0.0 );
  } else {
    gsl_odeiv_control_init( B_.__c, P_.__gsl_error_tol, 0.0, 1.0, 0.0 );
  }

  if ( B_.__e == 0 ){
    B_.__e = gsl_odeiv_evolve_alloc( 8 );
  } else {
    gsl_odeiv_evolve_reset( B_.__e );
  }

  B_.__sys.function = traub_psc_alpha_nestml_dynamics;
  B_.__sys.jacobian = NULL;
  B_.__sys.dimension = 8;
  B_.__sys.params = reinterpret_cast< void* >( this );
  B_.__step = nest::Time::get_resolution().get_ms();
  B_.__integration_step = nest::Time::get_resolution().get_ms();
}

void traub_psc_alpha_nestml::calibrate(){
  B_.logger_.init();
  
  
  V_.RefractoryCounts =nest::Time(nest::Time::ms((double) (P_.t_ref))).get_steps();
}

/* ----------------------------------------------------------------
* Update and spike handling functions
* ---------------------------------------------------------------- */

/*
 *
 */
void traub_psc_alpha_nestml::update(nest::Time const & origin,const long from, const long to){
  double __t = 0;

  for ( long lag = from ; lag < to ; ++lag ) {
    B_.spikeInh_grid_sum_ = get_spikeInh().get_value(lag);
    B_.spikeExc_grid_sum_ = get_spikeExc().get_value(lag);
    B_.I_stim_grid_sum_ = get_I_stim().get_value(lag);
      
    


    double U_old = S_.ode_state[State_::V_m];
    

    __t = 0;
    // numerical integration with adaptive step size control:
    // ------------------------------------------------------
    // gsl_odeiv_evolve_apply performs only a single numerical
    // integration step, starting from t and bounded by step;
    // the while-loop ensures integration over the whole simulation
    // step (0, step] if more than one integration step is needed due
    // to a small integration step size;
    // note that (t+IntegrationStep > step) leads to integration over
    // (t, step] and afterwards setting t to step, but it does not
    // enforce setting IntegrationStep to step-t; this is of advantage
    // for a consistent and efficient integration across subsequent
    // simulation intervals
    while ( __t < B_.__step )
    {
      const int status = gsl_odeiv_evolve_apply(B_.__e,
                                                B_.__c,
                                                B_.__s,
                                                &B_.__sys,              // system of ODE
                                                &__t,                   // from t
                                                B_.__step,              // to t <= step
                                                &B_.__integration_step, // integration step size
                                                S_.ode_state);          // neuronal state

      if ( status != GSL_SUCCESS ) {
        throw nest::GSLSolverFailure( get_name(), status );
      }
    }
    



    if (S_.r>0) {
      

      S_.r -= 1;
    }else if(S_.ode_state[State_::V_m]>P_.V_Tr&&U_old>P_.V_Tr) {
      

      S_.r = V_.RefractoryCounts;
      
      set_spiketime(nest::Time::step(origin.get_steps()+lag+1));
      nest::SpikeEvent se;
      nest::kernel().event_delivery_manager.send(*this, se, lag);
    } /* if end */

    

    S_.ode_state[State_::I_syn_in] += (B_.spikeInh_grid_sum_ / 1.0) * 0;
    

    S_.ode_state[State_::I_syn_in__d] += (B_.spikeInh_grid_sum_ / 1.0) * numerics::e / P_.tau_syn_in;
    

    S_.ode_state[State_::I_syn_ex] += (B_.spikeExc_grid_sum_ / 1.0) * 0;
    

    S_.ode_state[State_::I_syn_ex__d] += (B_.spikeExc_grid_sum_ / 1.0) * numerics::e / P_.tau_syn_ex;

    // voltage logging
    B_.logger_.record_data(origin.get_steps()+lag);
  }

}

// Do not move this function as inline to h-file. It depends on
// universal_data_logger_impl.h being included here.
void traub_psc_alpha_nestml::handle(nest::DataLoggingRequest& e){
  B_.logger_.handle(e);
}


void traub_psc_alpha_nestml::handle(nest::SpikeEvent &e){
  assert(e.get_delay_steps() > 0);
  
  const double weight = e.get_weight();
  const double multiplicity = e.get_multiplicity();
  
  if ( weight < 0.0 ){ // inhibitory
    get_spikeInh().
        add_value(e.get_rel_delivery_steps( nest::kernel().simulation_manager.get_slice_origin()),
                      
                       weight * multiplicity );
  }
  if ( weight >= 0.0 ){ // excitatory
    get_spikeExc().
        add_value(e.get_rel_delivery_steps( nest::kernel().simulation_manager.get_slice_origin()),
                       weight * multiplicity );
  }
}

void traub_psc_alpha_nestml::handle(nest::CurrentEvent& e){
  assert(e.get_delay_steps() > 0);

  const double current = e.get_current();		// we assume that in NEST, this returns a current in pA
  const double weight = e.get_weight();
  get_I_stim().add_value(
               e.get_rel_delivery_steps( nest::kernel().simulation_manager.get_slice_origin()),
               weight * current );
  
}
