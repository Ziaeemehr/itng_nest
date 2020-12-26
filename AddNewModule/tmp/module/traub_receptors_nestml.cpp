/*
*  traub_receptors_nestml.cpp
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
*  2020-03-26 06:36:34.787825
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

#include "traub_receptors_nestml.h"


/* ----------------------------------------------------------------
* Recordables map
* ---------------------------------------------------------------- */
nest::RecordablesMap<traub_receptors_nestml> traub_receptors_nestml::recordablesMap_;

namespace nest
{
  // Override the create() method with one call to RecordablesMap::insert_()
  // for each quantity to be recorded.
  template <> void RecordablesMap<traub_receptors_nestml>::create(){
  // use standard names whereever you can for consistency!  

  insert_("alpha_n_init", &traub_receptors_nestml::get_alpha_n_init);

  insert_("beta_n_init", &traub_receptors_nestml::get_beta_n_init);

  insert_("alpha_m_init", &traub_receptors_nestml::get_alpha_m_init);

  insert_("beta_m_init", &traub_receptors_nestml::get_beta_m_init);

  insert_("alpha_h_init", &traub_receptors_nestml::get_alpha_h_init);

  insert_("beta_h_init", &traub_receptors_nestml::get_beta_h_init);

  insert_("V_m", &traub_receptors_nestml::get_V_m);

  insert_("Act_n", &traub_receptors_nestml::get_Act_n);

  insert_("Act_m", &traub_receptors_nestml::get_Act_m);

  insert_("Inact_h", &traub_receptors_nestml::get_Inact_h);

  insert_("g_AMPA__d", &traub_receptors_nestml::get_g_AMPA__d);

  insert_("g_AMPA", &traub_receptors_nestml::get_g_AMPA);

  insert_("g_NMDA__d", &traub_receptors_nestml::get_g_NMDA__d);

  insert_("g_NMDA", &traub_receptors_nestml::get_g_NMDA);

  insert_("g_GABAA__d", &traub_receptors_nestml::get_g_GABAA__d);

  insert_("g_GABAA", &traub_receptors_nestml::get_g_GABAA);

  insert_("g_GABAB__d", &traub_receptors_nestml::get_g_GABAB__d);

  insert_("g_GABAB", &traub_receptors_nestml::get_g_GABAB);

  insert_("I_syn_ampa", &traub_receptors_nestml::get_I_syn_ampa);

  insert_("I_syn_nmda", &traub_receptors_nestml::get_I_syn_nmda);

  insert_("I_syn_gaba_a", &traub_receptors_nestml::get_I_syn_gaba_a);

  insert_("I_syn_gaba_b", &traub_receptors_nestml::get_I_syn_gaba_b);

  insert_("I_syn", &traub_receptors_nestml::get_I_syn);
  }
}

/* ----------------------------------------------------------------
 * Default constructors defining default parameters and state
 * Note: the implementation is empty. The initialization is of variables
 * is a part of the traub_receptors_nestml's constructor.
 * ---------------------------------------------------------------- */
traub_receptors_nestml::Parameters_::Parameters_(){}

traub_receptors_nestml::State_::State_(){}

/* ----------------------------------------------------------------
* Parameter and state extractions and manipulation functions
* ---------------------------------------------------------------- */

traub_receptors_nestml::Buffers_::Buffers_(traub_receptors_nestml &n):
  logger_(n), spike_inputs_( std::vector< nest::RingBuffer >( SUP_SPIKE_RECEPTOR - 1 ) ), __s( 0 ), __c( 0 ), __e( 0 ){
  // Initialization of the remaining members is deferred to
  // init_buffers_().
}

traub_receptors_nestml::Buffers_::Buffers_(const Buffers_ &, traub_receptors_nestml &n):
  logger_(n), spike_inputs_( std::vector< nest::RingBuffer >( SUP_SPIKE_RECEPTOR - 1 ) ), __s( 0 ), __c( 0 ), __e( 0 ){
  // Initialization of the remaining members is deferred to
  // init_buffers_().
}

/* ----------------------------------------------------------------
 * Default and copy constructor for node, and destructor
 * ---------------------------------------------------------------- */
traub_receptors_nestml::traub_receptors_nestml():Archiving_Node(), P_(), S_(), B_(*this)
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
  
  P_.AMPA_g_peak = 0.1*1.0; // as nS
  
  P_.AMPA_E_rev = 0.0*1.0; // as mV
  
  P_.AMPA_Tau_1 = 0.5*1.0; // as ms
  
  P_.AMPA_Tau_2 = 2.4*1.0; // as ms
  
  P_.NMDA_g_peak = 0.075*1.0; // as nS
  
  P_.NMDA_Tau_1 = 4.0*1.0; // as ms
  
  P_.NMDA_Tau_2 = 40.0*1.0; // as ms
  
  P_.NMDA_E_rev = 0.0*1.0; // as mV
  
  P_.NMDA_Vact = (-58.0*1.0); // as mV
  
  P_.NMDA_Sact = 2.5*1.0; // as mV
  
  P_.GABA_A_g_peak = 0.33*1.0; // as nS
  
  P_.GABA_A_Tau_1 = 1.0*1.0; // as ms
  
  P_.GABA_A_Tau_2 = 7.0*1.0; // as ms
  
  P_.GABA_A_E_rev = (-70.0*1.0); // as mV
  
  P_.GABA_B_g_peak = 0.0132*1.0; // as nS
  
  P_.GABA_B_Tau_1 = 60.0*1.0; // as ms
  
  P_.GABA_B_Tau_2 = 200.0*1.0; // as ms
  
  P_.GABA_B_E_rev = (-90.0*1.0); // as mV
  
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
  
  S_.ode_state[State_::g_AMPA__d] = 0.0*1.0 / 1.0; // as nS / ms
  
  S_.ode_state[State_::g_AMPA] = 0.0*1.0; // as nS
  
  S_.ode_state[State_::g_NMDA__d] = 0.0*1.0 / 1.0; // as nS / ms
  
  S_.ode_state[State_::g_NMDA] = 0.0*1.0; // as nS
  
  S_.ode_state[State_::g_GABAA__d] = 0.0*1.0 / 1.0; // as nS / ms
  
  S_.ode_state[State_::g_GABAA] = 0.0*1.0; // as nS
  
  S_.ode_state[State_::g_GABAB__d] = 0.0*1.0 / 1.0; // as nS / ms
  
  S_.ode_state[State_::g_GABAB] = 0.0*1.0; // as nS
}

traub_receptors_nestml::traub_receptors_nestml(const traub_receptors_nestml& __n):
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
  P_.AMPA_g_peak = __n.P_.AMPA_g_peak;
  P_.AMPA_E_rev = __n.P_.AMPA_E_rev;
  P_.AMPA_Tau_1 = __n.P_.AMPA_Tau_1;
  P_.AMPA_Tau_2 = __n.P_.AMPA_Tau_2;
  P_.NMDA_g_peak = __n.P_.NMDA_g_peak;
  P_.NMDA_Tau_1 = __n.P_.NMDA_Tau_1;
  P_.NMDA_Tau_2 = __n.P_.NMDA_Tau_2;
  P_.NMDA_E_rev = __n.P_.NMDA_E_rev;
  P_.NMDA_Vact = __n.P_.NMDA_Vact;
  P_.NMDA_Sact = __n.P_.NMDA_Sact;
  P_.GABA_A_g_peak = __n.P_.GABA_A_g_peak;
  P_.GABA_A_Tau_1 = __n.P_.GABA_A_Tau_1;
  P_.GABA_A_Tau_2 = __n.P_.GABA_A_Tau_2;
  P_.GABA_A_E_rev = __n.P_.GABA_A_E_rev;
  P_.GABA_B_g_peak = __n.P_.GABA_B_g_peak;
  P_.GABA_B_Tau_1 = __n.P_.GABA_B_Tau_1;
  P_.GABA_B_Tau_2 = __n.P_.GABA_B_Tau_2;
  P_.GABA_B_E_rev = __n.P_.GABA_B_E_rev;
  P_.I_e = __n.P_.I_e;
  
  S_.r = __n.S_.r;
  
  S_.ode_state[State_::V_m] = __n.S_.ode_state[State_::V_m];
  S_.ode_state[State_::Act_n] = __n.S_.ode_state[State_::Act_n];
  S_.ode_state[State_::Act_m] = __n.S_.ode_state[State_::Act_m];
  S_.ode_state[State_::Inact_h] = __n.S_.ode_state[State_::Inact_h];
  S_.ode_state[State_::g_AMPA__d] = __n.S_.ode_state[State_::g_AMPA__d];
  S_.ode_state[State_::g_AMPA] = __n.S_.ode_state[State_::g_AMPA];
  S_.ode_state[State_::g_NMDA__d] = __n.S_.ode_state[State_::g_NMDA__d];
  S_.ode_state[State_::g_NMDA] = __n.S_.ode_state[State_::g_NMDA];
  S_.ode_state[State_::g_GABAA__d] = __n.S_.ode_state[State_::g_GABAA__d];
  S_.ode_state[State_::g_GABAA] = __n.S_.ode_state[State_::g_GABAA];
  S_.ode_state[State_::g_GABAB__d] = __n.S_.ode_state[State_::g_GABAB__d];
  S_.ode_state[State_::g_GABAB] = __n.S_.ode_state[State_::g_GABAB];
  
  V_.AMPAInitialValue = __n.V_.AMPAInitialValue;
  V_.NMDAInitialValue = __n.V_.NMDAInitialValue;
  V_.GABA_AInitialValue = __n.V_.GABA_AInitialValue;
  V_.GABA_BInitialValue = __n.V_.GABA_BInitialValue;
  V_.RefractoryCounts = __n.V_.RefractoryCounts;
  
}

traub_receptors_nestml::~traub_receptors_nestml(){ 
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

void traub_receptors_nestml::init_state_(const Node& proto){
  const traub_receptors_nestml& pr = downcast<traub_receptors_nestml>(proto);
  S_ = pr.S_;
}



extern "C" inline int traub_receptors_nestml_dynamics(double, const double ode_state[], double f[], void* pnode){
  typedef traub_receptors_nestml::State_ State_;
  // get access to node so we can almost work as in a member function
  assert( pnode );
  const traub_receptors_nestml& node = *( reinterpret_cast< traub_receptors_nestml* >( pnode ) );

  // ode_state[] here is---and must be---the state vector supplied by the integrator,
  // not the state vector in the node, node.S_.ode_state[].
  
  f[State_::V_m] = (-(((node.get_g_Na() * ode_state[State_::Act_m] * ode_state[State_::Act_m] * ode_state[State_::Act_m] * ode_state[State_::Inact_h] * (ode_state[State_::V_m] - node.get_E_Na())) + (node.get_g_K() * ode_state[State_::Act_n] * ode_state[State_::Act_n] * ode_state[State_::Act_n] * ode_state[State_::Act_n] * (ode_state[State_::V_m] - node.get_E_K())) + (node.get_g_L() * (ode_state[State_::V_m] - node.get_E_L())))) + node.get_I_e() + node.B_.I_stim_grid_sum_ + ((-(ode_state[State_::g_AMPA]) * (ode_state[State_::V_m] - node.get_AMPA_E_rev())) + (-(ode_state[State_::g_NMDA]) * (ode_state[State_::V_m] - node.get_NMDA_E_rev()) / (1 + std::exp((node.get_NMDA_Vact() - ode_state[State_::V_m]) / node.get_NMDA_Sact()))) + (-(ode_state[State_::g_GABAA]) * (ode_state[State_::V_m] - node.get_GABA_A_E_rev())) + (-(ode_state[State_::g_GABAB]) * (ode_state[State_::V_m] - node.get_GABA_B_E_rev())))) / node.get_C_m();
  f[State_::Act_n] = ((0.032 * (ode_state[State_::V_m] / 1.0 + 52.0) / (1.0 - std::exp(-((ode_state[State_::V_m] / 1.0 + 52.0)) / 5.0))) * (1 - ode_state[State_::Act_n]) - (0.5 * std::exp(-((ode_state[State_::V_m] / 1.0 + 57.0)) / 40.0)) * ode_state[State_::Act_n]) / 1.0;
  f[State_::Act_m] = ((0.32 * (ode_state[State_::V_m] / 1.0 + 54.0) / (1.0 - std::exp(-((ode_state[State_::V_m] / 1.0 + 54.0)) / 4.0))) * (1 - ode_state[State_::Act_m]) - (0.28 * (ode_state[State_::V_m] / 1.0 + 27.0) / (std::exp((ode_state[State_::V_m] / 1.0 + 27.0) / 5.0) - 1.0)) * ode_state[State_::Act_m]) / 1.0;
  f[State_::Inact_h] = ((0.128 * std::exp(-((ode_state[State_::V_m] / 1.0 + 50.0)) / 18.0)) * (1 - ode_state[State_::Inact_h]) - (4.0 / (1.0 + std::exp(-((ode_state[State_::V_m] / 1.0 + 27.0)) / 5.0))) * ode_state[State_::Inact_h]) / 1.0;
  f[State_::g_AMPA__d] = -(ode_state[State_::g_AMPA__d]) / node.get_AMPA_Tau_1();
  f[State_::g_AMPA] = ode_state[State_::g_AMPA__d] - ode_state[State_::g_AMPA] / node.get_AMPA_Tau_2();
  f[State_::g_NMDA__d] = -(ode_state[State_::g_NMDA__d]) / node.get_NMDA_Tau_1();
  f[State_::g_NMDA] = ode_state[State_::g_NMDA__d] - ode_state[State_::g_NMDA] / node.get_NMDA_Tau_2();
  f[State_::g_GABAA__d] = -(ode_state[State_::g_GABAA__d]) / node.get_GABA_A_Tau_1();
  f[State_::g_GABAA] = ode_state[State_::g_GABAA__d] - ode_state[State_::g_GABAA] / node.get_GABA_A_Tau_2();
  f[State_::g_GABAB__d] = -(ode_state[State_::g_GABAB__d]) / node.get_GABA_B_Tau_1();
  f[State_::g_GABAB] = ode_state[State_::g_GABAB__d] - ode_state[State_::g_GABAB] / node.get_GABA_B_Tau_2();
  return GSL_SUCCESS;
}



void traub_receptors_nestml::init_buffers_(){
  get_AMPA().clear(); //includes resize
  get_NMDA().clear(); //includes resize
  get_GABA_A().clear(); //includes resize
  get_GABA_B().clear(); //includes resize
  get_I_stim().clear(); //includes resize
  
  B_.logger_.reset(); // includes resize
  Archiving_Node::clear_history();
  
  if ( B_.__s == 0 ){
    B_.__s = gsl_odeiv_step_alloc( gsl_odeiv_step_rkf45, 12 );
  } else {
    gsl_odeiv_step_reset( B_.__s );
  }

  if ( B_.__c == 0 ){
    B_.__c = gsl_odeiv_control_y_new( P_.__gsl_error_tol, 0.0 );
  } else {
    gsl_odeiv_control_init( B_.__c, P_.__gsl_error_tol, 0.0, 1.0, 0.0 );
  }

  if ( B_.__e == 0 ){
    B_.__e = gsl_odeiv_evolve_alloc( 12 );
  } else {
    gsl_odeiv_evolve_reset( B_.__e );
  }

  B_.__sys.function = traub_receptors_nestml_dynamics;
  B_.__sys.jacobian = NULL;
  B_.__sys.dimension = 12;
  B_.__sys.params = reinterpret_cast< void* >( this );
  B_.__step = nest::Time::get_resolution().get_ms();
  B_.__integration_step = nest::Time::get_resolution().get_ms();
}

void traub_receptors_nestml::calibrate(){
  B_.logger_.init();
  
  
  V_.AMPAInitialValue =compute_synapse_constant(P_.AMPA_Tau_1, P_.AMPA_Tau_2, P_.AMPA_g_peak);
  
  
  V_.NMDAInitialValue =compute_synapse_constant(P_.NMDA_Tau_1, P_.NMDA_Tau_2, P_.NMDA_g_peak);
  
  
  V_.GABA_AInitialValue =compute_synapse_constant(P_.GABA_A_Tau_1, P_.GABA_A_Tau_2, P_.GABA_A_g_peak);
  
  
  V_.GABA_BInitialValue =compute_synapse_constant(P_.GABA_B_Tau_1, P_.GABA_B_Tau_2, P_.GABA_B_g_peak);
  
  
  V_.RefractoryCounts =nest::Time(nest::Time::ms((double) (P_.t_ref))).get_steps();
}

/* ----------------------------------------------------------------
* Update and spike handling functions
* ---------------------------------------------------------------- */

/*
 *
 */
void traub_receptors_nestml::update(nest::Time const & origin,const long from, const long to){
  double __t = 0;

  for ( long lag = from ; lag < to ; ++lag ) {
    B_.AMPA_grid_sum_ = get_AMPA().get_value(lag);
    B_.NMDA_grid_sum_ = get_NMDA().get_value(lag);
    B_.GABA_A_grid_sum_ = get_GABA_A().get_value(lag);
    B_.GABA_B_grid_sum_ = get_GABA_B().get_value(lag);
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
    

    S_.ode_state[State_::g_AMPA__d] += V_.AMPAInitialValue * B_.AMPA_grid_sum_ / 1.0;
    

    S_.ode_state[State_::g_NMDA__d] += V_.NMDAInitialValue * B_.NMDA_grid_sum_ / 1.0;
    

    S_.ode_state[State_::g_GABAA__d] += V_.GABA_AInitialValue * B_.GABA_A_grid_sum_ / 1.0;
    

    S_.ode_state[State_::g_GABAB__d] += V_.GABA_BInitialValue * B_.GABA_B_grid_sum_ / 1.0;
    



    if (S_.r>0) {
      

      S_.r -= 1;
    }else if(S_.ode_state[State_::V_m]>P_.V_Tr&&U_old>P_.V_Tr) {
      

      S_.r = V_.RefractoryCounts;
      
      set_spiketime(nest::Time::step(origin.get_steps()+lag+1));
      nest::SpikeEvent se;
      nest::kernel().event_delivery_manager.send(*this, se, lag);
    } /* if end */


    // voltage logging
    B_.logger_.record_data(origin.get_steps()+lag);
  }

}

// Do not move this function as inline to h-file. It depends on
// universal_data_logger_impl.h being included here.
void traub_receptors_nestml::handle(nest::DataLoggingRequest& e){
  B_.logger_.handle(e);
}

//
double traub_receptors_nestml::compute_synapse_constant(double Tau_1, double Tau_2, double g_peak) const

{
  


  double exact_integration_adjustment = ((1 / Tau_2) - (1 / Tau_1)) * 1.0;
  


  double t_peak = (Tau_2 * Tau_1) * std::log(Tau_2 / Tau_1) / (Tau_2 - Tau_1) / 1.0;
  


  double normalisation_factor = 1 / (std::exp((-t_peak) / Tau_1) - std::exp((-t_peak) / Tau_2));
  


  return g_peak * normalisation_factor * exact_integration_adjustment;

}

void traub_receptors_nestml::handle(nest::SpikeEvent &e){
  assert(e.get_delay_steps() > 0);
  assert( e.get_rport() < static_cast< int >( B_.spike_inputs_.size() ) );

  B_.spike_inputs_[ e.get_rport() ].add_value(
    e.get_rel_delivery_steps( nest::kernel().simulation_manager.get_slice_origin() ),
    e.get_weight() * e.get_multiplicity() );
  
}

void traub_receptors_nestml::handle(nest::CurrentEvent& e){
  assert(e.get_delay_steps() > 0);

  const double current = e.get_current();		// we assume that in NEST, this returns a current in pA
  const double weight = e.get_weight();
  get_I_stim().add_value(
               e.get_rel_delivery_steps( nest::kernel().simulation_manager.get_slice_origin()),
               weight * current );
  
}
