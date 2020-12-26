
/*
*  olm_receptors_nestml.h
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
*  2020-01-25 18:13:23.657668
*/
#ifndef OLM_RECEPTORS_NESTML
#define OLM_RECEPTORS_NESTML

#include "config.h"


#ifdef HAVE_GSL

// External includes:
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

// forwards the declaration of the function
/**
 * Function computing right-hand side of ODE for GSL solver.
 * @note Must be declared here so we can befriend it in class.
 * @note Must have C-linkage for passing to GSL. Internally, it is
 *       a first-class C++ function, but cannot be a member function
 *       because of the C-linkage.
 * @note No point in declaring it inline, since it is called
 *       through a function pointer.
 * @param void* Pointer to model neuron instance.
 */
extern "C" inline int olm_receptors_nestml_dynamics( double, const double y[], double f[], void* pnode );


// Includes from nestkernel:
#include "archiving_node.h"
#include "connection.h"
#include "event.h"
#include "nest_types.h"
#include "ring_buffer.h"
#include "universal_data_logger.h"


// Includes from sli:
#include "dictdatum.h"

/* BeginDocumentation
  Name: olm_receptors_nestml.

  Description:  
    

  Parameters:
  The following parameters can be set in the status dictionary.
  t_ref [ms]  Refractory period 2.0
  g_Na [nS]  Sodium peak conductance
  g_K [nS]  Potassium peak conductance
  g_L [nS]  Leak conductance
  C_m [pF]  Membrane Capacitance
  E_Na [mV]  Sodium reversal potential
  E_K [mV]  Potassium reversal potentia
  E_L [mV]  Leak reversal Potential (aka resting potential)
  E_H [mV]  
  V_Tr [mV]  Spike Threshold
  AMPA_g_peak [nS]  Parameters for synapse of type AMPA, GABA_A, GABA_B and NMDA
 peak conductance
  AMPA_E_rev [mV]  reversal potential
  AMPA_Tau_1 [ms]  rise time
  AMPA_Tau_2 [ms]  decay time, Tau_1 < Tau_2
  NMDA_g_peak [nS]  peak conductance
  NMDA_Tau_1 [ms]  rise time
  NMDA_Tau_2 [ms]  decay time, Tau_1 < Tau_2
  NMDA_E_rev [mV]  reversal potential
  NMDA_Vact [mV]  inactive for V << Vact, inflection of sigmoid
  NMDA_Sact [mV]  scale of inactivation
  GABA_A_g_peak [nS]  peak conductance
  GABA_A_Tau_1 [ms]  rise time
  GABA_A_Tau_2 [ms]  decay time, Tau_1 < Tau_2
  GABA_A_E_rev [mV]  reversal potential
  GABA_B_g_peak [nS]  peak conductance
  GABA_B_Tau_1 [ms]  rise time
  GABA_B_Tau_2 [ms]  decay time, Tau_1 < Tau_2
  GABA_B_E_rev [mV]  reversal potential for intrinsic current
  I_e [pA]  constant external input current
  

  Dynamic state variables:
  r [integer]  number of steps in the current refractory phase
  

  Initial values:
  Inact_h [real]  Inactivation variable h for Na
  Act_n [real]  Activation variable n for K 
  V_m [mV]  Membrane potential
  

  References: Empty

  Sends: nest::SpikeEvent

  Receives: Spike, Current, DataLoggingRequest
*/
class olm_receptors_nestml : public nest::Archiving_Node{
public:
  /**
  * The constructor is only used to create the model prototype in the model manager.
  */
  olm_receptors_nestml();

  /**
  * The copy constructor is used to create model copies and instances of the model.
  * @node The copy constructor needs to initialize the parameters and the state.
  *       Initialization of buffers and interal variables is deferred to
  *       @c init_buffers_() and @c calibrate().
  */
  olm_receptors_nestml(const olm_receptors_nestml &);

  /**
  * Releases resources.
  */
  ~olm_receptors_nestml();

  /**
   * Import sets of overloaded virtual functions.
   * @see Technical Issues / Virtual Functions: Overriding, Overloading, and
   * Hiding
   */
  using nest::Node::handles_test_event;
  using nest::Node::handle;

  /**
  * Used to validate that we can send nest::SpikeEvent to desired target:port.
  */
  nest::port send_test_event(nest::Node& target, nest::rport receptor_type, nest::synindex, bool);

  /**
  * @defgroup mynest_handle Functions handling incoming events.
  * We tell nest that we can handle incoming events of various types by
  * defining @c handle() and @c connect_sender() for the given event.
  * @{
  */
  void handle(nest::SpikeEvent &);        //! accept spikes
  void handle(nest::CurrentEvent &);      //! accept input current
  void handle(nest::DataLoggingRequest &);//! allow recording with multimeter

  nest::port handles_test_event(nest::SpikeEvent&, nest::port);
  nest::port handles_test_event(nest::CurrentEvent&, nest::port);
  nest::port handles_test_event(nest::DataLoggingRequest&, nest::port);
  /** @} */

  // SLI communication functions:
  void get_status(DictionaryDatum &) const;
  void set_status(const DictionaryDatum &);

private:
  /**
     * Synapse types to connect to
     * @note Excluded upper and lower bounds are defined as INF_, SUP_.
     *       Excluding port 0 avoids accidental connections.
     */
    enum SynapseTypes
    {
      INF_SPIKE_RECEPTOR = 0,
      AMPA ,
      NMDA ,
      GABA_A ,
      GABA_B ,
      SUP_SPIKE_RECEPTOR
    };
  //! Reset parameters and state of neuron.

  //! Reset state of neuron.
  void init_state_(const Node& proto);

  //! Reset internal buffers of neuron.
  void init_buffers_();

  //! Initialize auxiliary quantities, leave parameters and state untouched.
  void calibrate();

  //! Take neuron through given time interval
  void update(nest::Time const &, const long, const long);

  // The next two classes need to be friends to access the State_ class/member
  friend class nest::RecordablesMap<olm_receptors_nestml>;
  friend class nest::UniversalDataLogger<olm_receptors_nestml>;

  /**
  * Free parameters of the neuron.
  *
  *
  *
  * These are the parameters that can be set by the user through @c SetStatus.
  * They are initialized from the model prototype when the node is created.
  * Parameters do not change during calls to @c update() and are not reset by
  * @c ResetNetwork.
  *
  * @note Parameters_ need neither copy constructor nor @c operator=(), since
  *       all its members are copied properly by the default copy constructor
  *       and assignment operator. Important:
  *       - If Parameters_ contained @c Time members, you need to define the
  *         assignment operator to recalibrate all members of type @c Time . You
  *         may also want to define the assignment operator.
  *       - If Parameters_ contained members that cannot copy themselves, such
  *         as C-style arrays, you need to define the copy constructor and
  *         assignment operator to copy those members.
  */
  struct Parameters_{
        
        

    //!  Refractory period 2.0
    double t_ref;

    //!  Sodium peak conductance
    double g_Na;

    //!  Potassium peak conductance
    double g_K;

    //!  Leak conductance
    double g_L;

    double g_H;

    double g_A;

    //!  Membrane Capacitance
    double C_m;

    //!  Sodium reversal potential
    double E_Na;

    //!  Potassium reversal potentia
    double E_K;

    //!  Leak reversal Potential (aka resting potential)
    double E_L;

    //!  
    double E_H;

    double E_A;

    //!  Spike Threshold
    double V_Tr;

    //!  Parameters for synapse of type AMPA, GABA_A, GABA_B and NMDA
    //!  peak conductance
    double AMPA_g_peak;

    //!  reversal potential
    double AMPA_E_rev;

    //!  rise time
    double AMPA_Tau_1;

    //!  decay time, Tau_1 < Tau_2
    double AMPA_Tau_2;

    //!  peak conductance
    double NMDA_g_peak;

    //!  rise time
    double NMDA_Tau_1;

    //!  decay time, Tau_1 < Tau_2
    double NMDA_Tau_2;

    //!  reversal potential
    double NMDA_E_rev;

    //!  inactive for V << Vact, inflection of sigmoid
    double NMDA_Vact;

    //!  scale of inactivation
    double NMDA_Sact;

    //!  peak conductance
    double GABA_A_g_peak;

    //!  rise time
    double GABA_A_Tau_1;

    //!  decay time, Tau_1 < Tau_2
    double GABA_A_Tau_2;

    //!  reversal potential
    double GABA_A_E_rev;

    //!  peak conductance
    double GABA_B_g_peak;

    //!  rise time
    double GABA_B_Tau_1;

    //!  decay time, Tau_1 < Tau_2
    double GABA_B_Tau_2;

    //!  reversal potential for intrinsic current
    double GABA_B_E_rev;

    //!  constant external input current
    double I_e;

    double __gsl_error_tol;
    /** Initialize parameters to their default values. */
    Parameters_();
  };

  /**
  * Dynamic state of the neuron.
  *
  *
  *
  * These are the state variables that are advanced in time by calls to
  * @c update(). In many models, some or all of them can be set by the user
  * through @c SetStatus. The state variables are initialized from the model
  * prototype when the node is created. State variables are reset by @c ResetNetwork.
  *
  * @note State_ need neither copy constructor nor @c operator=(), since
  *       all its members are copied properly by the default copy constructor
  *       and assignment operator. Important:
  *       - If State_ contained @c Time members, you need to define the
  *         assignment operator to recalibrate all members of type @c Time . You
  *         may also want to define the assignment operator.
  *       - If State_ contained members that cannot copy themselves, such
  *         as C-style arrays, you need to define the copy constructor and
  *         assignment operator to copy those members.
  */
  struct State_{
    //! Symbolic indices to the elements of the state vector y
    enum StateVecElems{
    //  Inactivation variable h for Na
      Inact_h,
      //  Activation variable n for K 
      Act_n,
      
      Act_r_h,
      
      Act_A,
      
      Act_B,
      //  Membrane potential
      V_m,
      
      g_AMPA__d,
      
      g_AMPA,
      
      g_NMDA__d,
      
      g_NMDA,
      
      g_GABAA__d,
      
      g_GABAA,
      
      g_GABAB__d,
      
      g_GABAB,
      STATE_VEC_SIZE
    };
    //! state vector, must be C-array for GSL solver
    double ode_state[STATE_VEC_SIZE];    

    //!  number of steps in the current refractory phase
    long r;

    double alpha_n_init;

    double beta_n_init;

    double alpha_m_init;

    double beta_m_init;

    double alpha_h_init;

    double beta_h_init;

    double r_h_inf_init;

    double Act_A_inf_init;

    double Act_B_inf_init;    

    double I_syn_ampa;

    double I_syn_nmda;

    double I_syn_gaba_a;

    double I_syn_gaba_b;

    double I_syn;

    double alpha_n;

    double beta_n;

    double alpha_m;

    double beta_m;

    double alpha_h;

    double beta_h;

    double Act_m_inf;

    double I_Na;

    double I_K;

    double I_L;

    double I_H;

    double I_A;

    double r_h_inf;

    double tau_r_h;

    double Act_A_inf;

    double Act_B_inf;

    double tau_A_A;

    double tau_B_A;
        State_();
  };

  /**
  * Internal variables of the neuron.
  *
  *
  *
  * These variables must be initialized by @c calibrate, which is called before
  * the first call to @c update() upon each call to @c Simulate.
  * @node Variables_ needs neither constructor, copy constructor or assignment operator,
  *       since it is initialized by @c calibrate(). If Variables_ has members that
  *       cannot destroy themselves, Variables_ will need a destructor.
  */
  struct Variables_ {    

    double AMPAInitialValue;
        

    double NMDAInitialValue;
        

    double GABA_AInitialValue;
        

    double GABA_BInitialValue;
        

    //!  refractory time in steps
    long RefractoryCounts;
    
  };

  /**
    * Buffers of the neuron.
    * Ususally buffers for incoming spikes and data logged for analog recorders.
    * Buffers must be initialized by @c init_buffers_(), which is called before
    * @c calibrate() on the first call to @c Simulate after the start of NEST,
    * ResetKernel or ResetNetwork.
    * @node Buffers_ needs neither constructor, copy constructor or assignment operator,
    *       since it is initialized by @c init_nodes_(). If Buffers_ has members that
    *       cannot destroy themselves, Buffers_ will need a destructor.
    */
  struct Buffers_ {
    Buffers_(olm_receptors_nestml &);
    Buffers_(const Buffers_ &, olm_receptors_nestml &);

    /** Logger for all analog data */
    nest::UniversalDataLogger<olm_receptors_nestml> logger_;
    
    std::vector<long> receptor_types_;
    /** buffers and sums up incoming spikes/currents */
    std::vector< nest::RingBuffer > spike_inputs_;

    
    inline nest::RingBuffer& get_AMPA() {  return spike_inputs_[AMPA - 1]; }
    double AMPA_grid_sum_;
    
    inline nest::RingBuffer& get_NMDA() {  return spike_inputs_[NMDA - 1]; }
    double NMDA_grid_sum_;
    
    inline nest::RingBuffer& get_GABA_A() {  return spike_inputs_[GABA_A - 1]; }
    double GABA_A_grid_sum_;
    
    inline nest::RingBuffer& get_GABA_B() {  return spike_inputs_[GABA_B - 1]; }
    double GABA_B_grid_sum_;
    
    //!< Buffer incoming pAs through delay, as sum
    nest::RingBuffer I_stim;
    inline nest::RingBuffer& get_I_stim() {return I_stim;}
    double I_stim_grid_sum_;
    /** GSL ODE stuff */
    gsl_odeiv_step* __s;    //!< stepping function
    gsl_odeiv_control* __c; //!< adaptive stepsize control function
    gsl_odeiv_evolve* __e;  //!< evolution function
    gsl_odeiv_system __sys; //!< struct describing system

    // IntergrationStep_ should be reset with the neuron on ResetNetwork,
    // but remain unchanged during calibration. Since it is initialized with
    // step_, and the resolution cannot change after nodes have been created,
    // it is safe to place both here.
    double __step;             //!< step size in ms
    double __integration_step; //!< current integration time step, updated by GSL
    };
  inline long get_r() const {
    return S_.r;
  }
  inline void set_r(const long __v) {
    S_.r = __v;
  }

  inline double get_alpha_n_init() const {
    return 0.018 * (S_.ode_state[State_::V_m] / 1.0 - 25.0) / (1.0 - std::exp((-(S_.ode_state[State_::V_m] / 1.0 - 25.0)) / 25.0));
  }

  inline double get_beta_n_init() const {
    return 0.0036 * (35.0 - S_.ode_state[State_::V_m] / 1.0) / (1.0 - std::exp((-(35.0 - S_.ode_state[State_::V_m] / 1.0)) / 12.0));
  }

  inline double get_alpha_m_init() const {
    return 0.1 * (S_.ode_state[State_::V_m] / 1.0 + 38.0) / (1.0 - std::exp((-0.1) * (S_.ode_state[State_::V_m] / 1.0 + 38.0)));
  }

  inline double get_beta_m_init() const {
    return 4.0 * std::exp((-(S_.ode_state[State_::V_m] / 1.0 + 65.0)) / 18.0);
  }

  inline double get_alpha_h_init() const {
    return 0.07 * std::exp((-(S_.ode_state[State_::V_m] / 1.0 + 63.0)) / 20.0);
  }

  inline double get_beta_h_init() const {
    return 1.0 / (std::exp((-(S_.ode_state[State_::V_m] / 1.0 + 33.0)) / 10.0) + 1.0);
  }

  inline double get_r_h_inf_init() const {
    return 1.0 / (1.0 + std::exp((S_.ode_state[State_::V_m] / 1.0 + 84.0) / 10.2));
  }

  inline double get_Act_A_inf_init() const {
    return 1.0 / (1.0 + std::exp((-(S_.ode_state[State_::V_m] / 1.0 + 14.0)) / 16.6));
  }

  inline double get_Act_B_inf_init() const {
    return 1.0 / (1.0 + std::exp((S_.ode_state[State_::V_m] / 1.0 + 71.0) / 7.3));
  }

  inline double get_Inact_h() const {
    return S_.ode_state[State_::Inact_h];
  }
  inline void set_Inact_h(const double __v) {
    S_.ode_state[State_::Inact_h] = __v;
  }

  inline double get_Act_n() const {
    return S_.ode_state[State_::Act_n];
  }
  inline void set_Act_n(const double __v) {
    S_.ode_state[State_::Act_n] = __v;
  }

  inline double get_Act_r_h() const {
    return S_.ode_state[State_::Act_r_h];
  }
  inline void set_Act_r_h(const double __v) {
    S_.ode_state[State_::Act_r_h] = __v;
  }

  inline double get_Act_A() const {
    return S_.ode_state[State_::Act_A];
  }
  inline void set_Act_A(const double __v) {
    S_.ode_state[State_::Act_A] = __v;
  }

  inline double get_Act_B() const {
    return S_.ode_state[State_::Act_B];
  }
  inline void set_Act_B(const double __v) {
    S_.ode_state[State_::Act_B] = __v;
  }

  inline double get_V_m() const {
    return S_.ode_state[State_::V_m];
  }
  inline void set_V_m(const double __v) {
    S_.ode_state[State_::V_m] = __v;
  }

  inline double get_g_AMPA__d() const {
    return S_.ode_state[State_::g_AMPA__d];
  }
  inline void set_g_AMPA__d(const double __v) {
    S_.ode_state[State_::g_AMPA__d] = __v;
  }

  inline double get_g_AMPA() const {
    return S_.ode_state[State_::g_AMPA];
  }
  inline void set_g_AMPA(const double __v) {
    S_.ode_state[State_::g_AMPA] = __v;
  }

  inline double get_g_NMDA__d() const {
    return S_.ode_state[State_::g_NMDA__d];
  }
  inline void set_g_NMDA__d(const double __v) {
    S_.ode_state[State_::g_NMDA__d] = __v;
  }

  inline double get_g_NMDA() const {
    return S_.ode_state[State_::g_NMDA];
  }
  inline void set_g_NMDA(const double __v) {
    S_.ode_state[State_::g_NMDA] = __v;
  }

  inline double get_g_GABAA__d() const {
    return S_.ode_state[State_::g_GABAA__d];
  }
  inline void set_g_GABAA__d(const double __v) {
    S_.ode_state[State_::g_GABAA__d] = __v;
  }

  inline double get_g_GABAA() const {
    return S_.ode_state[State_::g_GABAA];
  }
  inline void set_g_GABAA(const double __v) {
    S_.ode_state[State_::g_GABAA] = __v;
  }

  inline double get_g_GABAB__d() const {
    return S_.ode_state[State_::g_GABAB__d];
  }
  inline void set_g_GABAB__d(const double __v) {
    S_.ode_state[State_::g_GABAB__d] = __v;
  }

  inline double get_g_GABAB() const {
    return S_.ode_state[State_::g_GABAB];
  }
  inline void set_g_GABAB(const double __v) {
    S_.ode_state[State_::g_GABAB] = __v;
  }

  inline double get_t_ref() const {
    return P_.t_ref;
  }
  inline void set_t_ref(const double __v) {
    P_.t_ref = __v;
  }

  inline double get_g_Na() const {
    return P_.g_Na;
  }
  inline void set_g_Na(const double __v) {
    P_.g_Na = __v;
  }

  inline double get_g_K() const {
    return P_.g_K;
  }
  inline void set_g_K(const double __v) {
    P_.g_K = __v;
  }

  inline double get_g_L() const {
    return P_.g_L;
  }
  inline void set_g_L(const double __v) {
    P_.g_L = __v;
  }

  inline double get_g_H() const {
    return P_.g_H;
  }
  inline void set_g_H(const double __v) {
    P_.g_H = __v;
  }

  inline double get_g_A() const {
    return P_.g_A;
  }
  inline void set_g_A(const double __v) {
    P_.g_A = __v;
  }

  inline double get_C_m() const {
    return P_.C_m;
  }
  inline void set_C_m(const double __v) {
    P_.C_m = __v;
  }

  inline double get_E_Na() const {
    return P_.E_Na;
  }
  inline void set_E_Na(const double __v) {
    P_.E_Na = __v;
  }

  inline double get_E_K() const {
    return P_.E_K;
  }
  inline void set_E_K(const double __v) {
    P_.E_K = __v;
  }

  inline double get_E_L() const {
    return P_.E_L;
  }
  inline void set_E_L(const double __v) {
    P_.E_L = __v;
  }

  inline double get_E_H() const {
    return P_.E_H;
  }
  inline void set_E_H(const double __v) {
    P_.E_H = __v;
  }

  inline double get_E_A() const {
    return P_.E_A;
  }
  inline void set_E_A(const double __v) {
    P_.E_A = __v;
  }

  inline double get_V_Tr() const {
    return P_.V_Tr;
  }
  inline void set_V_Tr(const double __v) {
    P_.V_Tr = __v;
  }

  inline double get_AMPA_g_peak() const {
    return P_.AMPA_g_peak;
  }
  inline void set_AMPA_g_peak(const double __v) {
    P_.AMPA_g_peak = __v;
  }

  inline double get_AMPA_E_rev() const {
    return P_.AMPA_E_rev;
  }
  inline void set_AMPA_E_rev(const double __v) {
    P_.AMPA_E_rev = __v;
  }

  inline double get_AMPA_Tau_1() const {
    return P_.AMPA_Tau_1;
  }
  inline void set_AMPA_Tau_1(const double __v) {
    P_.AMPA_Tau_1 = __v;
  }

  inline double get_AMPA_Tau_2() const {
    return P_.AMPA_Tau_2;
  }
  inline void set_AMPA_Tau_2(const double __v) {
    P_.AMPA_Tau_2 = __v;
  }

  inline double get_NMDA_g_peak() const {
    return P_.NMDA_g_peak;
  }
  inline void set_NMDA_g_peak(const double __v) {
    P_.NMDA_g_peak = __v;
  }

  inline double get_NMDA_Tau_1() const {
    return P_.NMDA_Tau_1;
  }
  inline void set_NMDA_Tau_1(const double __v) {
    P_.NMDA_Tau_1 = __v;
  }

  inline double get_NMDA_Tau_2() const {
    return P_.NMDA_Tau_2;
  }
  inline void set_NMDA_Tau_2(const double __v) {
    P_.NMDA_Tau_2 = __v;
  }

  inline double get_NMDA_E_rev() const {
    return P_.NMDA_E_rev;
  }
  inline void set_NMDA_E_rev(const double __v) {
    P_.NMDA_E_rev = __v;
  }

  inline double get_NMDA_Vact() const {
    return P_.NMDA_Vact;
  }
  inline void set_NMDA_Vact(const double __v) {
    P_.NMDA_Vact = __v;
  }

  inline double get_NMDA_Sact() const {
    return P_.NMDA_Sact;
  }
  inline void set_NMDA_Sact(const double __v) {
    P_.NMDA_Sact = __v;
  }

  inline double get_GABA_A_g_peak() const {
    return P_.GABA_A_g_peak;
  }
  inline void set_GABA_A_g_peak(const double __v) {
    P_.GABA_A_g_peak = __v;
  }

  inline double get_GABA_A_Tau_1() const {
    return P_.GABA_A_Tau_1;
  }
  inline void set_GABA_A_Tau_1(const double __v) {
    P_.GABA_A_Tau_1 = __v;
  }

  inline double get_GABA_A_Tau_2() const {
    return P_.GABA_A_Tau_2;
  }
  inline void set_GABA_A_Tau_2(const double __v) {
    P_.GABA_A_Tau_2 = __v;
  }

  inline double get_GABA_A_E_rev() const {
    return P_.GABA_A_E_rev;
  }
  inline void set_GABA_A_E_rev(const double __v) {
    P_.GABA_A_E_rev = __v;
  }

  inline double get_GABA_B_g_peak() const {
    return P_.GABA_B_g_peak;
  }
  inline void set_GABA_B_g_peak(const double __v) {
    P_.GABA_B_g_peak = __v;
  }

  inline double get_GABA_B_Tau_1() const {
    return P_.GABA_B_Tau_1;
  }
  inline void set_GABA_B_Tau_1(const double __v) {
    P_.GABA_B_Tau_1 = __v;
  }

  inline double get_GABA_B_Tau_2() const {
    return P_.GABA_B_Tau_2;
  }
  inline void set_GABA_B_Tau_2(const double __v) {
    P_.GABA_B_Tau_2 = __v;
  }

  inline double get_GABA_B_E_rev() const {
    return P_.GABA_B_E_rev;
  }
  inline void set_GABA_B_E_rev(const double __v) {
    P_.GABA_B_E_rev = __v;
  }

  inline double get_I_e() const {
    return P_.I_e;
  }
  inline void set_I_e(const double __v) {
    P_.I_e = __v;
  }

  inline double get_AMPAInitialValue() const {
    return V_.AMPAInitialValue;
  }
  inline void set_AMPAInitialValue(const double __v) {
    V_.AMPAInitialValue = __v;
  }

  inline double get_NMDAInitialValue() const {
    return V_.NMDAInitialValue;
  }
  inline void set_NMDAInitialValue(const double __v) {
    V_.NMDAInitialValue = __v;
  }

  inline double get_GABA_AInitialValue() const {
    return V_.GABA_AInitialValue;
  }
  inline void set_GABA_AInitialValue(const double __v) {
    V_.GABA_AInitialValue = __v;
  }

  inline double get_GABA_BInitialValue() const {
    return V_.GABA_BInitialValue;
  }
  inline void set_GABA_BInitialValue(const double __v) {
    V_.GABA_BInitialValue = __v;
  }

  inline long get_RefractoryCounts() const {
    return V_.RefractoryCounts;
  }
  inline void set_RefractoryCounts(const long __v) {
    V_.RefractoryCounts = __v;
  }

  inline double get_I_syn_ampa() const {
    return (-S_.ode_state[State_::g_AMPA]) * (S_.ode_state[State_::V_m] - P_.AMPA_E_rev);
  }

  inline double get_I_syn_nmda() const {
    return (-S_.ode_state[State_::g_NMDA]) * (S_.ode_state[State_::V_m] - P_.NMDA_E_rev) / (1 + std::exp((P_.NMDA_Vact - S_.ode_state[State_::V_m]) / P_.NMDA_Sact));
  }

  inline double get_I_syn_gaba_a() const {
    return (-S_.ode_state[State_::g_GABAA]) * (S_.ode_state[State_::V_m] - P_.GABA_A_E_rev);
  }

  inline double get_I_syn_gaba_b() const {
    return (-S_.ode_state[State_::g_GABAB]) * (S_.ode_state[State_::V_m] - P_.GABA_B_E_rev);
  }

  inline double get_I_syn() const {
    return ((-S_.ode_state[State_::g_AMPA]) * (S_.ode_state[State_::V_m] - P_.AMPA_E_rev)) + ((-S_.ode_state[State_::g_NMDA]) * (S_.ode_state[State_::V_m] - P_.NMDA_E_rev) / (1 + std::exp((P_.NMDA_Vact - S_.ode_state[State_::V_m]) / P_.NMDA_Sact))) + ((-S_.ode_state[State_::g_GABAA]) * (S_.ode_state[State_::V_m] - P_.GABA_A_E_rev)) + ((-S_.ode_state[State_::g_GABAB]) * (S_.ode_state[State_::V_m] - P_.GABA_B_E_rev));
  }

  inline double get_alpha_n() const {
    return 0.018 * (S_.ode_state[State_::V_m] / 1.0 - 25.0) / (1.0 - std::exp((-(S_.ode_state[State_::V_m] / 1.0 - 25.0)) / 25.0));
  }

  inline double get_beta_n() const {
    return 0.0036 * (35.0 - S_.ode_state[State_::V_m] / 1.0) / (1.0 - std::exp((-(35.0 - S_.ode_state[State_::V_m] / 1.0)) / 12.0));
  }

  inline double get_alpha_m() const {
    return 0.1 * (S_.ode_state[State_::V_m] / 1.0 + 38.0) / (1.0 - std::exp((-0.1) * (S_.ode_state[State_::V_m] / 1.0 + 38.0)));
  }

  inline double get_beta_m() const {
    return 4.0 * std::exp((-(S_.ode_state[State_::V_m] / 1.0 + 65.0)) / 18.0);
  }

  inline double get_alpha_h() const {
    return 0.07 * std::exp((-(S_.ode_state[State_::V_m] / 1.0 + 63.0)) / 20.0);
  }

  inline double get_beta_h() const {
    return 1.0 / (std::exp((-(S_.ode_state[State_::V_m] / 1.0 + 33.0)) / 10.0) + 1.0);
  }

  inline double get_Act_m_inf() const {
    return (0.1 * (S_.ode_state[State_::V_m] / 1.0 + 38.0) / (1.0 - std::exp((-0.1) * (S_.ode_state[State_::V_m] / 1.0 + 38.0)))) / ((0.1 * (S_.ode_state[State_::V_m] / 1.0 + 38.0) / (1.0 - std::exp((-0.1) * (S_.ode_state[State_::V_m] / 1.0 + 38.0)))) + (4.0 * std::exp((-(S_.ode_state[State_::V_m] / 1.0 + 65.0)) / 18.0)));
  }

  inline double get_I_Na() const {
    return P_.g_Na * ((0.1 * (S_.ode_state[State_::V_m] / 1.0 + 38.0) / (1.0 - std::exp((-0.1) * (S_.ode_state[State_::V_m] / 1.0 + 38.0)))) / ((0.1 * (S_.ode_state[State_::V_m] / 1.0 + 38.0) / (1.0 - std::exp((-0.1) * (S_.ode_state[State_::V_m] / 1.0 + 38.0)))) + (4.0 * std::exp((-(S_.ode_state[State_::V_m] / 1.0 + 65.0)) / 18.0)))) * ((0.1 * (S_.ode_state[State_::V_m] / 1.0 + 38.0) / (1.0 - std::exp((-0.1) * (S_.ode_state[State_::V_m] / 1.0 + 38.0)))) / ((0.1 * (S_.ode_state[State_::V_m] / 1.0 + 38.0) / (1.0 - std::exp((-0.1) * (S_.ode_state[State_::V_m] / 1.0 + 38.0)))) + (4.0 * std::exp((-(S_.ode_state[State_::V_m] / 1.0 + 65.0)) / 18.0)))) * ((0.1 * (S_.ode_state[State_::V_m] / 1.0 + 38.0) / (1.0 - std::exp((-0.1) * (S_.ode_state[State_::V_m] / 1.0 + 38.0)))) / ((0.1 * (S_.ode_state[State_::V_m] / 1.0 + 38.0) / (1.0 - std::exp((-0.1) * (S_.ode_state[State_::V_m] / 1.0 + 38.0)))) + (4.0 * std::exp((-(S_.ode_state[State_::V_m] / 1.0 + 65.0)) / 18.0)))) * S_.ode_state[State_::Inact_h] * (S_.ode_state[State_::V_m] - P_.E_Na);
  }

  inline double get_I_K() const {
    return P_.g_K * S_.ode_state[State_::Act_n] * S_.ode_state[State_::Act_n] * S_.ode_state[State_::Act_n] * S_.ode_state[State_::Act_n] * (S_.ode_state[State_::V_m] - P_.E_K);
  }

  inline double get_I_L() const {
    return P_.g_L * (S_.ode_state[State_::V_m] - P_.E_L);
  }

  inline double get_I_H() const {
    return P_.g_H * S_.ode_state[State_::Act_r_h] * (S_.ode_state[State_::V_m] - P_.E_H);
  }

  inline double get_I_A() const {
    return P_.g_A * S_.ode_state[State_::Act_A] * S_.ode_state[State_::Act_B] * (S_.ode_state[State_::V_m] - P_.E_A);
  }

  inline double get_r_h_inf() const {
    return 1.0 / (1.0 + std::exp((S_.ode_state[State_::V_m] / 1.0 + 84.0) / 10.2));
  }

  inline double get_tau_r_h() const {
    return 1.0 / (std::exp((-14.59) - 0.086 * S_.ode_state[State_::V_m] / 1.0) + std::exp((-1.87) + 0.0701 * S_.ode_state[State_::V_m] / 1.0));
  }

  inline double get_Act_A_inf() const {
    return 1.0 / (1.0 + std::exp((-(S_.ode_state[State_::V_m] / 1.0 + 14.0)) / 16.6));
  }

  inline double get_Act_B_inf() const {
    return 1.0 / (1.0 + std::exp((S_.ode_state[State_::V_m] / 1.0 + 71.0) / 7.3));
  }

  inline double get_tau_A_A() const {
    return 5.0;
  }

  inline double get_tau_B_A() const {
    return 1.0 / (9e-06 / std::exp((S_.ode_state[State_::V_m] / 1.0 - 26.0) / 28.5) + 0.014 / (0.2 + std::exp((-(S_.ode_state[State_::V_m] / 1.0 + 70.0)) / 11.0)));
  }


  
  inline nest::RingBuffer& get_AMPA() {return B_.get_AMPA();};
  
  inline nest::RingBuffer& get_NMDA() {return B_.get_NMDA();};
  
  inline nest::RingBuffer& get_GABA_A() {return B_.get_GABA_A();};
  
  inline nest::RingBuffer& get_GABA_B() {return B_.get_GABA_B();};
  
  inline nest::RingBuffer& get_I_stim() {return B_.get_I_stim();};
  

  // Generate function header
  
  //
double compute_synapse_constant(double, double, double) const
;
  
  /**
  * @defgroup pif_members Member variables of neuron model.
  * Each model neuron should have precisely the following four data members,
  * which are one instance each of the parameters, state, buffers and variables
  * structures. Experience indicates that the state and variables member should
  * be next to each other to achieve good efficiency (caching).
  * @note Devices require one additional data member, an instance of the @c Device
  *       child class they belong to.
  * @{
  */
  Parameters_ P_;  //!< Free parameters.
  State_      S_;  //!< Dynamic state.
  Variables_  V_;  //!< Internal Variables
  Buffers_    B_;  //!< Buffers.

  //! Mapping of recordables names to access functions
  static nest::RecordablesMap<olm_receptors_nestml> recordablesMap_;

  friend int olm_receptors_nestml_dynamics( double, const double y[], double f[], void* pnode );
  
/** @} */
}; /* neuron olm_receptors_nestml */

inline nest::port olm_receptors_nestml::send_test_event(
    nest::Node& target, nest::rport receptor_type, nest::synindex, bool){
  // You should usually not change the code in this function.
  // It confirms that the target of connection @c c accepts @c nest::SpikeEvent on
  // the given @c receptor_type.
  nest::SpikeEvent e;
  e.set_sender(*this);
  return target.handles_test_event(e, receptor_type);
}

inline nest::port olm_receptors_nestml::handles_test_event(nest::SpikeEvent&, nest::port receptor_type){
  assert( B_.spike_inputs_.size() == 4 );

    if ( !( INF_SPIKE_RECEPTOR < receptor_type && receptor_type < SUP_SPIKE_RECEPTOR ) )
    {
      throw nest::UnknownReceptorType( receptor_type, get_name() );
      return 0;
    }
    else {
      return receptor_type - 1;
    }
}



inline nest::port olm_receptors_nestml::handles_test_event(
    nest::CurrentEvent&, nest::port receptor_type){
  // You should usually not change the code in this function.
  // It confirms to the connection management system that we are able
  // to handle @c CurrentEvent on port 0. You need to extend the function
  // if you want to differentiate between input ports.
  if (receptor_type != 0)
  throw nest::UnknownReceptorType(receptor_type, get_name());
  return 0;
}

inline nest::port olm_receptors_nestml::handles_test_event(
    nest::DataLoggingRequest& dlr, nest::port receptor_type){
  // You should usually not change the code in this function.
  // It confirms to the connection management system that we are able
  // to handle @c DataLoggingRequest on port 0.
  // The function also tells the built-in UniversalDataLogger that this node
  // is recorded from and that it thus needs to collect data during simulation.
  if (receptor_type != 0)
  throw nest::UnknownReceptorType(receptor_type, get_name());

  return B_.logger_.connect_logging_device(dlr, recordablesMap_);
}

// TODO call get_status on used or internal components
inline void olm_receptors_nestml::get_status(DictionaryDatum &__d) const{  
  def<double>(__d, "t_ref", get_t_ref());
      
  def<double>(__d, "g_Na", get_g_Na());
      
  def<double>(__d, "g_K", get_g_K());
      
  def<double>(__d, "g_L", get_g_L());
      
  def<double>(__d, "g_H", get_g_H());
      
  def<double>(__d, "g_A", get_g_A());
      
  def<double>(__d, "C_m", get_C_m());
      
  def<double>(__d, "E_Na", get_E_Na());
      
  def<double>(__d, "E_K", get_E_K());
      
  def<double>(__d, "E_L", get_E_L());
      
  def<double>(__d, "E_H", get_E_H());
      
  def<double>(__d, "E_A", get_E_A());
      
  def<double>(__d, "V_Tr", get_V_Tr());
      
  def<double>(__d, "AMPA_g_peak", get_AMPA_g_peak());
      
  def<double>(__d, "AMPA_E_rev", get_AMPA_E_rev());
      
  def<double>(__d, "AMPA_Tau_1", get_AMPA_Tau_1());
      
  def<double>(__d, "AMPA_Tau_2", get_AMPA_Tau_2());
      
  def<double>(__d, "NMDA_g_peak", get_NMDA_g_peak());
      
  def<double>(__d, "NMDA_Tau_1", get_NMDA_Tau_1());
      
  def<double>(__d, "NMDA_Tau_2", get_NMDA_Tau_2());
      
  def<double>(__d, "NMDA_E_rev", get_NMDA_E_rev());
      
  def<double>(__d, "NMDA_Vact", get_NMDA_Vact());
      
  def<double>(__d, "NMDA_Sact", get_NMDA_Sact());
      
  def<double>(__d, "GABA_A_g_peak", get_GABA_A_g_peak());
      
  def<double>(__d, "GABA_A_Tau_1", get_GABA_A_Tau_1());
      
  def<double>(__d, "GABA_A_Tau_2", get_GABA_A_Tau_2());
      
  def<double>(__d, "GABA_A_E_rev", get_GABA_A_E_rev());
      
  def<double>(__d, "GABA_B_g_peak", get_GABA_B_g_peak());
      
  def<double>(__d, "GABA_B_Tau_1", get_GABA_B_Tau_1());
      
  def<double>(__d, "GABA_B_Tau_2", get_GABA_B_Tau_2());
      
  def<double>(__d, "GABA_B_E_rev", get_GABA_B_E_rev());
      
  def<double>(__d, "I_e", get_I_e());
      
  def<long>(__d, "r", get_r());
      
  def<double>(__d, "alpha_n_init", get_alpha_n_init());
      
  def<double>(__d, "beta_n_init", get_beta_n_init());
      
  def<double>(__d, "alpha_m_init", get_alpha_m_init());
      
  def<double>(__d, "beta_m_init", get_beta_m_init());
      
  def<double>(__d, "alpha_h_init", get_alpha_h_init());
      
  def<double>(__d, "beta_h_init", get_beta_h_init());
      
  def<double>(__d, "r_h_inf_init", get_r_h_inf_init());
      
  def<double>(__d, "Act_A_inf_init", get_Act_A_inf_init());
      
  def<double>(__d, "Act_B_inf_init", get_Act_B_inf_init());
      
  def<double>(__d, "Inact_h", get_Inact_h());
      
  def<double>(__d, "Act_n", get_Act_n());
      
  def<double>(__d, "Act_r_h", get_Act_r_h());
      
  def<double>(__d, "Act_A", get_Act_A());
      
  def<double>(__d, "Act_B", get_Act_B());
      
  def<double>(__d, "V_m", get_V_m());
      
  def<double>(__d, "g_AMPA__d", get_g_AMPA__d());
      
  def<double>(__d, "g_AMPA", get_g_AMPA());
      
  def<double>(__d, "g_NMDA__d", get_g_NMDA__d());
      
  def<double>(__d, "g_NMDA", get_g_NMDA());
      
  def<double>(__d, "g_GABAA__d", get_g_GABAA__d());
      
  def<double>(__d, "g_GABAA", get_g_GABAA());
      
  def<double>(__d, "g_GABAB__d", get_g_GABAB__d());
      
  def<double>(__d, "g_GABAB", get_g_GABAB());
    DictionaryDatum __receptor_type = new Dictionary();
  ( *__receptor_type )[ "AMPA" ] = AMPA;
  ( *__receptor_type )[ "NMDA" ] = NMDA;
  ( *__receptor_type )[ "GABA_A" ] = GABA_A;
  ( *__receptor_type )[ "GABA_B" ] = GABA_B;
  
  ( *__d )[ "receptor_types" ] = __receptor_type;
  

  (*__d)[nest::names::recordables] = recordablesMap_.get_list();
  
  def< double >(__d, nest::names::gsl_error_tol, P_.__gsl_error_tol);
  if ( P_.__gsl_error_tol <= 0. ){
    throw nest::BadProperty( "The gsl_error_tol must be strictly positive." );
  }
  

}

inline void olm_receptors_nestml::set_status(const DictionaryDatum &__d){

  double tmp_t_ref = get_t_ref();
  updateValue<double>(__d, "t_ref", tmp_t_ref);


  double tmp_g_Na = get_g_Na();
  updateValue<double>(__d, "g_Na", tmp_g_Na);


  double tmp_g_K = get_g_K();
  updateValue<double>(__d, "g_K", tmp_g_K);


  double tmp_g_L = get_g_L();
  updateValue<double>(__d, "g_L", tmp_g_L);


  double tmp_g_H = get_g_H();
  updateValue<double>(__d, "g_H", tmp_g_H);


  double tmp_g_A = get_g_A();
  updateValue<double>(__d, "g_A", tmp_g_A);


  double tmp_C_m = get_C_m();
  updateValue<double>(__d, "C_m", tmp_C_m);


  double tmp_E_Na = get_E_Na();
  updateValue<double>(__d, "E_Na", tmp_E_Na);


  double tmp_E_K = get_E_K();
  updateValue<double>(__d, "E_K", tmp_E_K);


  double tmp_E_L = get_E_L();
  updateValue<double>(__d, "E_L", tmp_E_L);


  double tmp_E_H = get_E_H();
  updateValue<double>(__d, "E_H", tmp_E_H);


  double tmp_E_A = get_E_A();
  updateValue<double>(__d, "E_A", tmp_E_A);


  double tmp_V_Tr = get_V_Tr();
  updateValue<double>(__d, "V_Tr", tmp_V_Tr);


  double tmp_AMPA_g_peak = get_AMPA_g_peak();
  updateValue<double>(__d, "AMPA_g_peak", tmp_AMPA_g_peak);


  double tmp_AMPA_E_rev = get_AMPA_E_rev();
  updateValue<double>(__d, "AMPA_E_rev", tmp_AMPA_E_rev);


  double tmp_AMPA_Tau_1 = get_AMPA_Tau_1();
  updateValue<double>(__d, "AMPA_Tau_1", tmp_AMPA_Tau_1);


  double tmp_AMPA_Tau_2 = get_AMPA_Tau_2();
  updateValue<double>(__d, "AMPA_Tau_2", tmp_AMPA_Tau_2);


  double tmp_NMDA_g_peak = get_NMDA_g_peak();
  updateValue<double>(__d, "NMDA_g_peak", tmp_NMDA_g_peak);


  double tmp_NMDA_Tau_1 = get_NMDA_Tau_1();
  updateValue<double>(__d, "NMDA_Tau_1", tmp_NMDA_Tau_1);


  double tmp_NMDA_Tau_2 = get_NMDA_Tau_2();
  updateValue<double>(__d, "NMDA_Tau_2", tmp_NMDA_Tau_2);


  double tmp_NMDA_E_rev = get_NMDA_E_rev();
  updateValue<double>(__d, "NMDA_E_rev", tmp_NMDA_E_rev);


  double tmp_NMDA_Vact = get_NMDA_Vact();
  updateValue<double>(__d, "NMDA_Vact", tmp_NMDA_Vact);


  double tmp_NMDA_Sact = get_NMDA_Sact();
  updateValue<double>(__d, "NMDA_Sact", tmp_NMDA_Sact);


  double tmp_GABA_A_g_peak = get_GABA_A_g_peak();
  updateValue<double>(__d, "GABA_A_g_peak", tmp_GABA_A_g_peak);


  double tmp_GABA_A_Tau_1 = get_GABA_A_Tau_1();
  updateValue<double>(__d, "GABA_A_Tau_1", tmp_GABA_A_Tau_1);


  double tmp_GABA_A_Tau_2 = get_GABA_A_Tau_2();
  updateValue<double>(__d, "GABA_A_Tau_2", tmp_GABA_A_Tau_2);


  double tmp_GABA_A_E_rev = get_GABA_A_E_rev();
  updateValue<double>(__d, "GABA_A_E_rev", tmp_GABA_A_E_rev);


  double tmp_GABA_B_g_peak = get_GABA_B_g_peak();
  updateValue<double>(__d, "GABA_B_g_peak", tmp_GABA_B_g_peak);


  double tmp_GABA_B_Tau_1 = get_GABA_B_Tau_1();
  updateValue<double>(__d, "GABA_B_Tau_1", tmp_GABA_B_Tau_1);


  double tmp_GABA_B_Tau_2 = get_GABA_B_Tau_2();
  updateValue<double>(__d, "GABA_B_Tau_2", tmp_GABA_B_Tau_2);


  double tmp_GABA_B_E_rev = get_GABA_B_E_rev();
  updateValue<double>(__d, "GABA_B_E_rev", tmp_GABA_B_E_rev);


  double tmp_I_e = get_I_e();
  updateValue<double>(__d, "I_e", tmp_I_e);


  long tmp_r = get_r();
  updateValue<long>(__d, "r", tmp_r);

  
// ignores 'alpha_n_init' double' since it is an function and setter isn't defined

  
// ignores 'beta_n_init' double' since it is an function and setter isn't defined

  
// ignores 'alpha_m_init' double' since it is an function and setter isn't defined

  
// ignores 'beta_m_init' double' since it is an function and setter isn't defined

  
// ignores 'alpha_h_init' double' since it is an function and setter isn't defined

  
// ignores 'beta_h_init' double' since it is an function and setter isn't defined

  
// ignores 'r_h_inf_init' double' since it is an function and setter isn't defined

  
// ignores 'Act_A_inf_init' double' since it is an function and setter isn't defined

  
// ignores 'Act_B_inf_init' double' since it is an function and setter isn't defined

  

  double tmp_Inact_h = get_Inact_h();
  updateValue<double>(__d, "Inact_h", tmp_Inact_h);

  

  double tmp_Act_n = get_Act_n();
  updateValue<double>(__d, "Act_n", tmp_Act_n);

  

  double tmp_Act_r_h = get_Act_r_h();
  updateValue<double>(__d, "Act_r_h", tmp_Act_r_h);

  

  double tmp_Act_A = get_Act_A();
  updateValue<double>(__d, "Act_A", tmp_Act_A);

  

  double tmp_Act_B = get_Act_B();
  updateValue<double>(__d, "Act_B", tmp_Act_B);

  

  double tmp_V_m = get_V_m();
  updateValue<double>(__d, "V_m", tmp_V_m);

  

  double tmp_g_AMPA__d = get_g_AMPA__d();
  updateValue<double>(__d, "g_AMPA__d", tmp_g_AMPA__d);

  

  double tmp_g_AMPA = get_g_AMPA();
  updateValue<double>(__d, "g_AMPA", tmp_g_AMPA);

  

  double tmp_g_NMDA__d = get_g_NMDA__d();
  updateValue<double>(__d, "g_NMDA__d", tmp_g_NMDA__d);

  

  double tmp_g_NMDA = get_g_NMDA();
  updateValue<double>(__d, "g_NMDA", tmp_g_NMDA);

  

  double tmp_g_GABAA__d = get_g_GABAA__d();
  updateValue<double>(__d, "g_GABAA__d", tmp_g_GABAA__d);

  

  double tmp_g_GABAA = get_g_GABAA();
  updateValue<double>(__d, "g_GABAA", tmp_g_GABAA);

  

  double tmp_g_GABAB__d = get_g_GABAB__d();
  updateValue<double>(__d, "g_GABAB__d", tmp_g_GABAB__d);

  

  double tmp_g_GABAB = get_g_GABAB();
  updateValue<double>(__d, "g_GABAB", tmp_g_GABAB);

  // We now know that (ptmp, stmp) are consistent. We do not
  // write them back to (P_, S_) before we are also sure that
  // the properties to be set in the parent class are internally
  // consistent.
  Archiving_Node::set_status(__d);

  // if we get here, temporaries contain consistent set of properties


  set_t_ref(tmp_t_ref);



  set_g_Na(tmp_g_Na);



  set_g_K(tmp_g_K);



  set_g_L(tmp_g_L);



  set_g_H(tmp_g_H);



  set_g_A(tmp_g_A);



  set_C_m(tmp_C_m);



  set_E_Na(tmp_E_Na);



  set_E_K(tmp_E_K);



  set_E_L(tmp_E_L);



  set_E_H(tmp_E_H);



  set_E_A(tmp_E_A);



  set_V_Tr(tmp_V_Tr);



  set_AMPA_g_peak(tmp_AMPA_g_peak);



  set_AMPA_E_rev(tmp_AMPA_E_rev);



  set_AMPA_Tau_1(tmp_AMPA_Tau_1);



  set_AMPA_Tau_2(tmp_AMPA_Tau_2);



  set_NMDA_g_peak(tmp_NMDA_g_peak);



  set_NMDA_Tau_1(tmp_NMDA_Tau_1);



  set_NMDA_Tau_2(tmp_NMDA_Tau_2);



  set_NMDA_E_rev(tmp_NMDA_E_rev);



  set_NMDA_Vact(tmp_NMDA_Vact);



  set_NMDA_Sact(tmp_NMDA_Sact);



  set_GABA_A_g_peak(tmp_GABA_A_g_peak);



  set_GABA_A_Tau_1(tmp_GABA_A_Tau_1);



  set_GABA_A_Tau_2(tmp_GABA_A_Tau_2);



  set_GABA_A_E_rev(tmp_GABA_A_E_rev);



  set_GABA_B_g_peak(tmp_GABA_B_g_peak);



  set_GABA_B_Tau_1(tmp_GABA_B_Tau_1);



  set_GABA_B_Tau_2(tmp_GABA_B_Tau_2);



  set_GABA_B_E_rev(tmp_GABA_B_E_rev);



  set_I_e(tmp_I_e);



  set_r(tmp_r);


  // ignores 'alpha_n_init' double' since it is an function and setter isn't defined

  // ignores 'beta_n_init' double' since it is an function and setter isn't defined

  // ignores 'alpha_m_init' double' since it is an function and setter isn't defined

  // ignores 'beta_m_init' double' since it is an function and setter isn't defined

  // ignores 'alpha_h_init' double' since it is an function and setter isn't defined

  // ignores 'beta_h_init' double' since it is an function and setter isn't defined

  // ignores 'r_h_inf_init' double' since it is an function and setter isn't defined

  // ignores 'Act_A_inf_init' double' since it is an function and setter isn't defined

  // ignores 'Act_B_inf_init' double' since it is an function and setter isn't defined


  set_Inact_h(tmp_Inact_h);



  set_Act_n(tmp_Act_n);



  set_Act_r_h(tmp_Act_r_h);



  set_Act_A(tmp_Act_A);



  set_Act_B(tmp_Act_B);



  set_V_m(tmp_V_m);



  set_g_AMPA__d(tmp_g_AMPA__d);



  set_g_AMPA(tmp_g_AMPA);



  set_g_NMDA__d(tmp_g_NMDA__d);



  set_g_NMDA(tmp_g_NMDA);



  set_g_GABAA__d(tmp_g_GABAA__d);



  set_g_GABAA(tmp_g_GABAA);



  set_g_GABAB__d(tmp_g_GABAB__d);



  set_g_GABAB(tmp_g_GABAB);


  
  updateValue< double >(__d, nest::names::gsl_error_tol, P_.__gsl_error_tol);
  if ( P_.__gsl_error_tol <= 0. ){
    throw nest::BadProperty( "The gsl_error_tol must be strictly positive." );
  }
  
};

#endif /* #ifndef OLM_RECEPTORS_NESTML */
#endif /* HAVE GSL */