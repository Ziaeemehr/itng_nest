
/*
*  wb_cond_exp_nestml.h
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
*  2020-03-26 06:36:35.460951
*/
#ifndef WB_COND_EXP_NESTML
#define WB_COND_EXP_NESTML

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
extern "C" inline int wb_cond_exp_nestml_dynamics( double, const double y[], double f[], void* pnode );


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
  Name: wb_cond_exp_nestml.

  Description:  
    
  Name: wb_cond_exp - Wang buzsaki model

  Description:

  wb_cond_exp is an implementation of a modified Hodkin-Huxley model
  (1) Post-synaptic currents
   Incoming spike events induce a post-synaptic change of conductance modeled
   by an exponential function.

  (2) Spike Detection
   Spike detection is done by a combined threshold-and-local-maximum search: if
   there is a local maximum above a certain threshold of the membrane potential,
   it is considered a spike.

  References:

  Wang, X.J. and Buzsaki, G., (1996) Gamma oscillation by synaptic inhibition in a hippocampal interneuronal network model. Journal of neuroscience, 16(20), pp.6402-6413.


  SeeAlso: hh_cond_exp_traub



  Parameters:
  The following parameters can be set in the status dictionary.
  t_ref [ms]  Refractory period
  g_Na [nS]  Sodium peak conductance
  g_K [nS]  Potassium peak conductance
  g_L [nS]  Leak conductance
  C_m [pF]  Membrane Capacitance
  E_Na [mV]  Sodium reversal potential
  E_K [mV]  Potassium reversal potentia
  E_L [mV]  Leak reversal Potential (aka resting potential)
  V_Tr [mV]  Spike Threshold
  tau_syn_ex [ms]  Rise time of the excitatory synaptic alpha function i
  tau_syn_in [ms]  Rise time of the inhibitory synaptic alpha function
  E_ex [mV]  Excitatory synaptic reversal potential
  E_in [mV]  Inhibitory synaptic reversal potential
  I_e [pA]  constant external input current
  

  Dynamic state variables:
  r [integer]  number of steps in the current refractory phase
  

  Initial values:
  V_m [mV]  Membrane potential
  Inact_h [real]  Act_m real = alpha_m_init / ( alpha_m_init + beta_m_init )
  

  References: Empty

  Sends: nest::SpikeEvent

  Receives: Spike, Current, DataLoggingRequest
*/
class wb_cond_exp_nestml : public nest::Archiving_Node{
public:
  /**
  * The constructor is only used to create the model prototype in the model manager.
  */
  wb_cond_exp_nestml();

  /**
  * The copy constructor is used to create model copies and instances of the model.
  * @node The copy constructor needs to initialize the parameters and the state.
  *       Initialization of buffers and interal variables is deferred to
  *       @c init_buffers_() and @c calibrate().
  */
  wb_cond_exp_nestml(const wb_cond_exp_nestml &);

  /**
  * Releases resources.
  */
  ~wb_cond_exp_nestml();

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
  friend class nest::RecordablesMap<wb_cond_exp_nestml>;
  friend class nest::UniversalDataLogger<wb_cond_exp_nestml>;

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
        
        

    //!  Refractory period
    double t_ref;

    //!  Sodium peak conductance
    double g_Na;

    //!  Potassium peak conductance
    double g_K;

    //!  Leak conductance
    double g_L;

    //!  Membrane Capacitance
    double C_m;

    //!  Sodium reversal potential
    double E_Na;

    //!  Potassium reversal potentia
    double E_K;

    //!  Leak reversal Potential (aka resting potential)
    double E_L;

    //!  Spike Threshold
    double V_Tr;

    //!  Rise time of the excitatory synaptic alpha function i
    double tau_syn_ex;

    //!  Rise time of the inhibitory synaptic alpha function
    double tau_syn_in;

    //!  Excitatory synaptic reversal potential
    double E_ex;

    //!  Inhibitory synaptic reversal potential
    double E_in;

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
    //  Membrane potential
      V_m,
      
      Act_n,
      //  Act_m real = alpha_m_init / ( alpha_m_init + beta_m_init )
      Inact_h,
      
      g_in,
      
      g_ex,
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

    double I_syn_exc;

    double I_syn_inh;

    double alpha_n;

    double beta_n;

    double alpha_m;

    double beta_m;

    double alpha_h;

    double beta_h;

    //!  alias Act_m real = alpha_m / ( alpha_m + beta_m ) 
    //!  function I_Na  pA = g_Na * Act_m * Act_m * Act_m * Inact_h * ( V_m - E_Na )
    double I_Na;

    double I_K;

    double I_L;
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
    Buffers_(wb_cond_exp_nestml &);
    Buffers_(const Buffers_ &, wb_cond_exp_nestml &);

    /** Logger for all analog data */
    nest::UniversalDataLogger<wb_cond_exp_nestml> logger_;
    
    inline nest::RingBuffer& get_spikeInh() {return spikeInh;}
    //!< Buffer incoming nSs through delay, as sum
    nest::RingBuffer spikeInh;
    double spikeInh_grid_sum_;
    
    inline nest::RingBuffer& get_spikeExc() {return spikeExc;}
    //!< Buffer incoming nSs through delay, as sum
    nest::RingBuffer spikeExc;
    double spikeExc_grid_sum_;
    
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
    return (-0.05) / (1.0 * 1.0) * (S_.ode_state[State_::V_m] + 34.0*1.0) / (std::exp((-0.1) * (S_.ode_state[State_::V_m] + 34.0*1.0)) - 1.0);
  }

  inline double get_beta_n_init() const {
    return 0.625 / 1.0 * std::exp((-(S_.ode_state[State_::V_m] + 44.0*1.0)) / 80.0*1.0);
  }

  inline double get_alpha_m_init() const {
    return 0.1 / (1.0 * 1.0) * (S_.ode_state[State_::V_m] + 35.0*1.0) / (1.0 - std::exp((-0.1*1.0) * (S_.ode_state[State_::V_m] + 35.0*1.0)));
  }

  inline double get_beta_m_init() const {
    return 4.0 / (1.0) * std::exp((-(S_.ode_state[State_::V_m] + 60.0*1.0)) / 18.0*1.0);
  }

  inline double get_alpha_h_init() const {
    return 0.35 / 1.0 * std::exp((-(S_.ode_state[State_::V_m] + 58.0*1.0)) / 20.0*1.0);
  }

  inline double get_beta_h_init() const {
    return 5.0 / (std::exp((-0.1) / 1.0 * (S_.ode_state[State_::V_m] + 28.0*1.0)) + 1.0) / 1.0;
  }

  inline double get_V_m() const {
    return S_.ode_state[State_::V_m];
  }
  inline void set_V_m(const double __v) {
    S_.ode_state[State_::V_m] = __v;
  }

  inline double get_Act_n() const {
    return S_.ode_state[State_::Act_n];
  }
  inline void set_Act_n(const double __v) {
    S_.ode_state[State_::Act_n] = __v;
  }

  inline double get_Inact_h() const {
    return S_.ode_state[State_::Inact_h];
  }
  inline void set_Inact_h(const double __v) {
    S_.ode_state[State_::Inact_h] = __v;
  }

  inline double get_g_in() const {
    return S_.ode_state[State_::g_in];
  }
  inline void set_g_in(const double __v) {
    S_.ode_state[State_::g_in] = __v;
  }

  inline double get_g_ex() const {
    return S_.ode_state[State_::g_ex];
  }
  inline void set_g_ex(const double __v) {
    S_.ode_state[State_::g_ex] = __v;
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

  inline double get_V_Tr() const {
    return P_.V_Tr;
  }
  inline void set_V_Tr(const double __v) {
    P_.V_Tr = __v;
  }

  inline double get_tau_syn_ex() const {
    return P_.tau_syn_ex;
  }
  inline void set_tau_syn_ex(const double __v) {
    P_.tau_syn_ex = __v;
  }

  inline double get_tau_syn_in() const {
    return P_.tau_syn_in;
  }
  inline void set_tau_syn_in(const double __v) {
    P_.tau_syn_in = __v;
  }

  inline double get_E_ex() const {
    return P_.E_ex;
  }
  inline void set_E_ex(const double __v) {
    P_.E_ex = __v;
  }

  inline double get_E_in() const {
    return P_.E_in;
  }
  inline void set_E_in(const double __v) {
    P_.E_in = __v;
  }

  inline double get_I_e() const {
    return P_.I_e;
  }
  inline void set_I_e(const double __v) {
    P_.I_e = __v;
  }

  inline long get_RefractoryCounts() const {
    return V_.RefractoryCounts;
  }
  inline void set_RefractoryCounts(const long __v) {
    V_.RefractoryCounts = __v;
  }

  inline double get_I_syn_exc() const {
    return S_.ode_state[State_::g_ex] * (S_.ode_state[State_::V_m] - P_.E_ex);
  }

  inline double get_I_syn_inh() const {
    return S_.ode_state[State_::g_in] * (S_.ode_state[State_::V_m] - P_.E_in);
  }

  inline double get_alpha_n() const {
    return (-0.05) / (1.0 * 1.0) * (S_.ode_state[State_::V_m] + 34.0*1.0) / (std::exp((-0.1) * (S_.ode_state[State_::V_m] + 34.0*1.0)) - 1.0);
  }

  inline double get_beta_n() const {
    return 0.625 / 1.0 * std::exp((-(S_.ode_state[State_::V_m] + 44.0*1.0)) / 80.0*1.0);
  }

  inline double get_alpha_m() const {
    return 0.1 / (1.0 * 1.0) * (S_.ode_state[State_::V_m] + 35.0*1.0) / (1.0 - std::exp((-0.1*1.0) * (S_.ode_state[State_::V_m] + 35.0*1.0)));
  }

  inline double get_beta_m() const {
    return 4.0 / (1.0) * std::exp((-(S_.ode_state[State_::V_m] + 60.0*1.0)) / 18.0*1.0);
  }

  inline double get_alpha_h() const {
    return 0.35 / 1.0 * std::exp((-(S_.ode_state[State_::V_m] + 58.0*1.0)) / 20.0*1.0);
  }

  inline double get_beta_h() const {
    return 5.0 / (std::exp((-0.1) / 1.0 * (S_.ode_state[State_::V_m] + 28.0*1.0)) + 1.0) / 1.0;
  }

  inline double get_I_Na() const {
    return P_.g_Na * (0.1 / (1.0 * 1.0) * (S_.ode_state[State_::V_m] + 35.0*1.0) / (1.0 - std::exp((-0.1*1.0) * (S_.ode_state[State_::V_m] + 35.0*1.0)))) / ((0.1 / (1.0 * 1.0) * (S_.ode_state[State_::V_m] + 35.0*1.0) / (1.0 - std::exp((-0.1*1.0) * (S_.ode_state[State_::V_m] + 35.0*1.0)))) + (4.0 / (1.0) * std::exp((-(S_.ode_state[State_::V_m] + 60.0*1.0)) / 18.0*1.0))) * (0.1 / (1.0 * 1.0) * (S_.ode_state[State_::V_m] + 35.0*1.0) / (1.0 - std::exp((-0.1*1.0) * (S_.ode_state[State_::V_m] + 35.0*1.0)))) / ((0.1 / (1.0 * 1.0) * (S_.ode_state[State_::V_m] + 35.0*1.0) / (1.0 - std::exp((-0.1*1.0) * (S_.ode_state[State_::V_m] + 35.0*1.0)))) + (4.0 / (1.0) * std::exp((-(S_.ode_state[State_::V_m] + 60.0*1.0)) / 18.0*1.0))) * (0.1 / (1.0 * 1.0) * (S_.ode_state[State_::V_m] + 35.0*1.0) / (1.0 - std::exp((-0.1*1.0) * (S_.ode_state[State_::V_m] + 35.0*1.0)))) / ((0.1 / (1.0 * 1.0) * (S_.ode_state[State_::V_m] + 35.0*1.0) / (1.0 - std::exp((-0.1*1.0) * (S_.ode_state[State_::V_m] + 35.0*1.0)))) + (4.0 / (1.0) * std::exp((-(S_.ode_state[State_::V_m] + 60.0*1.0)) / 18.0*1.0))) * S_.ode_state[State_::Inact_h] * (S_.ode_state[State_::V_m] - P_.E_Na);
  }

  inline double get_I_K() const {
    return P_.g_K * S_.ode_state[State_::Act_n] * S_.ode_state[State_::Act_n] * S_.ode_state[State_::Act_n] * S_.ode_state[State_::Act_n] * (S_.ode_state[State_::V_m] - P_.E_K);
  }

  inline double get_I_L() const {
    return P_.g_L * (S_.ode_state[State_::V_m] - P_.E_L);
  }


  
  inline nest::RingBuffer& get_spikeInh() {return B_.get_spikeInh();};
  
  inline nest::RingBuffer& get_spikeExc() {return B_.get_spikeExc();};
  
  inline nest::RingBuffer& get_I_stim() {return B_.get_I_stim();};
  

  // Generate function header
  
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
  static nest::RecordablesMap<wb_cond_exp_nestml> recordablesMap_;

  friend int wb_cond_exp_nestml_dynamics( double, const double y[], double f[], void* pnode );
  
/** @} */
}; /* neuron wb_cond_exp_nestml */

inline nest::port wb_cond_exp_nestml::send_test_event(
    nest::Node& target, nest::rport receptor_type, nest::synindex, bool){
  // You should usually not change the code in this function.
  // It confirms that the target of connection @c c accepts @c nest::SpikeEvent on
  // the given @c receptor_type.
  nest::SpikeEvent e;
  e.set_sender(*this);
  return target.handles_test_event(e, receptor_type);
}

inline nest::port wb_cond_exp_nestml::handles_test_event(nest::SpikeEvent&, nest::port receptor_type){
  
    // You should usually not change the code in this function.
    // It confirms to the connection management system that we are able
    // to handle @c SpikeEvent on port 0. You need to extend the function
    // if you want to differentiate between input ports.
    if (receptor_type != 0)
      throw nest::UnknownReceptorType(receptor_type, get_name());
    return 0;
}



inline nest::port wb_cond_exp_nestml::handles_test_event(
    nest::CurrentEvent&, nest::port receptor_type){
  // You should usually not change the code in this function.
  // It confirms to the connection management system that we are able
  // to handle @c CurrentEvent on port 0. You need to extend the function
  // if you want to differentiate between input ports.
  if (receptor_type != 0)
  throw nest::UnknownReceptorType(receptor_type, get_name());
  return 0;
}

inline nest::port wb_cond_exp_nestml::handles_test_event(
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
inline void wb_cond_exp_nestml::get_status(DictionaryDatum &__d) const{  
  def<double>(__d, "t_ref", get_t_ref());
      
  def<double>(__d, "g_Na", get_g_Na());
      
  def<double>(__d, "g_K", get_g_K());
      
  def<double>(__d, "g_L", get_g_L());
      
  def<double>(__d, "C_m", get_C_m());
      
  def<double>(__d, "E_Na", get_E_Na());
      
  def<double>(__d, "E_K", get_E_K());
      
  def<double>(__d, "E_L", get_E_L());
      
  def<double>(__d, "V_Tr", get_V_Tr());
      
  def<double>(__d, "tau_syn_ex", get_tau_syn_ex());
      
  def<double>(__d, "tau_syn_in", get_tau_syn_in());
      
  def<double>(__d, "E_ex", get_E_ex());
      
  def<double>(__d, "E_in", get_E_in());
      
  def<double>(__d, "I_e", get_I_e());
      
  def<long>(__d, "r", get_r());
      
  def<double>(__d, "alpha_n_init", get_alpha_n_init());
      
  def<double>(__d, "beta_n_init", get_beta_n_init());
      
  def<double>(__d, "alpha_m_init", get_alpha_m_init());
      
  def<double>(__d, "beta_m_init", get_beta_m_init());
      
  def<double>(__d, "alpha_h_init", get_alpha_h_init());
      
  def<double>(__d, "beta_h_init", get_beta_h_init());
      
  def<double>(__d, "V_m", get_V_m());
      
  def<double>(__d, "Act_n", get_Act_n());
      
  def<double>(__d, "Inact_h", get_Inact_h());
      
  def<double>(__d, "g_in", get_g_in());
      
  def<double>(__d, "g_ex", get_g_ex());
    

  (*__d)[nest::names::recordables] = recordablesMap_.get_list();
  
  def< double >(__d, nest::names::gsl_error_tol, P_.__gsl_error_tol);
  if ( P_.__gsl_error_tol <= 0. ){
    throw nest::BadProperty( "The gsl_error_tol must be strictly positive." );
  }
  

}

inline void wb_cond_exp_nestml::set_status(const DictionaryDatum &__d){

  double tmp_t_ref = get_t_ref();
  updateValue<double>(__d, "t_ref", tmp_t_ref);


  double tmp_g_Na = get_g_Na();
  updateValue<double>(__d, "g_Na", tmp_g_Na);


  double tmp_g_K = get_g_K();
  updateValue<double>(__d, "g_K", tmp_g_K);


  double tmp_g_L = get_g_L();
  updateValue<double>(__d, "g_L", tmp_g_L);


  double tmp_C_m = get_C_m();
  updateValue<double>(__d, "C_m", tmp_C_m);


  double tmp_E_Na = get_E_Na();
  updateValue<double>(__d, "E_Na", tmp_E_Na);


  double tmp_E_K = get_E_K();
  updateValue<double>(__d, "E_K", tmp_E_K);


  double tmp_E_L = get_E_L();
  updateValue<double>(__d, "E_L", tmp_E_L);


  double tmp_V_Tr = get_V_Tr();
  updateValue<double>(__d, "V_Tr", tmp_V_Tr);


  double tmp_tau_syn_ex = get_tau_syn_ex();
  updateValue<double>(__d, "tau_syn_ex", tmp_tau_syn_ex);


  double tmp_tau_syn_in = get_tau_syn_in();
  updateValue<double>(__d, "tau_syn_in", tmp_tau_syn_in);


  double tmp_E_ex = get_E_ex();
  updateValue<double>(__d, "E_ex", tmp_E_ex);


  double tmp_E_in = get_E_in();
  updateValue<double>(__d, "E_in", tmp_E_in);


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

  

  double tmp_V_m = get_V_m();
  updateValue<double>(__d, "V_m", tmp_V_m);

  

  double tmp_Act_n = get_Act_n();
  updateValue<double>(__d, "Act_n", tmp_Act_n);

  

  double tmp_Inact_h = get_Inact_h();
  updateValue<double>(__d, "Inact_h", tmp_Inact_h);

  

  double tmp_g_in = get_g_in();
  updateValue<double>(__d, "g_in", tmp_g_in);

  

  double tmp_g_ex = get_g_ex();
  updateValue<double>(__d, "g_ex", tmp_g_ex);

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



  set_C_m(tmp_C_m);



  set_E_Na(tmp_E_Na);



  set_E_K(tmp_E_K);



  set_E_L(tmp_E_L);



  set_V_Tr(tmp_V_Tr);



  set_tau_syn_ex(tmp_tau_syn_ex);



  set_tau_syn_in(tmp_tau_syn_in);



  set_E_ex(tmp_E_ex);



  set_E_in(tmp_E_in);



  set_I_e(tmp_I_e);



  set_r(tmp_r);


  // ignores 'alpha_n_init' double' since it is an function and setter isn't defined

  // ignores 'beta_n_init' double' since it is an function and setter isn't defined

  // ignores 'alpha_m_init' double' since it is an function and setter isn't defined

  // ignores 'beta_m_init' double' since it is an function and setter isn't defined

  // ignores 'alpha_h_init' double' since it is an function and setter isn't defined

  // ignores 'beta_h_init' double' since it is an function and setter isn't defined


  set_V_m(tmp_V_m);



  set_Act_n(tmp_Act_n);



  set_Inact_h(tmp_Inact_h);



  set_g_in(tmp_g_in);



  set_g_ex(tmp_g_ex);


  
  updateValue< double >(__d, nest::names::gsl_error_tol, P_.__gsl_error_tol);
  if ( P_.__gsl_error_tol <= 0. ){
    throw nest::BadProperty( "The gsl_error_tol must be strictly positive." );
  }
  
};

#endif /* #ifndef WB_COND_EXP_NESTML */
#endif /* HAVE GSL */