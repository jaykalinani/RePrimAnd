#include "test_utils.h"
#include "con2prim_imhd.h"
#include "eos_thermal.h"
#include <boost/format.hpp>

using EOS_Toolkit::eos_thermal;
using EOS_Toolkit::real_t;
using EOS_Toolkit::sm_metric3;
using EOS_Toolkit::sm_tensor1;
using EOS_Toolkit::atmosphere;
using EOS_Toolkit::con2prim_mhd;
using EOS_Toolkit::prim_vars_mhd;
using EOS_Toolkit::cons_vars_mhd;



struct env_idealgas {
  real_t rho_max; 
  real_t eps_max; 
  real_t atmo_rho; 
  real_t atmo_eps;
  real_t c2p_strict; 
  real_t c2p_zmax;
  real_t c2p_acc;
};


struct env_hybrideos {
  real_t atmo_rho; 
  real_t c2p_strict; 
  real_t c2p_zmax;
  real_t c2p_acc;
};

class prim2con_mhd {
  public:

  eos_thermal eos;
  sm_metric3 g;
  atmosphere atmo;
  
  void setup_prim_cons(prim_vars_mhd& pv, 
                       cons_vars_mhd& cv) const;

  bool check_isfinite(const cons_vars_mhd& cv) const;
  bool check_isfinite(const prim_vars_mhd& pv) const;
  bool check_isfinite(const cons_vars_mhd& cv,
                      const prim_vars_mhd& pv) const;
};

