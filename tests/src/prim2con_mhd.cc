#define BOOST_TEST_MODULE con2prim
#include <boost/test/unit_test.hpp>
#include <limits>

#include "test_config.h"
#include "unitconv.h"
#include "test_con2prim_mhd.h"

#include "eos_thermal.h"
#include "eos_idealgas.h"
#include "eos_thermal_file.h"
#include "eos_hybrid.h"

using boost::format;
using boost::str;

using namespace std;
using namespace EOS_Toolkit;


void prim2con_mhd::setup_prim_cons(
        prim_vars_mhd& pv, cons_vars_mhd& cv) const
{
  pv            = prim_vars_mhd(rho, eps, ye, press0, 
                                vel0, wl0, E0, B0);
  cv.from_prim(pv, g);
}

bool test_con2prim_mhd::check_isfinite(
          const cons_vars_mhd& cv) const
{
  failcount hope("All conserved vars finite");

  hope.isfinite(cv.dens, "dens");
  hope.isfinite(cv.tau, "tau");
  hope.isfinite(cv.tracer_ye, "tracer_ye");
  hope(check_isfinite(cv.scon), "scon finite");
  hope(check_isfinite(cv.bcons), "bcons finite");

  return hope;
}

bool test_con2prim_mhd::check_isfinite(
          const prim_vars_mhd& pv) const
{
  failcount hope("All primitive vars finite");
  
  hope.isfinite(pv.rho, "rho"); 
  hope.isfinite(pv.eps, "eps"); 
  hope.isfinite(pv.ye, "ye");
  hope.isfinite(pv.press, "press");
  hope.isfinite(pv.w_lor, "w_lor");
  hope(check_isfinite(pv.vel), "vel finite");
  hope(check_isfinite(pv.E), "E finite");
  hope(check_isfinite(pv.B), "B finite");
  
  return hope;
}

bool test_con2prim_mhd::check_isfinite(
          const cons_vars_mhd& cv,
          const prim_vars_mhd& pv) const
{
  return check_isfinite(pv) && check_isfinite(cv);
}


test_con2prim_mhd make_env(const env_idealgas& e) 
{
  
  auto eos = make_eos_idealgas(1.0, e.eps_max, e.rho_max);

  const real_t atmo_ye{ 0.5 };
  const real_t atmo_cut{ e.atmo_rho * 1.01 };
  const real_t atmo_p{ 
    eos.at_rho_eps_ye(e.atmo_rho, e.atmo_eps, atmo_ye).press() 
  };

  atmosphere atmo{e.atmo_rho, e.atmo_eps, atmo_ye, atmo_p, atmo_cut};

  const real_t rho_strict{ e.c2p_strict * e.atmo_rho };
  const bool   ye_lenient{ false };
  const int    max_iter{ 30 };
  const real_t max_b{ 10. };
  con2prim_mhd cv2pv(eos, rho_strict, ye_lenient, e.c2p_zmax, max_b, 
                     atmo, e.c2p_acc, max_iter);

  sm_metric3 g;
  g.minkowski();
  
  return test_con2prim_mhd(eos, g, atmo, cv2pv);
};

}



