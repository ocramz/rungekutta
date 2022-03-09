-- Module for testing Runge-Kutta integration of ODEs
-- Copyright 2009, Uwe Hollerbach <uh@alumni.caltech.edu>
-- BSD License. See file 'LICENSE' for details.
-- $Id: testrk.hs,v 1.6 2009/04/19 18:29:55 uwe Exp $

module Main where
import Numeric.RungeKutta

scaler s y = s*y
summer y1 y2 = y1 + y2

-- Solve the test problem dy/dt = -t^3 y^3
--
-- which has the solution y = 1/sqrt(C + t^4/2)
-- where C = 1/y_0^2 - t_0^4/2

-- deriv t y = -((t*y)^3)
-- t0 = 0.0
-- y0 = 1.75
-- cc = 1.0/(y0^2) - 0.5*(t0^4)
-- exact t = 1.0/(sqrt (cc + 0.5*t^4))

-- Solve the test problem dy/dt = -0.4 y
--
-- which has the solution y = y_0*exp(-0.4*(t - t_0))

con = -0.4
deriv t y = con*y
t0 = 0.0
y0 = 1.75
exact t = y0*exp(con*(t - t0))

showst t y = (show y) ++ "\t" ++ (show (y - (exact t)))

gen_soln1 stepr t st =
  let (tn,stn) = stepr (t,st)
  in if t >= 5.0
        then putStrLn (showst t st)
        else gen_soln1 stepr tn stn

gen_soln2 stepr t st =
  let (tn,stn,en) = stepr (t,st)
  in if t >= 5.0
        then putStrLn (showst t st)
        else gen_soln2 stepr tn stn

do_case1 stepr n =
  let h = if n < 0 then 2^(-n) else 0.5^n
      sep = if n <= 4 then "\t\t" else "\t"
  in putStr ((show h) ++ sep) >> gen_soln1 (stepr h) t0 y0

do_case2 stepr n =
  let h = if n < 0 then 2^(-n) else 0.5^n
      sep = if n <= 4 then "\t\t" else "\t"
  in putStr ((show h) ++ sep) >> gen_soln2 (stepr h) t0 y0

do_stepper1 (stepr,stats) =
  putStrLn stats >>
  putStrLn "# step yf err" >>
  (mapM (do_case1 (stepr scaler summer deriv)) [(-2) .. 12]) >>
  putStrLn "# All done!\n"

do_stepper2 (stepr,stats) =
  putStrLn stats >>
  putStrLn "# step yf err" >>
  (mapM (do_case2 (stepr scaler summer deriv)) [(-2) .. 12]) >>
  putStrLn "# All done!\n"

main =
 do putStrLn "#### Non-Adaptive Solvers"
    mapM do_stepper1 [(rkfe, show_rkfe),
                      (rk3, show_rk3),
                      (rk4a, show_rk4a),
                      (rk4b, show_rk4b)]
    putStrLn "#### Adaptive Solvers"
    mapM do_stepper2 [(rkhe, show_rkhe),
                      (rkbs, show_rkbs),
                      (rkf45, show_rkf45),
                      (rkck, show_rkck),
                      (rkdp, show_rkdp),
                      (rkf78, show_rkf78),
                      (rkv65, show_rkv65)]
    putStrLn "#### Auxiliary Solvers: Error Estimators from Adaptives"
    mapM do_stepper1 [(rkhe_aux, show_rkhe_aux),
                      (rkbs_aux, show_rkbs_aux),
                      (rkf45_aux, show_rkf45_aux),
                      (rkck_aux, show_rkck_aux),
                      (rkdp_aux, show_rkdp_aux),
                      (rkf78_aux, show_rkf78_aux),
                      (rkv65_aux, show_rkv65_aux)]
