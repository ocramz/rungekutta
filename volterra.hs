-- Module for testing Runge-Kutta integration of ODEs
-- Copyright 2009, Uwe Hollerbach <uh@alumni.caltech.edu>
-- BSD License. See file 'LICENSE' for details.
-- $Id: volterra.hs,v 1.5 2009/04/19 19:03:11 uwe Exp $
--
-- Solve the Lotka-Volterra equation
--
-- dy1                   dy2
-- --- = y1*(a - b*y2)   --- = -y2*(c - d*y1)
--  dt                    dt
--
-- Represent the solution as a tuple (y1,y2), then we need three functions
-- to manipulate these: a summer, a scaler, and a derivative operator: all
-- very trivial functions.
--
-- Y1 is the prey population, Y2 the predator

module Main where
import System
import System.IO

import Numeric.RungeKutta

type LVState = (Double, Double)

scaler :: Double -> LVState -> LVState
scaler s (y1,y2) = (s*y1, s*y2)

summer :: LVState -> LVState -> LVState
summer (y11,y12) (y21,y22)= (y11 + y21, y12 + y22)

-- Parameters of the model

a = 3    -- prey population growth parameter
b = 2    -- prey decay due to presence of predator
c = 3    -- predator starvation parameter
d = 2    -- predator growth due to presence of prey

deriv :: Double -> LVState -> LVState
deriv t (y1,y2) = (y1*(a - b*y2), -y2*(c - d*y1))

shower = show_rk4b
stepper = rk4b scaler summer deriv

step = 0.01

show_st t (y1,y2) = (show t) ++ " " ++ (show y1) ++ " " ++ (show y2)

gen_soln t st =
  let (tn,stn) = stepper step (t,st)
  in putStrLn (show_st t st) >>
     if t > 40.00001
        then putStrLn "# All done!"
        else gen_soln tn stn

gen_start [] = return (2.0, 0.3)
gen_start (y1:y2:[]) = return (read y1, read y2)
gen_start _ = error "bad number of args specified!"

main = hPutStrLn stderr shower >> getArgs >>= gen_start >>= (gen_soln 0.0)
