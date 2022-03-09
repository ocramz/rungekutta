-- Module for testing Runge-Kutta integration of ODEs
-- Copyright 2009, Uwe Hollerbach <uh@alumni.caltech.edu>
-- BSD License. See file 'LICENSE' for details.
-- $Id: volterra2.hs,v 1.5 2009/04/19 19:03:11 uwe Exp $
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

oracle :: Double -> LVState -> Either Double Double
oracle h (y1,y2) =
  let e = (abs y1) + (abs y2)
  in if e < 0.5e-8
        then Right (1.414*h)		-- step too small, accept but grow
        else if e < 1.0e-8
             then Right h		-- step just right
             else Left (0.5*h)		-- step too large, reject and shrink

shower = show_rkdp
stepper = rkdp scaler summer deriv

show_st t (y1,y2) = (show t) ++ " " ++ (show y1) ++ " " ++ (show y2)

gen_soln h t st =
  let hf = 40 - t
      hs = if hf > h then h else hf
      (tn, stn, etn) = stepper hs (t,st)
  in if hs < 1 && hs == hf
        then putStrLn (show_st tn stn) >> putStrLn "# All done!"
        else case oracle h etn of
                  Left bad -> shs "de" bad >> gen_soln bad t st
                  Right good -> (if good > h
                                    then shs "in" good
                                    else shnil) >>
                                putStrLn (show_st t st) >>
                                gen_soln good tn stn
  where shs pre h = hPutStrLn stderr (pre ++ "creasing step to " ++ (show h))
        shnil = hPutStr stderr ""

-- initial step guaranteed to fail!
step = 100.0

gen_start [] = return (2.0, 0.3)
gen_start (y1:y2:[]) = return (read y1, read y2)
gen_start _ = error "bad number of args specified!"

main = hPutStrLn stderr shower >> getArgs >>= gen_start >>= (gen_soln step 0.0)
