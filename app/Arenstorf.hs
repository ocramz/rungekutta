-- Module for testing Runge-Kutta integration of ODEs
-- Copyright 2009, Uwe Hollerbach <uh@alumni.caltech.edu>
-- BSD License. See file 'LICENSE' for details.
-- $Id: arenstorf.hs,v 1.2 2009/04/19 18:29:55 uwe Exp $
--
-- Solve the restricted 3-body problem, seek Arenstorf orbits
--
--	x1'' = x1 + 2 x2' - mu'*(x1 + mu)/D1 - mu*(x1 - mu')/D2
--	x2'' = x2 - 2 x1' - mu'*x2/D1 - mu*x2/D2
--	D1 = ((x1 + mu)^2 + x2^2)^(3/2)
--	D2 = ((x1 - mu')^2 + x2^2)^(3/2)
--	mu' = 1 - mu
--
-- Rewrite the second derivatives in terms of new variables v1 and v2
-- to get a system of first-order equations
--
-- Represent the solution as a tuple (x1,v1,x2,v2), then we need three
-- functions to manipulate these: a summer, a scaler, and a derivative
-- operator: all very trivial functions.

module Main where

import System.IO (hPutStrLn, hPutStr, stderr)

import Numeric.RungeKutta

-- AState = (x1, v1, x2, v2)
type AState = (Double, Double, Double, Double)

scaler :: Double -> AState -> AState
scaler s (c1,c2,c3,c4) = (s*c1, s*c2, s*c3, s*c4)

summer :: AState -> AState -> AState
summer (a1,a2,a3,a4) (b1,b2,b3,b4) = (a1+b1, a2+b2, a3+b3, a4+b4)

deriv :: Double -> AState -> AState
deriv t (x1,v1,x2,v2) =
  let yp = x1 + mu
      ym = x1 - mu'
      d1 = sqrt (yp*yp + x2*x2)
      d1c = d1*d1*d1
      d2 = sqrt (ym*ym + x2*x2)
      d2c = d2*d2*d2
      f1 = x1 + 2*v2 - (mu'*yp/d1c) - (mu*ym/d2c)
      f2 = x2 - 2*v1 - (mu'*x2/d1c) - (mu*x2/d2c)
  in (v1,f1,v2,f2)

-- Step-size selector oracle

oracle :: Double -> AState -> Either Double Double
oracle h (x1,v1,x2,v2) =
  let e = (abs x1) + (abs v1) + (abs x2) + (abs v2)
  in if e < 0.5*tol
        then Right (1.414*h) -- step too small, accept but grow
        else if e < tol
             then Right h -- step just right
             else Left (0.5*h) -- step too large, reject and shrink

show_pos t (x1,v1,x2,v2) = (show t) ++ " " ++ (show x1) ++ " " ++ (show x2)
show_st  t (x1,v1,x2,v2) = (show t) ++ " " ++ (show x1) ++ " " ++ (show x2)
                                    ++ " " ++ (show v1) ++ " " ++ (show v2)

gen_soln h t st =
  let hf = tmax - t
      hs = if hf > h then h else hf
      (tn, stn, etn) = stepper hs (t,st)
  in if hs < 1 && hs == hf
        then putStrLn (show_st tn stn) >> putStrLn "# All done!"
        else case oracle h etn of
                  Left bad -> shs "de" bad >> gen_soln bad t st
                  Right good -> (if good > h
                                    then shs "in" good
                                    else shnil) >>
                                putStrLn (show_pos t st) >>
                                gen_soln good tn stn
  where shs pre h = hPutStrLn stderr (pre ++ "creasing step to " ++ (show h))
        shnil = hPutStr stderr ""

-- Parameters of the model, initial conditions, integration time

mu = 0.012277471
mu' = 1 - mu
-- from Hairer Norsett Wanner p 130: 3-loop orbit
tmax3 = 17.0652165601579625588917206249
start3 = (0.994, 0, 0, -2.00158510637908252240537862224)

-- from Hairer Norsett Wanner p 186: 2-loop orbit
tmax2 = 11.124340337266085134999734047
start2 = (0.994, 0, 0, -2.0317326295573368357302057924)

tmax = tmax2
gen_start = return start2

-- Numerical parameters: tolerance to which we want to work,
-- method we want to use, initial stepsize

tol = 1.0e-9
shower = show_rkf78
stepper = rkf78 scaler summer deriv
step = 100.0 -- initial step guaranteed to fail!

main = hPutStrLn stderr shower >> gen_start >>= (gen_soln step 0.0)
