{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE CPP #-}
module Numeric.QpOases where

import Control.Monad
import Control.Monad.IO.Class
import Control.Monad.Trans.Except
import qualified Data.Vector.Storable as VS
import Foreign
import Foreign.C.Types
import Numeric.LinearAlgebra
import Numeric.LinearAlgebra.Devel

data OasesError = OasesError deriving (Show)

#define HST_POSDEF           2

solveQP :: Matrix Double -> Vector Double -> Matrix Double
        -> Vector Double -> Vector Double
        -> ExceptT OasesError IO (Vector Double, Double)
solveQP h g a lbA ubA = liftIO $
  mat' h $ \hRow hCol hPtr ->
  vec' g $ \gSize gPtr ->
  mat' a $ \aRow aCol aPtr ->
  vec' lbA $ \lbASize lbAPtr ->
  vec' ubA $ \ubASize ubAPtr -> do 
    let nVar = gSize
    guard $ nVar == hRow
    guard $ nVar == hCol
    guard $ nVar == aCol
    let nConstr = aRow
    guard $ nConstr == lbASize
    guard $ nConstr == ubASize
    primalFP <- mallocForeignPtrArray (fromIntegral nVar)
    (obj, status) <- withForeignPtr primalFP $ \xPtr ->
      allocaArray (fromIntegral $ nVar + nConstr) $ \yPtr ->
      alloca $ \objPtr ->
      alloca $ \statusPtr ->
      with (5 * fromIntegral (nVar + nConstr)) $ \nWSRPtr -> do
        _ <- c_SQProblem_setup nVar nConstr HST_POSDEF
        _ <- c_SQProblem_init hPtr gPtr aPtr nullPtr nullPtr lbAPtr ubAPtr nWSRPtr nullPtr nullPtr xPtr yPtr objPtr statusPtr
        !obj <- peek objPtr
        !status <- peek statusPtr
        _ <- c_SQProblem_cleanup
        return $ (obj, status)
    let solutionVec = VS.unsafeFromForeignPtr0 primalFP (fromIntegral nVar)
    return (solutionVec, obj)
  where
    mat' m f = mat (cmat m) $ \church -> church $ \nrow ncol ptr -> f nrow ncol ptr
    vec' v f = vec v $ \church -> church $ \size ptr -> f size ptr

data Options

foreign import ccall "SQProblem_init"
    c_SQProblem_init
        :: Ptr Double -> Ptr Double -> Ptr Double -> Ptr Double -> Ptr Double -> Ptr Double -> Ptr Double -> Ptr Int -> Ptr Double -> Ptr Options -> Ptr Double -> Ptr Double -> Ptr Double -> Ptr Int -> IO CInt

foreign import ccall "SQProblem_setup"
    c_SQProblem_setup
        :: CInt -> CInt -> CInt -> IO CInt

foreign import ccall "SQProblem_cleanup"
    c_SQProblem_cleanup
        :: IO CInt
