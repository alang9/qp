{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE CPP #-}
module Numeric.QpOases where

import Control.Monad
import Control.Monad.IO.Class
import Control.Monad.Trans.Except
import qualified Data.Vector.Storable as VS
import Foreign hiding (void)
import Foreign.C.Types
import Numeric.LinearAlgebra
import Numeric.LinearAlgebra.Devel

data OasesError = OasesError deriving (Show)

data HessianType = HST_ZERO | HST_IDENTITY | HST_POSDEF | HST_POSDEF_NULLSPACE | HST_SEMIDEF | HST_INDEF | HST_UNKNOWN deriving (Enum, Show)

newtype SQProblem = SQProblem (ForeignPtr SQProblem) deriving (Eq, Show)

setupSQP :: Int -> Int -> HessianType -> IO SQProblem
setupSQP nV nC h = do
    ptr2 <- alloca $ \ptr -> do
        !_ <- sqproblem_setup (fromIntegral nV) (fromIntegral nC) (fromIntegral $ fromEnum h) ptr
        peek ptr
    SQProblem <$> newForeignPtr sqproblem_cleanup ptr2

withSQProblem :: SQProblem -> (Ptr SQProblem -> IO a) -> IO a
withSQProblem (SQProblem fPtr) f = withForeignPtr fPtr f

initSQP :: SQProblem -> Matrix Double -> Vector Double -> Matrix Double
        -> Vector Double -> Vector Double
        -> ExceptT OasesError IO (Vector Double, Double)
initSQP sqp h g a lbA ubA = liftIO $
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
      with (10 * fromIntegral (nVar + nConstr)) $ \nWSRPtr -> do
        _ <- withSQProblem sqp $ \ptr ->
               sqproblem_init ptr hPtr gPtr aPtr nullPtr nullPtr lbAPtr ubAPtr
               nWSRPtr nullPtr xPtr yPtr objPtr statusPtr
        !obj <- peek objPtr
        !status <- peek statusPtr
        return $ (obj, status)
    let solutionVec = VS.unsafeFromForeignPtr0 primalFP (fromIntegral nVar)
    return (solutionVec, obj)
  where
    mat' m f = mat (cmat m) $ \church -> church $ \nrow ncol ptr -> f nrow ncol ptr
    vec' v f = vec v $ \church -> church $ \size ptr -> f size ptr

data Options

foreign import ccall "sqproblem_init"
    sqproblem_init
        :: Ptr SQProblem -> Ptr Double -> Ptr Double -> Ptr Double -> Ptr Double
        -> Ptr Double -> Ptr Double -> Ptr Double -> Ptr Int -> Ptr Double
        -> Ptr Double -> Ptr Double -> Ptr Double -> Ptr Int -> IO CInt

foreign import ccall "sqproblem_setup"
    sqproblem_setup
        :: CInt -> CInt -> CInt -> Ptr (Ptr SQProblem)-> IO CInt

foreign import ccall "&sqproblem_cleanup"
    sqproblem_cleanup
        :: FunPtr (Ptr SQProblem -> IO ())