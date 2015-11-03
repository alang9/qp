{-# LANGUAGE RecordWildCards #-}

module Data.Matrix.CSC where

import Control.Arrow (second)
import Control.Exception
import Data.Maybe
import qualified Data.Vector.Storable as VS
import qualified Data.Vector as VB
import Foreign.C.Types
import Numeric.LinearAlgebra

data CSC a = CSC
    { cscVals :: VS.Vector a
    , cscCols :: VS.Vector CInt
    , cscRows :: VS.Vector CInt
    , cscNRows :: CInt
    , cscNCols :: CInt
    } deriving (Eq, Show)

fromMatrix :: (Element a, VS.Storable a, Eq a, Num a) => Matrix a -> CSC a
fromMatrix mat = CSC{..}
  where
    cscNRows = fromIntegral $ rows mat
    cscNCols = fromIntegral $ cols mat
    cscRows = VS.fromList $ concat irs
    cscVals = VS.fromList $ concat vals
    cscCols = VS.fromList $ scanl (\a b -> a + fromIntegral (length b)) 0 vals
    (irs, vals) = unzip $ map (go . VS.toList) (toColumns mat)
    go col = unzip $ catMaybes $ zipWith (\idx val -> if val == 0 then Nothing else Just (fromIntegral idx, val)) [0..] col

toMatrix :: (Element a, Num a, VS.Storable a) => CSC a -> Matrix a
toMatrix CSC {..} = fromColumns cols
  where
    cols = zipWith go rows vals
    go rs vs = VS.replicate (fromIntegral cscNRows) 0 VS.//
        (zipWith (\r v -> (fromIntegral r, v)) (VS.toList rs) (VS.toList vs))
    rows = split cscRows cscCols
    vals = split cscVals cscCols
    split v = snd . VS.foldr splitter (v, []) . VS.init
      where
        splitter idx (v, xs) = second (:xs) $ VS.splitAt (fromIntegral idx) v    

(===) :: (VS.Storable a) => CSC a -> CSC a -> CSC a
a === b = assert (cscNCols a == cscNCols b) $ CSC
    { cscNCols = cscNCols a
    , cscNRows = cscNRows a + cscNRows b
    , cscRows = rows
    , cscVals = vals
    , cscCols = cols
    }
  where
    cols = VS.zipWith (+) (cscCols a) (cscCols b)
    rows = VS.concat $ concat $
        zipWith (\a b -> [a, b]) (split (cscRows a) (cscCols a))
        (split (VS.map (+ cscNRows a) $ cscRows b) (cscCols b))
    vals = VS.concat $ concat $
        zipWith (\a b -> [a, b]) (split (cscVals a) (cscCols a))
        (split (cscVals b) (cscCols b))
    split v = snd . VS.foldr splitter (v, [])
      where
        splitter idx (v, xs) = second (:xs) $ VS.splitAt (fromIntegral idx) v

(|||) :: (VS.Storable a) => CSC a -> CSC a -> CSC a
a ||| b = assert (cscNRows a == cscNRows b) $ CSC
    { cscNCols = cscNCols a + cscNCols b
    , cscNRows = cscNRows a
    , cscRows = cscRows a VS.++ cscRows b
    , cscVals = cscVals a VS.++ cscVals b
    , cscCols = cscCols a VS.++ (VS.map (+ VS.last (cscCols a)) $ VS.tail (cscCols b))
    }

diagBlock :: (VS.Storable a) => CSC a -> CSC a -> CSC a
diagBlock a b = CSC
    { cscNCols = cscNCols a + cscNCols b
    , cscNRows = cscNRows a + cscNRows b
    , cscRows = cscRows a VS.++ VS.map (+ cscNRows a) (cscRows b)
    , cscVals = cscVals a VS.++ cscVals b
    , cscCols = cscCols a VS.++ (VS.map (+ VS.last (cscCols a)) $ VS.tail (cscCols b))
    }

zero :: (Integral a, VS.Storable b) => a -> a -> CSC b
zero rows cols = CSC
    { cscNCols = fromIntegral cols
    , cscNRows = fromIntegral rows
    , cscRows = VS.empty
    , cscVals = VS.empty
    , cscCols = VS.replicate (fromIntegral cols + 1) 0
    }

fromBlocks :: (VS.Storable a) => [[CSC a]] -> CSC a
fromBlocks = foldr1 (===) . map (foldr1 (|||))

scale :: (VS.Storable a, Num a) => a -> CSC a -> CSC a
scale mult csc = csc { cscVals = VS.map (* mult) (cscVals csc) }
