{-# LANGUAGE NoMonomorphismRestriction,DeriveDataTypeable,TemplateHaskell #-}
-- +Require ModML-Units
-- +Require ModML-Reactions
-- +Require typehash

module Rice2008
where
import qualified ModML.Units.UnitsDAEModel as U
import qualified ModML.Core.BasicDAEModel as B
import ModML.Units.UnitsDAEOpAliases
import qualified ModML.Reactions.Reactions as R
import qualified Data.Data as D
import qualified Data.TypeHash as D
import ModML.Units.SIUnits

uProbability = U.dimensionless
uProbabilityR = U.liftUnits uProbability
uDistance = uMicro $*$ uMetre
uDistanceR = U.liftUnits uDisplacement

R.declareNamedTaggedEntity [e|uProbabilityR|] "tmNNoXB" "Non-permissive tropomyosin not near cross-bridge"
R.declareNamedTaggedEntity [e|uProbabilityR|] "tmPNoXB" "Permissive tropomyosin not near cross-bridge"
R.declareNamedTaggedEntity [e|uProbabilityR|] "tmNXB" "Non-permissive tropomyosin near cross-bridge"
R.declareNamedTaggedEntity [e|uProbabilityR|] "tmPXB" "Permissive tropomyosin near cross-bridge"
R.declareNamedTaggedEntity [e|uProbabilityR|] "xbPreR" "Cross bridges (pre rotation)"
R.declareNamedTaggedEntity [e|uProbabilityR|] "xbPostR" "Cross bridges (post rotation)"
R.declareNamedTaggedEntity [e|uDistanceR|] "xXbPreR" "Average distortion of cross bridges (pre rotation)"
R.declareNamedTaggedEntity [e|uDistanceR|] "xXbPostR" "Average distortion of cross bridges (post rotation)"
R.declareNamedTaggedEntity [e|uProbabilityR|] "caTropH" "Troponin with calcium bound to the high-affinity regulatory site"
R.declareNamedTaggedEntity [e|uProbabilityR|] "caTropL" "Troponin with calcium bound to the low-affinity regulatory site"

R.declareRealVariable [e|uDistanceR|] "sarcomereLength" "Sarcomere length"

calciumBindingToSite site = do
  

reactionModel = do
  R.newAllCompartmentProcess 
      

unitsModel :: Monad m => U.ModelBuilderT m ()
unitsModel = do
  R.runReactionBuilderInUnitBuilder reactionModel

model = B.buildModel $ do
  U.unitsToCore uSecond unitsModel
