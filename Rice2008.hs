{-# LANGUAGE NoMonomorphismRestriction,DeriveDataTypeable,TemplateHaskell #-}
-- +Require ModML-Units
-- +Require ModML-Reactions
-- +Require typehash
-- +Require containers

module Rice2008
where
import qualified ModML.Units.UnitsDAEModel as U
import qualified ModML.Core.BasicDAEModel as B
import ModML.Units.UnitsDAEOpAliases
import qualified ModML.Reactions.Reactions as R
import qualified Data.Data as D
import qualified Data.TypeHash as D
import ModML.Units.SIUnits
import qualified Control.Monad as M
import Data.Maybe
import Data.List
import qualified Data.Map as M
import Data.Map ((!))
import qualified Debug.Trace as DT

type RExB m = U.ModelBuilderT m U.RealExpression
data Parameters m = Parameters {
      maxSarcomereLength :: RExB m,                           -- SL_{max}
      minSarcomereLength :: RExB m,                           -- SL_{min}
      thickFilamentLength :: RExB m,                          -- length_{thick}
      hbareLength :: RExB m,                                  -- length_{hbare}
      thinFilamentLength :: RExB m,                           -- length_{thin}
      temperature :: RExB m,                                  -- TmpC
      calciumOnTrop :: ReactionParameters m,                  -- k_{on} / Qk_{on}
      calciumOffTropL :: ReactionParameters m,                -- k_{offL} ...
      calciumOffTropH :: ReactionParameters m,                -- k_{offH} ...
      tropomyosinNToP :: ReactionParameters m,                -- k_{n_p} ...
      tropomyosinPToN :: ReactionParameters m,                -- k_{p_n} ...
      crossBridgeFormation :: ReactionParameters m,           -- f_{app} ...
      crossBridgeDissociation :: ReactionParameters m,        -- g_{app} ...
      crossBridgeRotation :: ReactionParameters m,            -- h_f...
      crossBridgeReverseRotation :: ReactionParameters m,     -- h_b...
      rotatedCrossBridgeDissociation :: ReactionParameters m, -- g_{xb}...
      permissiveHalfActivationConstant :: RExB m,             -- perm_{50}
      permissiveHillCoefficient :: RExB m,                    -- n_{perm}
      overlapModStrongToWeak :: RExB m,                       -- gslmod
      preRotStrainFactor :: RExB m,                           -- hfmdc
      strainEffectPositive :: RExB m,                         -- \sigma_p
      strainEffectNegative :: RExB m,                         -- \sigma_n
      meanStrain :: RExB m,                                   -- x_0
      strainScalingFactor :: RExB m,                          -- \phi
      restingSarcomereLength :: RExB m,                       -- SL_{rest}
      passiveTitinConstant :: RExB m,                         -- PCon_{titin}
      passiveTitinExponent :: RExB m,                         -- PExp_{titin}
      sarcomereLengthCollagen :: RExB m,                      -- SL_{collagen}
      passiveCollagenConstant :: RExB m,                      -- PCon_{collagen}
      passiveCollagenExponent :: RExB m,                      -- PExp_{collagen}
      normalisedMass :: RExB m,                               -- Mass
      normalisedViscosity :: RExB m,                          -- Viscosity
      constantAfterload :: RExB m,                            -- F^{constant}_{afterload}
      stiffness :: RExB m,                                    -- KSE
      initialSarcomereLength :: RExB m,
      initialtmNNoXB :: RExB m,
      initialtmPNoXB :: RExB m,
      initialtmNXB :: RExB m,
      initialtmPXB :: RExB m,
      initialXBPreR :: RExB m,
      initialXBPostR :: RExB m,
      initialxXBPreR :: RExB m,
      initialxXBPostR :: RExB m,
      initialCaTropH :: RExB m,
      initialCaTropL :: RExB m,
      calciumTransient :: TransientParameters m,
      modelContext :: ModelContext,
      contractionType :: ContractionType
    }
data ReactionParameters m = ReactionParameters {
      baseRate :: RExB m,
      otherMod :: Maybe (RExB m),
      speciesMod :: Maybe (RExB m),
      q10 :: RExB m
    }
data TransientParameters m = TransientParameters {
  transientStartTime :: RExB m,
  transientBase :: RExB m,
  transientAmplitude :: RExB m,
  transientTime1 :: RExB m,
  transientTime2 :: RExB m
}
data ModelContext = IsolatedCell | Trabeculae
data ContractionType = Isotonic |
                       Isometric -- i.e. constant length

-- Model units
U.declareBaseType "normalisedForce" "normalisedForceBase"
uProbability = U.dimensionless
uProbabilityR = U.liftUnits uProbability
uDistance = uMicro $*$ uMetre
uDistanceR = U.liftUnits uDistance
uConcentration = uMicro $*$ uMole $*$ uLitre $**$ (-1)
uConcentrationR = U.liftUnits uConcentration
uFlux = uConcentration $*$ uSecond $**$ (-1)
uFluxR = U.liftUnits uFlux
uNthOrderRate n = uConcentration $**$ (-n) $*$ uSecond $**$ (-1)
uNormalisedForce = M.liftM U.singletonUnit normalisedForceBase
uNormMass = uNormalisedForce $*$ uSecond$**$2 $*$ uDistance$**$(-1)
uNormViscosity = uNormMass $*$ uSecond$**$(-1)
uNormStiffness = uNormalisedForce $*$ (uDistance $**$(-1))
uForceIntegral = uNormalisedForce $*$ uSecond

R.declareNamedTaggedEntity [e|uProbabilityR|] "Non-permissive tropomyosin not near cross-bridge" "tmNNoXB"
R.declareNamedTaggedEntity [e|uProbabilityR|] "Permissive tropomyosin not near cross-bridge" "tmPNoXB"
R.declareNamedTaggedEntity [e|uProbabilityR|] "Non-permissive tropomyosin near cross-bridge" "tmNXB"
R.declareNamedTaggedEntity [e|uProbabilityR|] "Permissive tropomyosin near cross-bridge" "tmPXB"
R.declareNamedTaggedEntity [e|uProbabilityR|] "Cross bridges (pre rotation)" "xbPreR"
R.declareNamedTaggedEntity [e|uProbabilityR|] "Cross bridges (post rotation)" "xbPostR"
R.declareNamedTaggedEntity [e|uProbabilityR|] "Troponin with calcium bound to the high-affinity regulatory site" "caTropH"
R.declareNamedTaggedEntity [e|uProbabilityR|] "Troponin with calcium bound to the low-affinity regulatory site" "caTropL"
R.declareNamedTaggedEntity [e|uConcentrationR|] "Calcium^(2+) concentration" "calcium"

R.declareNamedTaggedCompartment "Cardiac Muscle Site" "cardiacMuscleSite"

U.declareRealVariable [e|uDistanceR|] "Sarcomere length" "sarcomereLength" -- SL
U.declareRealVariable [e|uDistanceR|] "Mean distortion pre-rotation" "meanDistortionPreR" -- xXBPreR
U.declareRealVariable [e|uDistanceR|] "Mean distortion post-rotation" "meanDistortionPostR" -- xXBPostR
U.declareRealVariable [e|U.dimensionless|] "Fraction of strongly bound crossbridges" "fractSBXB" -- Fract_{SBXB}
U.declareRealVariable [e|uForceIntegral|] "Integral of Force" "integralForce" -- Integral_{Force}

model = parameterisedModel defaultParameters

parameterisedModel p = B.buildModel $ do
  U.unitsToCore uSecond (unitsModel p)

unitsModel :: Monad m => Parameters m -> U.ModelBuilderT m ()
unitsModel p = do
  fromPre <- unitsModelBeforeReaction p
  (cem, _, ces) <-
    R.runReactionBuilderInUnitBuilder' (reactionModelWithCalciumTransient fromPre p)
  unitsModelAfterReaction p fromPre cem ces

unitsModelBeforeReaction :: Monad m => Parameters m -> U.ModelBuilderT m (RExB m)
unitsModelBeforeReaction _ = 
  U.mkNamedRealVariable uProbability
    "Amount of Permissive tropomyosin near cross bridge in Cardiac Muscle Site"

unitsModelAfterReaction ::
  Monad m => Parameters m -> RExB m -> M.Map R.CompartmentEntity U.RealExpression ->
    (R.CompartmentEntity, R.CompartmentEntity, R.CompartmentEntity, R.CompartmentEntity, R.CompartmentEntity) -> 
     U.ModelBuilderT m ()
unitsModelAfterReaction p permTropNorm cem (tmNXBMuscle, xbPreRMuscle, xbPostRMuscle, preRMuscle, postRMuscle) =
  do
    -- Normalisation equation...
    permTropNorm `U.newEq` (U.realConstant uProbability 1 .-. (return $ cem!tmNXBMuscle) .-.
                            (return $ cem!xbPreRMuscle) .-. (return $ cem!xbPostRMuscle))
    
    let preRVar = return $ cem!preRMuscle
    let postRVar = return $ cem!postRMuscle
      
    physicalModel p preRVar postRVar
    return ()

-- The physical model...
physicalModel :: Monad m => Parameters m -> RExB m -> RExB m -> U.ModelBuilderT m ()
physicalModel p preRVar postRVar = do
  meanDistortionModel p
  sarcomereLengthModel p preRVar postRVar

meanDistortionModel p = do
  fappT <- U.realCommonSubexpression $ standardRate (crossBridgeFormation p) (temperature p)
  gappT <- U.realCommonSubexpression $ xbPreRToPermissiveRate p
  hbT <- U.realCommonSubexpression $ standardRate (crossBridgeReverseRotation p) (temperature p)
  gxbT <- U.realCommonSubexpression $ xbPostRToPermissiveRate p
  hfT <- U.realCommonSubexpression $ crossBridgeRotationRate p
  preRTerm <- U.realCommonSubexpression $ fappT .*. (hbT .+. gxbT)
  postRTerm <- U.realCommonSubexpression $ fappT .*. hfT
  tot <- U.realCommonSubexpression $ preRTerm .+. postRTerm .+. gxbT .*. (hfT .+. gappT) .+. gappT .*. hbT
  xbDutyFracPreR <- U.realCommonSubexpression $ preRTerm ./. tot
  xbDutyFracPostR <- U.realCommonSubexpression $ postRTerm ./. tot
  U.newBoundaryEq {- when -} (U.realConstant U.boundUnits 0 .==. U.boundVariable)
                             (U.realVariable meanDistortionPreR) {- == -} (initialxXBPreR p)
  (U.derivative (U.realVariable meanDistortionPreR)) `U.newEq`
    (U.dConstant 0.5 .*. U.derivative (U.realVariable sarcomereLength) .+. (strainScalingFactor p ./. xbDutyFracPreR) .*. (fappT .*. (U.negateX (U.realVariable meanDistortionPreR)) .+. hbT .*. ((U.realVariable meanDistortionPostR) .-. meanStrain p .-. (U.realVariable meanDistortionPreR))))
  U.newBoundaryEq {- when -} (U.realConstant U.boundUnits 0 .==. U.boundVariable)
                             (U.realVariable meanDistortionPostR) {- == -} (initialxXBPostR p)
  (U.derivative (U.realVariable meanDistortionPostR)) `U.newEq`
    (U.dConstant 0.5 .*. U.derivative (U.realVariable sarcomereLength) .+. (strainScalingFactor p ./. xbDutyFracPostR) .*. (hfT .*. ((U.realVariable meanDistortionPreR) .+. meanStrain p .-. (U.realVariable meanDistortionPostR))))

sarcomereLengthModel p preRVar postRVar = do
  fapp <- U.realCommonSubexpression $ baseRate (crossBridgeFormation p)
  gapp <- U.realCommonSubexpression $ baseRate (crossBridgeDissociation p)
  hb <- U.realCommonSubexpression $ baseRate (crossBridgeReverseRotation p)
  gxb <- U.realCommonSubexpression $ baseRate (rotatedCrossBridgeDissociation p)
  hf <- U.realCommonSubexpression $ baseRate (crossBridgeRotation p)
  preRContrib <- U.realCommonSubexpression $ fapp .*. (hb .+. gxb)
  postRContrib <- U.realCommonSubexpression $ fapp .*. hf
  denom <- U.realCommonSubexpression $ preRContrib .+. postRContrib .+. gxb .*. hf .+. gapp .*. (hb .+. gxb)
  xbMaxPreR <- U.realCommonSubexpression $ preRContrib ./. denom
  xbMaxPostR <- U.realCommonSubexpression $ postRContrib ./. denom
  (U.realVariable fractSBXB) `U.newEq` ((preRVar .+. postRVar) ./. (xbMaxPreR .+. xbMaxPostR))
  U.newBoundaryEq {- when -} (U.realConstant U.boundUnits 0 .==. U.boundVariable)
                             (U.realVariable sarcomereLength) {- == -} (initialSarcomereLength p)
  (U.derivative $ U.realVariable sarcomereLength) `U.newEq`
    (((U.realVariable integralForce) .+.
      (initialSarcomereLength p .-. (U.realVariable sarcomereLength)) .*. normalisedViscosity p)
     ./. normalisedMass p)
  U.newBoundaryEq {- when -} (U.realConstant U.boundUnits 0 .==. U.boundVariable)
                             (U.realVariable integralForce) {- == -} (U.realConstant uForceIntegral 0)
    
  vpassiveForce <- U.mkNamedRealVariable uNormalisedForce "Passive Force"
  vpreloadForce <- U.mkNamedRealVariable uNormalisedForce "Preload Force"
  vafterloadForce <- U.mkNamedRealVariable uNormalisedForce "Afterload Force"
  vactiveForce <- U.mkNamedRealVariable uNormalisedForce "Active Force"
  vpassiveForce `U.newEq` (passiveForce p (U.realVariable sarcomereLength))
  vpreloadForce `U.newEq` (preloadForce p)
  vafterloadForce `U.newEq` (afterloadForce p)
  vactiveForce `U.newEq` (singleOverlapThick p (U.realVariable sarcomereLength) .*.
                          (U.realVariable meanDistortionPreR .*.
                           preRVar .+. U.realVariable meanDistortionPostR .*. postRVar) ./.
                          (meanStrain p .*. xbMaxPostR))

  
  (U.derivative $ U.realVariable integralForce) `U.newEq`
    U.negateX (vactiveForce .+. vpassiveForce .-. vpreloadForce .-. vafterloadForce)

reactionModelWithCalciumTransient fromPre p = do
  R.addEntityInstance
       (calcium `R.inCompartment` cardiacMuscleSite)
       (R.entityClamped (standardTransient (calciumTransient p) U.boundVariable))
  reactionModel fromPre p

-- The reaction model...
reactionModel normalisedPermissive p = do
  R.newAllCompartmentProcesses [
    calciumBindingToTroponinSite p caTropH,
    calciumDisassociatingTroponinSite p (calciumOffTropH p) caTropH,
    calciumBindingToTroponinSite p caTropL,
    calciumDisassociatingTroponinSite p (calciumOffTropL p) caTropL,
    nToPNotNearXB p, pToNNotNearXB p,
    nToPNearXB p, pToNNearXB p, pToXBPreR p,
    xbPreRToPermissive p, xbPostRToPermissive p,
    crossBridgePreToPost p, crossBridgePostToPre p
                  ]

  -- Now the entity instances where we have an initial value...
  let zeroProbFlux = U.realConstant (uNthOrderRate 0) 0
  R.addEntityInstance
       (tmNNoXB `R.inCompartment` cardiacMuscleSite)
       (R.entityFromProcesses (initialtmNNoXB p) zeroProbFlux)
  R.addEntityInstance
       (tmPNoXB `R.inCompartment` cardiacMuscleSite)
       (R.entityFromProcesses (initialtmPNoXB p) zeroProbFlux)
  R.addEntityInstance
       (tmNXB `R.inCompartment` cardiacMuscleSite)
       (R.entityFromProcesses (initialtmNXB p) zeroProbFlux)
  -- R.addEntityInstance
  --      (tmPXB `R.inCompartment` cardiacMuscleSite)
  --      (R.entityFromProcesses (initialtmPXB p) zeroProbFlux)
  R.addEntityInstance
       (xbPreR `R.inCompartment` cardiacMuscleSite)
       (R.entityFromProcesses (initialXBPreR p) zeroProbFlux)
  R.addEntityInstance
       (xbPostR `R.inCompartment` cardiacMuscleSite)
       (R.entityFromProcesses (initialXBPostR p) zeroProbFlux)
  -- Clamped so the relative amounts add to 1.0...
  R.addEntityInstance
       (tmPXB `R.inCompartment` cardiacMuscleSite)
       (R.entityClamped normalisedPermissive)

  tmNXBMuscle <- tmNXB `R.inCompartment` cardiacMuscleSite
  xbPreRMuscle <- xbPreR `R.inCompartment` cardiacMuscleSite
  xbPostRMuscle <- xbPostR `R.inCompartment` cardiacMuscleSite
  preRMuscle <- xbPreR `R.inCompartment` cardiacMuscleSite
  postRMuscle <- xbPostR `R.inCompartment` cardiacMuscleSite
  return (tmNXBMuscle, xbPreRMuscle, xbPostRMuscle, preRMuscle, postRMuscle)

standardRate :: Monad m => ReactionParameters m -> RExB m -> RExB m
standardRate (ReactionParameters{baseRate=baseRate,otherMod=otherMod,speciesMod=speciesMod,q10=q10}) temp = do
    let l = catMaybes [Just baseRate, otherMod, speciesMod]
    foldl' (.*.) (q10 .**. ((temp .-. U.realConstant uCelsius 37) ./. U.dConstant 10)) l

standardTransient :: Monad m => TransientParameters m -> RExB m -> RExB m
standardTransient p t = do
  timeRatio <- U.realCommonSubexpression ((transientTime1 p) ./. (transientTime2 p))
  t' <- U.realCommonSubexpression (t .-. (transientStartTime p))
  let dim1 = U.dConstant 1
  let dimm1 = U.dConstant (-1)
  let beta = timeRatio .**. (dimm1 ./. (timeRatio .-. dim1)) .-.
             timeRatio .**. (dimm1 ./. (dim1 .-. (transientTime2 p) ./. (transientTime1 p)))
  U.ifX (t .<=. (transientStartTime p))
    {- then -} (transientBase p)
    {- else -} $ ((transientAmplitude p .-. transientBase p) ./. beta) .*.
                   (U.expX (U.negateX (t' ./. (transientTime1 p))) .-.
                    U.expX (U.negateX (t' ./. (transientTime2 p))))
                   .+. (transientBase p)

-- Functions for sarcomere geometry...
singleOverlapNearestZ p x = U.minX (thickFilamentLength p) x ./. U.dConstant 2
singleOverlapNearestCentreLine p x = U.maxX (x ./. U.dConstant 2 .-.
                                            (x .-. thinFilamentLength p))
                                       (hbareLength p ./. U.dConstant 2)
lengthSingleOverlap p x = singleOverlapNearestZ p x .-. singleOverlapNearestCentreLine p x
singleOverlapThick p x = (U.dConstant 2 .*. lengthSingleOverlap p x) ./. (thickFilamentLength p .-. hbareLength p)
singleOverlapThin p x = lengthSingleOverlap p x ./. thinFilamentLength p
-- Functions for normalised passive force...
titinForce p x = U.ifX (x .>=. restingSarcomereLength p)
                 {- then -}
                   (passiveTitinConstant p .*.
                      (U.expX (passiveTitinExponent p .*.
                               (x .-. restingSarcomereLength p)) .-.
                       U.dConstant 1))
                   {- else -}
                   (U.negateX $ passiveTitinConstant p .*.
                      (U.expX (passiveTitinExponent p .*.
                               (restingSarcomereLength p .-. x)) .-.
                       U.dConstant 1))
collagenForce p x = U.ifX (x .>=. sarcomereLengthCollagen p)
                     {- then -}
                     (passiveCollagenConstant p .*.
                        (U.expX (passiveCollagenExponent p .*.
                                 (x .-. sarcomereLengthCollagen p)) .-.
                         U.dConstant 1))
                     {- else -}
                     (U.realConstant uNormalisedForce 0)

passiveForce p x =
  do
    vtitinForce <- U.mkNamedRealVariable uNormalisedForce "Titin Force"
    vtitinForce `U.newEq` (titinForce p x)
    vcollagenForce <- U.mkNamedRealVariable uNormalisedForce "Collagen Force"
    vcollagenForce `U.newEq` (collagenForce p x)
    case (modelContext p)
      of
        IsolatedCell -> vtitinForce
        Trabeculae -> vtitinForce .+. vcollagenForce

preloadForce p = passiveForce p (initialSarcomereLength p)
afterloadForce p = case (contractionType p)
                   of
                     Isometric -> stiffness p .*. (initialSarcomereLength p .-. U.realVariable sarcomereLength)
                     Isotonic -> constantAfterload p

calciumBindingToTroponinSite p site compartment = do
  cavar <- R.addEntity R.EssentialForProcess R.CantBeCreatedByProcess R.NotModifiedByProcess 0 (calcium `R.withCompartment` compartment)
  sitevar <- R.addEntity R.NotEssentialForProcess R.CanBeCreatedByProcess R.ModifiedByProcess 1 (site `R.withCompartment` compartment)
  let calciumTroponinBindingRateT = standardRate (calciumOnTrop p) (temperature p)
  R.rateEquation $ calciumTroponinBindingRateT .*. cavar .*. (U.realConstant uProbability 1 .-.  sitevar)

calciumDisassociatingTroponinSite p rp site compartment = do
  sitevar <- R.addEntity R.EssentialForProcess R.CantBeCreatedByProcess R.ModifiedByProcess (-1) (site `R.withCompartment` compartment)
  R.rateEquation $ (standardRate rp (temperature p)) .*. sitevar

tropRegulatory p catroph catropl = (U.dConstant 1 .-. singleOverlapThin p (U.realVariable sarcomereLength)) .*. catroph .+.
                                   singleOverlapThin p (U.realVariable sarcomereLength) .*. catropl
permissiveTotal p catroph catropl = (U.dConstant 1 ./. (U.dConstant 1 .+. (permissiveHalfActivationConstant p ./. tropRegulatory p catroph catropl).**. permissiveHillCoefficient p)) .**. U.dConstant 0.5
inversePermissiveTotal p catroph catropl = U.minX (U.dConstant 1 ./. permissiveTotal p catroph catropl) $ U.dConstant 100

nToP pent nent p c = do
  pvar <- R.addEntity R.NotEssentialForProcess R.CanBeCreatedByProcess R.ModifiedByProcess 1 (pent `R.withCompartment` c)
  nvar <- R.addEntity R.EssentialForProcess R.CantBeCreatedByProcess R.ModifiedByProcess (-1) (nent `R.withCompartment` c)
  catroph <- R.addEntity R.NotEssentialForProcess R.CantBeCreatedByProcess R.NotModifiedByProcess 0 (caTropH `R.withCompartment` c)
  catropl <- R.addEntity R.NotEssentialForProcess R.CantBeCreatedByProcess R.NotModifiedByProcess 0 (caTropL `R.withCompartment` c)
  R.rateEquation $ (standardRate ((tropomyosinNToP p) {otherMod = Just $ permissiveTotal p catroph catropl}) (temperature p)) .*. nvar

pToN pent nent p c = do
  pvar <- R.addEntity R.EssentialForProcess R.CantBeCreatedByProcess R.ModifiedByProcess (-1) (pent `R.withCompartment` c)
  nvar <- R.addEntity R.NotEssentialForProcess R.CanBeCreatedByProcess R.ModifiedByProcess 1 (nent `R.withCompartment` c)
  catroph <- R.addEntity R.NotEssentialForProcess R.CantBeCreatedByProcess R.NotModifiedByProcess 0 (caTropH `R.withCompartment` c)
  catropl <- R.addEntity R.NotEssentialForProcess R.CantBeCreatedByProcess R.NotModifiedByProcess 0 (caTropL `R.withCompartment` c)
  R.rateEquation $ (standardRate ((tropomyosinPToN p) {otherMod = Just $ inversePermissiveTotal p catroph catropl}) (temperature p)) .*. pvar

nToPNotNearXB = nToP tmPNoXB tmNNoXB
pToNNotNearXB = pToN tmPNoXB tmNNoXB
nToPNearXB = nToP tmPXB tmNXB
pToNNearXB = pToN tmPXB tmNXB

pToXBPreR p c = do
  pvar <- R.addEntity R.EssentialForProcess R.CantBeCreatedByProcess R.ModifiedByProcess (-1) (tmPXB `R.withCompartment` c)
  xbprervar <- R.addEntity R.NotEssentialForProcess R.CanBeCreatedByProcess R.ModifiedByProcess 1 (xbPreR `R.withCompartment` c)
  R.rateEquation $ (standardRate (crossBridgeFormation p) (temperature p)) .*. pvar

xbToPermissive rate fromXB p c = do
  xbvar <- R.addEntity R.EssentialForProcess R.CantBeCreatedByProcess R.ModifiedByProcess (-1) (fromXB `R.withCompartment` c)
  pvar <- R.addEntity R.NotEssentialForProcess R.CanBeCreatedByProcess R.ModifiedByProcess 1 (tmPXB `R.withCompartment` c)
  R.rateEquation $ rate .*. xbvar

xbPreRToPermissiveRate p =
    standardRate
      ((crossBridgeDissociation p){
         otherMod=Just (U.dConstant 1 .+.
                        (U.dConstant 1 .-. singleOverlapThick p (U.realVariable sarcomereLength)) .*.
                        (overlapModStrongToWeak p))})
      (temperature p)

xbPreRToPermissive p c = 
  xbToPermissive (xbPreRToPermissiveRate p) xbPreR p c

xbPostRToPermissiveRate p =
    let
        mod = U.ifX (U.realVariable meanDistortionPostR .<. meanStrain p)
                {- then -} (U.expX $ strainEffectPositive p .*.
                                 (((meanStrain p .-. U.realVariable meanDistortionPostR) ./. meanStrain p).**.U.dConstant 2))
                {- else -} (U.expX $ strainEffectNegative p .*.
                                 (((meanStrain p .-. U.realVariable meanDistortionPostR) ./. meanStrain p).**.U.dConstant 2))
    in
      standardRate ((rotatedCrossBridgeDissociation p){
                      otherMod=Just mod }) (temperature p)

xbPostRToPermissive p c =
      xbToPermissive (xbPostRToPermissiveRate p) xbPostR p c

crossBridgeRotationRate p = standardRate ((crossBridgeRotation p){otherMod = Just . U.expX $ (U.negateX (U.signX $ U.realVariable meanDistortionPreR)) .*. preRotStrainFactor p .*. (U.realVariable meanDistortionPreR ./. meanStrain p).**.U.dConstant 2}) (temperature p)
crossBridgePreToPost p c = do
  xbPreRvar <- R.addEntity R.EssentialForProcess R.CantBeCreatedByProcess R.ModifiedByProcess (-1) (xbPreR `R.withCompartment` c)
  xbPostRvar <- R.addEntity R.NotEssentialForProcess R.CanBeCreatedByProcess R.ModifiedByProcess 1 (xbPostR `R.withCompartment` c)
  R.rateEquation $ crossBridgeRotationRate p .*. xbPreRvar
crossBridgePostToPre p c = do
  xbPreRvar <- R.addEntity R.NotEssentialForProcess R.CanBeCreatedByProcess R.ModifiedByProcess 1 (xbPreR `R.withCompartment` c)
  xbPostRvar <- R.addEntity R.EssentialForProcess R.CantBeCreatedByProcess R.ModifiedByProcess (-1) (xbPostR `R.withCompartment` c)
  R.rateEquation $ standardRate (crossBridgeReverseRotation p) (temperature p) .*. xbPostRvar

defaultReactionParameters =
    ReactionParameters {
      baseRate=U.realConstant (uNthOrderRate 0) 0,
      otherMod=Nothing,
      speciesMod=Nothing,
      q10=U.dConstant 1
    }

defaultParameters =
  Parameters {
    maxSarcomereLength = U.realConstant uDistance 2.4,
    minSarcomereLength = U.realConstant uDistance 1.4,
    thickFilamentLength = U.realConstant uDistance 1.65,
    hbareLength = U.realConstant uDistance 0.1,
    thinFilamentLength = U.realConstant uDistance 1.2,
    temperature = U.realConstant uCelsius 22.5,
    calciumOnTrop = defaultReactionParameters { baseRate=U.realConstant (uNthOrderRate 1) 50,
                                                q10=U.dConstant 1.5},
    calciumOffTropL = defaultReactionParameters { baseRate=U.realConstant (uNthOrderRate 0) 250,
                                                  q10=U.dConstant 1.3},
    calciumOffTropH = defaultReactionParameters { baseRate=U.realConstant (uNthOrderRate 0) 25,
                                                  q10=U.dConstant 1.3},
    tropomyosinNToP = defaultReactionParameters { baseRate=U.realConstant (uNthOrderRate 0) 500,
                                                  q10=U.dConstant 1.6},
    tropomyosinPToN = defaultReactionParameters { baseRate=U.realConstant (uNthOrderRate 0) 50,
                                                  q10=U.dConstant 1.6 },
    crossBridgeFormation = defaultReactionParameters { baseRate=U.realConstant (uNthOrderRate 0) 500,
                                                       q10=U.dConstant 6.25 },
    crossBridgeDissociation = defaultReactionParameters { baseRate=U.realConstant (uNthOrderRate 0) 70,
                                                          q10=U.dConstant 2.5 },
    crossBridgeRotation = defaultReactionParameters { baseRate=U.realConstant (uNthOrderRate 0) 2000,
                                                      q10=U.dConstant 6.25 },
    crossBridgeReverseRotation = defaultReactionParameters { baseRate=U.realConstant (uNthOrderRate 0) 400,
                                                             q10=U.dConstant 6.25 },
    rotatedCrossBridgeDissociation = defaultReactionParameters { baseRate=U.realConstant (uNthOrderRate 0) 70,
                                                                 q10 = U.dConstant 6.25 },
    permissiveHalfActivationConstant = U.dConstant 0.5,
    permissiveHillCoefficient = U.dConstant 15,
    overlapModStrongToWeak = U.dConstant 6,
    preRotStrainFactor = U.dConstant 5,
    strainEffectPositive = U.dConstant 8,
    strainEffectNegative = U.dConstant 1,
    meanStrain = U.realConstant uDistance 0.007,
    strainScalingFactor = U.dConstant 2,
    restingSarcomereLength = U.realConstant uDistance 1.9,
    passiveTitinConstant = U.realConstant uNormalisedForce 0.002,
    passiveTitinExponent = U.dConstant 10,
    sarcomereLengthCollagen = U.realConstant uDistance 2.25,
    passiveCollagenConstant = U.realConstant uNormalisedForce 0.02,
    passiveCollagenExponent = U.dConstant 70,
    normalisedMass = U.realConstant uNormMass 0.00005,
    normalisedViscosity = U.realConstant uNormViscosity 0.003,
    constantAfterload = U.realConstant uNormalisedForce 0.001,
    stiffness = U.realConstant uNormStiffness 50.0,
    initialSarcomereLength = U.realConstant uDistance 2.3,
    initialtmNNoXB = U.realConstant uProbability 0.9999,
    initialtmPNoXB = U.realConstant uProbability 0.0001,
    initialtmNXB = U.realConstant uProbability 0.9998,
    initialtmPXB = U.realConstant uProbability 0.0000,
    initialXBPreR = U.realConstant uProbability 0.0001,
    initialXBPostR = U.realConstant uProbability 0.0001,
    initialxXBPreR = U.realConstant uDistance 0,
    initialxXBPostR = U.realConstant uDistance 0.007,
    initialCaTropH = U.realConstant uProbability 0.001,
    initialCaTropL = U.realConstant uProbability 0.001,
    calciumTransient = TransientParameters {
                         transientStartTime = U.realConstant uSecond 0.0,
                         transientBase = U.realConstant uConcentration 0.09,
                         transientAmplitude = U.realConstant uConcentration 1.45,
                         transientTime1 = U.realConstant uSecond 0.02,
                         transientTime2 = U.realConstant uSecond 0.11 },
    modelContext = Trabeculae,
    contractionType = Isometric
  }
