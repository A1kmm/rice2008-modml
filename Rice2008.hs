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
import qualified Control.Monad as M
import Data.Maybe
import Data.List

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
uNormStiffness = uNormMass $*$ uDistance
uForceIntegral = uNormalisedForce $*$ uSecond
type RExB m = U.ModelBuilderT m U.RealExpression

data ReactionParameters m = ReactionParameters {
      baseRate :: RExB m,
      otherMod :: Maybe (RExB m),
      speciesMod :: Maybe (RExB m),
      q10 :: RExB m
    }
standardRate :: Monad m => ReactionParameters m -> RExB m -> RExB m
standardRate (ReactionParameters{baseRate=baseRate,otherMod=otherMod,speciesMod=speciesMod,q10=q10}) temp = do
    let l = catMaybes [Just baseRate, otherMod, speciesMod]
    foldl' (.*.) (q10 .**. ((temp .-. U.realConstant uCelsius 37) ./. U.dConstant 10)) l

data TransientParameters m = TransientParameters {
  transientStartTime :: RExB m,
  transientBase :: RExB m,
  transientAmplitude :: RExB m,
  transientTime1 :: RExB m,
  transientTime2 :: RExB m
}
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
                 (U.expX (U.negateX (t' ./. (transientTime1 p))) .-. U.expX (U.negateX (t' ./. (transientTime2 p))))
                 .+. (transientBase p)

data ModelContext = IsolatedCell | Trabeculae

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
      xPosition :: RExB m,
      modelContext :: ModelContext
    }

defaultReactionParameters = ReactionParameters {baseRate=U.realConstant (uNthOrderRate 0) 0, otherMod=Nothing, speciesMod=Nothing,
                                                q10=U.dConstant 1 }

defaultParameters =
  Parameters {
    maxSarcomereLength = U.realConstant uDistance 2.4,
    minSarcomereLength = U.realConstant uDistance 1.4,
    thickFilamentLength = U.realConstant uDistance 1.65,
    hbareLength = U.realConstant uDistance 0.1,
    thinFilamentLength = U.realConstant uDistance 1.2,
    temperature = U.realConstant uCelsius 37,
    calciumOnTrop = defaultReactionParameters { baseRate=U.realConstant (uNthOrderRate 1) 50,
                                                q10=U.dConstant 1.5},
    calciumOffTropL = defaultReactionParameters { baseRate=U.realConstant (uNthOrderRate 0) 250,
                                                  q10=U.dConstant 1.3},
    calciumOffTropH = defaultReactionParameters { baseRate=U.realConstant (uNthOrderRate 0) 25,
                                                  q10=U.dConstant 1.3},
    tropomyosinNToP = defaultReactionParameters { baseRate=U.realConstant (uNthOrderRate 0) 50,
                                                  q10=U.dConstant 1.6},
    tropomyosinPToN = defaultReactionParameters { baseRate=U.realConstant (uNthOrderRate 0) 500,
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
                                                                 q10 = U.dConstant 2.5 },
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
    constantAfterload = U.realConstant uNormalisedForce 0.5,
    stiffness = U.realConstant uNormStiffness 100.5,
    initialSarcomereLength = U.realConstant uDistance 1.9,
    initialtmNNoXB = U.realConstant uProbability 0.99,
    initialtmPNoXB = U.realConstant uProbability 0.01,
    initialtmNXB = U.realConstant uProbability 0.97,
    initialtmPXB = U.realConstant uProbability 0.01,
    initialXBPreR = U.realConstant uProbability 0.01,
    initialXBPostR = U.realConstant uProbability 0.01,
    initialxXBPreR = U.realConstant uDistance 0,
    initialxXBPostR = U.realConstant uDistance 1,
    initialCaTropH = U.realConstant uProbability 0,
    initialCaTropL = U.realConstant uProbability 0,
    calciumTransient = TransientParameters {
                         transientStartTime = U.realConstant uSecond 0.0,
                         transientBase = U.realConstant uConcentration 0.09,
                         transientAmplitude = U.realConstant uConcentration 1.45,
                         transientTime1 = U.realConstant uSecond 0.02,
                         transientTime2 = U.realConstant uSecond 0.11 },
    xPosition = U.realConstant uDistance 0.5,
    modelContext = Trabeculae
  }

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
U.declareRealVariable [e|uDistanceR|] "Mean distortion, pre-rotation" "meanDistortionPreR" -- xXBPreR
U.declareRealVariable [e|uDistanceR|] "Mean distortion, post-rotation" "meanDistortionPostR" -- xXBPostR
U.declareRealVariable [e|uDistanceR|] "Integral of Force" "integralForce" -- Integral_{Force}

-- Functions for sarcomere geometry...
singleOverlapNearestZ p = U.minX (thickFilamentLength p) (xPosition p) ./. U.dConstant 2
singleOverlapNearestCentreLine p = U.maxX (xPosition p ./. U.dConstant 2 .-.
                                           (xPosition p .-. thinFilamentLength p))
                                          (hbareLength p ./. U.dConstant 2)
lengthSingleOverlap p = singleOverlapNearestZ p .-. singleOverlapNearestCentreLine p
singleOverlapThick p = (U.dConstant 2 .*. lengthSingleOverlap p) ./. (thickFilamentLength p .-. hbareLength p)
singleOverlapThin p = lengthSingleOverlap p ./. thinFilamentLength p
-- Functions for normalised passive force...
titinForce p = U.ifX (xPosition p .>=. restingSarcomereLength p)
                 {- then -}
                 (passiveTitinConstant p .*.
                    (U.expX (passiveTitinExponent p .*.
                             (xPosition p .-. restingSarcomereLength p)) .-.
                     U.dConstant 1))
                 {- else -}
                 (U.negateX $ passiveTitinConstant p .*.
                    (U.expX (passiveTitinExponent p .*.
                             (restingSarcomereLength p .-. xPosition p)) .-.
                     U.dConstant 1))
collagenForce p = U.ifX (xPosition p .>=. sarcomereLengthCollagen p)
                   {- then -}
                   (passiveCollagenConstant p .*.
                      (U.expX (passiveCollagenExponent p .*.
                               (xPosition p .-. sarcomereLengthCollagen p)) .-.
                       U.dConstant 1))
                   {- else -}
                   (U.realConstant uNormalisedForce 0)

passiveForce p = case (modelContext p)
                 of
                   IsolatedCell -> titinForce p
                   Trabeculae -> titinForce p .+. collagenForce p

calciumBindingToTroponinSite p site compartment = do
  cavar <- R.addEntity R.EssentialForProcess R.CantBeCreatedByProcess R.NotModifiedByProcess 0 (calcium `R.withCompartment` compartment)
  sitevar <- R.addEntity R.NotEssentialForProcess R.CanBeCreatedByProcess R.ModifiedByProcess 1 (site `R.withCompartment` compartment)
  let calciumTroponinBindingRateT = standardRate (calciumOnTrop p) (temperature p)
  R.rateEquation $ calciumTroponinBindingRateT .*. cavar .*. (U.realConstant uProbability 1 .-.  sitevar)

calciumDisassociatingTroponinSite p rp site compartment = do
  sitevar <- R.addEntity R.EssentialForProcess R.CantBeCreatedByProcess R.ModifiedByProcess (-1) (site `R.withCompartment` compartment)
  R.rateEquation $ (standardRate rp (temperature p)) .*. sitevar

tropRegulatory p catroph catropl = (U.dConstant 1 .-. singleOverlapThin p) .*. catroph .+.
                                   singleOverlapThin p .*. catropl
permissiveTotal p catroph catropl = (U.dConstant 1 ./. (U.dConstant 1 .+. (permissiveHalfActivationConstant p ./. tropRegulatory p catroph catropl).**. permissiveHillCoefficient)) .**. U.dConstant 0.5
inversePermissiveTotal p catroph catropl = U.minX (U.dConstant 1 ./. permissiveTotal p catroph catropl) $ U.dConstant 100

nToP pent nent p c = do
  pvar <- R.addEntity R.NotEssentialForProcess R.CanBeCreatedByProcess R.ModifiedByProcess 1 (pent `R.withCompartment` c)
  nvar <- R.addEntity R.EssentialForProcess R.CantBeCreatedByProcess R.ModifiedByProcess (-1) (nent `R.withCompartment` c)
  catroph <- R.addEntity R.NotEssentialForProcess R.CantBeCreatedByProcess R.NotModifiedByProcess 1 (caTropH `R.withCompartment` c)
  catropl <- R.addEntity R.NotEssentialForProcess R.CantBeCreatedByProcess R.NotModifiedByProcess 1 (caTropL `R.withCompartment` c)
  R.rateEquation $ (standardRate (tropomyosinNToP p $ {otherMod = permissiveTotal p catroph catropl}) (temperature p)) .*. nvar

pToN pent nent p c = do
  pvar <- R.addEntity R.EssentialForProcess R.CantBeCreatedByProcess R.ModifiedByProcess (-1) (pent `R.withCompartment` c)
  nvar <- R.addEntity R.NotEssentialForProcess R.CanBeCreatedByProcess R.ModifiedByProcess 1 (nent `R.withCompartment` c)
  catroph <- R.addEntity R.NotEssentialForProcess R.CantBeCreatedByProcess R.NotModifiedByProcess 1 (caTropH `R.withCompartment` c)
  catropl <- R.addEntity R.NotEssentialForProcess R.CantBeCreatedByProcess R.NotModifiedByProcess 1 (caTropL `R.withCompartment` c)
  R.rateEquation $ (standardRate (tropomyosinPToN p $ {otherMod = inversePermissiveTotal p catroph catropl}) (temperature p)) .*. pvar

nToPNotNearXB = nToP tmPNoXB tmNNoXB
pToNNotNearXB = nToP tmPNoXB tmNNoXB
nToPNearXB = nToP tmPXB tmNXB
pToNNearXB = nToP tmPXB tmNXB

pToXBPreR p c = do
  pvar <- R.addEntity R.EssentialForProcess R.CantBeCreatedByProcess R.ModifiedByProcess (-1) (tmPXB `R.withCompartment` c)
  xbprervar <- R.addEntity R.NotEssentialForProcess R.CanBeCreatedByProcess R.ModifiedByProcess (-1) (xbPreR `R.withCompartment` c)
  R.rateEquation $ (standardRate (crossBridgeFormation p) (temperature p)) .*. pvar

xbToPermissive rate fromXB p c = do
  xbvar <- R.addEntity R.EssentialForProcess R.CantBeCreatedByProcess R.ModifiedByProcess (-1) (fromXB `R.withCompartment` c)
  pvar <- R.addEntity R.NotEssentialForProcess R.CanBeCreatedByProcess R.ModifiedByProcess 1 (tmPXB `R.withCompartment` c)
  R.rateEquation $ rate .*. xbvar

xbPreRToPermissiveRate p =
    standardRate
      ((crossBridgeDissociation p){
         otherMod=Just (U.dConstant 1 .+.
                        (U.dConstant 1 .-. singleOverlapThick p) .*.
                        (overlapModStrongToWeak p))})
      (temperature p)

xbPreRToPermissive p c = 
  xbToPermissive (xbPreRToPermissiveRate p) xbPreR p c

xbPostRToPermissiveRate p =
    let
        mod = U.ifX (meanDistortionPostR .<. meanStrain p)
                {- then -} (U.expX $ strainEffectPositive p .*.
                                 (((meanStrain p .-. meanDistortionPostR) ./. meanStrain p).**.U.dConstant 2))
                {- else -} (U.expX $ strainEffectNegative p .*.
                                 (((meanStrain p .-. meanDistortionPostR) ./. meanStrain p).**.U.dConstant 2))
    in
      standardRate ((rotatedCrossBridgeDissociation p){
                      otherMod=Just mod }) (temperature p)

xbPostRToPermissive p c =
      xbToPermissive (xbPostRToPermissiveRate p) xbPostR p c

crossBridgeRotationRate p = standardRate ((crossBridgeRotation p){otherMod = (U.negateX (U.signX meanDistortionPreR)) .*. preRotStrainFactor p .*. (meanDistortionPreR ./. meanStrain p).**.U.dConstant 2}) (temperature p)
crossBridgePreToPost p c = do
  xbPreRvar <- R.addEntity R.EssentialForProcess R.CantBeCreatedByProcess R.ModifiedByProcess (-1) (xbPreR `R.withCompartment` c)
  xbPostRvar <- R.addEntity R.NotEssentialForProcess R.CanBeCreatedByProcess R.ModifiedByProcess 1 (xbPostR `R.withCompartment` c)
  R.rateEquation $ crossBridgeRotationRate p .*. xbPreRvar
crossBridgePostToPre p c = do
  xbPreRvar <- R.addEntity R.NotEssentialForProcess R.CanBeCreatedByProcess R.ModifiedByProcess 1 (xbPreR `R.withCompartment` c)
  xbPostRvar <- R.addEntity R.EssentialForProcess R.CantBeCreatedByProcess R.ModifiedByProcess (-1) (xbPostR `R.withCompartment` c)
  R.rateEquation $ standardRate (crossBridgeReverseRotationRate p) (temperature p) .*. xbPostRvar

reactionModel p = do
  R.newAllCompartmentProcess $ calciumBindingToTroponinSite p caTropH
  R.newAllCompartmentProcess $ calciumDisassociatingTroponinSite p (calciumOffTropH p) caTropH
  R.newAllCompartmentProcess $ calciumBindingToTroponinSite p caTropL
  R.newAllCompartmentProcess $ calciumDisassociatingTroponinSite p (calciumOffTropL p) caTropL
  R.newAllCompartmentProcess $ nToPNotNearXB p
  R.newAllCompartmentProcess $ pToNNotNearXB p
  R.newAllCompartmentProcess $ nToPNearXB p
  R.newAllCompartmentProcess $ pToNNearXB p
  R.newAllCompartmentProcess $ xbPreRToPermissive p
  R.newAllCompartmentProcess $ xbPostRToPermissive p
  R.newAllCompartmentProcess $ crossBridgePreToPost p
  R.newAllCompartmentProcess $ crossBridgePostToPre p

physicalModel p = do
  fappT <- U.realCommonSubexpressionX $ standardRate (crossBridgeFormation p) (temperature p)
  gappT <- U.realCommonSubexpressionX $ xbPreRToPermissiveRate p
  hbT <- U.realCommonSubexpressionX $ standardRate (crossBridgeReverseRotation p) (temperature p)
  gxbT <- U.realCommonSubexpressionX $ xbPostRToPermissiveRate p
  hfT <- U.realCommonSubexpressionX $ crossBridgeRotationRate p
  preRTerm <- U.realCommonSubexpressionX $ fappT .*. (hbT .+. gxbT)
  postRTerm <- U.realCommonSubexpressionX $ fappT .*. hfT
  tot <- U.realCommonSubexpressionX $ preRTerm .+. postRTerm .+. gxbT .*. (hfT .+. gappT) .+. gappT .*. hbT
  xbDutyFracPreR <- U.realCommonSubexpressionX $ preRTerm ./. tot
  xbDutyFracPostR <- U.realCommonSubexpressionX $ postRTerm ./. tot
  (U.derivative meanDistortionPreR) `U.newEq`
    (U.dConstant 0.5 .*. U.derivative sarcomereLength

unitsModel :: Monad m => U.ModelBuilderT m ()
unitsModel = do
  R.runReactionBuilderInUnitBuilder (reactionModel defaultParameters)
  

model = B.buildModel $ do
  U.unitsToCore uSecond unitsModel
  physicalModel