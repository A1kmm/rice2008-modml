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

U.declareBaseType "normalisedForce" "normalisedForceBase"
uNormalisedForce = M.liftM U.singletonUnit normalisedForceBase
uProbability = U.dimensionless
uProbabilityR = U.liftUnits uProbability
uDistance = uMicro $*$ uMetre
uDistanceR = U.liftUnits uDistance
uConcentration = uMicro $*$ uMol $*$ uLitre $**$ (-1)
uConcentrationR = U.liftUnits uConcentration
uFlux = uConcentration $*$ uSecond $**$ (-1)
uFluxR = U.liftUnits uFlux
uNthOrderRate n = uConcentration $**$ (-n) $*$ uSecond $**$ (-1)

type RExB = U.ModelBuilderT m U.RealExpression

data ReactionParameters = Parameters {
      baseRate :: RExB,
      otherMod :: Maybe RExB,
      speciesMod :: Maybe RExB,
      q10 :: RExB
    }
standardRate :: Monad m => ReactionParameters -> RExB -> RExB
standardRate (ReactionParameters{baseRate=baseRate,otherMod=otherMod,speciesMod=speciesMod,q10=q10}) temp = do
    l <- catMaybes [Just baseRate, otherMod, speciesMod]
    foldl' (.*.) (q10 .**. ((temp .-. U.realConstant uCelsius 37) ./. realConstant uDimensionless 10)) l

data TransientParameters = TransientParameters {
  transientStartTime :: RExB,
  transientBase :: RExB,
  transientAmplitude :: RExB,
  transientTime1 :: RExB,
  transientTime2 :: RExB
}
standardTransient :: TransientParameters -> RExB -> RExB
standardTransient p t = do
  timeRatio <- realCommonSubexpression ((transientTime1 p) ./. (transientTime2 p))
  t' <- realCommonSubexpression (t .-. (transientStartTime p))
  let dim1 = U.dConstant 1
  let dimm1 = U.dConstant (-1)
  let beta = timeRatio .**. (dimm1 ./. (timeRatio .-. dim1)) .-.
             timeRatio .**. (dimm1 ./. (dim1 .-. (transientTime2 p) ./. (transientTime1 p)))
  ifX (t .<=. (transientStartTime p))
    {- then -} (transientBase p)
    {- else -} $ ((transientAmplitude p .-. transientBase p) ./. beta) .*.
                 (U.expX (U.negateX (t' ./. (transientTime1 p))) .-. U.expX (U.negateX (t' ./. (transientTime2 p))))
                 .+. (transientBase p)

data Parameters = Parameters {
      maxSarcomereLength :: RExB,                           -- SL_{max}
      minSarcomereLength :: RExB,                           -- SL_{min}
      thickFilamentLength :: RExB,                          -- length_{thick}
      hbareLength :: RExB,                                  -- length_{hbare}
      thinFilamentLength :: RExB,                           -- length_{thin}
      temperature :: RExB,                                  -- TmpC
      calciumOnTrop :: ReactionParameters,                  -- k_{on} / Qk_{on}
      calciumOffTropL :: ReactionParameters,                -- k_{offL} ...
      calciumOffTropH :: ReactionParameters,                -- k_{offH} ...
      tropomyosinNToP :: ReactionParameters,                -- k_{n_p} ...
      tropomyosinPToN :: ReactionParameters,                -- k_{p_n} ...
      crossBridgeFormation :: ReactionParameters,           -- f_{app} ...
      crossBridgeDissociation :: ReactionParameters,        -- g_{app} ...
      crossBridgeRotation :: ReactionParameters,            -- h_f...
      crossBridgeReverseRotation :: ReactionParameters,     -- h_b...
      rotatedCrossBridgeDissociation :: ReactionParameters, -- g_{xb}...
      permissiveHalfActivationConstant :: RExB,             -- perm_{50}
      permissiveHillCoefficient :: RExB,                    -- n_{perm}
      overlapModStrongToWeak :: RExB,                       -- gslmod
      preRotStrainFactor :: RExB,                           -- hfmdc
      strainEffectPositive :: RExB,                         -- \sigma_p
      strainEffectNegative :: RExB,                         -- \sigma_n
      meanStrain :: RExB,                                   -- x_0
      strainScalingFactor :: RExB,                          -- \phi
      restingSarcomereLength :: RExB,                       -- SL_{rest}
      passiveTitinConstant :: RExB,                         -- PCon_{titin}
      passiveTitinExponent :: RExB,                         -- PExp_{titin}
      sarcomereLengthCollagen :: RExB,                      -- SL_{collagen}
      passiveCollagenConstant :: RExB,                      -- PCon_{collagen}
      passiveCollagenExponent :: RExB,                      -- PExp_{collagen}
      normalisedMass :: RExB,                               -- Mass
      normalisedViscosity :: RExB,                          -- Viscosity
      constantAfterload :: RExB,                            -- F^{constant}_{afterload}
      stiffness :: RExB,                                    -- KSE
      initialSarcomereLength :: RExB,
      initialtmNNoXB :: RExB,
      initialtmPNoXB :: RExB,
      initialtmNXB :: RExB,
      initialtmPXB :: RExB,
      initialXBPreR :: RExB,
      initialXBPostR :: RExB,
      initialxXbPreR :: RExB,
      initialxXbPostR :: RExB,
      initialCaTropH :: RExB,
      initialCaTropL :: RExB,
      calciumTransient :: TransientParameters,
      xPosition :: RExB
    }

defaultReactionParameters = ReactionParameters {baseRate=U.realConstant (uNthOrderRate 0) 0, otherMod=Nothing, speciesMod=Nothing,
                                                q10=U.dConstant 1 }

defaultParameters =
  Parameters {
    maxSarcomereLength = realConstant uDistance 2.4,
    minSarcomereLength = realConstant uDistance 1.4,
    thickFilamentLength = realConstant uDistance 1.65,
    hbareLength = realConstant uDistance 0.1,
    thinFilamentLength = realConstant uDistance 1.2,
    temperature = realConstant uCelsius 37,
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
    passive
  }

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
R.declareNamedTaggedEntity [e|uConcentrationR|] "calcium" "Calcium^(2+) concentration"

R.declareRealVariable [e|uDistanceR|] "sarcomereLength" "Sarcomere length"

calciumBindingToTroponinSite p site = do
  cavar <- R.addEntity R.EssentialForProcess R.CantBeCreatedByProcess R.NotModifiedByProcess 0 calcium
  sitevar <- R.addEntity R.EssentialForProcess R.CanBeCreatedByProcess R.ModifiedByProcess 1 site
  let calciumTroponinBindingRateT = standardRate p (calciumBinding p)
  R.rateEquation $ calciumTroponinBindingRateT .*. cavar .*. (U.realConstant uProbability 1 .-.  sitevar)

calciumDisassociatingTroponinSite p rp site = do
  sitevar <- R.addEntity R.EssentialForProcess R.CantBeCreatedByProcess R.ModifiedByProcess (-1) site
  R.rateEquation $ (standardRate baserate Nothing Nothing q10) .*. sitevar

nToPNotNearXB = do
  

reactionModel = do
  R.newAllCompartmentProcess (calciumBindingToTroponinSite caTropH)
  R.newAllCompartmentProcess (calciumDisassociatingTroponinSite caTropH calciumOffTropBaseRate calciumOffTropQ10)
  R.newAllCompartmentProcess (calciumBindingToTroponinSite caTropL)
  R.newAllCompartmentProcess (calciumDisassociatingTroponinSite caTropH calciumOffTropBaseRate calciumOffTropQ10)

unitsModel :: Monad m => U.ModelBuilderT m ()
unitsModel = do
  R.runReactionBuilderInUnitBuilder reactionModel

model = B.buildModel $ do
  U.unitsToCore uSecond unitsModel
