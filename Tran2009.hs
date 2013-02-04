{-# LANGUAGE NoMonomorphismRestriction,DeriveDataTypeable,TemplateHaskell #-}
-- +Require ModML-Units
-- +Require ModML-Reactions
-- +Require typehash
-- +Require containers
-- +Require mtl

module Tran2009
where
  
import qualified Rice2008 as Rice
import qualified ModML.Units.UnitsDAEModel as U
import qualified ModML.Core.BasicDAEModel as B
import ModML.Units.UnitsDAEOpAliases
import qualified ModML.Reactions.Reactions as R
import qualified Data.Data as D
import qualified Data.TypeHash as D
import ModML.Units.SIUnits
import qualified Control.Monad as M

R.declareNamedTaggedEntity [e|Rice.uConcentrationR|] "Magnesium-ATP" "mgATP"
R.declareNamedTaggedEntity [e|Rice.uConcentrationR|] "Magnesium-ADP" "mgADP"
R.declareNamedTaggedEntity [e|Rice.uConcentrationR|] "Inorganic Phosphate" "pinorganic"
R.declareNamedTaggedEntity [e|Rice.uConcentrationR|] "Proton" "proton"

data Parameters = Parameters { riceParameters :: Rice.Parameters, 
                               protonDissociationConstant :: Rice.RExB,      -- kdHCa
                               protonCooperativity :: Rice.RExB,             -- m
                               adpDissociationConstant :: Rice.RExB,         -- kdADP
                               physiologicalInorganicPhosphate :: Rice.RExB, -- Pi'
                               physiologicalMgADP :: Rice.RExB,              -- MgADP'
                               referenceProtonConcentration :: Rice.RExB     -- H'
                             }

defaultParameters = Parameters { riceParameters = Rice.defaultParameters,
                                 protonDissociationConstant = U.realConstant Rice.uConcentration 2E-5,
                                 protonCooperativity = U.dConstant 1,
                                 adpDissociationConstant = U.realConstant Rice.uConcentration 4E-3,
                                 physiologicalInorganicPhosphate = U.realConstant Rice.uConcentration 2,
                                 referenceProtonConcentration = U.realConstant Rice.uConcentration 2
                               }

model = parameterisedModel defaultParameters
parameterisedModel p = B.buildModel $ do
  U.unitsToCore uSecond (unitsModel p)
unitsModel p@(Parameters { riceParameters = pr }) = do
  fromPre <- unitsModelBeforeReaction p
  (cem, _, ces) <-
    R.runReactionBuilderInUnitBuilder' (do
                                           Rice.reactionModelWithCalciumTransient fromPre p
                                           reactionModel p
                                       )
  unitsModelAfterReaction p fromPre cem ces
unitsModelBeforeReaction (Parameters { riceParameters = p }) = Rice.unitsModelBeforeReaction p
unitsModelAfterReaction (Parameters { riceParameters = p }) fromPre ceMap fromReaction =
  Rice.unitsModelAfterReaction p fromPre ceMap fromReaction

reactionModel p = do
  -- Take out processes that have changed from the Rice et. al. model, so we can add replacement processes.
  R.removeAllCompartmentProcessesInvolving [[(Rice.xbPreR, -1), (Rice.tmPXB, 1)],
                                            [(Rice.xbPostR, -1), (Rice.xbPreR, 1)],
                                            [(Rice.xbPostR, -1), (Rice.tmPXB, 1)]]
  R.newAllCompartmentProcesses [xbPreRToPermissive p, crossBridgePostToPre p,
                                xbPostRToPermissive p]

xbPreRToPermissive (p@Parameters { riceParameters = rp }) c = do
  xbvar <- R.addEntity R.EssentialForProcess R.CantBeCreatedByProcess R.ModifiedByProcess (-1) (Rice.xbPreR `R.withCompartment` c)
  pvar <- R.addEntity R.NotEssentialForProcess R.CanBeCreatedByProcess R.ModifiedByProcess 1 (Rice.tmPXB `R.withCompartment` c)
  pivar <- R.addEntity R.EssentialForProcess R.CantBeCreatedByProcess R.NotModifiedByProcess 0 (pinorganic `R.withCompartment` c)
  R.rateEquation $ Rice.xbPreRToPermissiveRate rp .*. (pivar ./. physiologicalInorganicPhosphate p) .*. xbvar
  
crossBridgeReverseRotationRateProtonSensitive (p@Parameters { riceParameters = rp}) =
  Rice.standardRate
    ((Rice.crossBridgeReverseRotation rp){
        otherMod=Just (referenceProtonConcentration rp)})
    (Rice.temperature rp)

crossBridgePostToPre (p@Parameters { riceParameters = rp }) c = do
  xbPreRvar <- R.addEntity R.NotEssentialForProcess R.CanBeCreatedByProcess R.ModifiedByProcess 1 (xbPreR `R.withCompartment` c)
  xbPostRvar <- R.addEntity R.EssentialForProcess R.CantBeCreatedByProcess R.ModifiedByProcess (-1) (xbPostR `R.withCompartment` c)
  protonVar <- R.addEntity R.EssentialForProcess R.CantBeCreatedByProcess R.NotModifiedByProcess 0 (proton `R.withCompartment` c)
  mgADPVar <- R.addEntity R.EssentialForProcess R.CantBeCreatedByProcess R.NotModifiedByProcess 0 (mgADP `R.withCompartment` c)
  R.rateEquation $ Rice.standardRate (crossBridgeReverseRotationRateProtonSensitive p) (Rice.temperature rp) .*.
                     xbPostRvar .*. protonVar .*.
                     ((adpDissociationConstant p .+. physiologicalMgADP p) ./. physiologicalMgADP p) .*.
                     mgADPVar ./. (adpDissociationConstant p .+. mgADPVar)

xbPostRToPermissive (p@Parameters { riceParameters = rp}) c = do
  xbPreRvar <- R.addEntity R.EssentialForProcess R.CantBeCreatedByProcess R.ModifiedByProcess (-1) (Rice.xbPreR `R.withCompartment` c)
  pvar <- R.addEntity R.NotEssentialForProcess R.CanBeCreatedByProcess R.ModifiedByProcess 1 (Rice.tmPXB `R.withCompartment` c)
  mgADPVar <- R.addEntity R.EssentialForProcess R.CantBeCreatedByProcess R.NotModifiedByProcess 0 (mgADP `R.withCompartment` c)
  R.rateEquation $ (R.xbPostRToPermissiveRate rp) .*. xbPostRvar .*. (mgADPVar ./. physiologicalMgADP p)  .*. ((adpDissociationConstant p .+. physiologicalMgADP p) ./. (adpDissociationConstant p .+. mgADPVar))
  
permissiveToXBPostRBaseRate (p@Parameters { riceParameters = rp}) =
  adpDissociationConstant p .*. (Rice.baseRate (Rice.crossBridgeFormation rp)) .*.
  Rice.baseRate (Rice.crossBridgeRotation rp) .*.
  Rice.baseRate (Rice.rotatedCrossBridgeDissociation rp) .*. referenceProtonConcentration p ./.
  (Rice.baseRate (Rice.crossBridgeDissociation rp) .*.
   Rice.baseRate (Rice.crossBridgeReverseRotation rp) .*.
   U.expX (-(mgATPFreeEnergyChange p) ./. ())
  )