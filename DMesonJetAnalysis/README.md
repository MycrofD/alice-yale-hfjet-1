# How to run the analysis

## Detector response

Start with generating the detector response and reconstruction efficiency:

`./runResponseAnalysis.py LHC15i2response_Train1399.yaml `

When this is done, it should be repeated using results from the previous to apply the reconstruction efficiency:

`./runResponseAnalysis.py LHC15i2response_Train1399_efficiency.yaml`

## Reflection templates

Generate reflection templates without reconstruction efficiency:

`./runDataAnalysis.py LHC15i2reflections_Train1399.yaml`

`cd ../rawYieldUnc`

`./GenerateReflTemplate.py LHC15i2genreflections_Train1399.yaml`

and with reconstruction efficiency:

`./runDataAnalysis.py LHC15i2reflections_Train1399_efficiency.yaml`

`cd ../rawYieldUnc`

`./GenerateReflTemplate.py LHC15i2genreflections_Train1399_efficiency.yaml`

## Reconstruct data

Analyze data without reconstruction efficiency:

`./runDataAnalysis.py LHC10analysis_Train1116.yaml`

and with reconstruction efficiency:

`./runDataAnalysis.py LHC10analysis_Train1116_efficiency.yaml`

## Raw yield uncertainty

Go in `rawYieldUnc` folder:

`cd ../rawYieldUnc`.

Then execute:

`./ExtractDZeroJetRawYieldUncertainty.py ../DMesonJetAnalysis/LHC10analysis_Train1116.yaml --refl DoubleGaus -b`

`./ExtractDZeroJetRawYieldReflectionUncertainty.py ../DMesonJetAnalysis/LHC10analysis_Train1116.yaml -b` 

`./ExtractDZeroJetRawYieldUncertainty.py ../DMesonJetAnalysis/LHC10analysis_Train1116.yaml --no-refl -b`

`./ExtractDZeroJetRawYieldUncertainty.py ../DMesonJetAnalysis/LHC10analysis_Train1116_efficiency.yaml --refl DoubleGaus -b`

`./ExtractDZeroJetRawYieldReflectionUncertainty.py ../DMesonJetAnalysis/LHC10analysis_Train1116_efficiency.yaml -b`

`./ExtractDZeroJetRawYieldUncertainty.py ../DMesonJetAnalysis/LHC10analysis_Train1116_efficiency.yaml --no-refl -b`

## B Feed-Down

` ./MCSimulationSystematics.py FDCorrection_EvtGen.yaml `

`./MCSimulationSystematics.py FDCorrection_EvtGen_1416.yaml`

## Unfolding

`./runDataUnfolding.py Unfolding_LHC10_Train1116_LHC15i2_Train1399_efficiency.yaml -b`

Note: no need to run the above if you run the full systematics

### Full Systematics

##### Tracking Efficiency

 `./runDataUnfolding.py Unfolding_LHC10_Train1116_LHC15i2_Train1416_efficiency.yaml -b`

##### B Feed-Down and Raw Yield

 `./runDataUnfolding.py Unfolding_LHC10_Train1116_LHC15i2_Train1399_efficiency.yaml --fd-syst --ry-syst -b`

##### Reflections

`./runDataUnfolding.py Unfolding_LHC10_Train1116_LHC15i2_Train1399_efficiency.yaml --refl-ros 5 -b`

`./runDataUnfolding.py Unfolding_LHC10_Train1116_LHC15i2_Train1399_efficiency.yaml --refl-ros 15 -b`