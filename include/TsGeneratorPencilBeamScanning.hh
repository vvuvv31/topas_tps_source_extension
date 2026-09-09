#ifndef TsGeneratorPencilBeamScanning_hh
#define TsGeneratorPencilBeamScanning_hh

#include "TsVGenerator.hh"
#include "TsSourcePencilBeamScanning.hh"

class TsGeneratorPencilBeamScanning : public TsVGenerator
{
public:
	TsGeneratorPencilBeamScanning(TsParameterManager* pM, TsGeometryManager* gM,
		TsGeneratorManager* pgM, G4String sourceName);
	~TsGeneratorPencilBeamScanning();

	void ResolveParameters();
	void UpdateForNewRun(G4bool rebuiltSomeComponents);
	void GeneratePrimaries(G4Event* anEvent);

private:
	void SamplePrimary(const TsPBSPreparedSpot& prepared, G4Event* anEvent);

	TsSourcePencilBeamScanning* fPBS;
};

#endif
