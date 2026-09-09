#ifndef TsSourcePencilBeamScanning_hh
#define TsSourcePencilBeamScanning_hh

#include "TsSource.hh"
#include "TsPBSBeamModel.hh"
#include "TsPBSCoordinateModel.hh"
#include "TsPBSSpotPlan.hh"

#include <vector>

struct TsPBSPreparedSpot {
	TsPBSSpot spot;
	G4long histories;
	TsPBSBeamOptics optics;
	TsPBSSourceRay ray;
};

class TsSourcePencilBeamScanning : public TsSource
{
public:
	TsSourcePencilBeamScanning(TsParameterManager* pM, TsSourceManager* psM, G4String sourceName);
	~TsSourcePencilBeamScanning();

	void ResolveParameters();

	const std::vector<TsPBSPreparedSpot>& PreparedSpots() const { return fPreparedSpots; }
	G4long TotalHistories() const { return fTotalHistories; }

	// Map a run-global history / Geant4 event ID onto the owning spot.
	// Thread-safe: read-only after ResolveParameters. Returns nullptr if out of range.
	const TsPBSPreparedSpot* SpotForHistory(G4long historyIndex) const;

private:
	TsPBSSpotPlan fPlan;
	TsPBSBeamModel fBeamModel;
	std::vector<TsPBSPreparedSpot> fPreparedSpots;
	std::vector<G4long> fHistoryBegin;
	G4long fTotalHistories;
};

#endif
