// Particle Source for PencilBeamScanning

#include "TsSourcePencilBeamScanning.hh"

#include "TsParameterManager.hh"

#include "G4SystemOfUnits.hh"

#include <algorithm>
#include <cmath>
#include <initializer_list>
#include <sstream>
#include <stdexcept>

TsSourcePencilBeamScanning::TsSourcePencilBeamScanning(TsParameterManager* pM, TsSourceManager* psM, G4String sourceName)
: TsSource(pM, psM, sourceName), fTotalHistories(0)
{
	ResolveParameters();
}

TsSourcePencilBeamScanning::~TsSourcePencilBeamScanning()
{}

void TsSourcePencilBeamScanning::ResolveParameters()
{
	TsSource::ResolveParameters();
	fPreparedSpots.clear();
	fHistoryBegin.clear();
	fTotalHistories = 0;

	try {
		const G4String planFile = fPm->GetStringParameter(GetFullParmName("SpotPlanFile"));
		const G4String modelFile = fPm->GetStringParameter(GetFullParmName("BeamModelFile"));

		G4bool skipZero = true;
		if (fPm->ParameterExists(GetFullParmName("SkipZeroWeightSpots")))
			skipZero = fPm->GetBooleanParameter(GetFullParmName("SkipZeroWeightSpots"));

		G4bool interpolate = false;
		if (fPm->ParameterExists(GetFullParmName("InterpolateBeamModel")))
			interpolate = fPm->GetBooleanParameter(GetFullParmName("InterpolateBeamModel"));

		G4double historiesScale = 1.0;
		if (fPm->ParameterExists(GetFullParmName("HistoriesScale")))
			historiesScale = fPm->GetUnitlessParameter(GetFullParmName("HistoriesScale"));
		if (historiesScale < 0.)
			throw std::runtime_error("HistoriesScale must be >= 0");

		G4String weightMode = "Histories";
		if (fPm->ParameterExists(GetFullParmName("WeightMode")))
			weightMode = fPm->GetStringParameter(GetFullParmName("WeightMode"));
		G4String weightModeLower = weightMode;
		G4StrUtil::to_lower(weightModeLower);
		if (weightModeLower != "histories")
			throw std::runtime_error("only WeightMode=Histories is supported");

		G4String conventionName = "TPS";
		if (fPm->ParameterExists(GetFullParmName("SpotCoordinateConvention")))
			conventionName = fPm->GetStringParameter(GetFullParmName("SpotCoordinateConvention"));
		else if (fPm->ParameterExists(GetFullParmName("CoordinateConvention")))
			conventionName = fPm->GetStringParameter(GetFullParmName("CoordinateConvention"));
		const TsPBSSpotConvention convention = TsPBSCoordinateModel::ParseConvention(std::string(conventionName));

		auto firstExisting = [this](std::initializer_list<const char*> names) -> G4String {
			for (const char* name : names) {
				const G4String full = GetFullParmName(name);
				if (fPm->ParameterExists(full))
					return full;
			}
			return "";
		};

		const G4String vsadXName = firstExisting({"VirtualScanningMagneticX", "VirtualSADX"});
		const G4String vsadYName = firstExisting({"VirtualScanningMagneticY", "VirtualSADY"});
		if (vsadXName.empty() || vsadYName.empty())
			throw std::runtime_error("VirtualScanningMagneticX and VirtualScanningMagneticY are required");
		const G4double vsadX = fPm->GetDoubleParameter(vsadXName, "Length") / mm;
		const G4double vsadY = fPm->GetDoubleParameter(vsadYName, "Length") / mm;

		const G4String sadName = firstExisting({
			"VirtualSourceToIsocenterDistance", "SAD", "SourceToIsocenterDistance"});
		if (sadName.empty())
			throw std::runtime_error("VirtualSourceToIsocenterDistance is required");
		const G4double sad = fPm->GetDoubleParameter(sadName, "Length") / mm;

		fPlan.Load(std::string(planFile), skipZero);
		fBeamModel.Load(std::string(modelFile));

		G4int firstSpot = 0;
		G4int lastSpot = -1;
		if (fPm->ParameterExists(GetFullParmName("FirstSpot")))
			firstSpot = fPm->GetIntegerParameter(GetFullParmName("FirstSpot"));
		if (fPm->ParameterExists(GetFullParmName("LastSpot")))
			lastSpot = fPm->GetIntegerParameter(GetFullParmName("LastSpot"));
		if (firstSpot < 0)
			throw std::runtime_error("FirstSpot must be >= 0");
		if (static_cast<std::size_t>(firstSpot) >= fPlan.Size())
			throw std::runtime_error("FirstSpot is past the end of the loaded spot plan");
		std::size_t lastIndex = fPlan.Size() - 1;
		if (lastSpot >= 0)
			lastIndex = static_cast<std::size_t>(lastSpot);
		if (lastIndex >= fPlan.Size())
			throw std::runtime_error("LastSpot is past the end of the loaded spot plan");
		if (static_cast<std::size_t>(firstSpot) > lastIndex)
			throw std::runtime_error("FirstSpot is greater than LastSpot");

		for (std::size_t i = static_cast<std::size_t>(firstSpot); i <= lastIndex; ++i) {
			const TsPBSSpot& spot = fPlan.At(i);
			const G4long histories = static_cast<G4long>(std::llround(spot.weight * historiesScale));
			if (histories <= 0)
				continue;

			TsPBSPreparedSpot prepared;
			prepared.spot = spot;
			prepared.histories = histories;
			prepared.optics = fBeamModel.Lookup(spot.energyMeV, interpolate);
			prepared.ray = TsPBSCoordinateModel::Compute(
				spot.xIsoMm, spot.yIsoMm, vsadX, vsadY, sad, convention);
			fPreparedSpots.push_back(prepared);
			fTotalHistories += histories;
		}

		if (fPreparedSpots.empty() || fTotalHistories <= 0)
			throw std::runtime_error("all selected spots have zero histories");
		if (fTotalHistories > 1000000000L)
			throw std::runtime_error("total histories exceed 1e9");

		fHistoryBegin.resize(fPreparedSpots.size() + 1);
		fHistoryBegin[0] = 0;
		for (std::size_t i = 0; i < fPreparedSpots.size(); ++i)
			fHistoryBegin[i + 1] = fHistoryBegin[i] + fPreparedSpots[i].histories;

		fNumberOfHistoriesInRun = fTotalHistories;

		G4cout << "PencilBeamScanning source " << fSourceName
			<< ": " << fPreparedSpots.size() << " spots, "
			<< fTotalHistories << " histories"
			<< " (skipped " << fPlan.ZeroWeightCount() << " zero-weight rows)"
			<< G4endl;
		G4cout << "  plan: " << fPlan.FileName() << G4endl;
		G4cout << "  beam model: " << fBeamModel.FileName() << G4endl;
		G4cout << "  SpotCoordinateConvention = " << conventionName
			<< (fPm->ParameterExists(GetFullParmName("SpotCoordinateConvention")) ||
				fPm->ParameterExists(GetFullParmName("CoordinateConvention"))
				? "" : " (default)") << G4endl;
		G4cout << "  VirtualScanningMagnetic X/Y = " << vsadX << " / " << vsadY
			<< " mm, VirtualSourceToIsocenterDistance = " << sad << " mm" << G4endl;
		G4cout << "  History-to-spot map is event-ID based (MT-safe)." << G4endl;
	} catch (const std::exception& exc) {
		G4cerr << "Topas is exiting due to a serious error in source " << fSourceName << G4endl;
		G4cerr << exc.what() << G4endl;
		fPm->AbortSession(1);
	}
}

const TsPBSPreparedSpot* TsSourcePencilBeamScanning::SpotForHistory(G4long historyIndex) const
{
	if (historyIndex < 0 || fHistoryBegin.size() < 2)
		return nullptr;
	if (historyIndex >= fHistoryBegin.back())
		return nullptr;

	auto it = std::upper_bound(fHistoryBegin.begin(), fHistoryBegin.end(), historyIndex);
	const std::size_t spotIndex = static_cast<std::size_t>(std::distance(fHistoryBegin.begin(), it) - 1);
	return &fPreparedSpots[spotIndex];
}
