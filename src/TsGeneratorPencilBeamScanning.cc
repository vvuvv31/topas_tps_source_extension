// Particle Generator for PencilBeamScanning
//
// MT-safe: each Geant4 event ID is a unique history in [0, N). The owning
// spot is looked up from a read-only prefix table on the source. Workers
// never share mutable scheduler state.

#include "TsGeneratorPencilBeamScanning.hh"

#include "TsParameterManager.hh"
#include "TsSourcePencilBeamScanning.hh"

#include "G4Event.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#include <cmath>

TsGeneratorPencilBeamScanning::TsGeneratorPencilBeamScanning(TsParameterManager* pM, TsGeometryManager* gM,
	TsGeneratorManager* pgM, G4String sourceName)
: TsVGenerator(pM, gM, pgM, sourceName), fPBS(0)
{
	ResolveParameters();
}

TsGeneratorPencilBeamScanning::~TsGeneratorPencilBeamScanning()
{}

void TsGeneratorPencilBeamScanning::ResolveParameters()
{
	TsVGenerator::ResolveParameters();
	fPBS = dynamic_cast<TsSourcePencilBeamScanning*>(GetSource());
	if (!fPBS) {
		G4cerr << "Topas is exiting due to a serious error in source " << fSourceName << G4endl;
		G4cerr << "PencilBeamScanning generator could not find its matching source." << G4endl;
		fPm->AbortSession(1);
	}
}

void TsGeneratorPencilBeamScanning::UpdateForNewRun(G4bool rebuiltSomeComponents)
{
	TsVGenerator::UpdateForNewRun(rebuiltSomeComponents);
}

void TsGeneratorPencilBeamScanning::GeneratePrimaries(G4Event* anEvent)
{
	if (CurrentSourceHasGeneratedEnough())
		return;

	if (!fPBS)
		return;

	const G4long histIndex = anEvent->GetEventID();
	const TsPBSPreparedSpot* prepared = fPBS->SpotForHistory(histIndex);
	if (!prepared)
		return;

	SamplePrimary(*prepared, anEvent);
}

void TsGeneratorPencilBeamScanning::SamplePrimary(const TsPBSPreparedSpot& prepared, G4Event* anEvent)
{
	const TsPBSBeamOptics& optics = prepared.optics;
	const TsPBSSourceRay& ray = prepared.ray;

	const G4double ux = G4RandGauss::shoot();
	const G4double vx = G4RandGauss::shoot();
	const G4double uy = G4RandGauss::shoot();
	const G4double vy = G4RandGauss::shoot();

	const G4double dxMm = optics.sigmaXMm * ux;
	const G4double dyMm = optics.sigmaYMm * uy;
	const G4double xprime = optics.sigmaXpRad * (optics.corrX * ux + vx * std::sqrt(1.0 - optics.corrX * optics.corrX));
	const G4double yprime = optics.sigmaYpRad * (optics.corrY * uy + vy * std::sqrt(1.0 - optics.corrY * optics.corrY));

	double dirX = 0.;
	double dirY = 0.;
	double dirZ = 1.;
	TsPBSCoordinateModel::DirectionFromAngles(ray, ray.thetaXRad + xprime, ray.thetaYRad + yprime, dirX, dirY, dirZ);

	const G4double meanEnergy = prepared.spot.energyMeV * MeV;
	const G4double energySigma = meanEnergy * optics.energySpreadPercent / 100.;
	G4double energy = meanEnergy;
	if (energySigma > 0.) {
		do {
			energy = G4RandGauss::shoot(meanEnergy, energySigma);
		} while (energy <= 0.);
	}

	TsPrimaryParticle p;
	if (ray.convention == TsPBSSpotConvention::TPS) {
		p.posX = (ray.xSrcMm + dxMm) * mm;
		p.posY = ray.ySrcMm * mm;
		p.posZ = (ray.zSrcMm + dyMm) * mm;
	} else {
		p.posX = (ray.xSrcMm + dxMm) * mm;
		p.posY = (ray.ySrcMm + dyMm) * mm;
		p.posZ = ray.zSrcMm * mm;
	}
	p.dCos1 = dirX;
	p.dCos2 = dirY;
	p.dCos3 = dirZ;
	p.kEnergy = energy;
	p.weight = 1.;
	p.isNewHistory = true;
	SetParticleType(p);

	TransformPrimaryForComponent(&p);
	GenerateOnePrimary(anEvent, p);
	AddPrimariesToEvent(anEvent);
}
