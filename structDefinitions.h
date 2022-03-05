#ifndef STRUCTDEFINITIONS_H
#define STRUCTDEFINITIONS_H

typedef struct nLines_config
{
	int nLines_inputVectors, nLines_vwdSize, nLines_freeVolumeconfig, nLines_HBondAtoms;
} NLINES_CONFIG;

typedef struct distVar
{
	float maxDist, binSize_OOP, binSize_deg, binSize_dist;
	float binStart_dist, binEnd_dist, binStart_OOP, binEnd_OOP, binStart_deg, binEnd_deg;
	int nBins_dist, nBins_OOP, nBins_deg, size_degrees, size_oop;

	// Optional element to add, for easy access
	int nElements;
} DIST_VAR;

typedef struct orderParameter
{
	int atom1, atom2, atom3, atom4;
	float distance, theta_rad, theta_deg, orderParameter;
} ORDERPARAMETER;

typedef struct dumpinfo
{
	int timestep, nAtoms;
	float xlo, xhi, ylo, yhi, zlo, zhi;
} DUMPFILE_INFO;

typedef struct config
{
	int atom1, atom2;
	float radius;
} CONFIG;

typedef struct datafile_atoms
{
	int resNumber, ix, iy, iz;
	char resName[6], atomName[6], atomType2[6], molName[6];

	int id, molType, atomType;
	float charge, x, y, z;
} DATA_ATOMS;

typedef struct datafile_bonds
{
	int id, bondType, atom1, atom2, atom1Type, atom2Type;
	float x1, y1, z1, x2, y2, z2, xc, yc, zc;
} DATA_BONDS;

typedef struct datafile_angles
{
	int id, angleType, atom1, atom2, atom3;
} DATA_ANGLES;

typedef struct datafile_dihedrals
{
	int id, dihedralType, atom1, atom2, atom3, atom4;
} DATA_DIHEDRALS;

typedef struct datafile_impropers
{
	int id, improperType, atom1, atom2, atom3, atom4;
} DATA_IMPROPERS;

typedef struct datafileInfo
{
	int nAtoms, nBonds, nAngles, nDihedrals, nImpropers;
	int nAtomTypes, nBondTypes, nAngleTypes, nDihedralTypes, nImproperTypes;
} DATAFILE_INFO;

typedef struct distribution
{
	float binStart_OOP, binEnd_OOP, binStart_dist, binEnd_dist, binStart_deg, binEnd_deg;
	int count;
} DISTRIBUTION;

typedef struct freeVolumeVars
{
	float minProbeSize, maxProbeSize, delProbeSize, currentProbeSize;
	int nBins_probeSweep;
	float xLength, yLength, zLength;
	int nBins_dist_x, nBins_dist_y, nBins_dist_z;
	float delDistance;

	float binStart_dist, binEnd_dist, binSize_dist;
	int nBins_dist;
} FREEVOLUME_VARS;

typedef struct freeVolumeDistribution
{
	float binStart_dist, binEnd_dist;
	int nUnoccupied, nOccupied, atomType;
} FREEVOLUME_DISTRIBUTION;

typedef struct bounds
{
	float lo, hi;
} BOUNDS;

typedef struct HBondNetwork
{
	int bondPresent, bondAbsent;
} HBONDNETWORK;

typedef struct computeMSDvars
{
	float lowerBound, upperBound;
} MSD_VARS;

typedef struct allData_bondRDF_ACF
{
	float x1, y1, z1, x2, y2, z2, x3, y3, z3, distance;
} ALL_DATA_BONDRDF_ACF;

typedef struct openLogFiles_bondRDF_ACF
{
	int currentTrajCount, nLines;
	long long int nLinesTotal;
} LOGFILES_VARIABLES_BONDRDF_ACF;

#endif