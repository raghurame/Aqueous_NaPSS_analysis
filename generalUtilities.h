#ifndef GENERALUTILITIES_H
#define GENERALUTILITIES_H

float translatePeriodic (float coord, int image, float dimension);
float translatePeriodicDistance (float coord1, float coord2, float halfBoxDistance);
float findBondCenter (float coord1, int image1, float coord2, int image2, float xlo, float xhi);
float findConnectedAtom_periodicTranslation (float coord1, int image1, float coord2, int image2, float xlo, float xhi);
void setDistributionZero (DISTRIBUTION **rawArray, int arraySize);
int getIndex1d (int i, int j, int width);
long long int getIndex1d_from3d (int x, int xWidth, int y, int yWidth, int z, int zWidth);

#endif