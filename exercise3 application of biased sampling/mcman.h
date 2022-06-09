#ifndef _MCMAN_H_KDH_
#define _MCMAN_H_KDH_

#include "define.h"

class Cmcmanager{
private:
	int i,j,k;
	double loc[2];
	rng RN;

	double theta;
	double pdf;

	int Num_history,Num;
	int set;

	double qSum, qSqu;
	double q, qMean, qVar,qMeanSum,qMeanSqu, qMeanVar;

	double qSum2, qSqu2;
	double q2, qMean2, qVar2, qMeanSum2, qMeanSqu2, qMeanVar2;
	double Range_Ratio;
public:

	Cmcmanager();
	~Cmcmanager();

	int readINPUT(void);
	int mcrunNonBaiased(void);
	int mcrunStratification(void);
	int mcrunBaiasedLinear(void);
	int mcrunBaiasedExponential(void);
	int mcrunBaiasedConst(void);


};
#endif