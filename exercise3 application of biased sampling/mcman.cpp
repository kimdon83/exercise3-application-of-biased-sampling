#include "mcman.h"
#include "define.h"

Cmcmanager::Cmcmanager(){
	Num_history = 100; set = 100;

	q = 0.;	qMean = 0.;	qVar = 0.;
	qSum = 0.; qSqu = 0.; qMeanSqu = 0.; qMeanSum = 0.; qMeanVar = 0.;

	q2 = 0.;	qMean2 = 0.;	qVar2 = 0.;
	qSum2 = 0.; qSqu2 = 0.; qMeanSqu2 = 0.; qMeanSum2 = 0.; qMeanVar2 = 0.;

	Range_Ratio = 0.5;
}
Cmcmanager::~Cmcmanager(){


}

int Cmcmanager::readINPUT(void){
	cout << "------------------------------------------------------------------" << endl;
	cout << "|       Exercise #3. application of biased sampling              |" << endl;
	cout << "------------------------------------------------------------------" << endl;
	cout << "      input the number of history " << endl;
	//cin >> Num_history;
	Num_history = 100000;
	//cin >> set;
	set = 100;
	return 0;
}
int Cmcmanager::mcrunNonBaiased(void){
	cout << "1. Non-Baiased Sampling" << endl;

	for (j = 0; j < Num_history; j++){
		theta = 2 * _pi*RN.gnRN();
		q = 2 * _pi*sin(theta) / theta;
		qSum += q;
		qSqu += q*q;
	}
	qMean = qSum / Num_history;
	qVar = (-qMean*qMean + (qSqu) / (Num_history));
// cout the result
	cout << "Mean : " << qMean << endl;
	cout << "Variance : " << qVar << endl;
/////
	for (i = 0; i < set; i++){
		qMean = 0;
		qSum = 0;
		qSqu = 0.;
		for (j = 0; j < Num_history; j++){
			theta = 2 * _pi*RN.gnRN();
			q = 2 * _pi*sin(theta) / theta;
			qSum += q;
			qSqu += q*q;
		}
		qMean = qSum / Num_history;
		qMeanSum += qMean;
		qMeanSqu += qMean*qMean;
	}
	qMeanSum /= set;
	qMeanVar = ((qMeanSqu) / set - qMeanSum*qMeanSum);// / (set ;
	
// cout the result
	cout << "Variance of mean : " << qMeanVar<< endl;
	cout << endl;
	return 0;
}

int Cmcmanager::mcrunStratification(void){
	cout << "2. Stratified Sampling" << endl;
	q = 0.;	qMean = 0.;	qVar = 0.;
	qSum = 0.; qSqu = 0.; qMeanSqu = 0.; qMeanSum = 0.; qMeanVar = 0.;

	q2 = 0.;	qMean2 = 0.;	qVar2 = 0.;
	qSum2 = 0.; qSqu2 = 0.; qMeanSqu2 = 0.; qMeanSum2 = 0.; qMeanVar2 = 0.;
//q1 : theta (0:pi)	

	Num = Num_history*Range_Ratio;
	for (j = 0; j < Num; j++){
		theta = _pi*RN.gnRN();
		q = _pi*sin(theta) / theta;
		qSum += q;
		qSqu += q*q;
	}
	qMean = qSum / Num;
	qVar = (-qMean*qMean + (qSqu) / (Num));
	/////
	for (i = 0; i < set; i++){
		qMean = 0;
		qSum = 0;
		qSqu = 0.;
		for (j = 0; j < Num; j++){
			theta = _pi*RN.gnRN();
			q = _pi*sin(theta) / theta;
			qSum += q;
			qSqu += q*q;
		}
		qMean = qSum / Num;
		qMeanSum += qMean;
		qMeanSqu += qMean*qMean;
	}
	qMeanSum /= set;
	qMeanVar = ((qMeanSqu) / set - qMeanSum*qMeanSum);// / (set ;

//q2 : theta (pi:2pi)
	Num = Num_history*(1 - Range_Ratio);
	for (j = 0; j < Num; j++){
		theta = _pi*(1+RN.gnRN());
		q = _pi*sin(theta) / theta;
		qSum2 += q;
		qSqu2 += q*q;
	}
	qMean2 = qSum2 / Num;
	qVar2 = (-qMean2*qMean2 + (qSqu2) / (Num));
	/////
	for (i = 0; i < set; i++){
		qMean2 = 0;
		qSum2 = 0;
		qSqu2 = 0.;
		for (j = 0; j < Num; j++){
			theta = _pi*(1+RN.gnRN());
			q = _pi*sin(theta) / theta;
			qSum2 += q;
			qSqu2 += q*q;
		}
		qMean2 = qSum2 / Num;
		qMeanSum2 += qMean2;
		qMeanSqu2 += qMean2*qMean2;
	}
	qMeanSum2 /= set;
	qMeanVar2 = ((qMeanSqu2) / set - qMeanSum2*qMeanSum2);// / (set ;

	// cout the result
	cout << "Mean : " <<  qMean+qMean2 << endl;
	cout << "Variance : " << qVar+ qVar2<< endl;
	//cout << "Variance of mean : " << qMeanVar / Num_history / Range_Ratio + qMeanVar2 / Num_history / (1 - Range_Ratio) << endl;
	cout << "Variance of mean : " << qMeanVar+ qMeanVar2 << endl;
	cout << endl;
	return 0;
}

int Cmcmanager ::mcrunBaiasedLinear(void){
	cout << "3. Baiased Sampling Linear" << endl;
	q = 0.;	qMean = 0.;	qVar = 0.;
	qSum = 0.; qSqu = 0.; qMeanSqu = 0.; qMeanSum = 0.; qMeanVar = 0.;
	double a0, b0;
	double a1, b1;
	double a2, b2;
	double a3, b3;
	a0 = 1.13004;	b0 = -0.3469;
	a1 = 0.20013;	b1 = -0.01463;
	a2 = 0.502156;	b2 = -0.15415;
	a3 = 0.088932;	b3 = -0.0065;

	bool bReject;

	double rn1, rn2;	//random number or rejection technic
	for (j = 0; j < Num_history; j++){
		do{
			rn1 = 2*_pi*RN.gnRN();
			rn2 = 1.2 *RN.gnRN();
			if (rn1 <= _pi) bReject = (a0 + b0*rn1) <= rn2;
			else bReject = (a1 + b1*rn1) <= rn2;
		} while (bReject);
		theta = rn1;
		if (theta<_pi) q = sin(theta) / theta/(a2+b2*theta);
		else q = sin(theta) / theta / (a3 + b3*theta);

		qSum += q;
		qSqu += q*q;
	}
	qMean = qSum / Num_history;
	qVar = (-qMean*qMean + (qSqu) / (Num_history));
	// cout the result
	cout << "Mean : " << qMean << endl;
	cout << "Variance : " << qVar << endl;
	/////
	for (i = 0; i < set; i++){
		qMean = 0;
		qSum = 0;
		qSqu = 0.;
		for (j = 0; j < Num_history; j++){
			do{
				rn1 = 2 * _pi*RN.gnRN();
				rn2 = 1.2 *RN.gnRN();
				if (rn1 <= _pi) bReject = (a0 + b0*rn1) <= rn2;
				else bReject = (a1 + b1*rn1) <= rn2;
			} while (bReject);
			theta = rn1;
			if (theta<_pi) q = sin(theta) / theta / (a2 + b2*theta);
			else q = sin(theta) / theta / (a3 + b3*theta);

			qSum += q;
			qSqu += q*q;
		}
		qMean = qSum / Num_history;
		qMeanSum += qMean;
		qMeanSqu += qMean*qMean;
	}
	qMeanSum /= set;
	qMeanVar = ((qMeanSqu) / set - qMeanSum*qMeanSum);// / (set ;

	// cout the result
	cout << "Variance of mean : " << qMeanVar << endl;
	cout << endl;
	return 0;
}
int Cmcmanager::mcrunBaiasedExponential(void){
	cout << "4. Baiased Sampling Exponential" << endl;
	q = 0.;	qMean = 0.;	qVar = 0.;
	qSum = 0.; qSqu = 0.; qMeanSqu = 0.; qMeanSum = 0.; qMeanVar = 0.;
	double y0[4], A[4],R0[4];

	y0[0] = 1.79214;  A[0] = -0.73426; R0[0] = 0.2921;
	y0[1] = -0.80215; A[1] = 1; R0[1] = -0.01466;

	y0[2] = 0.791897;  A[2] = -0.32445; R0[2] = 0.129071;
	y0[3] = -0.354448; A[3] = 0.441872; R0[3] = -0.006477846;


	bool bReject;

	double rn1, rn2;	//random number or rejection technic
	for (j = 0; j < Num_history; j++){
		do{
			rn1 = 2 * _pi*RN.gnRN();
			rn2 = 1.1 *RN.gnRN();
			if (rn1 <= _pi) bReject = (y0[0] + A[0]*exp(R0[0]*rn1)) <= rn2;
			else bReject = (y0[1] + A[1] * exp(R0[1] * rn1)) <= rn2;
		} while (bReject);
		theta = rn1;
		if (theta<_pi) q = sin(theta) / theta / (y0[2] + A[2] * exp(R0[2] * theta));
		else           q = sin(theta) / theta / (y0[3] + A[3] * exp(R0[3] * theta));

		qSum += q;
		qSqu += q*q;
	}
	qMean = qSum / Num_history;
	qVar = (-qMean*qMean + (qSqu) / (Num_history));
	// cout the result
	cout << "Mean : " << qMean << endl;
	cout << "Variance : " << qVar << endl;
	/////
	for (i = 0; i < set; i++){
		qMean = 0;
		qSum = 0;
		qSqu = 0.;
		for (j = 0; j < Num_history; j++){
			do{
				rn1 = 2 * _pi*RN.gnRN();
				rn2 = 1.1 *RN.gnRN();
				if (rn1 <= _pi) bReject = (y0[0] + A[0] * exp(R0[0] * rn1)) <= rn2;
				else bReject = (y0[1] + A[1] * exp(R0[1] * rn1)) <= rn2;
			} while (bReject);
			theta = rn1;
			if (theta<_pi) q = sin(theta) / theta / (y0[2] + A[2] * exp(R0[2] * theta));
			else q = sin(theta) / theta / (y0[3] + A[3] * exp(R0[3] * theta));

			qSum += q;
			qSqu += q*q;
		}
		qMean = qSum / Num_history;
		qMeanSum += qMean;
		qMeanSqu += qMean*qMean;
	}
	qMeanSum /= set;
	qMeanVar = ((qMeanSqu) / set - qMeanSum*qMeanSum);// / (set ;

	// cout the result
	cout << "Variance of mean : " << qMeanVar << endl;
	cout << endl;
	return 0;
}
int Cmcmanager::mcrunBaiasedConst(void){
	cout << "5. Baiased Sampling Const" << endl;
	q = 0.;	qMean = 0.;	qVar = 0.;
	qSum = 0.; qSqu = 0.; qMeanSqu = 0.; qMeanSum = 0.; qMeanVar = 0.;
	double a0, b0;
	double a1, b1;
	double a2, b2;
	double a3, b3;
	a0 = 1.13004;	b0 = -0.3469;
	a1 = 0.20013;	b1 = -0.01463;
	a2 = 0.502156;	b2 = -0.15415;
	a3 = 0.088932;	b3 = -0.0065;

	bool bReject;

	double rn1, rn2;	//random number or rejection technic
	for (j = 0; j < Num_history; j++){
		do{
			rn1 = 2 * _pi*RN.gnRN();
			rn2 = 1.2 *RN.gnRN();
			if (rn1 <= _pi) bReject = (a0 + b0*rn1) <= rn2;
			else bReject = (a1 + b1*rn1) <= rn2;
		} while (bReject);
		theta = rn1;
		if (theta<_pi) q = sin(theta) / theta / (a2 + b2*theta);
		else q = sin(theta) / theta / (a3 + b3*theta);

		qSum += q;
		qSqu += q*q;
	}
	qMean = qSum / Num_history;
	qVar = (-qMean*qMean + (qSqu) / (Num_history));
	// cout the result
	cout << "Mean : " << qMean << endl;
	cout << "Variance : " << qVar << endl;
	/////
	for (i = 0; i < set; i++){
		qMean = 0;
		qSum = 0;
		qSqu = 0.;
		for (j = 0; j < Num_history; j++){
			do{
				rn1 = 2 * _pi*RN.gnRN();
				rn2 = 1.2 * _pi*RN.gnRN();
				if (rn1 <= _pi) bReject = (a0 + b0*rn1) <= rn2;
				else bReject = (a1 + b1*rn1) <= rn2;
			} while (bReject);
			theta = rn1;
			if (theta<_pi) q = sin(theta) / theta / (a2 + b2*theta);
			else q = sin(theta) / theta / (a3 + b3*theta);

			qSum += q;
			qSqu += q*q;
		}
		qMean = qSum / Num_history;
		qMeanSum += qMean;
		qMeanSqu += qMean*qMean;
	}
	qMeanSum /= set;
	qMeanVar = ((qMeanSqu) / set - qMeanSum*qMeanSum);// / (set ;

	// cout the result
	cout << "Variance of mean : " << qMeanVar << endl;
	cout << endl;
	return 0;
}