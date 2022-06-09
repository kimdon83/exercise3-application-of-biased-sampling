#include "define.h"


int main(void){
	Cmcmanager mcman;
	mcman.readINPUT();
	mcman.mcrunNonBaiased();
	mcman.mcrunStratification();
	mcman.mcrunBaiasedLinear();
	mcman.mcrunBaiasedExponential();
	mcman.mcrunBaiasedConst();
	
	system("pause");
}